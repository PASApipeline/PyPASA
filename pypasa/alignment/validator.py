"""
Alignment validator for PyPASA.

Port of validate_alignments_in_db.dbi plus the structural checks Perl runs
when it builds each alignment object (CDNA_alignment.pm). An alignment is
valid if it:
- is contiguous on the transcript,
- has consensus splice sites at every intron in one orientation (which also
  sets its spliced orientation),
- meets the minimum average percent identity,
- covers the minimum percent of the transcript's length.

Alignments are updated in place (valid, spliced_orient, percent_aligned).
"""

import logging
from typing import Dict, List, Optional

from pypasa.models import Alignment, Orientation

logger = logging.getLogger(__name__)


class AlignmentValidator:
    """Validate transcript alignments against Perl PASA's criteria."""
    
    def __init__(
        self,
        genome: Dict[str, str],
        transcript_lengths: Dict[str, int],
        min_percent_id: float = 95.0,
        min_percent_aligned: float = 75.0,
        num_bp_perfect_splice: int = 0,
    ):
        """
        Args:
            genome: scaffold -> sequence (for splice sites)
            transcript_lengths: transcript accession -> length (for percent aligned)
            min_percent_id: Minimum average percent identity
            min_percent_aligned: Minimum percent of transcript aligned
            num_bp_perfect_splice: Perl's NUM_BP_PERFECT_SPLICE_BOUNDARY; only 0 is
                supported so far
        """
        if num_bp_perfect_splice:
            raise NotImplementedError("NUM_BP_PERFECT_SPLICE_BOUNDARY > 0 is not ported yet")
        self.genome = genome
        self.transcript_lengths = transcript_lengths
        self.min_percent_id = min_percent_id
        self.min_percent_aligned = min_percent_aligned
        
        self.stats = {
            'total': 0,
            'passed': 0,
            'failed_percent_id': 0,
            'failed_percent_aligned': 0,
            'failed_splice': 0,
            'failed_contiguity': 0,
        }
    
    def validate(self, alignments: List[Alignment]) -> Dict[str, int]:
        """
        Validate alignments in place.
        
        Returns:
            Dictionary of validation statistics
        """
        for alignment in alignments:
            self._validate_one(alignment)
        
        logger.info(
            f"Validation complete: {self.stats['passed']}/{self.stats['total']} passed"
        )
        return self.stats
    
    def _validate_one(self, alignment: Alignment) -> None:
        self.stats['total'] += 1
        genome_seq = self.genome[alignment.scaffold]
        
        # Structural checks Perl runs when building the alignment object
        # (CDNA_alignment.pm); either one invalidates the alignment.
        structural_errors = []
        if not self._is_contiguous(alignment):
            structural_errors.append("incontiguous alignment")
        if len(alignment.segments) > 1:
            spliced_orient = self._identify_splice_junctions(alignment, genome_seq)
            if spliced_orient:
                alignment.spliced_orient = spliced_orient
            else:
                structural_errors.append("splice site validations failed")
        
        is_valid, reason = self._validate_alignment(alignment)
        if structural_errors:
            is_valid, reason = False, "; ".join(structural_errors)
        alignment.valid = is_valid
        
        if is_valid:
            self.stats['passed'] += 1
        else:
            if 'percent_id' in reason:
                self.stats['failed_percent_id'] += 1
            elif 'aligned' in reason:
                self.stats['failed_percent_aligned'] += 1
            elif 'splice' in reason:
                self.stats['failed_splice'] += 1
            elif 'incontiguous' in reason:
                self.stats['failed_contiguity'] += 1
            logger.debug(f"Alignment {alignment.align_acc} failed: {reason}")

    def _validate_alignment(self, alignment: Alignment) -> tuple:
        """
        Check if a single alignment passes quality filters.
        
        Args:
            alignment: Alignment to validate
            
        Returns:
            Tuple of (is_valid: bool, reason: str)
        """
        # Check percent identity
        if alignment.avg_percent_id is not None:
            if alignment.avg_percent_id < self.min_percent_id:
                return (
                    False,
                    f"percent_id {alignment.avg_percent_id:.2f} < {self.min_percent_id}",
                )
        
        # Check percent aligned using ACTUAL transcript length from FASTA
        # This is critical - we must use the full cDNA length, not just the aligned span
        transcript_length = self._get_transcript_length(alignment.cdna_acc)
        if transcript_length > 0:
            percent_aligned = (alignment.aligned_length / transcript_length) * 100.0
            
            # Store percent_aligned for later use
            alignment.percent_aligned = percent_aligned
            
            if percent_aligned < self.min_percent_aligned:
                return (
                    False,
                    f"percent_aligned {percent_aligned:.2f} < {self.min_percent_aligned}",
                )
        
        return (True, "passed")
    
    def _get_transcript_length(self, cdna_acc: str) -> int:
        """Full transcript length, or 0 (with a warning) if the accession is unknown."""
        length = self.transcript_lengths.get(cdna_acc, 0)
        if not length:
            logger.warning(f"Transcript {cdna_acc} not found among transcript lengths")
        return length
    
    def get_validation_report(self) -> str:
        """
        Generate a human-readable validation report.
        
        Returns:
            Formatted report string
        """
        total = self.stats['total']
        if total == 0:
            return "No alignments validated"
        
        report = []
        report.append("=" * 60)
        report.append("ALIGNMENT VALIDATION REPORT")
        report.append("=" * 60)
        report.append(f"Total alignments:     {total:>10,}")
        report.append(f"Passed:              {self.stats['passed']:>10,} ({100.0 * self.stats['passed'] / total:>6.2f}%)")
        report.append(f"Failed:              {total - self.stats['passed']:>10,} ({100.0 * (total - self.stats['passed']) / total:>6.2f}%)")
        report.append("")
        report.append("Failure reasons:")
        report.append(f"  Low percent identity:    {self.stats['failed_percent_id']:>10,}")
        report.append(f"  Low percent aligned:     {self.stats['failed_percent_aligned']:>10,}")
        report.append(f"  Splice boundary issues:  {self.stats['failed_splice']:>10,}")
        report.append(f"  Incontiguous:            {self.stats['failed_contiguity']:>10,}")
        report.append("=" * 60)
        
        return "\n".join(report)

    # Consensus (donor, acceptor) dinucleotides on the genome's + strand, per
    # transcribed orientation. AT-AC is allowed, as in Perl's default
    # ($ALLOW_ATAC_splice_pairs=1, extended AT-AC donor check off).
    CONSENSUS_SPLICE_PAIRS = {
        Orientation.PLUS: {("GT", "AG"), ("GC", "AG"), ("AT", "AC")},
        Orientation.MINUS: {("CT", "AC"), ("CT", "GC"), ("GT", "AT")},
    }
    
    def _identify_splice_junctions(self, alignment: Alignment, genome_seq: str) -> Optional[Orientation]:
        """
        Determine transcribed orientation from splice sites.
        
        Port of CDNA_alignment::identify_splice_junctions: try the aligned
        orientation, then the opposite; the first in which every intron has a
        consensus splice pair is the spliced orientation.
        
        Returns:
            Orientation.PLUS or Orientation.MINUS, or None if neither
            orientation validates (Perl: "Splice site validations failed")
        """
        aligned = Orientation(alignment.orient)
        if aligned == Orientation.UNKNOWN:
            aligned = Orientation.PLUS
        opposite = Orientation.MINUS if aligned == Orientation.PLUS else Orientation.PLUS
        
        pairs = []
        for prev, curr in zip(alignment.segments, alignment.segments[1:]):
            donor = genome_seq[prev.genome_rend:prev.genome_rend + 2].upper()
            acceptor = genome_seq[curr.genome_lend - 3:curr.genome_lend - 1].upper()
            pairs.append((donor, acceptor))
        
        for orient in (aligned, opposite):
            if all(pair in self.CONSENSUS_SPLICE_PAIRS[orient] for pair in pairs):
                return orient
        return None
    
    @staticmethod
    def _is_contiguous(alignment: Alignment) -> bool:
        """
        Check that consecutive segments are adjacent on the transcript.
        
        Port of CDNA_alignment::verify_contiguity, using the aligned
        orientation: on '+' each segment starts 1 after the previous one ends
        in transcript coords; on '-' (segments still in genome order) each
        ends 1 before the previous one starts.
        """
        minus = Orientation(alignment.orient) == Orientation.MINUS
        for prev, curr in zip(alignment.segments, alignment.segments[1:]):
            if minus:
                if curr.transcript_rend != prev.transcript_lend - 1:
                    return False
            elif curr.transcript_lend != prev.transcript_rend + 1:
                return False
        return True
