"""
Alignment validator for PyPASA.

This module validates imported alignments based on quality criteria:
- Minimum percent identity
- Minimum alignment length
- Minimum percent of transcript aligned
- Perfect splice boundaries (optional)

Corresponds to validate_alignments_in_db.dbi in Perl PASA.
"""

import logging
from typing import List, Dict, Optional

from pypasa.storage.models import Alignment, Orientation
from pypasa.storage.datastore import DataStore

logger = logging.getLogger(__name__)


class AlignmentValidator:
    """
    Validate transcript alignments based on quality criteria.
    
    Quality filters:
    - min_percent_id: Minimum average percent identity (default 95.0)
    - min_percent_aligned: Minimum % of transcript aligned (default 75.0)
    - min_avg_per_id: Alternative min percent identity threshold
    - require_perfect_splice: Require perfect dinucleotide splice boundaries
    """
    
    def __init__(
        self,
        datastore: DataStore,
        min_percent_id: float = 95.0,
        min_percent_aligned: float = 75.0,
        num_bp_perfect_splice: int = 0,
        genome_fasta: str = None,
        transcript_fasta: str = None,
    ):
        """
        Initialize validator.
        
        Args:
            datastore: DataStore for querying/updating alignments
            min_percent_id: Minimum average percent identity
            min_percent_aligned: Minimum percent of transcript aligned
            num_bp_perfect_splice: Number of bp required for perfect splice boundary (0=disabled)
            genome_fasta: Path to genome FASTA file (for splice site detection)
            transcript_fasta: Path to transcript FASTA file (for calculating percent_aligned)
        """
        self.datastore = datastore
        self.min_percent_id = min_percent_id
        self.min_percent_aligned = min_percent_aligned
        self.num_bp_perfect_splice = num_bp_perfect_splice
        self.genome_fasta = genome_fasta
        self.transcript_fasta = transcript_fasta
        self._genome_seqs = {}  # Cache for genome sequences
        self._transcript_lengths = {}  # Cache for transcript lengths
        
        self.stats = {
            'total': 0,
            'passed': 0,
            'failed_percent_id': 0,
            'failed_percent_aligned': 0,
            'failed_splice': 0,
        }
    
    def validate_alignments(self, scaffold: str = None) -> Dict[str, int]:
        """
        Validate alignments and mark them as valid/invalid.
        
        Args:
            scaffold: If provided, only validate alignments on this scaffold.
                     If None, validate all scaffolds.
        
        Returns:
            Dictionary of validation statistics
        """
        if scaffold:
            scaffolds = [scaffold]
        else:
            scaffolds = self.datastore.get_scaffolds()
        
        logger.info(f"Validating alignments on {len(scaffolds)} scaffold(s)")
        
        for scaffold in scaffolds:
            self._validate_scaffold(scaffold)
        
        logger.info(
            f"Validation complete: {self.stats['passed']}/{self.stats['total']} passed "
            f"({100.0 * self.stats['passed'] / self.stats['total']:.1f}%)"
        )
        
        return self.stats
    
    def _validate_scaffold(self, scaffold: str) -> None:
        """
        Validate all alignments on a single scaffold.
        
        Args:
            scaffold: Scaffold name
        """
        # Load alignments for this scaffold
        alignments = self.datastore.get_alignments_by_scaffold(scaffold, valid_only=False)
        
        if not alignments:
            return
        
        logger.debug(f"Validating {len(alignments)} alignments on {scaffold}")
        
        validity_updates = {}
        spliced_orient_updates = {}
        percent_aligned_updates = {}
        
        # Load genome sequence if needed for splice detection
        genome_seq = None
        if self.genome_fasta:
            genome_seq = self._get_genome_sequence(scaffold)
        
        for alignment in alignments:
            self.stats['total'] += 1
            
            # Detect transcribed orientation from splice sites
            if genome_seq and len(alignment.segments) > 1:
                spliced_orient = self._detect_spliced_orientation(alignment, genome_seq)
                if spliced_orient:
                    spliced_orient_updates[alignment.align_id] = spliced_orient
            
            # Validate alignment
            is_valid, reason = self._validate_alignment(alignment)
            
            validity_updates[alignment.align_id] = is_valid
            
            # Store percent_aligned if calculated
            if alignment.percent_aligned is not None:
                percent_aligned_updates[alignment.align_id] = alignment.percent_aligned
            
            if is_valid:
                self.stats['passed'] += 1
            else:
                # Track failure reason
                if 'percent_id' in reason:
                    self.stats['failed_percent_id'] += 1
                elif 'aligned' in reason:
                    self.stats['failed_percent_aligned'] += 1
                elif 'splice' in reason:
                    self.stats['failed_splice'] += 1
                
                logger.debug(
                    f"Alignment {alignment.cdna_acc} on {scaffold} failed: {reason}"
                )
        
        # Bulk update validity, spliced orientation, and percent_aligned in database
        if validity_updates:
            self.datastore.update_alignment_validity(validity_updates)
        if spliced_orient_updates:
            self.datastore.update_spliced_orientations(spliced_orient_updates)
        if percent_aligned_updates:
            self.datastore.update_percent_aligned(percent_aligned_updates)
    
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
        
        # TODO: Check splice boundaries if required
        # This would require genome sequence access
        if self.num_bp_perfect_splice > 0:
            # Placeholder - would need genome sequence
            pass
        
        return (True, "passed")
    
    def _get_transcript_length(self, cdna_acc: str) -> int:
        """
        Get the FULL length of a transcript from the FASTA file.
        
        This is critical for correct percent_aligned calculation.
        The max aligned position does NOT equal the transcript length!
        
        Args:
            cdna_acc: Transcript accession
            
        Returns:
            Full transcript length in bases, or 0 if not found
        """
        # Check cache first
        if cdna_acc in self._transcript_lengths:
            return self._transcript_lengths[cdna_acc]
        
        # Load from FASTA file
        if not self.transcript_fasta:
            logger.warning(
                f"No transcript FASTA provided - cannot calculate percent_aligned for {cdna_acc}"
            )
            return 0
        
        try:
            with open(self.transcript_fasta, 'r') as f:
                in_target = False
                seq_length = 0
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Check if we found our sequence
                        if in_target:
                            self._transcript_lengths[cdna_acc] = seq_length
                            return seq_length
                        
                        # Check if this header matches our accession
                        header = line[1:].split()[0]
                        if cdna_acc in header or header in cdna_acc:
                            in_target = True
                            seq_length = 0
                    elif in_target:
                        seq_length += len(line)
                
                # Check last sequence
                if in_target:
                    self._transcript_lengths[cdna_acc] = seq_length
                    return seq_length
        except Exception as e:
            logger.warning(f"Could not load transcript length for {cdna_acc}: {e}")
        
        logger.warning(f"Transcript {cdna_acc} not found in {self.transcript_fasta}")
        return 0
    
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
        report.append("=" * 60)
        
        return "\n".join(report)
    
    def _get_genome_sequence(self, scaffold: str) -> str:
        """
        Load genome sequence for a scaffold.
        
        Args:
            scaffold: Scaffold name
            
        Returns:
            Genome sequence as string
        """
        if scaffold in self._genome_seqs:
            return self._genome_seqs[scaffold]
        
        if not self.genome_fasta:
            return None
        
        # Simple FASTA parser
        try:
            with open(self.genome_fasta, 'r') as f:
                current_scaffold = None
                seq_lines = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous scaffold
                        if current_scaffold == scaffold and seq_lines:
                            self._genome_seqs[scaffold] = ''.join(seq_lines)
                            return self._genome_seqs[scaffold]
                        
                        # Start new scaffold
                        current_scaffold = line[1:].split()[0]
                        seq_lines = []
                    elif current_scaffold == scaffold:
                        seq_lines.append(line)
                
                # Check last scaffold
                if current_scaffold == scaffold and seq_lines:
                    self._genome_seqs[scaffold] = ''.join(seq_lines)
                    return self._genome_seqs[scaffold]
        except Exception as e:
            logger.warning(f"Could not load genome sequence for {scaffold}: {e}")
        
        return None
    
    def _detect_spliced_orientation(self, alignment: Alignment, genome_seq: str) -> Optional[Orientation]:
        """
        Detect transcribed orientation from splice dinucleotides.
        
        Checks intron boundaries for canonical splice sites:
        - GT-AG or GC-AG = forward (+) transcription
        - CT-AC or CT-GC = reverse (-) transcription
        
        Args:
            alignment: Alignment with multiple segments (introns)
            genome_seq: Genome sequence
            
        Returns:
            Orientation.PLUS, Orientation.MINUS, or None if ambiguous
        """
        if len(alignment.segments) < 2:
            return None  # No introns
        
        plus_votes = 0
        minus_votes = 0
        
        # Check each intron
        for i in range(len(alignment.segments) - 1):
            donor_end = alignment.segments[i].genome_rend
            acceptor_start = alignment.segments[i + 1].genome_lend
            
            # Extract splice dinucleotides (1-based to 0-based)
            try:
                donor = genome_seq[donor_end:donor_end + 2].upper()
                acceptor = genome_seq[acceptor_start - 3:acceptor_start - 1].upper()
                
                # Check for forward strand signals (GT-AG, GC-AG)
                if donor in ['GT', 'GC'] and acceptor == 'AG':
                    plus_votes += 1
                # Check for reverse strand signals (CT-AC, CT-GC)
                elif donor == 'CT' and acceptor in ['AC', 'GC']:
                    minus_votes += 1
            except IndexError:
                continue
        
        # Determine orientation based on votes
        if plus_votes > 0 and minus_votes == 0:
            return Orientation.PLUS
        elif minus_votes > 0 and plus_votes == 0:
            return Orientation.MINUS
        else:
            return None  # Ambiguous or no clear signal
