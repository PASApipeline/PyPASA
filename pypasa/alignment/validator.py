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
from typing import List, Dict

from pypasa.storage.models import Alignment
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
    ):
        """
        Initialize validator.
        
        Args:
            datastore: DataStore for querying/updating alignments
            min_percent_id: Minimum average percent identity
            min_percent_aligned: Minimum percent of transcript aligned
            num_bp_perfect_splice: Number of bp required for perfect splice boundary (0=disabled)
        """
        self.datastore = datastore
        self.min_percent_id = min_percent_id
        self.min_percent_aligned = min_percent_aligned
        self.num_bp_perfect_splice = num_bp_perfect_splice
        
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
        
        for alignment in alignments:
            self.stats['total'] += 1
            
            # Validate alignment
            is_valid, reason = self._validate_alignment(alignment)
            
            validity_updates[alignment.align_id] = is_valid
            
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
        
        # Bulk update validity in database
        if validity_updates:
            self.datastore.update_alignment_validity(validity_updates)
    
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
        
        # Check percent aligned (if we have transcript length)
        # Note: This requires transcript length info which we'll calculate
        # For now, we'll compute it from the segments
        transcript_span = self._get_transcript_span(alignment)
        if transcript_span > 0:
            percent_aligned = (alignment.aligned_length / transcript_span) * 100.0
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
    
    def _get_transcript_span(self, alignment: Alignment) -> int:
        """
        Calculate the span of the transcript covered by this alignment.
        
        Args:
            alignment: Alignment object
            
        Returns:
            Transcript span in bases
        """
        if not alignment.segments:
            return 0
        
        min_pos = min(s.transcript_lend for s in alignment.segments)
        max_pos = max(s.transcript_rend for s in alignment.segments)
        
        return max_pos - min_pos + 1
    
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
