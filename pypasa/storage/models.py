"""
Data models for PyPASA.

These models represent alignments, assemblies, and related structures,
replacing the database-centric approach of the original Perl implementation
with Python dataclasses backed by SQLite for persistence.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from enum import Enum


class Orientation(str, Enum):
    """Strand orientation."""
    PLUS = "+"
    MINUS = "-"
    UNKNOWN = "?"


@dataclass
class AlignmentSegment:
    """
    Represents a single segment (exon) of an alignment.
    
    Corresponds to rows in the 'alignment' table in original PASA.
    """
    genome_lend: int  # Genomic left end (1-based)
    genome_rend: int  # Genomic right end (1-based)
    transcript_lend: int  # Transcript left end (1-based)
    transcript_rend: int  # Transcript right end (1-based)
    percent_id: float  # Percent identity
    
    def __post_init__(self) -> None:
        """Validate coordinates."""
        if self.genome_lend > self.genome_rend:
            raise ValueError(f"Invalid genomic coords: {self.genome_lend} > {self.genome_rend}")
        if self.transcript_lend > self.transcript_rend:
            raise ValueError(
                f"Invalid transcript coords: {self.transcript_lend} > {self.transcript_rend}"
            )
    
    @property
    def genome_length(self) -> int:
        """Length of this segment on genome."""
        return self.genome_rend - self.genome_lend + 1
    
    @property
    def transcript_length(self) -> int:
        """Length of this segment on transcript."""
        return self.transcript_rend - self.transcript_lend + 1


@dataclass
class Alignment:
    """
    Represents a complete transcript alignment to the genome.
    
    This is a multi-exon alignment consisting of one or more segments.
    Replaces the 'align_link' and related tables in original PASA.
    """
    align_id: int
    cdna_acc: str  # Transcript accession
    scaffold: str  # Genome scaffold/chromosome
    segments: List[AlignmentSegment]
    orient: Orientation
    spliced_orient: Optional[Orientation] = None
    
    # Validation status
    valid: bool = False
    
    # Clustering
    cluster_id: Optional[int] = None
    
    # Quality metrics
    avg_percent_id: Optional[float] = None
    percent_aligned: Optional[float] = None
    
    # Full-length status
    is_fli: bool = False  # Full-length insert
    
    def __post_init__(self) -> None:
        """Calculate derived metrics."""
        if not self.segments:
            raise ValueError("Alignment must have at least one segment")
        
        # Sort segments by genome position
        self.segments.sort(key=lambda s: s.genome_lend)
        
        # Calculate average percent identity
        if self.avg_percent_id is None:
            total_length = sum(s.genome_length for s in self.segments)
            weighted_id = sum(s.percent_id * s.genome_length for s in self.segments)
            self.avg_percent_id = weighted_id / total_length if total_length > 0 else 0.0
    
    @property
    def genome_lend(self) -> int:
        """Leftmost genomic coordinate."""
        return min(s.genome_lend for s in self.segments)
    
    @property
    def genome_rend(self) -> int:
        """Rightmost genomic coordinate."""
        return max(s.genome_rend for s in self.segments)
    
    @property
    def effective_orient(self) -> Orientation:
        """
        Get the effective orientation for clustering/assembly.
        
        Uses spliced_orient (transcribed orientation from splice sites) if available,
        otherwise falls back to aligned orient.
        
        This matches Perl PASA's force_flexorient behavior where alignments are
        grouped by transcribed orientation, not aligned orientation.
        """
        return self.spliced_orient if self.spliced_orient else self.orient
    
    @property
    def genome_span(self) -> int:
        """Span on genome (including introns)."""
        return self.genome_rend - self.genome_lend + 1
    
    @property
    def num_segments(self) -> int:
        """Number of exons/segments."""
        return len(self.segments)
    
    @property
    def aligned_length(self) -> int:
        """Total aligned length (sum of exons on transcript, not genome)."""
        return sum(s.transcript_length for s in self.segments)
    
    def overlaps(self, other: "Alignment", min_overlap: int = 1) -> bool:
        """
        Check if this alignment overlaps another on the same scaffold.
        
        Args:
            other: Another alignment
            min_overlap: Minimum overlap required (default 1bp)
            
        Returns:
            True if alignments overlap
        """
        if self.scaffold != other.scaffold:
            return False
        
        # Check strand if both are oriented
        if (self.orient != Orientation.UNKNOWN and 
            other.orient != Orientation.UNKNOWN and
            self.orient != other.orient):
            return False
        
        # Check coordinate overlap
        overlap = (
            min(self.genome_rend, other.genome_rend) - 
            max(self.genome_lend, other.genome_lend) + 1
        )
        
        return overlap >= min_overlap
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "align_id": self.align_id,
            "cdna_acc": self.cdna_acc,
            "scaffold": self.scaffold,
            "segments": [
                {
                    "genome_lend": s.genome_lend,
                    "genome_rend": s.genome_rend,
                    "transcript_lend": s.transcript_lend,
                    "transcript_rend": s.transcript_rend,
                    "percent_id": s.percent_id,
                }
                for s in self.segments
            ],
            "orient": self.orient.value,
            "spliced_orient": self.spliced_orient.value if self.spliced_orient else None,
            "valid": self.valid,
            "cluster_id": self.cluster_id,
            "avg_percent_id": self.avg_percent_id,
            "percent_aligned": self.percent_aligned,
            "is_fli": self.is_fli,
        }


@dataclass
class Assembly:
    """
    Represents an assembled transcript structure from clustered alignments.
    
    This is the output of the PASA assembly algorithm.
    Replaces the 'clusters' and 'subcluster_link' tables.
    """
    assembly_id: int
    cluster_id: int
    scaffold: str
    lend: int
    rend: int
    orient: Orientation
    
    # Assembly structure
    segments: List[AlignmentSegment]
    
    # Component alignments
    component_align_ids: List[int] = field(default_factory=list)
    
    # Sequence data
    sequence: Optional[str] = None
    
    # Full-length status
    is_fli: bool = False
    
    @property
    def span(self) -> int:
        """Genomic span of assembly."""
        return self.rend - self.lend + 1
    
    @property
    def num_segments(self) -> int:
        """Number of exons."""
        return len(self.segments)
    
    def to_gff3_lines(self) -> List[str]:
        """
        Generate GFF3 format lines for this assembly.
        
        Returns:
            List of GFF3 lines (without trailing newlines)
        """
        lines = []
        
        # Gene/mRNA line
        attrs = f"ID=asmbl_{self.assembly_id};Target={self.assembly_id}"
        lines.append(
            f"{self.scaffold}\tPASA\tmRNA\t{self.lend}\t{self.rend}\t.\t"
            f"{self.orient.value}\t.\t{attrs}"
        )
        
        # Exon lines
        for i, seg in enumerate(self.segments, 1):
            attrs = f"ID=asmbl_{self.assembly_id}.exon{i};Parent=asmbl_{self.assembly_id}"
            lines.append(
                f"{self.scaffold}\tPASA\texon\t{seg.genome_lend}\t{seg.genome_rend}\t.\t"
                f"{self.orient.value}\t.\t{attrs}"
            )
        
        return lines
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "assembly_id": self.assembly_id,
            "cluster_id": self.cluster_id,
            "scaffold": self.scaffold,
            "lend": self.lend,
            "rend": self.rend,
            "orient": self.orient.value,
            "segments": [
                {
                    "genome_lend": s.genome_lend,
                    "genome_rend": s.genome_rend,
                    "transcript_lend": s.transcript_lend,
                    "transcript_rend": s.transcript_rend,
                    "percent_id": s.percent_id,
                }
                for s in self.segments
            ],
            "component_align_ids": self.component_align_ids,
            "sequence": self.sequence,
            "is_fli": self.is_fli,
        }
