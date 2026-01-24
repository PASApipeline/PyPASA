"""
GFF3 alignment importer for PyPASA.

This module parses GFF3 format alignment files (from GMAP, BLAT, minimap2)
and imports them into the DataStore.

GFF3 Format:
    scaffold  aligner  cDNA_match  start  end  %id  strand  .  ID=align_X;Target=acc start end [strand]
    
Multiple lines with same ID represent a multi-exon alignment.
"""

import re
import logging
from pathlib import Path
from typing import List, Dict, Iterator, Tuple
from collections import defaultdict

from pypasa.storage.models import Alignment, AlignmentSegment, Orientation
from pypasa.storage.datastore import DataStore

logger = logging.getLogger(__name__)


class GFF3AlignmentImporter:
    """
    Import transcript alignments from GFF3 format.
    
    Handles multi-exon alignments where multiple GFF3 lines
    share the same alignment ID.
    """
    
    def __init__(self, datastore: DataStore):
        """
        Initialize importer.
        
        Args:
            datastore: DataStore for storing alignments
        """
        self.datastore = datastore
        self.next_align_id = 1
        
    def import_gff3(self, gff3_path: str, aligner: str = "gmap") -> int:
        """
        Import alignments from GFF3 file.
        
        Args:
            gff3_path: Path to GFF3 file
            aligner: Name of aligner (gmap, blat, minimap2)
            
        Returns:
            Number of alignments imported
        """
        logger.info(f"Importing alignments from {gff3_path}")
        
        # Parse GFF3 into alignment groups
        alignment_groups = self._parse_gff3_file(gff3_path)
        
        # Convert to Alignment objects
        alignments = []
        for align_id, (cdna_acc, scaffold, orient, segments) in alignment_groups.items():
            try:
                alignment = Alignment(
                    align_id=self.next_align_id,
                    cdna_acc=cdna_acc,
                    scaffold=scaffold,
                    segments=segments,
                    orient=orient,
                    spliced_orient=None,  # Determined later
                    valid=False,  # Validation happens next
                )
                alignments.append(alignment)
                self.next_align_id += 1
            except Exception as e:
                logger.warning(f"Skipping alignment {align_id}: {e}")
                continue
        
        # Bulk insert
        if alignments:
            self.datastore.add_alignments(alignments)
            logger.info(f"Imported {len(alignments)} alignments from {gff3_path}")
        else:
            logger.warning("No valid alignments found in file")
        
        return len(alignments)
    
    def _parse_gff3_file(
        self, gff3_path: str
    ) -> Dict[str, Tuple[str, str, Orientation, List[AlignmentSegment]]]:
        """
        Parse GFF3 file and group segments by alignment ID.
        
        Args:
            gff3_path: Path to GFF3 file
            
        Returns:
            Dictionary mapping alignment_id -> (cdna_acc, scaffold, orient, segments)
        """
        alignment_groups: Dict[str, List[Dict]] = defaultdict(list)
        
        with open(gff3_path) as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                # Parse GFF3 line
                try:
                    entry = self._parse_gff3_line(line)
                    if entry:
                        alignment_groups[entry['align_id']].append(entry)
                except Exception as e:
                    logger.debug(f"Skipping malformed line: {line[:100]} - {e}")
                    continue
        
        # Convert grouped entries to alignments
        alignments = {}
        for align_id, entries in alignment_groups.items():
            try:
                alignment_data = self._merge_alignment_segments(entries)
                alignments[align_id] = alignment_data
            except Exception as e:
                logger.warning(f"Failed to merge alignment {align_id}: {e}")
                continue
        
        return alignments
    
    def _parse_gff3_line(self, line: str) -> Dict:
        """
        Parse a single GFF3 line.
        
        GFF3 format:
        seqid source type start end score strand phase attributes
        
        Returns:
            Dictionary with parsed fields, or None if invalid
        """
        fields = line.split('\t')
        if len(fields) != 9:
            return None
        
        scaffold = fields[0]
        aligner = fields[1]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        score = fields[5]  # Percent identity
        strand = fields[6]
        phase = fields[7]
        attributes = fields[8]
        
        # Parse attributes
        attr_dict = {}
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key.strip()] = value.strip()
        
        # Extract ID and Target
        align_id = attr_dict.get('ID')
        target_str = attr_dict.get('Target')
        
        if not (align_id and target_str):
            return None
        
        # Parse Target: "cdna_acc t_start t_end [t_strand]"
        target_parts = target_str.split()
        if len(target_parts) < 3:
            return None
        
        cdna_acc = target_parts[0]
        t_start = int(target_parts[1])
        t_end = int(target_parts[2])
        t_strand = target_parts[3] if len(target_parts) > 3 else '+'
        
        # Parse percent identity from score
        try:
            percent_id = float(score)
        except ValueError:
            percent_id = 100.0  # Default if not provided
        
        # Determine orientation
        if strand == '+':
            orient = Orientation.PLUS
        elif strand == '-':
            orient = Orientation.MINUS
        else:
            orient = Orientation.UNKNOWN
        
        return {
            'align_id': align_id,
            'cdna_acc': cdna_acc,
            'scaffold': scaffold,
            'genome_lend': min(start, end),
            'genome_rend': max(start, end),
            'transcript_lend': min(t_start, t_end),
            'transcript_rend': max(t_start, t_end),
            'orient': orient,
            'percent_id': percent_id,
        }
    
    def _merge_alignment_segments(
        self, entries: List[Dict]
    ) -> Tuple[str, str, Orientation, List[AlignmentSegment]]:
        """
        Merge multiple GFF3 entries (exons) into a single alignment.
        
        Args:
            entries: List of parsed GFF3 entries for same alignment
            
        Returns:
            Tuple of (cdna_acc, scaffold, orient, segments)
        """
        if not entries:
            raise ValueError("No entries to merge")
        
        # All entries should have same cdna_acc, scaffold, orient
        cdna_acc = entries[0]['cdna_acc']
        scaffold = entries[0]['scaffold']
        orient = entries[0]['orient']
        
        # Create segments
        segments = []
        for entry in entries:
            segment = AlignmentSegment(
                genome_lend=entry['genome_lend'],
                genome_rend=entry['genome_rend'],
                transcript_lend=entry['transcript_lend'],
                transcript_rend=entry['transcript_rend'],
                percent_id=entry['percent_id'],
            )
            segments.append(segment)
        
        # Sort segments by genome position
        segments.sort(key=lambda s: s.genome_lend)
        
        return (cdna_acc, scaffold, orient, segments)
