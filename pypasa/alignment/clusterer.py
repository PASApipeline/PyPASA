"""
Alignment clusterer for PyPASA.

This module clusters overlapping transcript alignments based on genomic position.
Uses Perl PASA's expanding window algorithm (Overlap_piler.pm) to create maximal
overlap piles. 

By default, matches Perl PASA's default clustering behavior (reassign_clusters_via_valid_align_coords.dbi):
- Orientation-agnostic: alignments cluster by coordinate overlap only
- Single-exon alignments (unknown spliced_orient) cluster with any overlapping alignment
- Mixed-orientation clusters are allowed at this stage
- Each scaffold is treated as completely distinct (gi|XXX and gi|revXXX are separate)

This creates CLUSTERS (first-pass grouping). After assembly, SUBCLUSTERS will be
created to identify alternative splice variants (≥50% overlap requirement), where
strand-specificity IS enforced (subcluster_builder.dbi).

Corresponds to:
- Overlap_piler.pm (cluster creation)
- reassign_clusters_via_valid_align_coords.dbi (default clustering mode)
- assign_clusters_by_stringent_alignment_overlap.dbi (optional strand-aware mode)
- subcluster_builder.dbi (subclustering - implemented after assembly)
"""

import logging
from typing import List, Dict, Set, Tuple
from collections import defaultdict

from pypasa.models import Alignment, Orientation

logger = logging.getLogger(__name__)


class AlignmentClusterer:
    """
    Cluster transcript alignments by genomic overlap using expanding window algorithm.
    
    This implements Perl PASA's Overlap_piler.pm algorithm:
    1. Sort alignments by left coordinate
    2. Process left-to-right with expanding window
    3. If alignment overlaps window: expand window, add to cluster
    4. If no overlap: close cluster, start new window
    
    This creates transitive clusters where alignments may not directly overlap
    but are connected through intermediate alignments.
    """
    
    def __init__(
        self,
        min_overlap_percent: float = 0.0,
        respect_strand: bool = False,
    ):
        """
        Initialize clusterer.
        
        Args:
            min_overlap_percent: Minimum percentage overlap for clustering (0-100)
                                If 0, any overlap triggers clustering (default mode)
                                Currently not used - reserved for future stringent mode
            respect_strand: If False (default), orientation-agnostic clustering - matches Perl's
                           default mode (reassign_clusters_via_valid_align_coords.dbi).
                           If True, strand-specific clustering - matches Perl's stringent mode
                           (assign_clusters_by_stringent_alignment_overlap.dbi).
        """
        self.min_overlap_percent = min_overlap_percent
        self.respect_strand = respect_strand
        
        self.next_cluster_id = 1
        self.stats = {
            'total_alignments': 0,
            'clustered_alignments': 0,
            'total_clusters': 0,
            'singleton_clusters': 0,
            'multi_alignment_clusters': 0,
        }
    
    def cluster_alignments(self, alignments: List[Alignment]) -> Dict[str, int]:
        """
        Cluster valid alignments by overlap, setting cluster_id in place.
        
        Cluster ids are assigned in (scaffold name, position) order.
        
        Returns:
            Dictionary of clustering statistics
        """
        by_scaffold: Dict[str, List[Alignment]] = defaultdict(list)
        for aln in alignments:
            if aln.valid:
                by_scaffold[aln.scaffold].append(aln)
        
        logger.info(f"Clustering alignments on {len(by_scaffold)} scaffold(s)")
        
        # Treat each scaffold independently (matching Perl PASA clustering behavior)
        # Note: gi|XXX and gi|revXXX are kept separate during clustering
        for scaffold in sorted(by_scaffold):
            self._cluster_scaffold(by_scaffold[scaffold])
        
        logger.info(
            f"Clustering complete: {self.stats['clustered_alignments']} alignments "
            f"in {self.stats['total_clusters']} clusters"
        )
        
        return self.stats
    
    def _cluster_scaffold(self, all_alignments: List[Alignment]) -> None:
        """
        Cluster the valid alignments of one scaffold.
        
        Uses expanding window algorithm (Perl's Overlap_piler.pm) to create
        maximal overlap piles.
        """
        
        self.stats['total_alignments'] += len(all_alignments)
        
        # Separate by strand if respect_strand is True
        if self.respect_strand:
            # Group by effective orientation (transcribed orientation with fallback to aligned)
            # This matches Perl PASA's behavior where clustering uses spliced_orient when available
            by_strand = defaultdict(list)
            for aln in all_alignments:
                by_strand[aln.effective_orient].append(aln)
            
            cluster_assignments = {}
            for strand, strand_alignments in by_strand.items():
                strand_clusters = self._expanding_window_cluster(strand_alignments)
                
                # Assign cluster IDs
                for cluster_members in strand_clusters:
                    cluster_id = self.next_cluster_id
                    self.next_cluster_id += 1
                    
                    for align_id in cluster_members:
                        cluster_assignments[align_id] = cluster_id
                    
                    # Update stats
                    self.stats['total_clusters'] += 1
                    if len(cluster_members) == 1:
                        self.stats['singleton_clusters'] += 1
                    else:
                        self.stats['multi_alignment_clusters'] += 1
        else:
            # Don't separate by strand
            clusters = self._expanding_window_cluster(all_alignments)
            
            cluster_assignments = {}
            for cluster_members in clusters:
                cluster_id = self.next_cluster_id
                self.next_cluster_id += 1
                
                for align_id in cluster_members:
                    cluster_assignments[align_id] = cluster_id
                
                # Update stats
                self.stats['total_clusters'] += 1
                if len(cluster_members) == 1:
                    self.stats['singleton_clusters'] += 1
                else:
                    self.stats['multi_alignment_clusters'] += 1
        
        self.stats['clustered_alignments'] += len(cluster_assignments)
        
        for aln in all_alignments:
            aln.cluster_id = cluster_assignments[aln.align_id]
    
    def _expanding_window_cluster(self, alignments: List[Alignment]) -> List[List[int]]:
        """
        Cluster alignments using expanding window algorithm (Perl's Overlap_piler.pm).
        
        Algorithm:
        1. Sort alignments by left coordinate (genome_lend)
        2. Initialize window with first alignment
        3. For each subsequent alignment:
           - If it overlaps the current window: expand window, add to cluster
           - Else: close current cluster, start new window
        
        This creates transitive clusters where alignments may not directly overlap
        but are connected through intermediate alignments.
        
        Args:
            alignments: List of alignments to cluster (same scaffold/strand)
            
        Returns:
            List of clusters, where each cluster is a list of align_ids
        """
        if not alignments:
            return []
        
        # Sort by left coordinate
        sorted_alns = sorted(alignments, key=lambda a: a.genome_lend)
        
        clusters = []
        current_cluster = [sorted_alns[0].align_id]
        
        # Initialize window with first alignment coordinates
        window_left = sorted_alns[0].genome_lend
        window_right = sorted_alns[0].genome_rend
        
        # Process remaining alignments
        for i in range(1, len(sorted_alns)):
            aln = sorted_alns[i]
            lend = aln.genome_lend
            rend = aln.genome_rend
            
            # Check if alignment overlaps current window
            if window_left <= rend and window_right >= lend:
                # Overlap - expand window and add to current cluster
                window_left = min(window_left, lend)
                window_right = max(window_right, rend)
                current_cluster.append(aln.align_id)
            else:
                # No overlap - close current cluster and start new one
                clusters.append(current_cluster)
                
                # Start new cluster
                current_cluster = [aln.align_id]
                window_left = lend
                window_right = rend
        
        # Don't forget the final cluster
        clusters.append(current_cluster)
        
        return clusters
    
    def get_clustering_report(self) -> str:
        """
        Generate a human-readable clustering report.
        
        Returns:
            Formatted report string
        """
        report = []
        report.append("=" * 60)
        report.append("ALIGNMENT CLUSTERING REPORT")
        report.append("=" * 60)
        
        total = self.stats['total_alignments']
        clustered = self.stats['clustered_alignments']
        clusters = self.stats['total_clusters']
        singletons = self.stats['singleton_clusters']
        multi = self.stats['multi_alignment_clusters']
        
        report.append(f"Total alignments:            {total:>10,}")
        report.append(f"Clustered alignments:        {clustered:>10,}")
        report.append(f"Total clusters:               {clusters:>10,}")
        report.append("")
        report.append("Cluster composition:")
        if clusters > 0:
            report.append(f"  Singleton clusters:             {singletons:>5,} ({singletons/clusters*100:>6.2f}%)")
            report.append(f"  Multi-alignment:            {multi:>5,} ({multi/clusters*100:>6.2f}%)")
        else:
            report.append("  Singleton clusters:                 0")
            report.append("  Multi-alignment:                0")
        report.append("")
        report.append("Clustering algorithm: Expanding window (Perl's Overlap_piler.pm)")
        report.append(f"Strand-specific: {self.respect_strand}")
        report.append("")
        report.append("Note: This creates CLUSTERS (first-pass grouping).")
        report.append("      SUBCLUSTERS (alternative splice variants) will be")
        report.append("      created after assembly step (future implementation).")
        report.append("=" * 60)
        
        return "\n".join(report)
