"""
Alignment clusterer for PyPASA.

This module clusters overlapping transcript alignments based on genomic position.
Alignments that overlap on the same scaffold (and strand) are grouped into clusters.

Multiple clustering strategies are supported:
- Default: Any overlap (even 1 bp)
- Stringent: Minimum percentage overlap required
- Gene-based: Cluster by overlap with annotated genes

Corresponds to populate_cdna_clusters.dbi and assign_clusters_by_* in Perl PASA.
"""

import logging
from typing import List, Dict, Set
from collections import defaultdict
from intervaltree import IntervalTree, Interval

from pypasa.storage.models import Alignment, Orientation
from pypasa.storage.datastore import DataStore

logger = logging.getLogger(__name__)


class AlignmentClusterer:
    """
    Cluster transcript alignments by genomic overlap.
    
    Clustering strategies:
    - default: Any overlap triggers clustering
    - stringent: Minimum percentage overlap required
    - gene_overlap: Cluster by overlapping gene annotations (future)
    """
    
    def __init__(
        self,
        datastore: DataStore,
        min_overlap_percent: float = 0.0,
        respect_strand: bool = True,
    ):
        """
        Initialize clusterer.
        
        Args:
            datastore: DataStore for querying/updating alignments
            min_overlap_percent: Minimum percentage overlap for clustering (0-100)
                                If 0, any overlap triggers clustering (default mode)
                                If >0, requires substantial overlap (stringent mode)
            respect_strand: If True, only cluster same-strand alignments
        """
        self.datastore = datastore
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
    
    def cluster_alignments(self, scaffold: str = None) -> Dict[str, int]:
        """
        Cluster alignments by overlap.
        
        Args:
            scaffold: If provided, only cluster alignments on this scaffold.
                     If None, cluster all scaffolds.
        
        Returns:
            Dictionary of clustering statistics
        """
        if scaffold:
            scaffolds = [scaffold]
        else:
            scaffolds = self.datastore.get_scaffolds()
        
        logger.info(f"Clustering alignments on {len(scaffolds)} scaffold(s)")
        
        # Group scaffolds by normalized name to handle gi|XXX vs gi|revXXX
        scaffold_groups = self._group_scaffolds_by_normalized_name(scaffolds)
        logger.info(f"Grouped into {len(scaffold_groups)} normalized scaffold groups")
        
        for normalized_name, scaffold_list in scaffold_groups.items():
            self._cluster_scaffold_group(normalized_name, scaffold_list)
        
        logger.info(
            f"Clustering complete: {self.stats['clustered_alignments']} alignments "
            f"in {self.stats['total_clusters']} clusters"
        )
        
        return self.stats
    
    def _group_scaffolds_by_normalized_name(self, scaffolds: List[str]) -> Dict[str, List[str]]:
        """
        Group scaffolds by normalized name to handle gi|XXX vs gi|revXXX.
        
        This ensures that alignments on both forward and reverse versions
        of the same scaffold are clustered together (matching Perl PASA behavior).
        
        Args:
            scaffolds: List of scaffold names
            
        Returns:
            Dictionary mapping normalized_name -> [list of scaffold variants]
        """
        from collections import defaultdict
        groups = defaultdict(list)
        
        for scaffold in scaffolds:
            # Normalize: gi|revXXX -> gi|XXX
            normalized = scaffold.replace("|rev", "|", 1) if scaffold.startswith("gi|rev") else scaffold
            groups[normalized].append(scaffold)
        
        return dict(groups)
    
    def _cluster_scaffold_group(self, normalized_name: str, scaffold_list: List[str]) -> None:
        """
        Cluster all valid alignments across a group of scaffold variants.
        
        This handles cases where the same transcript aligns to both gi|XXX and gi|revXXX,
        ensuring they are clustered together as they represent the same physical location.
        
        Args:
            normalized_name: Normalized scaffold name
            scaffold_list: List of scaffold variants (e.g., ['gi|68713', 'gi|rev68713'])
        """
        # Load alignments from ALL scaffold variants
        all_alignments = []
        for scaffold in scaffold_list:
            alignments = self.datastore.get_alignments_by_scaffold(scaffold, valid_only=True)
            all_alignments.extend(alignments)
        
        if not all_alignments:
            return
        
        if len(scaffold_list) > 1:
            logger.debug(
                f"Clustering {len(all_alignments)} alignments across "
                f"{len(scaffold_list)} scaffold variants: {scaffold_list}"
            )
        else:
            logger.debug(f"Clustering {len(all_alignments)} alignments on {scaffold_list[0]}")
        
        self.stats['total_alignments'] += len(all_alignments)
        
        # Build interval trees for fast overlap detection
        # Separate by strand if respect_strand is True
        if self.respect_strand:
            trees = self._build_strand_specific_trees(all_alignments)
        else:
            trees = {None: self._build_interval_tree(all_alignments)}
        
        # Find connected components (clusters)
        cluster_assignments = {}
        
        for strand, tree in trees.items():
            # Get alignments for this strand
            if strand is not None:
                strand_alignments = [a for a in all_alignments if a.orient == strand]
            else:
                strand_alignments = all_alignments
            
            # Find overlapping groups using union-find
            clusters = self._find_overlap_clusters(strand_alignments)
            
            # Assign cluster IDs
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
        
        # Bulk update cluster assignments in database
        if cluster_assignments:
            self.datastore.update_alignment_clusters(cluster_assignments)
            logger.debug(
                f"Assigned {len(cluster_assignments)} alignments to "
                f"{len(set(cluster_assignments.values()))} clusters on {scaffold}"
            )
    
    def _build_interval_tree(self, alignments: List[Alignment]) -> IntervalTree:
        """
        Build an interval tree for fast overlap queries.
        
        Args:
            alignments: List of alignments
            
        Returns:
            IntervalTree with alignment IDs as data
        """
        tree = IntervalTree()
        
        for alignment in alignments:
            # Add interval with alignment ID as data
            tree.addi(alignment.genome_lend, alignment.genome_rend + 1, alignment.align_id)
        
        return tree
    
    def _build_strand_specific_trees(
        self, alignments: List[Alignment]
    ) -> Dict[Orientation, IntervalTree]:
        """
        Build separate interval trees for each strand.
        
        Args:
            alignments: List of alignments
            
        Returns:
            Dictionary mapping Orientation -> IntervalTree
        """
        trees = {
            Orientation.PLUS: IntervalTree(),
            Orientation.MINUS: IntervalTree(),
            Orientation.UNKNOWN: IntervalTree(),
        }
        
        for alignment in alignments:
            tree = trees[alignment.orient]
            tree.addi(alignment.genome_lend, alignment.genome_rend + 1, alignment.align_id)
        
        return trees
    
    def _find_overlap_clusters(self, alignments: List[Alignment]) -> List[Set[int]]:
        """
        Find clusters of overlapping alignments using union-find.
        
        Args:
            alignments: List of alignments to cluster
            
        Returns:
            List of sets, each set containing align_ids that form a cluster
        """
        # Build parent mapping for union-find
        parent = {a.align_id: a.align_id for a in alignments}
        
        def find(x):
            """Find root of x with path compression."""
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        def union(x, y):
            """Union two sets."""
            root_x = find(x)
            root_y = find(y)
            if root_x != root_y:
                parent[root_x] = root_y
        
        # Check each pair of alignments for overlap
        alignment_dict = {a.align_id: a for a in alignments}
        
        for i, aln1 in enumerate(alignments):
            for aln2 in alignments[i + 1:]:
                if self._should_cluster(aln1, aln2):
                    union(aln1.align_id, aln2.align_id)
        
        # Group alignments by root
        clusters = defaultdict(set)
        for align_id in parent:
            root = find(align_id)
            clusters[root].add(align_id)
        
        return list(clusters.values())
    
    def _should_cluster(self, aln1: Alignment, aln2: Alignment) -> bool:
        """
        Determine if two alignments should be clustered together.
        
        Args:
            aln1: First alignment
            aln2: Second alignment
            
        Returns:
            True if alignments should be clustered
        """
        # Must be on same scaffold
        if aln1.scaffold != aln2.scaffold:
            return False
        
        # Check strand if required
        if self.respect_strand:
            # Skip if both have known strands that differ
            if (aln1.orient != Orientation.UNKNOWN and 
                aln2.orient != Orientation.UNKNOWN and
                aln1.orient != aln2.orient):
                return False
        
        # Calculate overlap
        overlap_start = max(aln1.genome_lend, aln2.genome_lend)
        overlap_end = min(aln1.genome_rend, aln2.genome_rend)
        overlap = overlap_end - overlap_start + 1
        
        if overlap < 1:
            return False  # No overlap
        
        # Check minimum overlap percentage if required
        if self.min_overlap_percent > 0:
            # Calculate overlap as percentage of shorter alignment
            min_length = min(aln1.genome_span, aln2.genome_span)
            overlap_percent = (overlap / min_length) * 100.0
            
            if overlap_percent < self.min_overlap_percent:
                return False
        
        return True
    
    def get_clustering_report(self) -> str:
        """
        Generate a human-readable clustering report.
        
        Returns:
            Formatted report string
        """
        total_aligns = self.stats['total_alignments']
        total_clusters = self.stats['total_clusters']
        
        if total_aligns == 0:
            return "No alignments clustered"
        
        report = []
        report.append("=" * 60)
        report.append("ALIGNMENT CLUSTERING REPORT")
        report.append("=" * 60)
        report.append(f"Total alignments:        {total_aligns:>10,}")
        report.append(f"Clustered alignments:    {self.stats['clustered_alignments']:>10,}")
        report.append(f"Total clusters:          {total_clusters:>10,}")
        report.append("")
        report.append("Cluster composition:")
        report.append(f"  Singleton clusters:    {self.stats['singleton_clusters']:>10,} ({100.0 * self.stats['singleton_clusters'] / total_clusters:>6.2f}%)")
        report.append(f"  Multi-alignment:       {self.stats['multi_alignment_clusters']:>10,} ({100.0 * self.stats['multi_alignment_clusters'] / total_clusters:>6.2f}%)")
        report.append("")
        
        if self.min_overlap_percent > 0:
            report.append(f"Clustering mode: Stringent (min {self.min_overlap_percent}% overlap)")
        else:
            report.append("Clustering mode: Default (any overlap)")
        
        report.append(f"Strand-specific: {self.respect_strand}")
        report.append("=" * 60)
        
        return "\n".join(report)
