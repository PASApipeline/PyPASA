"""
Cluster deduplicator module.

This module implements post-clustering deduplication to ensure only one
alignment per cDNA per cluster, matching Perl PASA's 
ensure_single_valid_alignment_per_cdna_per_cluster.pl behavior.
"""

import logging
from typing import Dict
from collections import defaultdict

from ..storage.datastore import DataStore

logger = logging.getLogger(__name__)


def deduplicate_clusters(datastore: DataStore, score_method: str = "weighted") -> Dict[str, int]:
    """
    Ensure only one valid alignment per cDNA per cluster.
    
    When a transcript has multiple alignments within the same cluster,
    keep only the best-scoring one and invalidate the others.
    
    This matches Perl PASA's ensure_single_valid_alignment_per_cdna_per_cluster.pl
    
    Args:
        datastore: DataStore instance
        score_method: Scoring method ('weighted', 'identity', or 'coverage')
    
    Returns:
        Dictionary with statistics about the deduplication process
    """
    logger.info("Deduplicating: ensuring single alignment per cDNA per cluster")
    
    stats = {
        'total_valid': 0,
        'kept': 0,
        'invalidated': 0,
        'clusters_processed': 0,
        'clusters_with_duplicates': 0
    }
    
    # Get all clusters
    clusters = datastore.get_clusters()
    logger.info(f"Processing {len(clusters)} clusters")
    
    alignments_to_invalidate = []
    
    for cluster in clusters:
        cluster_id = cluster.cluster_id
        
        # Get all valid alignments in this cluster
        alignments = [
            aln for aln in datastore.get_alignments_by_cluster(cluster_id)
            if aln.valid
        ]
        
        stats['total_valid'] += len(alignments)
        stats['clusters_processed'] += 1
        
        if not alignments:
            continue
        
        # Group by cDNA
        by_cdna = defaultdict(list)
        for aln in alignments:
            by_cdna[aln.cdna_acc].append(aln)
        
        # For each cDNA with multiple alignments, keep only the best
        has_duplicates = False
        for cdna_acc, aln_list in by_cdna.items():
            if len(aln_list) == 1:
                stats['kept'] += 1
                continue
            
            has_duplicates = True
            
            # Score each alignment
            scored = []
            for aln in aln_list:
                if score_method == "weighted":
                    # Balance identity and coverage
                    score = aln.avg_percent_id * aln.percent_aligned / 100.0
                elif score_method == "identity":
                    score = aln.avg_percent_id
                elif score_method == "coverage":
                    score = aln.percent_aligned
                else:
                    raise ValueError(f"Unknown score_method: {score_method}")
                
                scored.append((score, aln))
            
            # Sort by score descending
            scored.sort(reverse=True, key=lambda x: x[0])
            
            # Keep the best
            best_score, best_aln = scored[0]
            stats['kept'] += 1
            
            logger.debug(
                f"Cluster {cluster_id}, cDNA {cdna_acc}: "
                f"Keeping alignment {best_aln.align_id} (score={best_score:.2f})"
            )
            
            # Invalidate the rest
            for score, aln in scored[1:]:
                logger.debug(
                    f"  Invalidating alignment {aln.align_id} (score={score:.2f})"
                )
                alignments_to_invalidate.append(aln.align_id)
                stats['invalidated'] += 1
        
        if has_duplicates:
            stats['clusters_with_duplicates'] += 1
    
    # Apply invalidations
    if alignments_to_invalidate:
        logger.info(f"Invalidating {len(alignments_to_invalidate)} duplicate alignments")
        datastore.invalidate_alignments(alignments_to_invalidate)
    
    logger.info(
        f"Deduplication complete: {stats['kept']:,} kept, "
        f"{stats['invalidated']:,} invalidated"
    )
    
    return stats


def get_deduplication_report(stats: Dict[str, int]) -> str:
    """
    Generate a human-readable report of the deduplication process.
    
    Args:
        stats: Statistics dictionary from deduplicate_clusters()
    
    Returns:
        Formatted report string
    """
    lines = [
        "=" * 60,
        "CLUSTER DEDUPLICATION REPORT",
        "=" * 60,
        f"Total valid alignments:       {stats['total_valid']:,}",
        f"Kept (unique per cluster):   {stats['kept']:,} ({stats['kept']/stats['total_valid']*100:5.2f}%)",
        f"Invalidated (duplicates):    {stats['invalidated']:,} ({stats['invalidated']/stats['total_valid']*100:5.2f}%)",
        "",
        f"Clusters processed:          {stats['clusters_processed']:,}",
        f"Clusters with duplicates:    {stats['clusters_with_duplicates']:,} "
        f"({stats['clusters_with_duplicates']/stats['clusters_processed']*100:5.2f}%)",
        "=" * 60,
    ]
    
    return "\n".join(lines)
