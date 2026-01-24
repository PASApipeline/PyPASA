"""
Best alignment selection module.

This module implements the best-alignment selection logic to match Perl PASA behavior.
For each cDNA transcript, only the single best alignment globally is retained.
All scaffolds (including gi|XXX and gi|revXXX) are treated as completely distinct.
"""

import logging
from typing import Dict, List, Tuple
from collections import defaultdict

from ..storage.datastore import DataStore
from ..storage.models import Alignment

logger = logging.getLogger(__name__)


def select_best_alignments(
    datastore: DataStore,
    score_method: str = "weighted"
) -> Dict[str, int]:
    """
    Select the SINGLE best alignment for each cDNA GLOBALLY.
    
    This matches Perl PASA behavior with -N 1: when a transcript aligns to 
    multiple locations (different scaffolds or positions), only the single
    best-scoring alignment across ALL locations is kept. All others are 
    marked as invalid.
    
    Args:
        datastore: DataStore instance
        score_method: Scoring method ('weighted', 'identity', or 'coverage')
            - 'weighted': percent_id * percent_aligned (default, best balance)
            - 'identity': percent_id only (prioritize accuracy)
            - 'coverage': percent_aligned only (prioritize coverage)
    
    Returns:
        Dictionary with statistics about the selection process
    """
    logger.info("Selecting SINGLE best alignment per cDNA (globally)")
    
    stats = {
        'total_valid': 0,
        'kept': 0,
        'discarded': 0,
        'cdnas_processed': 0,
        'cdnas_with_multiple': 0
    }
    
    scaffolds = datastore.get_scaffolds()
    logger.info(f"Processing {len(scaffolds)} scaffold(s)")
    
    # Load ALL valid alignments (across all scaffolds)
    all_alignments = []
    for scaffold in scaffolds:
        alignments = datastore.get_alignments_by_scaffold(scaffold)
        all_alignments.extend([a for a in alignments if a.valid])
    
    stats['total_valid'] = len(all_alignments)
    logger.info(f"Loaded {stats['total_valid']:,} valid alignments")
    
    # Group by cDNA ONLY (not by scaffold) - GLOBAL best selection
    by_cdna = defaultdict(list)
    for aln in all_alignments:
        by_cdna[aln.cdna_acc].append(aln)
    
    stats['cdnas_processed'] = len(by_cdna)
    logger.info(f"Found {stats['cdnas_processed']:,} unique cDNAs")
    
    # For each cDNA, select the SINGLE best alignment across ALL scaffolds
    alignments_to_invalidate = []
    
    for cdna_acc, aln_list in by_cdna.items():
        if len(aln_list) == 1:
            # Only one alignment for this cDNA on this scaffold, keep it
            stats['kept'] += 1
            continue
        
        # Multiple alignments - need to select the SINGLE best
        stats['cdnas_with_multiple'] += 1
        # This happens when:
        # 1. Same cDNA aligns to multiple scaffolds (treated as distinct)
        # 2. Same cDNA aligns to multiple locations on same scaffold
        # 3. Same cDNA aligns to both gi|XXX and gi|revXXX (treated as distinct scaffolds)
        
        # Calculate scores
        scored_alignments = []
        for aln in aln_list:
            score = _calculate_score(aln, score_method)
            scored_alignments.append((score, aln))
        
        # Sort by score (descending) and keep ONLY the best
        scored_alignments.sort(reverse=True, key=lambda x: x[0])
        best_aln = scored_alignments[0][1]
        
        stats['kept'] += 1
        
        logger.debug(
            f"cDNA {cdna_acc}: Keeping BEST alignment {best_aln.align_id} "
            f"at {best_aln.scaffold}:{best_aln.genome_lend}-{best_aln.genome_rend} "
            f"(score={scored_alignments[0][0]:.2f}, discarding {len(scored_alignments)-1} others)"
        )
        
        # Mark all others as invalid
        for score, aln in scored_alignments[1:]:
            alignments_to_invalidate.append(aln.align_id)
            stats['discarded'] += 1
    
    stats['cdnas_processed'] = len(by_cdna)
    
    # Batch invalidate the discarded alignments
    if alignments_to_invalidate:
        datastore.invalidate_alignments(alignments_to_invalidate)
    
    logger.info(
        f"Best alignment selection complete: "
        f"{stats['kept']:,} kept, {stats['discarded']:,} discarded "
        f"({stats['discarded']*100//stats['total_valid'] if stats['total_valid'] > 0 else 0}%)"
    )
    
    return stats


def _calculate_score(alignment: Alignment, method: str) -> float:
    """
    Calculate alignment score based on the specified method.
    
    Args:
        alignment: Alignment to score
        method: Scoring method name
    
    Returns:
        Score value (higher is better)
    """
    percent_id = alignment.avg_percent_id or 0.0
    percent_aligned = alignment.percent_aligned or 0.0
    
    if method == "weighted":
        # Weighted combination (default): percent_id * percent_aligned
        # This balances both accuracy and coverage
        return percent_id * percent_aligned / 100.0
    
    elif method == "identity":
        # Prioritize alignment accuracy
        # Use percent_aligned as tiebreaker
        return percent_id + (percent_aligned / 10000.0)
    
    elif method == "coverage":
        # Prioritize alignment coverage
        # Use percent_id as tiebreaker
        return percent_aligned + (percent_id / 10000.0)
    
    else:
        raise ValueError(f"Unknown scoring method: {method}")


def get_selection_report(stats: Dict[str, int]) -> str:
    """
    Generate a human-readable report of the selection process.
    
    Args:
        stats: Statistics dictionary from select_best_alignments
    
    Returns:
        Formatted report string
    """
    lines = [
        "=" * 60,
        "BEST ALIGNMENT SELECTION REPORT",
        "=" * 60,
        f"Total valid alignments:      {stats['total_valid']:>10,}",
        f"Alignments kept:             {stats['kept']:>10,} ({stats['kept']*100//stats['total_valid'] if stats['total_valid'] > 0 else 0:>3}%)",
        f"Alignments discarded:        {stats['discarded']:>10,} ({stats['discarded']*100//stats['total_valid'] if stats['total_valid'] > 0 else 0:>3}%)",
        "",
        f"cDNAs processed:             {stats['cdnas_processed']:>10,}",
        f"cDNAs with multiple aligns:  {stats['cdnas_with_multiple']:>10,}",
        "=" * 60,
    ]
    
    return "\n".join(lines)
