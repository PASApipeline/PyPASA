#!/usr/bin/env python3
"""
Compare PyPASA clustering results with Perl PASA output.
"""

import gzip
import logging
import sqlite3
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_perl_clusters(perl_file):
    """Load cluster assignments from Perl PASA output."""
    clusters = {}  # cdna_acc -> cluster_id
    current_cluster = None
    
    with gzip.open(perl_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('// cluster:'):
                current_cluster = int(line.split(':')[1].strip())
            elif line and not line.startswith('//'):
                # Parse cdna accession (first field)
                cdna_acc = line.split(',')[0]
                clusters[cdna_acc] = current_cluster
    
    return clusters

def load_python_clusters(db_file):
    """Load cluster assignments from PyPASA database."""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Get valid alignments with clusters
    cursor.execute("""
        SELECT align_id, cdna_acc, cluster_id, scaffold, lend, rend, orient
        FROM alignments 
        WHERE valid = 1 AND cluster_id IS NOT NULL
        ORDER BY cluster_id
    """)
    
    alignments = []
    clusters_by_cdna = {}  # cdna_acc -> cluster_id
    clusters_by_id = defaultdict(list)  # cluster_id -> [alignments]
    
    for row in cursor.fetchall():
        align_id, cdna_acc, cluster_id, scaffold, lend, rend, orient = row
        alignments.append({
            'align_id': align_id,
            'cdna_acc': cdna_acc,
            'cluster_id': cluster_id,
            'scaffold': scaffold,
            'lend': lend,
            'rend': rend,
            'orient': orient
        })
        clusters_by_cdna[cdna_acc] = cluster_id
        clusters_by_id[cluster_id].append(cdna_acc)
    
    conn.close()
    
    return clusters_by_cdna, clusters_by_id, alignments

def analyze_clustering_differences(perl_clusters, python_clusters_by_cdna, python_clusters_by_id):
    """Analyze differences between Perl and Python clustering."""
    
    logger.info("=" * 80)
    logger.info("CLUSTERING COMPARISON ANALYSIS")
    logger.info("=" * 80)
    
    # Count clusters
    perl_cluster_ids = set(perl_clusters.values())
    python_cluster_ids = set(python_clusters_by_cdna.values())
    
    logger.info(f"\nCluster counts:")
    logger.info(f"  Perl PASA:  {len(perl_cluster_ids):,} clusters")
    logger.info(f"  PyPASA:     {len(python_cluster_ids):,} clusters")
    logger.info(f"  Difference: {len(python_cluster_ids) - len(perl_cluster_ids):,}")
    
    # Find cdnas in both
    perl_cdnas = set(perl_clusters.keys())
    python_cdnas = set(python_clusters_by_cdna.keys())
    
    common_cdnas = perl_cdnas & python_cdnas
    perl_only = perl_cdnas - python_cdnas
    python_only = python_cdnas - perl_cdnas
    
    logger.info(f"\nAlignment counts:")
    logger.info(f"  Perl PASA:  {len(perl_cdnas):,} alignments")
    logger.info(f"  PyPASA:     {len(python_cdnas):,} alignments")
    logger.info(f"  Common:     {len(common_cdnas):,} alignments")
    logger.info(f"  Perl only:  {len(perl_only):,} alignments")
    logger.info(f"  Python only: {len(python_only):,} alignments")
    
    # Analyze cluster size distributions
    logger.info(f"\nCluster size distributions:")
    
    perl_cluster_sizes = defaultdict(int)
    for cluster_id in perl_cluster_ids:
        size = sum(1 for cdna, cid in perl_clusters.items() if cid == cluster_id)
        perl_cluster_sizes[size] += 1
    
    python_cluster_sizes = defaultdict(int)
    for cluster_id, cdnas in python_clusters_by_id.items():
        size = len(cdnas)
        python_cluster_sizes[size] += 1
    
    logger.info(f"\nPerl PASA cluster sizes:")
    for size in sorted(perl_cluster_sizes.keys())[:10]:  # Top 10 sizes
        logger.info(f"  Size {size:3d}: {perl_cluster_sizes[size]:,} clusters")
    
    logger.info(f"\nPyPASA cluster sizes:")
    for size in sorted(python_cluster_sizes.keys())[:10]:  # Top 10 sizes
        logger.info(f"  Size {size:3d}: {python_cluster_sizes[size]:,} clusters")
    
    # Check if PyPASA is over-clustering (splitting Perl clusters)
    logger.info(f"\n" + "=" * 80)
    logger.info("Checking if PyPASA split Perl clusters...")
    logger.info("=" * 80)
    
    # Group Python alignments by Perl cluster
    perl_to_python = defaultdict(set)  # perl_cluster_id -> set of python_cluster_ids
    for cdna in common_cdnas:
        perl_cid = perl_clusters[cdna]
        python_cid = python_clusters_by_cdna[cdna]
        perl_to_python[perl_cid].add(python_cid)
    
    split_clusters = {pcid: pycids for pcid, pycids in perl_to_python.items() if len(pycids) > 1}
    
    logger.info(f"Found {len(split_clusters):,} Perl clusters split into multiple Python clusters")
    
    if split_clusters:
        logger.info(f"\nExample split clusters (showing first 5):")
        for i, (perl_cid, python_cids) in enumerate(list(split_clusters.items())[:5]):
            cdnas_in_perl = [cdna for cdna, cid in perl_clusters.items() if cid == perl_cid]
            logger.info(f"\n  Perl cluster {perl_cid}: {len(cdnas_in_perl)} alignments")
            logger.info(f"    → Split into {len(python_cids)} Python clusters: {sorted(python_cids)}")
            
            # Show distribution
            for pycid in sorted(python_cids):
                count = sum(1 for cdna in cdnas_in_perl if cdna in python_clusters_by_cdna 
                           and python_clusters_by_cdna[cdna] == pycid)
                logger.info(f"      Python cluster {pycid}: {count} alignments")
    
    # Check if PyPASA is under-clustering (merging Perl clusters)
    logger.info(f"\n" + "=" * 80)
    logger.info("Checking if PyPASA merged Perl clusters...")
    logger.info("=" * 80)
    
    # Group Perl alignments by Python cluster
    python_to_perl = defaultdict(set)  # python_cluster_id -> set of perl_cluster_ids
    for cdna in common_cdnas:
        perl_cid = perl_clusters[cdna]
        python_cid = python_clusters_by_cdna[cdna]
        python_to_perl[python_cid].add(perl_cid)
    
    merged_clusters = {pycid: pcids for pycid, pcids in python_to_perl.items() if len(pcids) > 1}
    
    logger.info(f"Found {len(merged_clusters):,} Python clusters containing multiple Perl clusters")
    
    if merged_clusters:
        logger.info(f"\nExample merged clusters (showing first 5):")
        for i, (python_cid, perl_cids) in enumerate(list(merged_clusters.items())[:5]):
            cdnas_in_python = python_clusters_by_id[python_cid]
            logger.info(f"\n  Python cluster {python_cid}: {len(cdnas_in_python)} alignments")
            logger.info(f"    ← Merged from {len(perl_cids)} Perl clusters: {sorted(perl_cids)}")
            
            # Show distribution
            for pcid in sorted(perl_cids):
                count = sum(1 for cdna in cdnas_in_python if cdna in perl_clusters 
                           and perl_clusters[cdna] == pcid)
                logger.info(f"      Perl cluster {pcid}: {count} alignments")

def main():
    # Paths
    base_dir = Path(__file__).parent.parent
    perl_file = base_dir.parent / "PASApipeline" / "sample_data" / "clusters_of_valid_alignments.txt.gz"
    python_db = base_dir / "test_output" / "pypasa_test.sqlite"
    
    logger.info(f"Loading Perl PASA clusters from: {perl_file}")
    perl_clusters = load_perl_clusters(perl_file)
    
    logger.info(f"Loading PyPASA clusters from: {python_db}")
    python_clusters_by_cdna, python_clusters_by_id, alignments = load_python_clusters(python_db)
    
    # Analyze differences
    analyze_clustering_differences(perl_clusters, python_clusters_by_cdna, python_clusters_by_id)

if __name__ == '__main__':
    main()
