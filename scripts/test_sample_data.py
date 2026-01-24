#!/usr/bin/env python3
"""
Quick test script for PyPASA using sample_data.

This script demonstrates the PyPASA pipeline on the provided sample data
and allows comparison with the Perl PASA implementation.
"""

import sys
import logging
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pypasa.storage.datastore import DataStore
from pypasa.alignment.importer import GFF3AlignmentImporter
from pypasa.alignment.validator import AlignmentValidator
from pypasa.alignment.best_selector import select_best_alignments, get_selection_report
from pypasa.alignment.clusterer import AlignmentClusterer
from pypasa.assembly.assembler import PASAAssembler

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Run PyPASA test on sample data."""
    
    # Paths (relative to PyPASA directory)
    base_dir = Path(__file__).parent.parent
    sample_data = base_dir.parent / "PASApipeline" / "sample_data"
    
    # Input files
    alignments_gff3 = sample_data / "custom_alignments.gff3"
    genome_fasta = sample_data / "genome_sample.fasta.gz"
    
    # Output
    output_dir = base_dir / "test_output"
    output_dir.mkdir(exist_ok=True)
    
    database = output_dir / "pypasa_test.sqlite"
    
    logger.info("=" * 70)
    logger.info("PyPASA Test on Sample Data")
    logger.info("=" * 70)
    logger.info(f"Alignments: {alignments_gff3}")
    logger.info(f"Genome: {genome_fasta}")
    logger.info(f"Database: {database}")
    logger.info(f"Output: {output_dir}")
    logger.info("")
    
    try:
        # Initialize datastore
        logger.info("Initializing database...")
        datastore = DataStore(str(database))
        
        # Step 1: Import alignments
        logger.info("\nStep 1/3: Importing alignments from GFF3")
        importer = GFF3AlignmentImporter(datastore)
        num_imported = importer.import_gff3(str(alignments_gff3), aligner="blat")
        logger.info(f"✓ Imported {num_imported:,} alignments")
        
        # Step 2: Validate alignments
        logger.info("\nStep 2/3: Validating alignments")
        validator = AlignmentValidator(
            datastore,
            min_percent_id=95.0,  # Match Perl PASA config
            min_percent_aligned=75.0,  # Match Perl PASA config
        )
        val_stats = validator.validate_alignments()
        logger.info(f"✓ Validated alignments: {val_stats['passed']:,}/{val_stats['total']:,} passed")
        
        # Step 3: Select SINGLE best alignment per cDNA (matching Perl PASA -N 1)
        logger.info("\nStep 3/4: Selecting best alignment per cDNA (globally)")
        best_stats = select_best_alignments(datastore, score_method="weighted")
        logger.info(f"✓ Kept {best_stats['kept']:,} best alignments, discarded {best_stats['discarded']:,}")
        
        # Step 4: Cluster alignments
        logger.info("\nStep 4/4: Clustering alignments")
        clusterer = AlignmentClusterer(
            datastore,
            min_overlap_percent=0.0,  # Default mode
        )
        clust_stats = clusterer.cluster_alignments()
        logger.info(f"✓ Created {clust_stats['total_clusters']:,} clusters")
        
        # Print detailed reports
        print("\n" + "=" * 70)
        print("DETAILED REPORTS")
        print("=" * 70)
        
        print("\n" + validator.get_validation_report())
        print("\n" + get_selection_report(best_stats))
        print("\n" + clusterer.get_clustering_report())
        
        # Database stats
        stats = datastore.get_stats()
        print("\n" + "=" * 60)
        print("FINAL DATABASE STATISTICS")
        print("=" * 60)
        print(f"Total alignments:    {stats['total_alignments']:>10,}")
        print(f"Valid alignments:    {stats['valid_alignments']:>10,}")
        print(f"Clusters:            {stats['clusters']:>10,}")
        print(f"Scaffolds:           {stats['scaffolds']:>10,}")
        print("=" * 60)
        
        # Test: Load a few alignments to verify they're stored correctly
        print("\n" + "=" * 60)
        print("SAMPLE ALIGNMENTS (first scaffold)")
        print("=" * 60)
        scaffolds = datastore.get_scaffolds()
        if scaffolds:
            first_scaffold = scaffolds[0]
            sample_alns = datastore.get_alignments_by_scaffold(first_scaffold, valid_only=True)
            print(f"Scaffold: {first_scaffold}")
            print(f"Valid alignments on this scaffold: {len(sample_alns)}")
            
            if sample_alns:
                for aln in sample_alns[:3]:  # Show first 3
                    print(f"\n  Alignment {aln.align_id}:")
                    print(f"    cDNA: {aln.cdna_acc}")
                    print(f"    Position: {aln.genome_lend:,}-{aln.genome_rend:,} ({aln.orient.value})")
                    print(f"    Segments: {aln.num_segments}")
                    print(f"    %ID: {aln.avg_percent_id:.2f}")
                    print(f"    Cluster: {aln.cluster_id}")
        
        datastore.close()
        
        logger.info("\n" + "=" * 70)
        logger.info("✓ TEST COMPLETED SUCCESSFULLY")
        logger.info("=" * 70)
        logger.info(f"\nResults saved to: {output_dir}")
        logger.info(f"Database: {database}")
        logger.info("\nNext steps:")
        logger.info("  1. Compare these results with Perl PASA output")
        logger.info("  2. Implement assembly step (currently placeholder)")
        logger.info("  3. Add GFF3/GTF output generation")
        
    except Exception as e:
        logger.error(f"\n✗ TEST FAILED: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
