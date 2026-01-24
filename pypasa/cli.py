"""
Command-line interface for PyPASA.

This module provides the main CLI entry point for running the
PyPASA alignment assembly pipeline.
"""

import click
import logging
import sys
from pathlib import Path

from pypasa.storage.datastore import DataStore
from pypasa.alignment.importer import GFF3AlignmentImporter
from pypasa.alignment.validator import AlignmentValidator
from pypasa.alignment.clusterer import AlignmentClusterer
from pypasa.assembly.assembler import PASAAssembler


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.pass_context
def cli(ctx, verbose):
    """PyPASA - Python implementation of PASA alignment assembly pipeline."""
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    setup_logging(verbose)


@cli.command()
@click.option('--genome', '-g', required=True, type=click.Path(exists=True),
              help='Genome FASTA file')
@click.option('--transcripts', '-t', required=True, type=click.Path(exists=True),
              help='Transcript FASTA file')
@click.option('--alignments', '-a', required=True, type=click.Path(exists=True),
              help='Alignments GFF3 file')
@click.option('--database', '-d', required=True, type=click.Path(),
              help='SQLite database path')
@click.option('--output-dir', '-o', default='./pasa_output', type=click.Path(),
              help='Output directory')
@click.option('--min-percent-id', default=95.0, type=float,
              help='Minimum percent identity')
@click.option('--min-percent-aligned', default=75.0, type=float,
              help='Minimum percent aligned')
@click.option('--min-overlap-percent', default=0.0, type=float,
              help='Minimum overlap percent for clustering (0=any overlap)')
@click.option('--pasa-home', type=click.Path(exists=True),
              help='PASA installation directory (for C++ assembler)')
@click.pass_context
def align_assemble(
    ctx,
    genome,
    transcripts,
    alignments,
    database,
    output_dir,
    min_percent_id,
    min_percent_aligned,
    min_overlap_percent,
    pasa_home,
):
    """Run complete alignment assembly pipeline."""
    
    logger = logging.getLogger(__name__)
    logger.info("=" * 70)
    logger.info("PyPASA Alignment Assembly Pipeline")
    logger.info("=" * 70)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    try:
        # Initialize datastore
        logger.info(f"Initializing database: {database}")
        datastore = DataStore(database)
        
        # Step 1: Import alignments
        logger.info("Step 1/4: Importing alignments")
        importer = GFF3AlignmentImporter(datastore)
        num_imported = importer.import_gff3(alignments, aligner="gmap")
        logger.info(f"  Imported {num_imported} alignments")
        
        # Step 2: Validate alignments
        logger.info("Step 2/4: Validating alignments")
        validator = AlignmentValidator(
            datastore,
            min_percent_id=min_percent_id,
            min_percent_aligned=min_percent_aligned,
            genome_fasta=genome,
            transcript_fasta=transcripts,
        )
        val_stats = validator.validate_alignments()
        logger.info(f"  {val_stats['passed']}/{val_stats['total']} alignments passed validation")
        
        # Step 3: Cluster alignments
        # NOTE: Following Perl PASA workflow - cluster ALL valid alignments.
        # A cDNA can have multiple alignments in different clusters.
        # No pre-clustering selection is performed (matches Perl PASA's actual behavior).
        logger.info("Step 3/4: Clustering alignments")
        logger.info("  (Allowing multiple alignments per cDNA - matching Perl PASA)")
        clusterer = AlignmentClusterer(
            datastore,
            min_overlap_percent=min_overlap_percent,
        )
        clust_stats = clusterer.cluster_alignments()
        logger.info(f"  Created {clust_stats['total_clusters']} clusters")
        
        # Step 4: Assemble clusters
        logger.info("Step 4/4: Assembling clusters")
        if pasa_home:
            assembler = PASAAssembler(datastore, pasa_home, genome)
            asm_stats = assembler.assemble_all_clusters(output_dir=output_dir)
            logger.info(f"  Generated {asm_stats['assemblies_generated']} assemblies")
        else:
            logger.warning("  Skipping assembly: --pasa-home not provided")
            logger.info("  To run assembly, provide path to PASA installation")
        
        # Print reports
        logger.info("")
        logger.info("=" * 70)
        logger.info("PIPELINE COMPLETE")
        logger.info("=" * 70)
        print("\n" + validator.get_validation_report())
        print("\n" + clusterer.get_clustering_report())
        if pasa_home:
            print("\n" + assembler.get_assembly_report())
        
        # Print database stats
        stats = datastore.get_stats()
        print("\n" + "=" * 60)
        print("DATABASE STATISTICS")
        print("=" * 60)
        print(f"Total alignments:    {stats['total_alignments']:>10,}")
        print(f"Valid alignments:    {stats['valid_alignments']:>10,}")
        print(f"Clusters:            {stats['clusters']:>10,}")
        print(f"Assemblies:          {stats['assemblies']:>10,}")
        print(f"Scaffolds:           {stats['scaffolds']:>10,}")
        print("=" * 60)
        
        datastore.close()
        
        logger.info(f"\nOutput saved to: {output_path}")
        logger.info(f"Database saved to: {database}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


@cli.command()
@click.option('--database', '-d', required=True, type=click.Path(exists=True),
              help='SQLite database path')
def stats(database):
    """Show database statistics."""
    datastore = DataStore(database)
    stats = datastore.get_stats()
    
    print("=" * 60)
    print("PyPASA DATABASE STATISTICS")
    print("=" * 60)
    print(f"Total alignments:    {stats['total_alignments']:>10,}")
    print(f"Valid alignments:    {stats['valid_alignments']:>10,}")
    print(f"Clusters:            {stats['clusters']:>10,}")
    print(f"Assemblies:          {stats['assemblies']:>10,}")
    print(f"Scaffolds:           {stats['scaffolds']:>10,}")
    print("=" * 60)
    
    datastore.close()


def main():
    """Entry point for CLI."""
    cli(obj={})


if __name__ == '__main__':
    main()
