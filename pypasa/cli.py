"""
Command-line interface for PyPASA.

This module provides the main CLI entry point for running the
PyPASA alignment assembly pipeline.
"""

import click
import logging
import sys

from pypasa.assembly.pasa_cpp import PasaCppBackend
from pypasa.fasta import read_fasta
from pypasa.output import write_all
from pypasa.pipeline import PipelineConfig, run_pipeline
from pypasa.run import PipelineRun


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
              help='Transcript FASTA file, as aligned (seqclean output if trimmed)')
@click.option('--alignments', '-a', required=True, multiple=True, type=click.Path(exists=True),
              help='Alignments GFF3 file (repeatable)')
@click.option('--untrimmed-transcripts', '-u', type=click.Path(exists=True),
              help='Transcripts before seqclean trimming (enables polyA-site analysis with --cln)')
@click.option('--cln', type=click.Path(exists=True),
              help="seqclean's .cln report (enables polyA-site analysis with -u)")
@click.option('--fl-accs', '-f', type=click.Path(exists=True),
              help='File of full-length cDNA accessions')
@click.option('--min-percent-id', default=95.0, type=float,
              help='Minimum average percent identity')
@click.option('--min-percent-aligned', default=75.0, type=float,
              help='Minimum percent of transcript aligned')
@click.option('--pasa-bin', type=click.Path(exists=True),
              help='C++ pasa assembler binary (PASApipeline/bin/pasa); assembly is skipped without it')
@click.option('--cpu', default=1, type=int,
              help='Clusters assembled concurrently')
@click.option('--checkpoint', '-o', type=click.Path(),
              help='Write the finished run here (gzipped JSON)')
@click.option('--out-prefix', type=click.Path(),
              help='Write assembly GFF3/GTF/BED/FASTA/described files as <prefix>.<suffix>')
def align_assemble(
    genome,
    transcripts,
    alignments,
    untrimmed_transcripts,
    cln,
    fl_accs,
    min_percent_id,
    min_percent_aligned,
    pasa_bin,
    cpu,
    checkpoint,
    out_prefix,
):
    """Run the alignment assembly pipeline."""
    logger = logging.getLogger(__name__)
    if bool(untrimmed_transcripts) != bool(cln):
        sys.exit("--untrimmed-transcripts and --cln go together")

    config = PipelineConfig(
        genome_fasta=genome,
        transcripts_fasta=transcripts,
        alignment_gff3s=list(alignments),
        fl_accs_file=fl_accs,
        untrimmed_transcripts_fasta=untrimmed_transcripts,
        cln_file=cln,
        min_percent_id=min_percent_id,
        min_percent_aligned=min_percent_aligned,
        backend=PasaCppBackend(pasa_bin) if pasa_bin else None,
        assembly_workers=cpu,
    )
    if not pasa_bin:
        logger.warning("No --pasa-bin given: skipping assembly")

    result = run_pipeline(config)
    for stage, stats in result.stats.items():
        print(f"{stage:<12} {result.seconds[stage]:>7.1f}s  {stats}")

    if checkpoint:
        result.run.save(checkpoint)
        logger.info(f"Run saved to: {checkpoint}")
    if out_prefix and result.run.assemblies:
        run_name = out_prefix.rstrip("/").split("/")[-1]
        for path in write_all(result.run, out_prefix, run_name, read_fasta(genome)):
            logger.info(f"Wrote {path}")


@cli.command()
@click.argument('checkpoint', type=click.Path(exists=True))
def stats(checkpoint):
    """Summarize a saved run."""
    run = PipelineRun.load(checkpoint)
    valid = [a for a in run.alignments if a.valid]
    print(f"Completed stages:    {', '.join(run.completed_stages)}")
    print(f"Total alignments:    {len(run.alignments):>10,}")
    print(f"Valid alignments:    {len(valid):>10,}")
    print(f"PolyA sites:         {len(run.polya_sites):>10,}")
    print(f"Clusters:            {len({a.cluster_id for a in valid if a.cluster_id is not None}):>10,}")
    print(f"Assemblies:          {len(run.assemblies):>10,}")
    print(f"Scaffolds:           {len({a.scaffold for a in run.alignments}):>10,}")


def main():
    """Entry point for CLI."""
    cli(obj={})


if __name__ == '__main__':
    main()
