#!/usr/bin/env python3
"""
Baseline run of the PyPASA pipeline on the minimap2 sample data, all in memory.

Usage: baseline_minimap2.py [CHECKPOINT.json.gz] [--no-assembly]

This is the canonical regression check referenced in HANDOFF.md. It expects
PASApipeline to be checked out next to PyPASA, with minimap2.spliced_alignments.gff3
and genome_sample.fasta in sample_data (untracked run products; see HANDOFF.md).

Expected (as of 2026-10-25): imported 14490, valid 12251, clusters 617, assemblies 817.
That matches a minimap2-only Perl run (no custom_alignments.gff3) alignment for
alignment: validation, spliced orientation, polyA sites (476), cluster membership
and assemblies. Save a checkpoint and run scripts/compare_with_perl_db.py on it
to check.
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pypasa.assembly.pasa_cpp import PasaCppBackend
from pypasa.pipeline import PipelineConfig, run_pipeline

logging.basicConfig(level=logging.WARNING)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoint", nargs="?", help="write the run here (gzipped JSON)")
    ap.add_argument("--no-assembly", action="store_true")
    ap.add_argument("--cpu", type=int, default=4, help="clusters assembled concurrently")
    args = ap.parse_args()

    pasa_home = Path(__file__).resolve().parents[2] / "PASApipeline"
    sample_data = pasa_home / "sample_data"
    config = PipelineConfig(  # sqlite.confs/alignAssembly.config; Launch_PASA_pipeline.pl -T -u -f
        genome_fasta=str(sample_data / "genome_sample.fasta"),
        transcripts_fasta=str(sample_data / "all_transcripts.fasta.clean"),
        alignment_gff3s=[str(sample_data / "minimap2.spliced_alignments.gff3")],
        fl_accs_file=str(sample_data / "FL_accs.txt"),
        untrimmed_transcripts_fasta=str(sample_data / "all_transcripts.fasta"),
        cln_file=str(sample_data / "all_transcripts.fasta.cln"),
        min_percent_id=95.0,
        min_percent_aligned=75.0,
        backend=None if args.no_assembly else PasaCppBackend(str(pasa_home / "bin" / "pasa")),
        assembly_workers=args.cpu,
    )
    result = run_pipeline(config)
    run = result.run

    print(
        f"imported {result.stats['import']['imported']} "
        f"valid {result.stats['validation']['passed']} "
        f"clusters {result.stats['clustering']['total_clusters']} "
        f"assemblies {len(run.assemblies)}"
    )
    print("seconds: " + ", ".join(f"{k} {v:.1f}" for k, v in result.seconds.items()))
    if args.checkpoint:
        run.save(args.checkpoint)


if __name__ == "__main__":
    main()
