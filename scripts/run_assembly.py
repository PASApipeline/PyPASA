#!/usr/bin/env python3
"""
Rerun only the assembly stage on a saved run, e.g. to try another backend.

Usage: run_assembly.py IN.json.gz OUT.json.gz [--backend pasa_cpp] [--pasa-bin PATH]

IN must be a checkpoint that has completed clustering (baseline_minimap2.py).
Any assemblies in it are replaced.
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pypasa.assembly.assembler import PASAAssembler
from pypasa.assembly.pasa_cpp import PasaCppBackend
from pypasa.run import PipelineRun

logging.basicConfig(level=logging.WARNING)


def main():
    default_bin = Path(__file__).resolve().parents[2] / "PASApipeline" / "bin" / "pasa"
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoint_in")
    ap.add_argument("checkpoint_out")
    ap.add_argument("--backend", choices=["pasa_cpp"], default="pasa_cpp")
    ap.add_argument("--pasa-bin", default=str(default_bin))
    args = ap.parse_args()

    run = PipelineRun.load(args.checkpoint_in)
    if "clustering" not in run.completed_stages:
        sys.exit(f"{args.checkpoint_in} hasn't completed clustering")

    assembler = PASAAssembler(PasaCppBackend(args.pasa_bin))
    run.assemblies = assembler.assemble_all_clusters(run.alignments)
    if "assembly" not in run.completed_stages:
        run.completed_stages.append("assembly")
    print(assembler.get_assembly_report())
    run.save(args.checkpoint_out)


if __name__ == "__main__":
    main()
