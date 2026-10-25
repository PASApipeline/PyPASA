#!/usr/bin/env python3
"""
Baseline run of PyPASA (import -> validate -> cluster) on the minimap2 sample data.

This is the canonical regression check referenced in HANDOFF.md. It expects
PASApipeline to be checked out next to PyPASA, and the Perl sample pipeline to
have been run (sample_data/runMe.SQLite.sh) so that minimap2.spliced_alignments.gff3
and genome_sample.fasta exist (they are not tracked in the PASApipeline repo).

Expected (as of 2026-10-25): imported 14490, valid 12465, clusters 615
Perl reference:              valid 12251 (sample_mydb_pasa.sqlite.valid_minimap2_alignments.gff3)

Takes ~4 minutes, mostly in validation.
"""

import sys
import logging
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pypasa.storage.datastore import DataStore
from pypasa.alignment.importer import GFF3AlignmentImporter
from pypasa.alignment.validator import AlignmentValidator
from pypasa.alignment.clusterer import AlignmentClusterer

logging.basicConfig(level=logging.WARNING)


def main():
    base_dir = Path(__file__).parent.parent
    sample_data = base_dir.parent / "PASApipeline" / "sample_data"
    database = Path(sys.argv[1]) if len(sys.argv) > 1 else base_dir / "test_output" / "baseline_minimap2.sqlite"

    database.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("", "-wal", "-shm"):
        Path(str(database) + suffix).unlink(missing_ok=True)

    datastore = DataStore(str(database))

    num_imported = GFF3AlignmentImporter(datastore).import_gff3(
        str(sample_data / "minimap2.spliced_alignments.gff3"), aligner="minimap2"
    )
    val_stats = AlignmentValidator(
        datastore,
        min_percent_id=95.0,  # sqlite.confs/alignAssembly.config
        min_percent_aligned=75.0,
        genome_fasta=str(sample_data / "genome_sample.fasta"),
        transcript_fasta=str(sample_data / "all_transcripts.fasta.clean"),
    ).validate_alignments()
    clust_stats = AlignmentClusterer(datastore, respect_strand=False).cluster_alignments()

    print(f"imported {num_imported} valid {val_stats['passed']} clusters {clust_stats['total_clusters']}")
    datastore.close()


if __name__ == "__main__":
    main()
