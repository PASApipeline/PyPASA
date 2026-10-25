"""
The PyPASA alignment-assembly pipeline, run entirely in memory.

Stages, in Perl PASA's order (Launch_PASA_pipeline.pl):
  import       GFF3 alignments -> Alignment objects
  fl_marking   flag full-length cDNAs (Perl: -f FL_accs.txt)
  validation   validate_alignments_in_db.dbi
  polya        polyA_site_transcript_mapper.dbi (when polyA-trimmed, Perl: -T -u)
  clustering   reassign_clusters_via_valid_align_coords.dbi
  assembly     assemble_clusters.dbi + assembly_db_loader.dbi

Each stage reads and updates the PipelineRun. Sequence files are loaded once
and shared across stages.
"""

import logging
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from pypasa.alignment.clusterer import AlignmentClusterer
from pypasa.alignment.importer import GFF3AlignmentImporter
from pypasa.alignment.polya import PolyASiteMapper, parse_cln
from pypasa.alignment.validator import AlignmentValidator
from pypasa.assembly.assembler import PASAAssembler
from pypasa.assembly.backend import AssemblerBackend
from pypasa.fasta import read_fasta, read_fasta_lengths
from pypasa.run import PipelineRun

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    genome_fasta: str
    transcripts_fasta: str  # transcripts as aligned (seqclean output when trimmed)
    alignment_gff3s: List[str]
    fl_accs_file: Optional[str] = None
    # polyA stage runs when both are given
    untrimmed_transcripts_fasta: Optional[str] = None
    cln_file: Optional[str] = None
    min_percent_id: float = 95.0  # Perl MIN_AVG_PER_ID
    min_percent_aligned: float = 75.0  # Perl MIN_PERCENT_ALIGNED
    # assembly runs when a backend is given
    backend: Optional[AssemblerBackend] = None
    assembly_workers: int = 1


@dataclass
class PipelineResult:
    run: PipelineRun
    stats: Dict[str, dict] = field(default_factory=dict)
    seconds: Dict[str, float] = field(default_factory=dict)


def run_pipeline(config: PipelineConfig) -> PipelineResult:
    result = PipelineResult(run=PipelineRun())
    run = result.run

    def stage(name, fn):
        t0 = time.perf_counter()
        result.stats[name] = fn() or {}
        result.seconds[name] = time.perf_counter() - t0
        run.completed_stages.append(name)
        logger.info(f"{name}: {result.stats[name]} ({result.seconds[name]:.1f}s)")

    genome = read_fasta(config.genome_fasta)

    def do_import():
        importer = GFF3AlignmentImporter()
        for path in config.alignment_gff3s:
            run.alignments.extend(importer.import_gff3(path))
        return {"imported": len(run.alignments)}
    stage("import", do_import)

    if config.fl_accs_file:
        def do_fl():
            with open(config.fl_accs_file) as fh:
                fl_accs = set(fh.read().split())
            n = 0
            for aln in run.alignments:
                if aln.cdna_acc in fl_accs:
                    aln.is_fli = True
                    n += 1
            return {"fl_alignments": n}
        stage("fl_marking", do_fl)

    stage("validation", lambda: AlignmentValidator(
        genome,
        read_fasta_lengths(config.transcripts_fasta),
        min_percent_id=config.min_percent_id,
        min_percent_aligned=config.min_percent_aligned,
    ).validate(run.alignments))

    if config.untrimmed_transcripts_fasta and config.cln_file:
        def do_polya():
            mapper = PolyASiteMapper(
                parse_cln(config.cln_file), genome, read_fasta(config.untrimmed_transcripts_fasta)
            )
            run.polya_sites = mapper.run(run.alignments)
            return mapper.stats
        stage("polya", do_polya)

    stage("clustering", lambda: AlignmentClusterer(respect_strand=False).cluster_alignments(run.alignments))

    if config.backend:
        def do_assembly():
            assembler = PASAAssembler(config.backend, workers=config.assembly_workers)
            run.assemblies = assembler.assemble_all_clusters(run.alignments)
            return assembler.stats
        stage("assembly", do_assembly)

    return result
