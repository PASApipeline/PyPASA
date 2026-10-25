# PyPASA handoff

Last updated 2026-10-25. This is the current reference for the port; everything in
`docs/history/` is an older lab notebook, worth reading for *why* but not for
current behavior or numbers.

## What this is

A Python port of the PASA alignment-assembly pipeline (`../PASApipeline`, Perl).
It is a **modernization**, not a line-by-line translation: PASA was designed 20+
years ago around a MySQL/SQLite database and much weaker machines. PyPASA keeps
PASA's results and changes how they're computed.

## Goals

1. **Same results as Perl PASA.** Every stage is checked against a Perl run,
   alignment by alignment, and must match. Any deviation is deliberate,
   documented below, and kept to labels/ordering or to places where Perl is
   effectively random.
2. **In memory, no database.** Stages are plain functions over in-memory
   objects. Persistence is an optional checkpoint file, not the backbone.
3. **Modular.** Each stage can be run, replaced and compared on its own. The
   assembler is a pluggable backend (C++ `pasa` now; LRAA is under
   consideration as the assembly step, see "Direction under discussion").
4. **Fast on today's hardware.** Load inputs once, parallelize where results
   don't depend on order.
5. **A modern results browser later.** PASA's web portal will be rebuilt from
   scratch on top of the run checkpoint/outputs, not ported.

The port is expected to go beyond alignment assembly (annotation comparison and
updates, alternative splicing, gene models); see "Direction under discussion".

## Status

On the PASApipeline sample data with minimap2 alignments, PyPASA matches a
minimap2-only Perl run exactly at every stage:

| Stage | Perl script(s) | PyPASA | Perl / PyPASA |
|---|---|---|---|
| Import | `import_spliced_alignments.dbi` | `alignment/importer.py` | 14,490 alignments |
| FL marking | `upload_transcript_data.dbi -f` | `pipeline.py` | |
| Validation | `validate_alignments_in_db.dbi` + `CDNA_alignment.pm` checks | `alignment/validator.py` | 12,251 valid; same set, same spliced orientations |
| PolyA sites | `polyA_site_transcript_mapper.dbi` | `alignment/polya.py` | 476 identical sites; sets strand of 154 single-exon alignments |
| Clustering | `reassign_clusters_via_valid_align_coords.dbi` (`Overlap_piler.pm`) | `alignment/clusterer.py` | 617 identical clusters |
| Assembly | `assemble_clusters.dbi`, `PASA_alignment_assembler.pm`, `assembly_db_loader.dbi` | `assembly/` | 817 identical assemblies (structure, orientation, members, `asmbl_N` names) |
| Subclusters | `subcluster_builder.dbi` + `subcluster_loader.dbi` | `assembly/subclusters.py` | 665 identical subclusters |
| Outputs | `PASA_transcripts_and_assemblies_to_GFF3.dbi`, `populate_mysql_assembly_*`, `describe_alignment_assemblies_cgi_convert.dbi` | `output.py` | GFF3, FASTA byte-identical; BED, GTF, described.txt identical up to the deviations below |

Speed on the same (heavily loaded) host: the whole pipeline takes ~1.5 s
(Perl: ~21 min to the same point). Assembly takes ~1.6 s with one worker and
~0.7 s with `--cpu 4`, using `pasa`'s batch mode (`-B`), which assembles many
strand passes in one process. With a `pasa` built without `-B` it falls back to
one process per pass (1,234 on the sample data, ~37 s).

## Running it

Expected layout (scripts use relative paths):

```
PASA_port/
  PASApipeline/   git@github.com:PASApipeline/PASApipeline.git  (cc1f8d7; `make` builds bin/pasa;
                  batch mode -B is on local branch pasa-batch-mode, not yet pushed)
  PyPASA/         git@github.com:PASApipeline/PyPASA.git
```

`pip install -e 'PyPASA[dev]'` (Python >= 3.9), then:

```bash
pypasa align-assemble -g genome.fa -t transcripts.fasta.clean -a alignments.gff3 \
    -u transcripts.fasta --cln transcripts.fasta.cln -f FL_accs.txt \
    --pasa-bin PASApipeline/bin/pasa --cpu 4 -o run.json.gz --out-prefix out/mydb
pypasa stats run.json.gz
```

Regression check on the sample data, against the Perl ground truth:

```bash
python3 PyPASA/scripts/baseline_minimap2.py run.json.gz --out-prefix out/mm2only_pasa.sqlite
#  imported 14490 valid 12251 clusters 617 assemblies 817
python3 PyPASA/scripts/compare_with_perl_db.py run.json.gz perl_ref/mm2only_pasa.sqlite
python3 PyPASA/scripts/compare_outputs_with_perl.py out/mm2only_pasa.sqlite perl_ref/mm2only_pasa.sqlite
python3 -m pytest PyPASA/tests
```

`scripts/run_assembly.py IN.json.gz OUT.json.gz` reruns only assembly and
subclustering on a saved run (e.g. to try another backend).

The sample data needs `minimap2.spliced_alignments.gff3` and
`genome_sample.fasta`, which are run products, not tracked in PASApipeline:
run the Perl reference below once and copy them from its output directory into
`PASApipeline/sample_data/`.

### Perl ground truth

`scripts/perl_reference/run_perl_mm2only.sh PASA_PORT_DIR OUT_DIR` builds it
(~25 min). Use this, not the stock `sample_data/runMe.SQLite.sh`: that also
imports `custom_alignments.gff3` into the same database, and those alignments are
clustered and assembled together with the minimap2 ones, so its assemblies can't
be compared with a minimap2-only run. Perl validation parameters come from
`sqlite.confs/alignAssembly.config`: `MIN_PERCENT_ALIGNED=75`,
`MIN_AVG_PER_ID=95`, `NUM_BP_PERFECT_SPLICE_BOUNDARY=0`,
`subcluster_builder.dbi -m=50`.

Host notes (Broad methods hosts): no system `g++`, so build PASApipeline with
gcc 13 from `/broad/software/free/Linux/redhat_7_x86_64/pkgs/gcc_13.2.0` and link
`slclust` with `-static-libstdc++ -static-libgcc`. The system perl lacks
`threads.pm` and `~/perl5`'s DBD::SQLite is built for a different perl, so use a
conda perl with `PERL5LIB` cleared (see the script header).

## Design rules (settled)

1. **Every scaffold is distinct.** `gi|68711` and `gi|rev68711` are never grouped.
2. **Initial clustering ignores strand** (Perl's default
   `reassign_clusters_via_valid_align_coords.dbi`). Strand separation happens in
   assembly (one pass per strand, as `force_flexorient('+')`/`('-')`) and in
   subclustering.
3. **All valid alignments are clustered**, including several per cDNA.
4. **Assembly membership is per alignment** (GFF3 `ID` = Perl `align_acc`), not
   per cDNA.
5. **Keep the C++ `pasa` assembler** rather than rewriting it in Python; other
   backends implement `assembly/backend.py`.
6. **Input order matters to `pasa`.** Alignments go in sorted by
   (lend, rend, import order), as Perl's indexed query returns them.

## Deliberate deviations from Perl

- **Subcluster numbering** is by cluster, then first assembly. Perl's follows
  `slclust`'s traversal and unordered SQL results. Only labels change (GTF
  `gene_id`, described.txt column 2).
- **BED** lines are in assembly order; Perl writes them in hash order.
- **Unspliced-assembly strand ties** (equal +/- aligned members, no spliced
  orientation) go to `+`; Perl breaks them by hash order, i.e. arbitrarily.

## Not ported yet

- Subsampling of clusters with > 5,000 alignments (`assemble_clusters.dbi -n`,
  `select_aligns_to_assemble`; partly random in Perl). PyPASA assembles all and
  warns.
- `NUM_BP_PERFECT_SPLICE_BOUNDARY > 0` (raises `NotImplementedError`).
- Clustering alternatives: `--stringent_alignment_overlap`, `--gene_overlap`;
  `--INVALIDATE_SINGLE_EXON_ESTS`; `--transcribed_is_aligned_orient`
  (strand-specific RNA-seq).
- Inputs: running aligners (PyPASA takes GFF3), `--trans_gtf`, TDN accessions,
  multiple aligners combined (`-N`, plus the `custom_alignments.gff3` sample).
- Downstream: `alignment_assembly_to_gene_models.dbi`,
  `build_comprehensive_transcriptome.dbi`, TransDecoder, annotation comparison,
  alt-splicing analysis, valid/failed alignment exports.

## LRAA as the assembler

LRAA (Long Read Alignment Assembler) has been evaluated as a replacement for the
`pasa` assembly step:

- `../LRAA_chaining/REPORT.md`: stock and patched LRAA chaining (patched reached
  ~96% identical assemblies; see the LRAA report for a correction to its numbers).
- `/broad/hptmp/bhaas/PASAPIPELINE_PAPER/LRAA_pasa_mode/REPORT.md`: an LRAA branch
  (`pasa_mode`, local only) whose `--pasa_mode` runs a port of PASA's chaining
  inside LRAA, from one whole-genome BAM. It reproduces Perl's 817 assemblies
  (structures and members) on the sample data, in ~28 s.

## Direction under discussion (paused here, 2026-10-25)

This commit is a checkpoint before a larger refactor, to be planned with the PASA
author. The idea being considered:

- **LRAA does alignment assembly.** LRAA already covers everything through
  assembly and also gives a read-to-isoform tracking file and expression values.
  PyPASA would use it (native chaining or `--pasa_mode`, undecided) instead of
  porting more of PASA's assembly machinery. The current `pasa` backend stays as
  the reference implementation and regression fixture.
- **PyPASA's effort moves downstream:** annotation comparison and updates
  (`cDNA_annotation_comparer.dbi`, `cDNA_annotation_updater.dbi`), alternative
  splicing (`classify_alt_splice_isoforms*.dbi`, `report_alt_splicing_findings.dbi`),
  gene models and the comprehensive transcriptome
  (`alignment_assembly_to_gene_models.dbi`, `build_comprehensive_transcriptome.dbi`).

Points raised so far:

- Exact parity with Perl's *assemblies* would no longer be a goal if LRAA defines
  them (unless `--pasa_mode`). Downstream modules should still be verified
  against Perl first by feeding them the assemblies we already match (the
  `sample_data/orig_annotations_sample.gff3` annotation gives the ground truth),
  then switched to LRAA.
- Downstream modules should read one assembly-set model (assemblies, members with
  FL status, gene grouping, optional expression) that any assembler fills.
- PASA's pre-assembly steps (validation, FL marking, polyA sites) don't exist in
  LRAA as such; keeping PyPASA's as a filter that writes LRAA's BAM is one option.
- Expression values are meaningful for long reads, less so for Trinity/EST/cDNA
  inputs; which inputs PyPASA centers on is open.
- LRAA's gene grouping might replace subclusters.

Other pending items:

1. `pasa -B` batch mode is on PASApipeline local branch `pasa-batch-mode`
   (cbdd7d2), not pushed.
2. The "not ported" items above, if still wanted after the refactor plan.
3. A second dataset (e.g. `custom_alignments.gff3`, or minimap2 + custom combined).

## Code map

- `pypasa/pipeline.py`: stage order, config, timing. Start here.
- `pypasa/run.py`: in-memory run state and the checkpoint format.
- `pypasa/models.py`: `Alignment`, `Assembly`, segments.
- `pypasa/alignment/`: importer, validator, polyA, clusterer.
- `pypasa/assembly/`: `backend.py` (interface), `pasa_cpp.py`,
  `cluster_assembler.py` (Perl's per-cluster logic), `assembler.py` (driver),
  `subclusters.py`.
- `pypasa/output.py`: GFF3/GTF/BED/FASTA/described writers.
- `scripts/`: baseline run, stage rerun, Perl comparisons, Perl reference run.
