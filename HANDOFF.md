# PyPASA handoff (2026-10-25)

Read this first. It overrides the other markdown files in this repo where they disagree.

## What this is

A Python port of the PASA alignment-assembly pipeline (`../PASApipeline`, Perl). The goal is to match Perl behavior; Perl PASA is ground truth. Work started 2026-01-24 with a single agent session and stopped after import, validation, and clustering. Assembly, subclustering, and output writers remain.

Expected layout (both repos side by side; scripts use relative paths):

```
PASA_port/
  PASApipeline/   git@github.com:PASApipeline/PASApipeline.git  (was at cc1f8d7)
  PyPASA/         git@github.com:PASApipeline/PyPASA.git        (main)
```

## State of the code

| Stage | File | State |
|---|---|---|
| GFF3 import | `pypasa/alignment/importer.py` | Done |
| Validation (%id, %aligned, splice orient) | `pypasa/alignment/validator.py` | Done; ~1.7% more valid than Perl (open issue below) |
| Clustering | `pypasa/alignment/clusterer.py` | Done; expanding-window port of `PerlLib/Overlap_piler.pm` |
| Assembly | `pypasa/assembly/assembler.py` | **Placeholder.** `_create_placeholder_assembly` (line 209) emits the longest alignment; `_call_pasa_binary` (250) and `_parse_assembly_output` (287) are not wired up |
| Subclustering | — | Not started |
| GFF3/GTF/BED/FASTA writers | — | Not started (no `exporter.py`, despite ARCHITECTURE.md) |
| Unit tests | `tests/` | Empty |
| Storage | `pypasa/storage/{models,datastore}.py` | SQLite, works |

`pypasa/alignment/best_selector.py` and `cluster_deduplicator.py` are vestigial. Nothing calls them. They came from an early "keep best alignment per cDNA" idea that was dropped once it was clear Perl clusters *all* valid alignments.

## Design rules that are settled (don't relitigate)

From ARCHITECTURE.md "Critical Design Rules", confirmed against the code:

1. **Every scaffold is distinct.** `gi|68711` and `gi|rev68711` are not grouped. Grouping them was tried and gave 361 clusters vs Perl's 611.
2. **Initial clustering is orientation-agnostic** (`respect_strand=False`), matching Perl's default `reassign_clusters_via_valid_align_coords.dbi`. Mixed-strand clusters are expected at this stage. Strand separation happens later: in assembly (Perl runs `force_flexorient('+')`, then `('-')`) and in subclustering.
3. **All valid alignments are clustered**, including multiple alignments per cDNA.
4. **Keep the C++ `pasa` assembler**; don't rewrite it in Python.

## Verified baseline

Run on 2026-10-25 on the current `main`:

```bash
python3 scripts/baseline_minimap2.py      # ~4 min, mostly validation
# imported 14490 valid 12465 clusters 615
```

| | PyPASA | Perl |
|---|---|---|
| minimap2 alignments imported | 14,490 | 14,490 |
| valid | 12,465 | 12,251 (IDs in `sample_mydb_pasa.sqlite.valid_minimap2_alignments.gff3`) |
| clusters | 615 | 611 per ARCHITECTURE.md* |

\*The 611 comparison was made by loading Perl's valid set into PyPASA first, so it isolates clustering from validation differences. The Perl DB it came from no longer exists (see below), so re-confirm it.

**The canonical dataset is minimap2** (`minimap2.spliced_alignments.gff3`). The older docs and `scripts/test_sample_data.py` / `compare_with_perl.py` use `custom_alignments.gff3` (22,784 alignments, Perl 605 clusters). That's a different dataset with different reference numbers, so don't mix them.

## Getting Perl ground truth on a new machine (do this first)

The Perl reference DB was at `/tmp/sample_mydb_pasa.sqlite` and is gone. Most debug snippets in the docs query it. A fresh PASApipeline clone also **lacks** `minimap2.spliced_alignments.gff3`, the uncompressed `genome_sample.fasta`, and every `sample_mydb_pasa.sqlite.*` output, because those are untracked run products.

**Easiest path:** copy the original machine's `PASApipeline/sample_data/` over the fresh clone's (rsync). That brings across the minimap2 alignments and all exported Perl outputs from the reference run. It does *not* bring the SQLite DB, so regenerate that if you need cluster-level comparisons.

**Otherwise, rerun Perl PASA:**

```bash
git clone --recurse-submodules git@github.com:PASApipeline/PASApipeline.git   # TransDecoder submodule is needed
cd PASApipeline && make
cd sample_data
# 1. Restrict aligners to minimap2 (the original reference run had this uncommitted local edit):
#    __run_sample_pipeline.pl:  my $ALIGNERS = "minimap2";   (was "blat,gmap,minimap2")
# 2. Optionally point DATABASE= in sqlite.confs/alignAssembly.config somewhere other than /tmp
./runMe.SQLite.sh          # requires minimap2 (and samtools, seqclean deps) on PATH
gunzip -k genome_sample.fasta.gz   # if genome_sample.fasta wasn't produced; the baseline reads it uncompressed
```

Note that the sample run imports *both* minimap2 alignments and `custom_alignments.gff3` into the same Perl DB (`--IMPORT_CUSTOM_ALIGNMENTS_GFF3`), so filter Perl tables by the minimap2 alignments when comparing.

Perl validation params (`sqlite.confs/alignAssembly.config`): `MIN_PERCENT_ALIGNED=75`, `MIN_AVG_PER_ID=95`, `NUM_BP_PERFECT_SPLICE_BOUNDARY=0`, `subcluster_builder.dbi -m=50`.

Then `pip install -e .` in PyPASA (Python ≥3.9, deps in `pyproject.toml`) and run the baseline above.

## Open issues, in priority order

1. **Integrate the C++ assembler.** This is the main next step.
   - Binary: `PASApipeline/bin/pasa` (built from `pasa_cpp/`; note: *not* `pasa-plugins/cdbtools/` as ARCHITECTURE.md says).
   - Input format examples: `pasa_cpp/pasa_cpp_sample_input*`, one alignment per line: `acc,orient,lend-rend,lend-rend,...`. See `pasa_cpp/README.pasa`.
   - Perl driver to mirror: `PerlLib/CDNA/PASA_alignment_assembler.pm`. Lines 115–118 run per strand via `force_flexorient` (lines ~317–332: use `spliced_orient` if +/-, else the given strand). `pasa_cpp_assemblies` (line 205) writes the input, calls the binary, and parses the output.
   - Target: assemblies matching `sample_mydb_pasa.sqlite.pasa_assemblies.gff3` / `.assemblies.fasta` exactly.
2. **Subclustering** after assembly: port `scripts/subcluster_builder.dbi` (≥50% overlap, `spliced_orient` must match exactly, so `'?'` only pairs with `'?'`). Background: `docs/history/SINGLE_EXON_CLUSTERING_ANALYSIS.md`.
3. **Validation gap (12,465 vs 12,251).** VALIDATION_FIX.md got 12,303 (the "99.58% agreement" ARCHITECTURE.md still quotes). ALIGNED_LENGTH_FIX.md then switched `aligned_length` to transcript coordinates, which moved the count to 12,465, *further* from Perl. That switch does match Perl's formula: `PerlLib/CDNA/CDNA_alignment.pm:600-637` sums `abs(mrend-mlend)+1` (transcript coords) per segment, divides by `cdna_length`, and length-weights `avg_per_id` the same way. So the gap is somewhere else. Candidates (the transcript FASTA is already ruled out: Perl and the baseline both use `all_transcripts.fasta.clean`): how `per_id` is parsed per segment, and the other rejection rules in `validate_alignments_in_db.dbi` (around lines 330–390). Start by diffing the valid ID sets to get the 214 and looking at a few.
4. **Output writers** (GFF3/GTF/BED/FASTA) once assemblies exist.
5. **Unit tests**: `tests/` is empty, and pytest `addopts` in `pyproject.toml` assume pytest-cov.
6. **Performance.** The docs claim "<10 seconds"; the current code takes ~4 min on 14.5k alignments, mostly in validation. Not a priority until correctness is done.

## Stale or wrong claims in the other docs

- QUICKREF.md: `--stringent-clustering` and `--no-strand-specific` CLI flags **don't exist**. The CLI has `--min-overlap-percent` (unused by the clusterer). `respect_strand` is only reachable from Python. Its "563 vs 303 clusters" numbers are from an old state.
- STATUS.md, PROGRESS.md, COMPARISON_RESULTS.md, CLUSTERING_SUMMARY.md, CLUSTERING_ANALYSIS.md, UPDATED_COMPARISON.md, EXPANDING_WINDOW_RESULTS.md, SPLICE_ORIENTATION_RESULTS.md: a chronological lab notebook. Cluster counts in them (1,034 → 758 → 710 → 713 → 358 → 615) track a sequence of fixes, and several conclusions were later reversed (best-alignment selection, strand-specific clustering, gi|rev grouping). Read them for *why*, not for current behavior.
- `PASApipeline/sample_data/clusters_of_valid_alignments.txt.gz` is the only Perl output a fresh clone contains, and `scripts/compare_with_perl.py` reads it. It was committed in **2015** from an older run. The "303 clusters" in STATUS/QUICKREF/CLUSTERING_ROOT_CAUSE come from it. **303 is not a target**, and `compare_with_perl.py` is stale with it.
- `docs/history/CLUSTERING_ROOT_CAUSE.md` concludes strand-specific clustering is "correct". That was **reversed**; see rule 2.
- ARCHITECTURE.md: `exporter.py` "exists" (it doesn't), the C++ path is wrong (above), the "99.58%" figure is stale (above). Otherwise it's the best design reference.
- `docs/history/` holds four analysis docs copied from the parent `PASA_port/` directory, which isn't under version control.

## Doc map

- **Current design:** this file → ARCHITECTURE.md ("Critical Design Rules") → the module docstring in `clusterer.py`
- **History/rationale:** everything else, plus `docs/history/`
- **Key Perl references:** `scripts/validate_alignments_in_db.dbi`, `scripts/reassign_clusters_via_valid_align_coords.dbi`, `PerlLib/Overlap_piler.pm`, `PerlLib/CDNA/PASA_alignment_assembler.pm`, `scripts/subcluster_builder.dbi`, `Launch_PASA_pipeline.pl`
