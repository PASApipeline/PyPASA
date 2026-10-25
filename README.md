# PyPASA

A Python port of the [PASA](https://github.com/PASApipeline/PASApipeline)
alignment-assembly pipeline. It reproduces Perl PASA's results exactly on its
sample data, from validated alignments through assemblies, subclusters and
output files, and runs in memory with no database.

See [HANDOFF.md](HANDOFF.md) for goals, status, design rules, how to check
against Perl, and what's left.

## Install

```bash
pip install -e '.[dev]'     # Python >= 3.9
```

Assembly uses PASA's C++ `pasa` binary: run `make` in a PASApipeline checkout and
point `--pasa-bin` at `PASApipeline/bin/pasa`.

## Usage

```bash
pypasa align-assemble \
    -g genome.fasta \
    -t transcripts.fasta.clean \        # transcripts as aligned (seqclean output)
    -a alignments.gff3 \                # repeatable
    -u transcripts.fasta --cln transcripts.fasta.cln \   # optional: polyA-site analysis
    -f FL_accs.txt \                    # optional: full-length cDNAs
    --pasa-bin PASApipeline/bin/pasa --cpu 4 \
    -o run.json.gz \                    # optional: save the run
    --out-prefix out/mydb               # optional: GFF3/GTF/BED/FASTA/described outputs

pypasa stats run.json.gz
```

## Layout

```
pypasa/
  pipeline.py        stage order and config
  run.py             in-memory run state, checkpoint save/load
  models.py          Alignment, Assembly
  alignment/         import, validation, polyA sites, clustering
  assembly/          backend interface, C++ pasa backend, per-cluster assembly, subclusters
  output.py          output writers
scripts/             sample-data baseline, Perl comparisons, Perl reference run
tests/
docs/history/        earlier development notes
```
