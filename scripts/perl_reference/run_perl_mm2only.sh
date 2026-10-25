#!/bin/bash
# Build the Perl PASA ground truth that PyPASA is compared against: the
# PASApipeline sample data, minimap2 alignments ONLY (no custom_alignments.gff3,
# which the stock runMe.SQLite.sh also loads and assembles together with them),
# alignment assembly only. Takes ~25 min on a loaded host.
#
# Usage: run_perl_mm2only.sh PASA_PORT_DIR OUT_DIR
#   PASA_PORT_DIR holds PASApipeline/ (built with `make`) and PyPASA/ side by side.
#   OUT_DIR gets mm2only_pasa.sqlite plus mm2only_pasa.sqlite.* outputs.
#
# Needs on PATH (or set PERL_BIN_DIR / MINIMAP2_DIR): a perl with DBI, DBD::SQLite,
# URI::Escape and threads; minimap2; samtools. On the Broad hosts the system perl
# lacks threads.pm and ~/perl5's DBD::SQLite is built for another perl, so we used
# a conda perl (micromamba create -c conda-forge -c bioconda perl perl-dbi
# perl-dbd-sqlite perl-uri) with PERL5LIB cleared, and minimap2 2.28.
set -euo pipefail
PORT=$(cd "$1" && pwd); OUT=$2
mkdir -p "$OUT" && OUT=$(cd "$OUT" && pwd)
SD=$PORT/PASApipeline/sample_data
export PERL5LIB= PERL_LOCAL_LIB_ROOT=
export PATH=${PERL_BIN_DIR:+$PERL_BIN_DIR:}${MINIMAP2_DIR:+$MINIMAP2_DIR:}$PATH

cd "$OUT"
for f in all_transcripts.fasta all_transcripts.fasta.clean all_transcripts.fasta.cln FL_accs.txt tdn.accs; do
    ln -sf "$SD/$f" .
done
[ -s genome_sample.fasta ] || gunzip -c "$SD/genome_sample.fasta.gz" > genome_sample.fasta
sed "s#^DATABASE=.*#DATABASE=$OUT/mm2only_pasa.sqlite#" "$SD/sqlite.confs/alignAssembly.config" > alignAssembly.config

"$PORT/PASApipeline/Launch_PASA_pipeline.pl" -c alignAssembly.config -C -r -R \
    -g genome_sample.fasta -t all_transcripts.fasta.clean -T -u all_transcripts.fasta \
    -f FL_accs.txt --TDN tdn.accs --ALIGNERS minimap2 --CPU 4 -N 1
