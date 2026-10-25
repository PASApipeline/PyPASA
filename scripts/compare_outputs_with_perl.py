#!/usr/bin/env python3
"""
Compare PyPASA's assembly output files with Perl PASA's, file by file.

Usage: compare_outputs_with_perl.py PYPASA_PREFIX PERL_PREFIX

Prefixes are as in <prefix>.pasa_assemblies.gff3 (pypasa.output.write_all).
GFF3 and FASTA must be byte-identical. BED is compared as a set of lines
(Perl writes it in hash order). GTF and described.txt are compared after
mapping PyPASA subcluster ids onto Perl's through the assemblies they contain,
since subcluster numbering differs by design (see pypasa/output.py).
Exits 1 on any difference.
"""

import re
import sys


def read(path):
    with open(path) as fh:
        return fh.read().splitlines(keepends=True)


def subcluster_by_asmbl(described_lines):
    return {f[2]: f[1] for f in (l.split("\t") for l in described_lines[1:])}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    py, perl = sys.argv[1], sys.argv[2]
    ok = True

    def report(name, same, detail=""):
        nonlocal ok
        ok &= same
        print(f"{name}: {'identical' if same else 'DIFFERS'}{detail}")

    for suffix in ("pasa_assemblies.gff3", "assemblies.fasta"):
        report(suffix, read(f"{py}.{suffix}") == read(f"{perl}.{suffix}"))

    report("pasa_assemblies.bed",
           sorted(read(f"{py}.pasa_assemblies.bed")) == sorted(read(f"{perl}.pasa_assemblies.bed")),
           " (as sorted lines)")

    py_desc = read(f"{py}.pasa_assemblies_described.txt")
    perl_desc = read(f"{perl}.pasa_assemblies_described.txt")
    py_sub, perl_sub = subcluster_by_asmbl(py_desc), subcluster_by_asmbl(perl_desc)
    to_perl = {}
    for asmbl, sid in py_sub.items():
        to_perl.setdefault(sid, set()).add(perl_sub.get(asmbl))
    consistent = all(len(v) == 1 for v in to_perl.values()) and len(to_perl) == len(set(perl_sub.values()))
    report("subcluster partition", consistent)
    to_perl = {k: next(iter(v)) for k, v in to_perl.items()}

    mapped_desc = [py_desc[0]] + [
        "\t".join([f[0], to_perl[f[1]]] + f[2:]) for f in (l.split("\t") for l in py_desc[1:])
    ]
    report("pasa_assemblies_described.txt", mapped_desc == perl_desc, " (subcluster ids mapped)")

    mapped_gtf = [
        re.sub(r'gene_id "PASA_cluster_(\d+)"', lambda m: f'gene_id "PASA_cluster_{to_perl[m.group(1)]}"', l)
        for l in read(f"{py}.pasa_assemblies.gtf")
    ]
    report("pasa_assemblies.gtf", mapped_gtf == read(f"{perl}.pasa_assemblies.gtf"), " (subcluster ids mapped)")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
