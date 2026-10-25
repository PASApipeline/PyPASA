#!/usr/bin/env python3
"""
Compare a PyPASA run checkpoint against a Perl PASA SQLite database, stage by stage.

Usage: compare_with_perl_db.py PYPASA_CHECKPOINT PERL_DB [--prog minimap2] [--examples N]

PYPASA_CHECKPOINT is written by e.g. scripts/baseline_minimap2.py.

Alignments are matched by alignment accession (GFF3 ID = Perl align_link.align_acc),
restricted to Perl rows with prog=PROG. Reports:
  1. validation: valid/invalid agreement
  2. spliced orientation: agreement among alignments both sides call valid
  3. clusters: identical cluster member sets
  4. assemblies: identical (contig, exon list) structures, then orientation and
     membership for structures found on both sides

Use a Perl DB built from the same alignments only (e.g. perl_ref_mm2only), since
Perl clusters and assembles every imported alignment set together.
"""

import argparse
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pypasa.run import PipelineRun


def ro(path):
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True)


def load_perl(path, prog):
    c = ro(path)
    aligns = {}
    for acc, valid, so, cid, contig in c.execute(
        """SELECT al.align_acc, al.validate, al.spliced_orient, al.cluster_id, cl.annotdb_asmbl_id
           FROM align_link al LEFT JOIN clusters cl ON cl.cluster_id = al.cluster_id
           WHERE al.prog = ?""",
        (prog,),
    ):
        aligns[acc] = dict(valid=bool(valid), spliced=so or "?", cluster=cid, contig=contig)

    segs = defaultdict(list)
    for align_id, lend, rend in c.execute(
        """SELECT a.align_id, a.lend, a.rend FROM alignment a
           JOIN align_link al ON al.align_id = a.align_id WHERE al.prog = 'assembler'"""
    ):
        segs[align_id].append((min(lend, rend), max(lend, rend)))
    members = defaultdict(set)
    for asmbl, acc in c.execute("SELECT asmbl_acc, cdna_acc FROM asmbl_link"):
        members[asmbl].add(acc)
    asmbls = []
    for align_id, acc, aorient, sorient, contig in c.execute(
        """SELECT al.align_id, al.align_acc, al.aligned_orient, al.spliced_orient, cl.annotdb_asmbl_id
           FROM align_link al JOIN clusters cl ON cl.cluster_id = al.cluster_id
           WHERE al.prog = 'assembler'"""
    ):
        asmbls.append(
            dict(name=acc, contig=contig, segs=tuple(sorted(segs[align_id])),
                 orient=aorient, spliced=sorient or "?", members=members[acc])
        )
    polya = set(c.execute(
        """SELECT al.align_acc, p.genomicCoord, p.transcribedOrient, p.trimType, p.alignEndOffset
           FROM transcriptPolyA p JOIN align_link al ON al.align_id = p.align_id"""
    ))
    subclusters = defaultdict(set)
    for sid, acc in c.execute("SELECT subcluster_id, cdna_acc FROM subcluster_link"):
        subclusters[sid].add(acc)
    return aligns, asmbls, polya, list(subclusters.values())


def load_pypasa(path):
    run = PipelineRun.load(path)
    aligns, by_id = {}, {}
    for a in run.alignments:
        aligns[a.align_acc] = dict(
            valid=a.valid, spliced=a.spliced_orient.value if a.spliced_orient else "?",
            cluster=a.cluster_id, contig=a.scaffold,
        )
        by_id[a.align_id] = a.align_acc
    asmbls = [
        dict(name=f"asmbl_{x.assembly_id}", contig=x.scaffold,
             segs=tuple(sorted((s.genome_lend, s.genome_rend) for s in x.segments)),
             orient=x.orient.value, spliced=x.spliced_orient.value,
             members=set(x.component_accs or [by_id[i] for i in x.component_align_ids]))
        for x in run.assemblies
    ]
    polya = {
        (by_id[p.align_id], p.genomic_coord, p.transcribed_orient, p.trim_type, p.align_end_offset)
        for p in run.polya_sites
    }
    subclusters = [{f"asmbl_{i}" for i in s.assembly_ids} for s in run.subclusters]
    return aligns, asmbls, polya, subclusters


def section(title):
    print(f"\n== {title} ==")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pypasa_checkpoint")
    ap.add_argument("perl_db")
    ap.add_argument("--prog", default="minimap2")
    ap.add_argument("--examples", type=int, default=5)
    args = ap.parse_args()
    n_ex = args.examples

    perl_al, perl_asm, perl_polya, perl_sub = load_perl(args.perl_db, args.prog)
    py_al, py_asm, py_polya, py_sub = load_pypasa(args.pypasa_checkpoint)
    if None in py_al:
        raise SystemExit("PyPASA run has no align_acc values; re-import with the current importer")

    section("alignments")
    shared = perl_al.keys() & py_al.keys()
    print(f"perl {len(perl_al)}  pypasa {len(py_al)}  shared {len(shared)}")

    section("validation")
    conf = Counter((perl_al[a]["valid"], py_al[a]["valid"]) for a in shared)
    print(f"both valid {conf[True, True]}  both invalid {conf[False, False]}  "
          f"perl-only valid {conf[True, False]}  pypasa-only valid {conf[False, True]}")
    for label, want in (("perl-only valid", (True, False)), ("pypasa-only valid", (False, True))):
        ex = sorted(a for a in shared if (perl_al[a]["valid"], py_al[a]["valid"]) == want)[:n_ex]
        if ex:
            print(f"  e.g. {label}: {ex}")

    section("spliced orientation (both valid)")
    both = [a for a in shared if perl_al[a]["valid"] and py_al[a]["valid"]]
    so = Counter((perl_al[a]["spliced"], py_al[a]["spliced"]) for a in both)
    agree = sum(n for (p, q), n in so.items() if p == q)
    print(f"agree {agree}/{len(both)}")
    for (p, q), n in sorted(so.items(), key=lambda kv: -kv[1]):
        if p != q:
            ex = sorted(a for a in both if (perl_al[a]["spliced"], py_al[a]["spliced"]) == (p, q))[:n_ex]
            print(f"  perl {p} / pypasa {q}: {n}  e.g. {ex}")

    section("polyA sites")
    print(f"perl {len(perl_polya)}  pypasa {len(py_polya)}  identical {len(perl_polya & py_polya)}")

    section("clusters (valid alignments)")
    def partition(al):
        groups = defaultdict(set)
        for a, d in al.items():
            if d["valid"] and d["cluster"] is not None:
                groups[d["cluster"]].add(a)
        return {frozenset(g) for g in groups.values()}
    pc, qc = partition({a: perl_al[a] for a in shared}), partition({a: py_al[a] for a in shared})
    print(f"perl {len(partition(perl_al))}  pypasa {len(partition(py_al))}  "
          f"(over shared alignments: perl {len(pc)}, pypasa {len(qc)}, identical {len(pc & qc)})")

    section("assemblies")
    def key(x):
        return (x["contig"], x["segs"])
    pk, qk = Counter(map(key, perl_asm)), Counter(map(key, py_asm))
    common = pk & qk
    print(f"perl {len(perl_asm)}  pypasa {len(py_asm)}  identical structures {sum(common.values())}  "
          f"perl-only {sum((pk - qk).values())}  pypasa-only {sum((qk - pk).values())}")
    p_by, q_by = defaultdict(list), defaultdict(list)
    for x in perl_asm:
        p_by[key(x)].append(x)
    for x in py_asm:
        q_by[key(x)].append(x)
    orient_mismatch, member_mismatch = [], []
    for k in common:
        if len(p_by[k]) != 1 or len(q_by[k]) != 1:
            continue
        p, q = p_by[k][0], q_by[k][0]
        if (p["orient"], p["spliced"]) != (q["orient"], q["spliced"]):
            orient_mismatch.append((p["name"], q["name"], p["orient"], p["spliced"], q["orient"], q["spliced"]))
        if p["members"] != q["members"]:
            member_mismatch.append((p["name"], q["name"], len(p["members"]), len(q["members"])))
    print(f"among identical structures: orientation differs {len(orient_mismatch)}, "
          f"membership differs {len(member_mismatch)}")
    for row in orient_mismatch[:n_ex]:
        print(f"  orient perl {row[0]} ({row[2]}/{row[3]})  pypasa {row[1]} ({row[4]}/{row[5]})")
    for row in member_mismatch[:n_ex]:
        print(f"  members perl {row[0]} n={row[2]}  pypasa {row[1]} n={row[3]}")
    for label, diff, by in (("perl-only", pk - qk, p_by), ("pypasa-only", qk - pk, q_by)):
        for k in list(diff)[:n_ex]:
            x = by[k][0]
            print(f"  {label} {x['name']} {x['contig']} {x['orient']} "
                  f"{','.join(f'{a}-{b}' for a, b in x['segs'])}")


    section("subclusters")
    # compare as sets of assembly structures; assembly names differ between runs
    p_key = {x["name"]: key(x) for x in perl_asm}
    q_key = {x["name"]: key(x) for x in py_asm}
    ps = Counter(frozenset(p_key[a] for a in sub) for sub in perl_sub)
    qs = Counter(frozenset(q_key[a] for a in sub) for sub in py_sub)
    print(f"perl {sum(ps.values())}  pypasa {sum(qs.values())}  identical {sum((ps & qs).values())}")


if __name__ == "__main__":
    main()
