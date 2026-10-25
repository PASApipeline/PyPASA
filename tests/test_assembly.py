from pathlib import Path

import pytest

from pypasa.assembly.backend import AssemblerBackend, RawAssembly
from pypasa.assembly.cluster_assembler import (
    AssemblyError,
    ClusterAlignment,
    assemble_cluster,
    force_flexorient,
)
from pypasa.assembly.assembler import remap_transcript_coords
from pypasa.assembly.pasa_cpp import PasaCppBackend, format_pasa_input, parse_pasa_output

PASA_HOME = Path(__file__).resolve().parents[2] / "PASApipeline"
PASA_BIN = PASA_HOME / "bin" / "pasa"


def aln(acc, segs, aligned="+", spliced="?", fli=False, align_id=None):
    return ClusterAlignment(acc, align_id or hash(acc) % 10000, tuple(segs), aligned, spliced, fli)


class ScriptedBackend(AssemblerBackend):
    """Returns canned assemblies per strand pass and records what it was given."""

    name = "scripted"

    def __init__(self, plus, minus):
        self.by_pass = [plus, minus]
        self.calls = []

    def assemble(self, inputs):
        self.calls.append(list(inputs))
        return self.by_pass[len(self.calls) - 1]


def test_force_flexorient_keeps_spliced_and_forces_unspliced():
    alns = [aln("a", [(1, 10), (20, 30)], spliced="-"), aln("b", [(5, 15)])]
    assert [i.orient for i in force_flexorient(alns, "+")] == ["-", "+"]
    assert [i.orient for i in force_flexorient(alns, "-")] == ["-", "-"]


def test_greedy_selection_skips_fully_covered_assemblies():
    alns = [aln("a", [(1, 100)]), aln("b", [(50, 150)]), aln("c", [(120, 200)])]
    backend = ScriptedBackend(
        plus=[RawAssembly([(1, 150)], ["a", "b"]), RawAssembly([(50, 200)], ["b", "c"])],
        minus=[RawAssembly([(1, 100)], ["a"])],
    )
    out = assemble_cluster(alns, backend)
    # both 2-member assemblies rank above the singleton; ties keep reverse order
    assert [t.accs for t in out] == [["b", "c"], ["a", "b"]]
    assert len(backend.calls) == 2


def test_orientation_from_spliced_members_and_majority_vote():
    alns = [
        aln("a", [(1, 10), (20, 30)], aligned="-", spliced="+"),
        aln("b", [(25, 40)], aligned="-"),
        aln("c", [(100, 120)], aligned="-"),
        aln("d", [(110, 130)], aligned="+"),
    ]
    backend = ScriptedBackend(
        plus=[RawAssembly([(1, 10), (20, 40)], ["a", "b"]), RawAssembly([(100, 130)], ["c", "d"])],
        minus=[],
    )
    out = {tuple(t.accs): t for t in assemble_cluster(alns, backend)}
    assert (out[("a", "b")].spliced_orient, out[("a", "b")].aligned_orient) == ("+", "+")
    # unspliced 1-1 tie goes to '+'
    assert (out[("c", "d")].spliced_orient, out[("c", "d")].aligned_orient) == ("?", "+")


def test_conflicting_spliced_orient_and_missing_coverage_raise():
    alns = [aln("a", [(1, 10), (20, 30)], spliced="+"), aln("b", [(1, 10), (20, 30)], spliced="-")]
    with pytest.raises(AssemblyError):
        assemble_cluster(alns, ScriptedBackend([RawAssembly([(1, 30)], ["a", "b"])], []))
    with pytest.raises(AssemblyError):
        assemble_cluster(alns, ScriptedBackend([RawAssembly([(1, 30)], ["a"])], []))


def test_remap_transcript_coords_minus_strand():
    segs = remap_transcript_coords([(1, 10), (21, 25)], "-", "+")
    assert [(s.transcript_lend, s.transcript_rend) for s in segs] == [(6, 15), (1, 5)]


def test_format_parse_roundtrip():
    alns = [aln("x", [(1, 10), (20, 30)], spliced="+")]
    assert format_pasa_input(force_flexorient(alns, "+")) == "x,+,1-10,20-30\n"
    line = "assembly: (0) contains alignments: [x,y] with structure [assembly_0,-,30-20,1-10] score: (2)"
    (raw,) = parse_pasa_output(line)
    assert raw.accs == ["x", "y"] and raw.segments == [(20, 30), (1, 10)]


@pytest.mark.skipif(not PASA_BIN.exists(), reason="PASApipeline/bin/pasa not built")
def test_pasa_cpp_sample_input():
    sample = PASA_HOME / "pasa_cpp" / "pasa_cpp_sample_input"
    alns = []
    for i, line in enumerate(sample.read_text().split()):
        acc, orient, *coords = line.split(",")
        segs = [tuple(int(x) for x in c.split("-")) for c in coords]
        alns.append(aln(acc, segs, aligned=orient, spliced=orient, align_id=i))
    out = assemble_cluster(alns, PasaCppBackend(str(PASA_BIN)))
    # pasa reports 3 assemblies for this cluster (see `pasa pasa_cpp_sample_input`)
    assert len(out) == 3
    assert {acc for t in out for acc in t.accs} == {a.acc for a in alns}
    assert all(t.spliced_orient == "-" for t in out)
