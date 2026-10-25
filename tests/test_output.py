import io

from pypasa.models import AlignmentSegment, Assembly, Orientation
from pypasa.output import alignment_text, spliced_sequence

# exon 1..10, intron 11..30 (GT...AG), exon 31..40
GENOME = "ACGTACGTAC" + "GT" + "N" * 16 + "AG" + "TTTTTGGGGG"


def two_exon(orient, spliced):
    tcoords = [(1, 10), (11, 20)] if orient == "+" else [(11, 20), (1, 10)]
    segs = [AlignmentSegment(g1, g2, t1, t2) for (g1, g2), (t1, t2) in zip([(1, 10), (31, 40)], tcoords)]
    return Assembly(1, 1, "chr", 1, 40, Orientation(orient), segs, spliced_orient=Orientation(spliced))


def test_alignment_text_marks_consensus_and_nonconsensus_sites():
    assert alignment_text(two_exon("+", "+"), GENOME) == \
        "orient(a+/s+) align: 1(1)-10(10)>GT....AG<31(11)-40(20)"
    bad = GENOME[:10] + "GA" + GENOME[12:]  # donor GA: not consensus on either strand
    assert alignment_text(two_exon("+", "?"), bad) == \
        "orient(a+/s?) align: 1(1)-10(10)>ga....AG<31(11)-40(20)"


def test_single_exon_keeps_stored_orientation():
    asm = Assembly(1, 1, "chr", 1, 10, Orientation.MINUS,
                   [AlignmentSegment(1, 10, 1, 10)], spliced_orient=Orientation.UNKNOWN)
    assert alignment_text(asm, GENOME) == "orient(a-/s?) align: 1(10)-10(1)"


def test_spliced_sequence_reverse_complements_minus():
    assert spliced_sequence(two_exon("+", "+"), GENOME) == "ACGTACGTACTTTTTGGGGG"
    assert spliced_sequence(two_exon("-", "-"), GENOME) == "CCCCCAAAAAGTACGTACGT"
