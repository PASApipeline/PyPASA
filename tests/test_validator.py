from pypasa.alignment.validator import AlignmentValidator
from pypasa.models import Alignment, AlignmentSegment, Orientation


def make(orient, segs):
    """segs: (genome_lend, genome_rend, transcript_lend, transcript_rend)"""
    return Alignment(
        align_id=1, cdna_acc="t", scaffold="chr", orient=Orientation(orient),
        segments=[AlignmentSegment(g1, g2, t1, t2, 100.0) for g1, g2, t1, t2 in segs],
    )


def genome_with_intron(donor, acceptor):
    # exon 1..10, intron 11..30, exon 31..40
    return "A" * 10 + donor + "N" * 16 + acceptor + "A" * 10


TWO_EXON_PLUS = [(1, 10, 1, 10), (31, 40, 11, 20)]


def test_splice_orientation_aligned_first_then_opposite():
    v = AlignmentValidator.__new__(AlignmentValidator)
    plus = make("+", TWO_EXON_PLUS)
    assert v._identify_splice_junctions(plus, genome_with_intron("GT", "AG")) == Orientation.PLUS
    assert v._identify_splice_junctions(plus, genome_with_intron("AT", "AC")) == Orientation.PLUS
    # reverse-strand consensus on a '+'-aligned transcript flips to '-'
    assert v._identify_splice_junctions(plus, genome_with_intron("CT", "AC")) == Orientation.MINUS
    assert v._identify_splice_junctions(plus, genome_with_intron("GT", "AT")) == Orientation.MINUS
    assert v._identify_splice_junctions(plus, genome_with_intron("GA", "AG")) is None


def test_splice_requires_every_intron_in_one_orientation():
    v = AlignmentValidator.__new__(AlignmentValidator)
    seq = "A" * 10 + "GT" + "N" * 16 + "AG" + "A" * 10 + "CT" + "N" * 16 + "AC" + "A" * 10
    aln = make("+", [(1, 10, 1, 10), (31, 40, 11, 20), (61, 70, 21, 30)])
    assert v._identify_splice_junctions(aln, seq) is None


def test_contiguity():
    assert AlignmentValidator._is_contiguous(make("+", TWO_EXON_PLUS))
    assert not AlignmentValidator._is_contiguous(make("+", [(1, 10, 1, 10), (31, 40, 13, 22)]))
    # '-' strand: genome-ordered segments run backwards on the transcript
    assert AlignmentValidator._is_contiguous(make("-", [(1, 10, 11, 20), (31, 40, 1, 10)]))
    assert not AlignmentValidator._is_contiguous(make("-", [(1, 10, 11, 20), (31, 40, 1, 9)]))
