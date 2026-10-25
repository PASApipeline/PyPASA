from pypasa.alignment.polya import PolyATrim, evaluate_site, infer_site, parse_cln
from pypasa.models import Alignment, AlignmentSegment, Orientation


def test_parse_cln_keeps_single_ended_trims(tmp_path):
    cln = tmp_path / "x.cln"
    cln.write_text(
        "a\t0.00\t1\t100\t120\ttrimpoly[+0, -20]\t\n"
        "b\t0.00\t15\t120\t120\ttrimpoly[+14, -0]\t\n"
        "c\t0.00\t1\t100\t100\t\t\n"
        "d\t0.00\t5\t90\t100\ttrimpoly[+4, -0]; trimpoly[+0, -10]\t\n"
    )
    trims = parse_cln(str(cln))
    assert set(trims) == {"a", "b"}
    assert (trims["a"].type, trims["a"].trim, trims["a"].untrimmed_length) == ("rend", 20, 120)
    assert (trims["b"].type, trims["b"].trim) == ("lend", 14)


def single_exon(orient, glend, grend, tlend, trend):
    return Alignment(1, "t", "chr", [AlignmentSegment(glend, grend, tlend, trend, 100.0)],
                     Orientation(orient))


def test_infer_site_strands():
    trim = PolyATrim("t", "rend", 20, 120)  # trimmed length 100
    site = infer_site(single_exon("+", 1001, 1099, 1, 99), trim)
    assert (site.transcribed_orient, site.genomic_coord, site.align_end_offset) == ("+", 1100, 1)
    site = infer_site(single_exon("-", 1001, 1100, 1, 100), trim)
    assert (site.transcribed_orient, site.genomic_coord, site.align_end_offset) == ("-", 1001, 0)
    polyT = PolyATrim("t", "lend", 20, 120)
    site = infer_site(single_exon("+", 1001, 1100, 1, 100), polyT)
    assert (site.transcribed_orient, site.genomic_coord) == ("-", 1001)


def test_evaluate_site_rejects_internal_priming():
    trim = PolyATrim("t", "rend", 20, 120)
    cdna = "C" * 100 + "A" * 20
    site = infer_site(single_exon("+", 11, 110, 1, 100), trim)
    good_genome = "C" * 110 + "G" * 50
    assert evaluate_site(site, trim, good_genome, cdna) is None
    site = infer_site(single_exon("+", 11, 110, 1, 100), trim)
    primed_genome = "C" * 110 + "A" * 50
    assert evaluate_site(site, trim, primed_genome, cdna) == "potential internal genomic priming"
