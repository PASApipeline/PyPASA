from pypasa.assembly.subclusters import build_subclusters
from pypasa.models import AlignmentSegment, Assembly, Orientation


def asm(aid, lend, rend, spliced="+", cluster=1):
    return Assembly(aid, cluster, "chr", lend, rend, Orientation.PLUS,
                    [AlignmentSegment(lend, rend, 1, rend - lend + 1)],
                    spliced_orient=Orientation(spliced))


def groups(assemblies):
    return [s.assembly_ids for s in build_subclusters(assemblies)]


def test_overlap_threshold_uses_either_span():
    # 1-100 vs 51-1000: 50 bp in common is 50% of the first, so linked
    assert groups([asm(1, 1, 100), asm(2, 51, 1000)]) == [[1, 2]]
    # 1-100 vs 52-1000: 49% / 5%, not linked
    assert groups([asm(1, 1, 100), asm(2, 52, 1000)]) == [[1], [2]]


def test_single_linkage_orientation_and_cluster_boundaries():
    a = [asm(1, 1, 100), asm(2, 40, 200), asm(3, 120, 260),  # chain 1-2-3; 1 and 3 do not overlap
         asm(4, 1, 100, spliced="?"),                         # same span, '?' vs '+'
         asm(5, 1, 100, cluster=2)]                           # other cluster
    assert groups(a) == [[1, 2, 3], [4], [5]]


def test_touching_spans_are_not_linked():
    # Perl requires lend < other.rend strictly
    assert groups([asm(1, 1, 100), asm(2, 100, 101)]) == [[1], [2]]
