"""
Group a cluster's assemblies into subclusters (candidate alternative isoforms of one gene).

Port of subcluster_builder.dbi (+ subcluster_loader.dbi, which just loads its
printed output). Within each cluster, two assemblies are linked when they:
- have the same spliced orientation ('?' pairs only with '?'),
- overlap (strictly: Perl's lend < other.rend and rend > other.lend), and
- share >= min_percent_overlap of either one's genomic span.
Subclusters are the single-linkage groups; unlinked assemblies are singletons.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List

from pypasa.models import Assembly

MIN_PERCENT_OVERLAP = 50.0  # subcluster_builder.dbi -m (sqlite.confs/alignAssembly.config)


@dataclass
class Subcluster:
    subcluster_id: int
    cluster_id: int
    assembly_ids: List[int]


def _linked(a: Assembly, b: Assembly, min_percent_overlap: float) -> bool:
    if a.spliced_orient != b.spliced_orient:
        return False
    if not (a.lend < b.rend and a.rend > b.lend):
        return False
    common = min(a.rend, b.rend) - max(a.lend, b.lend) + 1
    return (
        100.0 * common / a.span >= min_percent_overlap
        or 100.0 * common / b.span >= min_percent_overlap
    )


def build_subclusters(
    assemblies: List[Assembly], min_percent_overlap: float = MIN_PERCENT_OVERLAP
) -> List[Subcluster]:
    """
    Subcluster assemblies within each cluster.

    Returns:
        Subclusters ordered by cluster, then by their first assembly id;
        assembly ids within a subcluster are sorted.
    """
    by_cluster: Dict[int, List[Assembly]] = defaultdict(list)
    for asm in assemblies:
        by_cluster[asm.cluster_id].append(asm)

    subclusters = []
    for cluster_id in sorted(by_cluster):
        members = sorted(by_cluster[cluster_id], key=lambda a: a.assembly_id)
        parent = list(range(len(members)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                if _linked(members[i], members[j], min_percent_overlap):
                    parent[find(j)] = find(i)

        groups: Dict[int, List[int]] = defaultdict(list)
        for i, asm in enumerate(members):
            groups[find(i)].append(asm.assembly_id)
        for ids in sorted(groups.values()):
            subclusters.append(Subcluster(len(subclusters) + 1, cluster_id, sorted(ids)))

    return subclusters
