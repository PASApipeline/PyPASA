"""
Assemble every cluster of valid alignments.

Corresponds to assemble_clusters.dbi plus assembly_db_loader.dbi in Perl PASA.
The per-cluster logic is in cluster_assembler.py and the assembly algorithm
itself is a pluggable backend (backend.py); this module groups alignments
into clusters and turns results into Assembly objects.

Assemblies are numbered asmbl_1, asmbl_2, ... in (scaffold, cluster_id,
reported) order, as assembly_db_loader.dbi does per scaffold file.
"""

import logging
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List

from pypasa.assembly.backend import AssemblerBackend
from pypasa.assembly.cluster_assembler import (
    AssembledTranscript,
    ClusterAlignment,
    assemble_cluster,
)
from pypasa.models import Alignment, AlignmentSegment, Assembly, Orientation

logger = logging.getLogger(__name__)

# assemble_clusters.dbi -n default. Perl samples this many per oversized
# cluster (select_aligns_to_assemble); that sampling isn't ported yet.
MAX_ALIGNS_PER_CLUSTER = 5000


def to_cluster_alignment(aln: Alignment) -> ClusterAlignment:
    spliced = aln.spliced_orient.value if aln.spliced_orient else "?"
    return ClusterAlignment(
        acc=aln.align_acc or aln.cdna_acc,
        align_id=aln.align_id,
        segments=tuple((s.genome_lend, s.genome_rend) for s in aln.segments),
        aligned_orient=Orientation(aln.orient).value,
        spliced_orient=spliced,
        is_fli=aln.is_fli,
    )


def remap_transcript_coords(
    segments: List[tuple], spliced_orient: str, aligned_orient: str
) -> List[AlignmentSegment]:
    """
    Assign transcript coords 1..len along the transcribed strand.

    Port of CDNA_alignment::remap_cdna_segment_coords (spliced orient wins,
    else aligned orient). Returned segments stay in genome order.
    """
    orient = spliced_orient if spliced_orient in ("+", "-") else aligned_orient
    ordered = list(reversed(segments)) if orient == "-" else list(segments)
    mcoords: Dict[tuple, tuple] = {}
    pos = 0
    for lend, rend in ordered:
        length = rend - lend + 1
        mcoords[(lend, rend)] = (pos + 1, pos + length)
        pos += length
    return [
        AlignmentSegment(
            genome_lend=lend,
            genome_rend=rend,
            transcript_lend=mcoords[(lend, rend)][0],
            transcript_rend=mcoords[(lend, rend)][1],
        )
        for lend, rend in segments
    ]


def group_clusters(alignments: List[Alignment]) -> Dict[str, Dict[int, List[Alignment]]]:
    """
    Valid, clustered alignments as scaffold -> cluster_id -> alignments.

    Scaffolds and clusters are in sorted order. Within a cluster, alignments
    are ordered by (lend, rend, align_id): the order Perl's
    assemble_clusters.dbi query returns via its (cluster_id, lend, rend) index,
    which the pasa binary is sensitive to.
    """
    groups: Dict[str, Dict[int, List[Alignment]]] = defaultdict(lambda: defaultdict(list))
    for aln in alignments:
        if aln.valid and aln.cluster_id is not None:
            groups[aln.scaffold][aln.cluster_id].append(aln)
    return {
        scaffold: {
            cid: sorted(groups[scaffold][cid], key=lambda a: (a.genome_lend, a.genome_rend, a.align_id))
            for cid in sorted(groups[scaffold])
        }
        for scaffold in sorted(groups)
    }


class PASAAssembler:
    """Runs cluster assembly with a given backend."""

    def __init__(self, backend: AssemblerBackend, workers: int = 1):
        """
        Args:
            backend: the assembler backend
            workers: clusters assembled concurrently. Threads suffice because
                the pasa_cpp backend's time is spent in subprocesses.
        """
        self.backend = backend
        self.workers = workers
        self.next_assembly_id = 1
        self.stats = {
            "clusters_processed": 0,
            "assemblies_generated": 0,
        }

    def assemble_all_clusters(self, alignments: List[Alignment]) -> List[Assembly]:
        clusters_by_scaffold = group_clusters(alignments)
        logger.info(
            f"Assembling clusters on {len(clusters_by_scaffold)} scaffolds with {self.backend.name}"
        )

        jobs = []
        for scaffold, clusters in clusters_by_scaffold.items():
            for cluster_id, members in clusters.items():
                if len(members) > MAX_ALIGNS_PER_CLUSTER:
                    logger.warning(
                        f"cluster {cluster_id} has {len(members)} alignments "
                        f"(> {MAX_ALIGNS_PER_CLUSTER}); Perl would subsample, assembling all"
                    )
                jobs.append((scaffold, cluster_id, [to_cluster_alignment(a) for a in members]))

        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            results = pool.map(lambda job: assemble_cluster(job[2], self.backend), jobs)
            # map() yields in job order, so assembly ids are the same for any worker count
            assemblies = []
            for (scaffold, cluster_id, _), transcripts in zip(jobs, results):
                assemblies.extend(self._to_assembly(t, cluster_id, scaffold) for t in transcripts)
                self.stats["clusters_processed"] += 1

        self.stats["assemblies_generated"] = len(assemblies)
        logger.info(
            f"Assembly complete: {len(assemblies)} assemblies "
            f"from {self.stats['clusters_processed']} clusters"
        )
        return assemblies

    def _to_assembly(self, t: AssembledTranscript, cluster_id: int, scaffold: str) -> Assembly:
        assembly = Assembly(
            assembly_id=self.next_assembly_id,
            cluster_id=cluster_id,
            scaffold=scaffold,
            lend=t.segments[0][0],
            rend=t.segments[-1][1],
            orient=Orientation(t.aligned_orient),
            spliced_orient=Orientation(t.spliced_orient),
            segments=remap_transcript_coords(t.segments, t.spliced_orient, t.aligned_orient),
            component_align_ids=t.align_ids,
            component_accs=t.accs,
            is_fli=t.is_fli,
        )
        self.next_assembly_id += 1
        return assembly

    def get_assembly_report(self) -> str:
        report = [
            "=" * 60,
            f"PASA ASSEMBLY REPORT ({self.backend.name})",
            "=" * 60,
            f"Clusters processed:      {self.stats['clusters_processed']:>10,}",
            f"Assemblies generated:    {self.stats['assemblies_generated']:>10,}",
        ]
        if self.stats["clusters_processed"] > 0:
            avg = self.stats["assemblies_generated"] / self.stats["clusters_processed"]
            report.append(f"Avg assemblies/cluster:  {avg:>10.2f}")
        report.append("=" * 60)
        return "\n".join(report)
