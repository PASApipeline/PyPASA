"""
Assemble one cluster of valid alignments, independent of the assembler backend.

Port of assemble_alignments() and force_flexorient() in
PerlLib/CDNA/PASA_alignment_assembler.pm:

1. Sort alignments by lend (stable, so input order breaks ties).
2. Run the backend twice: once with unoriented ('?' spliced_orient) alignments
   forced to '+', once forced to '-'. Alignments with a spliced orientation
   keep it in both passes.
3. Order all assemblies ('+' pass first, then '-') by number of contained
   alignments, descending; ties keep reverse input order, matching Perl's
   `reverse sort {...}` over a stable sort.
4. Greedily report each assembly that contains at least one accession not yet
   covered, stopping once every accession is covered.
5. Set each reported assembly's spliced orientation from its members, its
   aligned orientation from the spliced one (or a majority vote of members'
   aligned orientations when unspliced), and its FL status.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

from pypasa.assembly.backend import AssemblerBackend, AssemblyInput, RawAssembly, Segment


@dataclass
class ClusterAlignment:
    """The fields of a valid alignment that assembly needs."""
    acc: str
    align_id: int
    segments: Tuple[Segment, ...]  # sorted by lend
    aligned_orient: str  # '+' or '-'
    spliced_orient: str  # '+', '-', or '?'
    is_fli: bool = False

    @property
    def lend(self) -> int:
        return self.segments[0][0]


@dataclass
class AssembledTranscript:
    """A reported assembly with PASA's orientation and FL calls applied."""
    segments: List[Segment]
    accs: List[str]
    spliced_orient: str
    aligned_orient: str
    is_fli: bool
    align_ids: List[int] = field(default_factory=list)


class AssemblyError(RuntimeError):
    pass


def force_flexorient(alignments: Sequence[ClusterAlignment], orient: str) -> List[AssemblyInput]:
    """Fix each alignment's orientation for one strand pass."""
    return [
        AssemblyInput(
            acc=a.acc,
            orient=a.spliced_orient if a.spliced_orient in ("+", "-") else orient,
            segments=a.segments,
        )
        for a in alignments
    ]


def assemble_cluster(
    alignments: Sequence[ClusterAlignment], backend: AssemblerBackend
) -> List[AssembledTranscript]:
    """Assemble one cluster; returns the reported (non-redundant) assemblies."""
    if not alignments:
        return []

    alignments = sorted(alignments, key=lambda a: a.lend)

    spliced = {a.acc: a.spliced_orient for a in alignments}
    aligned = {a.acc: a.aligned_orient for a in alignments}
    fli = {a.acc: a.is_fli for a in alignments}
    align_ids: Dict[str, List[int]] = {}
    for a in alignments:
        align_ids.setdefault(a.acc, []).append(a.align_id)

    raw: List[RawAssembly] = []
    for orient in ("+", "-"):
        raw.extend(backend.assemble(force_flexorient(alignments, orient)))

    ranked = list(reversed(sorted(raw, key=lambda r: len(r.accs))))

    covered = {acc: False for acc in spliced}
    reported = []
    for asm in ranked:
        have_unseen = False
        spliced_orient = "?"
        is_fli = False
        orient_counts: Counter = Counter()
        for acc in asm.accs:
            if not covered[acc]:
                have_unseen = True
            covered[acc] = True
            curr = spliced[acc]
            if curr != "?":
                if spliced_orient != "?" and spliced_orient != curr:
                    raise AssemblyError(
                        f"conflicting spliced orientations in assembly "
                        f"(have {spliced_orient}, {acc} has {curr})"
                    )
                spliced_orient = curr
            is_fli = is_fli or fli[acc]
            orient_counts[aligned[acc]] += 1

        if spliced_orient == "?":
            # Perl breaks count ties by hash-key order, i.e. arbitrarily;
            # here a tie goes to '+'.
            aligned_orient = max(("+", "-"), key=lambda o: orient_counts[o])
        else:
            aligned_orient = spliced_orient

        if have_unseen:
            reported.append(
                AssembledTranscript(
                    segments=list(asm.segments),
                    accs=list(asm.accs),
                    spliced_orient=spliced_orient,
                    aligned_orient=aligned_orient,
                    is_fli=is_fli,
                    align_ids=[i for acc in asm.accs for i in align_ids[acc]],
                )
            )
        if all(covered.values()):
            break
    else:
        missing = [acc for acc, seen in covered.items() if not seen]
        if missing:
            raise AssemblyError(f"assemblies don't cover all alignments; missing {missing[:5]}")

    return reported
