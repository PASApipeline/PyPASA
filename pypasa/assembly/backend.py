"""
Pluggable assembler backends.

A backend does one thing: given alignments that all carry a fixed orientation
('+' or '-'), return the maximal compatible assemblies and which input
accessions each one contains. Everything PASA-specific around that (running
both strand passes, choosing the reported subset, assigning orientation and
FL status) lives in cluster_assembler.py, so backends stay interchangeable.

Implementations:
- pasa_cpp.PasaCppBackend: the C++ `pasa` binary, as used by Perl PASA.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Sequence, Tuple

Segment = Tuple[int, int]  # (lend, rend), 1-based inclusive genome coords, lend <= rend


@dataclass(frozen=True)
class AssemblyInput:
    """One alignment as handed to a backend for a single strand pass."""
    acc: str
    orient: str  # '+' or '-' (fixed for this pass; see force_flexorient)
    segments: Tuple[Segment, ...]  # sorted by lend


@dataclass
class RawAssembly:
    """One backend-produced assembly, before PASA's selection step."""
    segments: List[Segment]  # sorted by lend
    accs: List[str]  # contained input accessions, in backend output order


class AssemblerBackend(ABC):
    """Assembles same-orientation alignments into maximal compatible assemblies."""

    name: str = "abstract"

    @abstractmethod
    def assemble(self, inputs: Sequence[AssemblyInput]) -> List[RawAssembly]:
        """
        Assemble one strand pass of one cluster.

        Args:
            inputs: alignments in layout order (sorted by lend)

        Returns:
            All assemblies the backend produces. Every input accession must be
            contained in at least one assembly.
        """

    def assemble_many(self, batches: Sequence[Sequence[AssemblyInput]]) -> List[List[RawAssembly]]:
        """
        Assemble many independent strand passes; results are in input order.

        Backends that can do this more cheaply than one assemble() per pass
        (e.g. one process for all of them) override it.
        """
        return [self.assemble(batch) for batch in batches]
