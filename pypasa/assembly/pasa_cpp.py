"""
Backend wrapping the C++ `pasa` assembler (PASApipeline/pasa_cpp).

Mirrors pasa_cpp_assemblies() in PerlLib/CDNA/PASA_alignment_assembler.pm:
write one alignment per line as `acc,orient,lend-rend,...`, run the binary,
and parse its `assembly: (N) contains alignments: [...] with structure [...]`
lines.

pasa builds with batch mode (-B) take many strand passes in one run, with
records separated by `//` lines; assemble_many() uses it when the binary
supports it, which avoids starting a process per pass.
"""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence

from pypasa.assembly.backend import AssemblerBackend, AssemblyInput, RawAssembly

_ASSEMBLY_LINE = re.compile(
    r"assembly:\s\(\d+\)\scontains\salignments:\s\[([^\]]+)\]\swith\sstructure\s\[([^\]]+)\]"
)


def format_pasa_input(inputs: Sequence[AssemblyInput]) -> str:
    """Render inputs in the pasa binary's text format."""
    lines = []
    for inp in inputs:
        if "," in inp.acc:
            raise ValueError(f"accession {inp.acc!r} contains a comma, which pasa input forbids")
        coords = ",".join(f"{lend}-{rend}" for lend, rend in inp.segments)
        lines.append(f"{inp.acc},{inp.orient},{coords}")
    return "\n".join(lines) + "\n"


def parse_pasa_output(text: str) -> List[RawAssembly]:
    """Parse assembly lines from pasa binary stdout."""
    assemblies = []
    for line in text.splitlines():
        m = _ASSEMBLY_LINE.search(line)
        if not m:
            continue
        acclist, structure = m.groups()
        fields = structure.split(",")
        # fields: assembly_name, orient, lend-rend, ...
        segments = []
        for coordset in fields[2:]:
            lend, rend = sorted(int(x) for x in coordset.split("-"))
            segments.append((lend, rend))
        assemblies.append(RawAssembly(segments=segments, accs=acclist.split(",")))
    return assemblies


class PasaCppBackend(AssemblerBackend):
    """Runs the C++ pasa binary once per strand pass."""

    name = "pasa_cpp"

    def __init__(self, pasa_bin: Optional[str] = None, tmpdir: Optional[str] = None):
        """
        Args:
            pasa_bin: path to the pasa binary; defaults to `pasa` on PATH
            tmpdir: where per-call input files go; defaults to $TMPDIR
        """
        resolved = pasa_bin or shutil.which("pasa")
        if not resolved or not os.access(resolved, os.X_OK):
            raise FileNotFoundError(
                f"pasa binary not found or not executable: {pasa_bin or 'pasa on PATH'}. "
                "Build it with 'make' in PASApipeline (it lands in PASApipeline/bin/pasa)."
            )
        self.pasa_bin = str(resolved)
        self.tmpdir = tmpdir
        self.supports_batch = self._supports_batch()

    def _supports_batch(self) -> bool:
        # An older pasa given '//' lines would loop forever, so check its usage text.
        usage = subprocess.run([self.pasa_bin], capture_output=True, text=True).stderr
        return "-B batch mode" in usage

    def _run(self, text: str, extra_args: List[str]) -> str:
        with tempfile.NamedTemporaryFile(
            mode="w", prefix="pasa.", suffix=".in", dir=self.tmpdir, delete=False
        ) as fh:
            fh.write(text)
            input_path = fh.name
        try:
            result = subprocess.run(
                [self.pasa_bin, input_path] + extra_args, capture_output=True, text=True
            )
            if result.returncode != 0:
                killer = Path(input_path).with_name("pasa_killer.input")
                shutil.copy(input_path, killer)
                raise RuntimeError(
                    f"pasa exited {result.returncode}; input saved to {killer}\n{result.stderr}"
                )
            return result.stdout
        finally:
            os.unlink(input_path)

    def assemble(self, inputs: Sequence[AssemblyInput]) -> List[RawAssembly]:
        return parse_pasa_output(self._run(format_pasa_input(inputs), []))

    def assemble_many(self, batches: Sequence[Sequence[AssemblyInput]]) -> List[List[RawAssembly]]:
        if not self.supports_batch:
            return super().assemble_many(batches)
        if any(len(b) == 0 for b in batches):
            raise ValueError("empty strand pass")
        text = "//\n".join(format_pasa_input(b) for b in batches)
        records = re.split(r"^//\n", self._run(text, ["-B"]), flags=re.M)[1:]
        if len(records) != len(batches):
            raise RuntimeError(f"pasa -B returned {len(records)} records for {len(batches)} inputs")
        return [parse_pasa_output(r) for r in records]
