"""
In-memory state of a PyPASA run, passed from stage to stage.

Stages read and update these objects directly; nothing goes through a
database. save()/load() write an optional gzipped JSON checkpoint, used to
resume or inspect a run and to compare against Perl PASA.
"""

import gzip
import json
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from typing import Dict, List

from pypasa.alignment.polya import PolyASite
from pypasa.assembly.subclusters import Subcluster
from pypasa.models import Alignment, Assembly, alignment_from_dict, assembly_from_dict

CHECKPOINT_VERSION = 1


@dataclass
class PipelineRun:
    alignments: List[Alignment] = field(default_factory=list)
    polya_sites: List[PolyASite] = field(default_factory=list)
    assemblies: List[Assembly] = field(default_factory=list)
    subclusters: List[Subcluster] = field(default_factory=list)
    completed_stages: List[str] = field(default_factory=list)

    def valid_alignments_by_scaffold(self) -> Dict[str, List[Alignment]]:
        """Valid alignments grouped by scaffold, scaffolds in sorted order, alignments in align_id order."""
        groups: Dict[str, List[Alignment]] = defaultdict(list)
        for aln in self.alignments:
            if aln.valid:
                groups[aln.scaffold].append(aln)
        return {s: groups[s] for s in sorted(groups)}

    def save(self, path: str) -> None:
        payload = {
            "version": CHECKPOINT_VERSION,
            "completed_stages": self.completed_stages,
            "alignments": [a.to_dict() for a in self.alignments],
            "polya_sites": [asdict(p) for p in self.polya_sites],
            "assemblies": [a.to_dict() for a in self.assemblies],
            "subclusters": [asdict(s) for s in self.subclusters],
        }
        with gzip.open(path, "wt") as fh:
            json.dump(payload, fh)

    @classmethod
    def load(cls, path: str) -> "PipelineRun":
        with gzip.open(path, "rt") as fh:
            payload = json.load(fh)
        if payload.get("version") != CHECKPOINT_VERSION:
            raise ValueError(f"{path}: checkpoint version {payload.get('version')}, expected {CHECKPOINT_VERSION}")
        return cls(
            alignments=[alignment_from_dict(d) for d in payload["alignments"]],
            polya_sites=[PolyASite(**d) for d in payload["polya_sites"]],
            assemblies=[assembly_from_dict(d) for d in payload["assemblies"]],
            subclusters=[Subcluster(**d) for d in payload.get("subclusters", [])],
            completed_stages=payload["completed_stages"],
        )
