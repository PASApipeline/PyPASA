"""Whole-file FASTA loading. Sequences are keyed by the first word of the header."""

from typing import Dict, List


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    acc = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                seqs[acc] = []
            elif acc is not None:
                seqs[acc].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_fasta_lengths(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    acc = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                lengths[acc] = 0
            elif acc is not None:
                lengths[acc] += len(line.strip())
    return lengths
