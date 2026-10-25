"""
PolyA-site analysis: assign transcribed orientation to single-exon alignments.

Port of scripts/polyA_site_transcript_mapper.dbi. Runs after validation and
before clustering/assembly (Launch_PASA_pipeline.pl runs it when transcripts
were polyA-trimmed by seqclean, i.e. -T -u).

For each valid alignment of a transcript that seqclean trimmed a poly-A (3')
or poly-T (5') tail from:
1. Infer the polyA site and transcribed strand from the aligned orientation and
   which end was trimmed. Skip if the alignment ends more than 2 bp short of
   the trim point.
2. Accept the site if the trimmed tail is >= 80% A and the 10 bp of genome just
   past the site are < 80% A (otherwise likely internal priming).
3. Record the site, and set spliced_orient on single-exon alignments whose
   spliced_orient is still unknown.
"""

import logging
import re
from dataclasses import dataclass
from typing import Dict, List, Optional

from pypasa.models import Alignment, Orientation

logger = logging.getLogger(__name__)

FLEX_SIZE = 2  # bp of unaligned sequence allowed next to the polyA site
WINDOW = 10  # bp of context on each side of the site
MIN_PERCENT_A = 80.0

_TRIMPOLY = re.compile(r"trimpoly\[\+(\d+),\s-(\d+)\]")
_COMPLEMENT = str.maketrans(
    "acgtrymkswhbvdnxACGTRYMKSWHBVDNX", "tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX"
)


@dataclass
class PolyATrim:
    acc: str
    type: str  # 'lend' (poly-T trimmed from 5') or 'rend' (poly-A trimmed from 3')
    trim: int
    untrimmed_length: int


@dataclass
class PolyASite:
    align_id: int
    genomic_coord: int
    transcribed_orient: str
    trim_type: str  # 'lend_PolyT' or 'rend_PolyA'
    align_end_offset: int


def parse_cln(path: str) -> Dict[str, PolyATrim]:
    """
    Read seqclean's .cln report for single-ended polyA/T trims.

    Transcripts trimmed at both ends are dropped, as in Perl.
    """
    trims: Dict[str, PolyATrim] = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            acc = re.sub(r"\s", "", fields[0])
            seqlength = int(fields[4])
            for m in _TRIMPOLY.finditer(line):
                lend_trim, rend_trim = int(m.group(1)), int(m.group(2))
                if lend_trim and not rend_trim:
                    end_type, size = "lend", lend_trim
                elif rend_trim and not lend_trim:
                    end_type, size = "rend", rend_trim
                else:
                    continue
                if acc in trims:
                    del trims[acc]  # Perl: "Eventually, restore this part."
                else:
                    trims[acc] = PolyATrim(acc, end_type, size, seqlength)
    return trims


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(_COMPLEMENT)


def percent_a(seq: str) -> float:
    return 100.0 * seq.upper().count("A") / len(seq)


def _perl_substr(s: str, start: int, length: int) -> str:
    if start < 0:
        start += len(s)
        if start < 0:
            return ""
    return s[start:start + length]


def _strip_window(s: str) -> str:
    """Perl's s/^\\w{10}//: drop the first WINDOW chars only if there are that many."""
    return s[WINDOW:] if re.match(rf"\w{{{WINDOW}}}", s) else s


def infer_site(aln: Alignment, trim: PolyATrim) -> PolyASite:
    """Putative polyA site and transcribed strand for one alignment."""
    lend, rend = aln.genome_lend, aln.genome_rend
    first, last = aln.segments[0], aln.segments[-1]
    # transcript coords at the alignment's genomic left and right ends
    if Orientation(aln.orient) == Orientation.MINUS:
        mlend, mrend = first.transcript_rend, last.transcript_lend
    else:
        mlend, mrend = first.transcript_lend, last.transcript_rend
    trimmed_length = trim.untrimmed_length - trim.trim

    plus = Orientation(aln.orient) != Orientation.MINUS
    if plus and trim.type == "lend":
        diff = mlend - 1
        site, transcribed = lend - diff, "-"
    elif plus and trim.type == "rend":
        diff = trimmed_length - mrend
        site, transcribed = rend + diff, "+"
    elif trim.type == "lend":
        diff = mrend - 1
        site, transcribed = rend + diff, "+"
    else:
        diff = trimmed_length - mlend
        site, transcribed = lend - diff, "-"

    trim_type = "lend_PolyT" if trim.type == "lend" else "rend_PolyA"
    return PolyASite(aln.align_id, site, transcribed, trim_type, diff)


def evaluate_site(site: PolyASite, trim: PolyATrim, genome_seq: str, cdna_seq: str) -> Optional[str]:
    """
    Check a putative site against genome and untrimmed transcript sequence.

    Returns None if the site is accepted (adjusting site.genomic_coord to take in
    a leading A, as Perl does), else the reason it was ignored.
    """
    pos = site.genomic_coord - 1 if site.transcribed_orient == "-" else site.genomic_coord
    region = _perl_substr(genome_seq, pos - WINDOW, 2 * WINDOW)
    if site.transcribed_orient == "-":
        region = reverse_complement(region)

    if trim.type == "lend":
        tail = reverse_complement(cdna_seq[: trim.trim + WINDOW])
    else:
        tail = _perl_substr(cdna_seq, len(cdna_seq) - trim.trim - WINDOW, trim.trim + WINDOW)

    trimmed_part = _strip_window(tail)
    genomic_part = _strip_window(region)
    if not trimmed_part or not genomic_part:
        return "missing end of sequence"
    if percent_a(genomic_part) >= MIN_PERCENT_A:
        return "potential internal genomic priming"
    if percent_a(trimmed_part) < MIN_PERCENT_A:
        return "trimmed region not A-rich"

    if genomic_part[0].upper() == "A" and trimmed_part[0].upper() == "A":
        site.genomic_coord += 1 if site.transcribed_orient == "+" else -1
    return None


class PolyASiteMapper:
    """Runs the polyA-site stage over in-memory alignments."""

    def __init__(self, trims: Dict[str, PolyATrim], genome: Dict[str, str], untrimmed_transcripts: Dict[str, str]):
        """
        Args:
            trims: from parse_cln() on seqclean's .cln report
            genome: scaffold -> sequence
            untrimmed_transcripts: accession -> sequence before seqclean trimming
        """
        self.trims = trims
        self.genome = genome
        self.untrimmed_transcripts = untrimmed_transcripts
        self.stats = {"candidates": 0, "too_far": 0, "ignored": 0, "sites": 0, "orients_set": 0}

    def run(self, alignments: List[Alignment]) -> List[PolyASite]:
        """
        Find polyA sites on valid alignments, and set spliced_orient on valid
        single-exon alignments that don't have one yet.

        Returns:
            Accepted polyA sites
        """
        sites = []
        for aln in alignments:
            if not aln.valid:
                continue
            trim = self.trims.get(aln.cdna_acc)
            if trim is None:
                continue
            self.stats["candidates"] += 1
            site = infer_site(aln, trim)
            if site.align_end_offset > FLEX_SIZE:
                self.stats["too_far"] += 1
                continue
            cdna_seq = self.untrimmed_transcripts.get(aln.cdna_acc)
            if not cdna_seq:
                raise KeyError(f"{aln.cdna_acc} not among untrimmed transcripts")
            reason = evaluate_site(site, trim, self.genome[aln.scaffold], cdna_seq)
            if reason:
                self.stats["ignored"] += 1
                logger.debug(f"{aln.align_acc}: ignored, {reason}")
                continue
            sites.append(site)
            if aln.num_segments == 1 and aln.spliced_orient in (None, Orientation.UNKNOWN):
                aln.spliced_orient = Orientation(site.transcribed_orient)
                self.stats["orients_set"] += 1

        self.stats["sites"] = len(sites)
        logger.info(f"PolyA sites: {self.stats}")
        return sites
