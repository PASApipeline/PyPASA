"""
Write a run's assemblies in Perl PASA's output formats.

  <name>.pasa_assemblies.gff3          PASA_transcripts_and_assemblies_to_GFF3.dbi -a
  <name>.pasa_assemblies.gtf           ... -a -T
  <name>.pasa_assemblies.bed           ... -a -B
  <name>.assemblies.fasta              populate_mysql_assembly_sequence_field.dbi
  <name>.pasa_assemblies_described.txt describe_alignment_assemblies_cgi_convert.dbi

Content matches Perl line for line. Two deliberate differences, both in
labels or ordering only:
- subcluster ids (GTF gene_id, described.txt) follow PyPASA's numbering
  (cluster, then first assembly), not slclust's traversal order;
- BED lines are in assembly order rather than Perl's hash order.
Perl's assembly row ids (align_<id>) are reproduced as
number-of-imported-alignments + assembly number.
"""

import re
from typing import Dict, List, TextIO, Tuple

from pypasa.alignment.validator import AlignmentValidator
from pypasa.models import Assembly, Orientation
from pypasa.run import PipelineRun

_COMPLEMENT = str.maketrans("ACGTacgtyrkmYRKM", "TGCAtgcarymkRYMK")


def _row_id(run: PipelineRun, asm: Assembly) -> int:
    return len(run.alignments) + asm.assembly_id


def _acc(asm: Assembly) -> str:
    return f"asmbl_{asm.assembly_id}"


def _subcluster_of(run: PipelineRun) -> Dict[int, int]:
    return {aid: s.subcluster_id for s in run.subclusters for aid in s.assembly_ids}


def write_gff3(run: PipelineRun, fh: TextIO, run_name: str) -> None:
    for asm in run.assemblies:
        # spliced orient falls back to aligned orient, so Target strand is always '+' here
        spliced = asm.spliced_orient if asm.spliced_orient != Orientation.UNKNOWN else asm.orient
        rel_orient = "+" if spliced == asm.orient else "-"
        for seg in asm.segments:
            fh.write("\t".join([
                asm.scaffold, f"assembler-{run_name}", "cDNA_match",
                str(seg.genome_lend), str(seg.genome_rend), ".", asm.orient.value, ".",
                f"ID=align_{_row_id(run, asm)};Target={_acc(asm)} "
                f"{seg.transcript_lend} {seg.transcript_rend} {rel_orient}",
            ]) + "\n")


def write_gtf(run: PipelineRun, fh: TextIO) -> None:
    subcluster = _subcluster_of(run)
    for asm in run.assemblies:
        attrs = (f'gene_id "PASA_cluster_{subcluster[asm.assembly_id]}"; '
                 f'transcript_id "align_id:{_row_id(run, asm)}|{_acc(asm)}";')
        rows = [("transcript", asm.lend, asm.rend)]
        rows += [("exon", s.genome_lend, s.genome_rend) for s in asm.segments]
        for feature, lend, rend in rows:
            fh.write("\t".join([
                asm.scaffold, "assembler", feature, str(lend), str(rend), ".",
                asm.orient.value, ".", attrs,
            ]) + "\n")
        fh.write("\n")  # Perl separates transcripts with a blank line


def write_bed(run: PipelineRun, fh: TextIO) -> None:
    for asm in run.assemblies:
        start = asm.lend - 1
        sizes = ",".join(str(s.genome_rend - s.genome_lend + 1) for s in asm.segments)
        starts = ",".join(str(s.genome_lend - 1 - start) for s in asm.segments)
        fh.write("\t".join([
            asm.scaffold, str(start), str(asm.rend), f"ID={_acc(asm)}", "0", asm.orient.value,
            str(start), str(asm.rend), "0", str(len(asm.segments)), sizes, starts,
        ]) + "\n")


def spliced_sequence(asm: Assembly, genome_seq: str) -> str:
    seq = "".join(genome_seq[s.genome_lend - 1:s.genome_rend] for s in asm.segments)
    if asm.orient == Orientation.MINUS:
        seq = seq[::-1].translate(_COMPLEMENT)
    return seq.upper()


def write_fasta(run: PipelineRun, fh: TextIO, genome: Dict[str, str]) -> None:
    # Perl orders by scaffold only; ties keep assembly order
    for asm in sorted(run.assemblies, key=lambda a: a.scaffold):
        seq = re.sub(r"(\w{60})", r"\1\n", spliced_sequence(asm, genome[asm.scaffold]))
        fh.write(f">{_acc(asm)}\n{seq}\n")


def _junction_flags(asm: Assembly, genome_seq: str) -> Tuple[str, List[Tuple[str, str, bool, bool]]]:
    """
    Re-run Perl's splice check on an assembly (CDNA_alignment.pm
    identify_splice_junctions + validate_splice_junctions).

    Returns:
        the computed spliced orientation ('?' if none validates), and per
        intron (donor chars, acceptor chars, donor is consensus, acceptor is consensus)
    """
    introns = []
    for prev, curr in zip(asm.segments, asm.segments[1:]):
        donor = genome_seq[prev.genome_rend:prev.genome_rend + 2].upper()
        acceptor = genome_seq[curr.genome_lend - 3:curr.genome_lend - 1].upper()
        introns.append((donor, acceptor))

    pairs = AlignmentValidator.CONSENSUS_SPLICE_PAIRS
    aligned = asm.orient if asm.orient != Orientation.UNKNOWN else Orientation.PLUS
    opposite = Orientation.MINUS if aligned == Orientation.PLUS else Orientation.PLUS
    computed = "?"
    if not introns:  # Perl only checks splice sites on multi-segment alignments
        return computed, []
    for orient in (aligned, opposite):
        if all(p in pairs[orient] for p in introns):
            computed, flag_orient = orient.value, orient
            break
    else:
        # fewest nonconsensus introns; ties go to '+' (Perl sorts ('+', '-') by error length)
        bad = {o: sum(p not in pairs[o] for p in introns) for o in (Orientation.PLUS, Orientation.MINUS)}
        flag_orient = min((Orientation.PLUS, Orientation.MINUS), key=lambda o: bad[o])

    flags = []
    for donor, acceptor in introns:
        if (donor, acceptor) in pairs[flag_orient]:
            flags.append((donor, acceptor, True, True))
        else:
            flags.append((
                donor, acceptor,
                any(d == donor for d, _ in pairs[flag_orient]),
                any(a == acceptor for _, a in pairs[flag_orient]),
            ))
    return computed, flags


def alignment_text(asm: Assembly, genome_seq: str) -> str:
    """Perl's CDNA_alignment::toToken, e.g. 'orient(a+/s+) align: 1(1)-10(10)>GT....AG<31(11)-40(20)'."""
    computed, flags = _junction_flags(asm, genome_seq)
    spliced = computed if computed != "?" else asm.spliced_orient.value
    parts = []
    for i, seg in enumerate(asm.segments):
        token = ""
        if i > 0:
            _, acceptor, _, ok = flags[i - 1]
            token += (acceptor.upper() if ok else acceptor.lower()) + "<"
        # transcript coords at the segment's genomic left/right ends
        if asm.orient == Orientation.MINUS:
            mleft, mright = seg.transcript_rend, seg.transcript_lend
        else:
            mleft, mright = seg.transcript_lend, seg.transcript_rend
        token += f"{seg.genome_lend}({mleft})-{seg.genome_rend}({mright})"
        if i < len(asm.segments) - 1:
            donor, _, ok, _ = flags[i]
            token += ">" + (donor.upper() if ok else donor.lower())
        parts.append(token)
    return f"orient(a{asm.orient.value}/s{spliced}) align: " + "....".join(parts)


def write_described(run: PipelineRun, fh: TextIO, genome: Dict[str, str]) -> None:
    subcluster = _subcluster_of(run)
    cdna_of = {a.align_acc: a.cdna_acc for a in run.alignments}
    fh.write("\t".join(["#scaffold", "subcluster_id", "asmbl_acc", "cdna_accs", "alignment_description"]) + "\n")
    for asm in run.assemblies:
        # Perl lists members in align_acc order (asmbl_link primary key), as cDNA accessions
        cdnas = ",".join(cdna_of[acc] for acc in sorted(asm.component_accs))
        fh.write("\t".join([
            asm.scaffold, str(subcluster[asm.assembly_id]), _acc(asm), cdnas,
            alignment_text(asm, genome[asm.scaffold]),
        ]) + "\n")


def write_all(run: PipelineRun, prefix: str, run_name: str, genome: Dict[str, str]) -> List[str]:
    """Write all assembly outputs as <prefix>.<suffix>; returns the paths."""
    writers = [
        ("pasa_assemblies.gff3", lambda fh: write_gff3(run, fh, run_name)),
        ("pasa_assemblies.gtf", lambda fh: write_gtf(run, fh)),
        ("pasa_assemblies.bed", lambda fh: write_bed(run, fh)),
        ("assemblies.fasta", lambda fh: write_fasta(run, fh, genome)),
        ("pasa_assemblies_described.txt", lambda fh: write_described(run, fh, genome)),
    ]
    paths = []
    for suffix, write in writers:
        path = f"{prefix}.{suffix}"
        with open(path, "w") as fh:
            write(fh)
        paths.append(path)
    return paths
