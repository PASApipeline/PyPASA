# Transcribed Orientation Issue

## Problem

PyPASA is currently clustering alignments by **aligned orientation** (which strand the alignment is on), but Perl PASA clusters by **transcribed orientation** (determined by splice dinucleotides).

## Background

There are TWO types of orientation for an alignment:

1. **Aligned Orientation** (`aligned_orient` in Perl, `orient` in PyPASA)
   - Which strand of the genome the alignment maps to (+ or -)
   - Simply the direction of the alignment

2. **Spliced/Transcribed Orientation** (`spliced_orient` in both)
   - The actual transcription direction based on splice site dinucleotides
   - GT-AG = forward strand (+)
   - CT-AC = reverse strand (-)
   - Can differ from aligned orientation!

## Why They Can Differ

Sometimes an alignment appears on the + strand, but the splice dinucleotides (CT-AC) indicate it's actually transcribed from the - strand. This happens when:
- The aligner maps to the wrong strand
- The transcript is from the opposite strand but happens to align

## Perl PASA Behavior

Found in `PASA_alignment_assembler.pm` lines 317-332:

```perl
sub force_flexorient {
    my $self = shift;
    my $orient = shift;
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    ## Fix orientations of fli and multi-segment alignments:
    for (my $i = 0; $i < $num_alignments; $i++) {
        my $alignment = $alignments_aref->[$i];
        my $num_segments = $alignment->get_num_segments();
        my $spliced_orientation = $alignment->get_spliced_orientation();
        
        if ($spliced_orientation =~ /^[+-]$/) { # set specifically
            $alignment->{fixed_orient} = $spliced_orientation;  ## USES SPLICED ORIENTATION
        } else {
            print "Setting $i to flex orient: $orient.\n" if $SEE;
            $alignment->{fixed_orient} = $orient;  ## Falls back to aligned orientation
        }
    }
}
```

**Key Point**: Perl PASA:
1. Checks if `spliced_orientation` is set (+ or -)
2. If YES: Uses `spliced_orientation` for clustering
3. If NO: Falls back to `aligned_orientation`

The assembler then runs separately for + and - orientations (lines 115-118):

```perl
$self->force_flexorient('+');
my @assemblies = $self->pasa_cpp_assemblies('+');
$self->force_flexorient('-');
push (@assemblies, $self->pasa_cpp_assemblies('-'));
```

## PyPASA Current Behavior

In `pypasa/alignment/clusterer.py` line 161:

```python
for aln in all_alignments:
    by_strand[aln.orient].append(aln)  # ← WRONG! Using aligned orientation
```

We're grouping by `aln.orient` (aligned orientation) instead of transcribed orientation.

## Schema Comparison

**Perl PASA** (`cdna_alignment_sqliteschema` lines 143-159):
```sql
CREATE TABLE align_link (
  ...
  aligned_orient TEXT DEFAULT NULL CHECK(aligned_orient IN ('', '+','-')),
  spliced_orient TEXT DEFAULT NULL CHECK(spliced_orient IN ('', '+','-', '?')),
  ...
)
```

**PyPASA** (`pypasa/storage/models.py` lines 55-67):
```python
@dataclass
class Alignment:
    ...
    orient: Orientation  # Aligned orientation
    spliced_orient: Optional[Orientation] = None  # Transcribed orientation
    ...
```

We have the field but aren't setting or using it!

## Required Fixes

### 1. Set spliced_orient During Validation

In `pypasa/alignment/validator.py`, add splice dinucleotide detection:

```python
def _determine_spliced_orient(self, alignment: Alignment, genome_seq: str) -> Optional[Orientation]:
    """
    Determine transcribed orientation from splice dinucleotides.
    
    Check if introns have canonical splice sites:
    - GT-AG (+ strand) = forward transcription
    - CT-AC (- strand) = reverse transcription
    
    Returns:
        Orientation.PLUS, Orientation.MINUS, or None if ambiguous
    """
    # For each intron (gap between segments):
    # - Extract donor/acceptor dinucleotides
    # - Check for GT-AG (forward) or CT-AC (reverse)
    # - If all introns agree: return that orientation
    # - If mixed or no introns: return None
```

### 2. Update Clustering to Use Transcribed Orientation

In `pypasa/alignment/clusterer.py` line 159-161:

```python
# OLD (WRONG):
for aln in all_alignments:
    by_strand[aln.orient].append(aln)

# NEW (CORRECT):
for aln in all_alignments:
    # Use transcribed orient if available, else fall back to aligned orient
    cluster_orient = aln.spliced_orient if aln.spliced_orient else aln.orient
    by_strand[cluster_orient].append(aln)
```

### 3. Add Helper Method to Alignment Model

In `pypasa/storage/models.py`:

```python
@property
def effective_orient(self) -> Orientation:
    """
    Get the effective orientation for clustering/assembly.
    
    Uses spliced_orient if determined from splice sites,
    otherwise falls back to aligned orient.
    
    This matches Perl PASA's force_flexorient behavior.
    """
    return self.spliced_orient if self.spliced_orient else self.orient
```

## Impact

This bug could explain:
1. **Why PyPASA has 758 clusters vs Perl's 605** - we're splitting some genes into separate + and - clusters when they should be together
2. **Assembly errors** - alignments from the same gene but flagged as different strands won't assemble together

## Priority

**HIGH** - This is a fundamental algorithmic difference that affects:
- Clustering accuracy
- Assembly correctness
- Gene structure prediction
- Downstream annotation quality

## References

- Perl PASA: `PerlLib/CDNA/PASA_alignment_assembler.pm` lines 115-118, 317-332
- Perl PASA: `PerlLib/CDNA/CDNA_alignment.pm` get_spliced_orientation()
- Perl Schema: `schema/cdna_alignment_sqliteschema` align_link table
