# Single-Exon Alignment Clustering Analysis

## The Question

How do single-exon alignments (with unknown spliced orientation '?') get integrated into clusters and subclusters?

## Background

Single-exon alignments:
- Have only ONE alignment segment (no introns)
- Cannot determine transcribed orientation from splice sites
- Get `spliced_orient = '?'` in Perl or `spliced_orient = NULL` in PyPASA
- Still have an `orient` value (+ or -) from the alignment direction

## Perl PASA Behavior

### Initial Clustering (Default Mode - Orientation-Agnostic)

**Script**: `reassign_clusters_via_valid_align_coords.dbi`

```perl
# Query does NOT filter by orientation
my $query = "select al.align_id from align_link al, clusters c "
    . " where c.cluster_id = al.cluster_id and c.annotdb_asmbl_id = ? 
    . " and al.validate = 1";

# Uses Overlap_assembler which only checks coordinates
```

**Behavior for single-exon alignments**:
- Clustered with ANY overlapping alignment (regardless of orientation)
- Single-exon '+' orient can cluster with multi-exon '-' spliced_orient
- Single-exon '-' orient can cluster with multi-exon '+' spliced_orient
- **Result**: Mixed-orientation clusters

**Example**: Perl cluster 81 has:
- 112 alignments with '-' spliced orientation
- 6 alignments with '+' spliced orientation  
- **20 alignments with '?' spliced orientation** (single-exon)

### Subclustering (Orientation-Specific)

**Script**: `subcluster_builder.dbi`

```perl
# Line 166: Strict orientation matching
unless ($i_spliced_orientation eq $j_spliced_orientation) {
    next; #must have same spliced orientation
}
```

**Behavior for single-exon alignments**:
- Uses exact string comparison: `'?' eq '?'` → TRUE, `'?' eq '+'` → FALSE
- Single-exon alignments ONLY cluster with other single-exon alignments
- They are **isolated** from multi-exon alignments
- Form separate subclusters containing only '?' orientation alignments

**Impact**:
- Single-exon alignments get their own subclusters
- Cannot be assembled with multi-exon alignments that show clear strand
- May miss legitimate isoforms if a gene has both single-exon and multi-exon transcripts

## PyPASA Behavior

### Initial Clustering (Orientation-Aware)

**Code**: `pypasa/storage/models.py` and `pypasa/alignment/clusterer.py`

```python
@property
def effective_orient(self) -> Orientation:
    """Uses spliced_orient if available, otherwise falls back to orient."""
    return self.spliced_orient if self.spliced_orient else self.orient

# Clustering groups by (scaffold, effective_orient)
key = (aln.scaffold, aln.effective_orient)
```

**Behavior for single-exon alignments**:
- `spliced_orient = NULL` → falls back to `orient` value
- Single-exon '+' orient → `effective_orient = '+'`
- Single-exon '-' orient → `effective_orient = '-'`
- Clustered with multi-exon alignments of same strand

**Example**: From PyPASA test data:
- Cluster 2: 22 single-exon '+' orient + 85 multi-exon '+' spliced + 14 flipped '-'/'+' 
- Cluster 37: 13 single-exon '-' orient (all single-exon)
- Cluster 38: 25 single-exon '-' orient (all single-exon)

### Subclustering (Not Yet Implemented)

PyPASA doesn't currently have subclustering, but if/when implemented:

**Option A - Match Perl's exact behavior**:
```python
# Treat NULL/None as distinct from '+' and '-'
if aln1.spliced_orient != aln2.spliced_orient:
    continue  # Don't subcluster
```
- Single-exon alignments form separate subclusters
- Consistent with Perl

**Option B - Use effective_orient**:
```python
# Use fallback logic for single-exon
if aln1.effective_orient != aln2.effective_orient:
    continue  # Don't subcluster
```
- Single-exon alignments subcluster with same-strand multi-exon alignments
- More biologically intuitive

## Comparison Table

| Aspect | Perl Default Clustering | Perl Subclustering | PyPASA Clustering |
|--------|------------------------|-------------------|-------------------|
| Single-exon '+' with multi-exon '+' | ✓ Together | ✗ Separate | ✓ Together |
| Single-exon '-' with multi-exon '-' | ✓ Together | ✗ Separate | ✓ Together |
| Single-exon '+' with multi-exon '-' | ✓ Together | ✗ Separate | ✗ Separate |
| Single-exon '?' with single-exon '?' | ✓ Together | ✓ Together | N/A (uses orient) |

## Biological Implications

### Perl's Approach (Exact Match in Subclustering)

**Pros**:
- Conservative: doesn't assume single-exon orientation is correct
- Avoids mixing if single-exon alignment is mismapped

**Cons**:
- Single-exon transcripts can't be assembled with multi-exon isoforms
- Misses legitimate alternative isoforms
- Creates separate subclusters for single-exon that may represent noise

### PyPASA's Current Approach (Use Aligned Orient)

**Pros**:
- Single-exon transcripts cluster with same-strand multi-exon transcripts
- Better for genes with both single-exon and multi-exon isoforms
- Uses available information (aligned orientation)

**Cons**:
- Assumes aligned orientation is biologically correct
- Could mix incorrect mappings if single-exon is on wrong strand

## Recommendations

### For Initial Clustering:

**PyPASA should offer both modes**:
1. **Orientation-agnostic mode** (match Perl default) - for backward compatibility
2. **Strand-specific mode** (current) - for better biological accuracy

Configuration option:
```python
--clustering-mode [default|strand-specific]
```

### For Subclustering (when implemented):

**Consider the biology**:
- If single-exon alignments are high-quality and trusted → use `effective_orient`
- If single-exon alignments are noisy or uncertain → use exact match ('?' separate)

Possible configuration:
```python
--single-exon-handling [strict|permissive]
# strict: '?' only matches '?' (Perl behavior)
# permissive: use effective_orient fallback
```

## Current State Summary

1. **Perl default clustering**: Single-exon alignments cluster with everything (orientation-agnostic)
2. **Perl subclustering**: Single-exon alignments isolated from multi-exon (exact '?' matching)
3. **PyPASA clustering**: Single-exon alignments cluster by their aligned orientation
4. **Result**: PyPASA has more clusters (713 vs 303) because it enforces strand-specificity from the start
