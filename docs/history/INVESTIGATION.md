# Investigation: PyPASA vs Perl PASA Clustering Differences

## Current Results (after matching validation thresholds)

| Metric | PyPASA | Perl PASA | Ratio |
|--------|--------|-----------|-------|
| **Imported** | 22,784 | unknown | - |
| **Valid** | 22,772 | unknown | - |
| **After Best-Selection** | 11,500 | 5,811 | 1.98x |
| **Clusters** | 553 | 303 | 1.82x |

**Validation Thresholds (now matching)**:
- MIN_AVG_PER_ID: 95%
- MIN_PERCENT_ALIGNED: 75%

## Key Perl PASA Findings

### 1. Clustering Algorithm
- **File**: `PerlLib/CDNA/Overlap_assembler.pm`
- **Method**: Simple expanding window algorithm
- **Steps**:
  1. Sort alignments by left coordinate (`lend`)
  2. Maintain expanding window `(exp_left, exp_right)`
  3. Check overlap: `exp_left <= rend && exp_right >= lend`
  4. Expand window to include overlapping alignments
  5. Create new cluster when no overlap

**This is essentially what PyPASA does!**

### 2. Best Alignment Selection
- **Finding**: Perl does `best_selector.py` type filtering BUT it's implicit
- **Location**: Not explicitly in code we've found so far
- **Hypothesis**: Either happens during:
  1. Alignment import (implicit filtering)
  2. Or through the `-N` parameter for GMAP/BLAT runs (but NOT custom alignments!)

### 3. Validation Thresholds (Sample Data)
From `sample_data/sqlite.confs/alignAssembly.config`:
```
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
```
✅ **PyPASA now matches these!**

### 4. Database Schema Observation
- Perl: ONE `align_link` record per alignment → multiple `alignment` segments
- PyPASA: ONE `Alignment` object → multiple `AlignmentSegment` objects
- ✅ **Structure matches!**

## Why PyPASA Has 2x More Alignments

### Theory 1: Scaffold Normalization (CONFIRMED ISSUE)
- PyPASA normalizes: `gi|revXXX` → `gi|XXX`
- But still treats them as **separate scaffolds** for clustering!
- **Problem**: We normalize for best-selection but NOT for clustering!

**Example**:
```
c1000_g0_i1 aligns to:
  - gi|68713 at 3407-4848 (+)     → cluster A on gi|68713
  - gi|rev68713 at 102102-103543 (-) → cluster B on gi|rev68713
```

These should be in ONE cluster (same physical location), but PyPASA creates TWO clusters!

### Theory 2: Different GFF3 Parsing
- Perl may group/merge alignments differently during GFF3 import
- **Need to check**: Are we creating more `align_id` records than Perl?

### Theory 3: Alignment Span Calculation
Looking at Perl code in `reassign_clusters_via_valid_align_coords.dbi`:
```perl
sub get_alignment_span {
    my ($dbproc, $align_id) = @_;
    my $query = "select lend, rend from alignment where align_id = $align_id";
    my @results = &do_sql_2D($dbproc, $query);
    my @coords;
    foreach my $result (@results) {
        push (@coords, @$result);
    }
    @coords = sort {$a<=>$b} @coords;
    my $min = shift @coords;
    my $max = pop @coords;
    return ($min, $max);
}
```

This gets ALL segment coordinates and takes min/max. PyPASA does the same in `Alignment.genome_lend` and `Alignment.genome_rend`.

## Critical Bug Found!

### The Clustering Scaffold Problem

**PyPASA clusters PER SCAFFOLD**:
```python
for scaffold in scaffolds:
    alignments = datastore.get_alignments_by_scaffold(scaffold)
    # cluster only these alignments
```

**Problem**: `gi|68713` and `gi|rev68713` are treated as DIFFERENT scaffolds!

This means:
1. Alignments on `gi|68713` cluster together
2. Alignments on `gi|rev68713` cluster SEPARATELY
3. Even if they're at the SAME physical location!

**Solution**: Before clustering, we need to:
1. Group scaffolds by normalized name
2. Load alignments from ALL scaffold variants together
3. Cluster the combined set

## Action Plan

### Immediate Fix (Priority 1)
✅ **DONE**: Match validation thresholds (95%/75%)
🔧 **TODO**: Fix clustering to handle scaffold normalization properly

### Code Change Needed
In `pypasa/alignment/clusterer.py`:
```python
# WRONG (current):
for scaffold in scaffolds:
    alignments = get_alignments_by_scaffold(scaffold)
    cluster(alignments)

# RIGHT (needed):
scaffold_groups = group_scaffolds_by_normalized_name(scaffolds)
for scaffold_group in scaffold_groups:
    all_alignments = []
    for scaffold in scaffold_group:
        all_alignments.extend(get_alignments_by_scaffold(scaffold))
    cluster(all_alignments)
```

### Verification Steps
1. Fix clustering scaffold grouping
2. Re-run test
3. Compare results:
   - Should get ~50% fewer clusters (from 553 → ~280-300)
   - Should match Perl's 303 clusters closely!

## Additional Investigation Needed

1. **Double-check alignment count**: 
   - Run: `python -c "import sqlite3; conn = sqlite3.connect('test_output/pypasa_test.sqlite'); print(conn.execute('SELECT COUNT(DISTINCT cdna_acc) FROM alignments WHERE valid=1').fetchone())"`
   - This should tell us unique cDNAs

2. **Check Perl database**:
   - If we can access the Perl PASA database, check:
     - `SELECT COUNT(*) FROM align_link WHERE validate=1`
     - `SELECT COUNT(DISTINCT align_acc) FROM align_link WHERE validate=1`

3. **Compare specific alignments**:
   - Pick a cDNA that has 2 alignments in PyPASA
   - Check if Perl has 1 or 2

## Expected Outcome After Fix

| Metric | Current | After Fix | Perl PASA |
|--------|---------|-----------|-----------|
| **Alignments** | 11,500 | 11,500 | 5,811 |
| **Clusters** | 553 | ~300 | 303 |

The cluster count should match Perl closely after fixing the scaffold grouping!

The alignment count difference (11,500 vs 5,811) is still concerning and needs further investigation, but it won't affect clustering correctness.
