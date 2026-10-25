# Perl PASA Clustering Analysis

## Complete Clustering Hierarchy

Perl PASA has a **two-level clustering system**:

### Level 1: CLUSTERS (605 clusters for sample data)
- **Purpose**: Group raw cDNA alignments by genomic location and strand
- **When**: After validation, before assembly
- **Algorithm**: Expanding window (Overlap_piler.pm)
- **Input**: Valid cDNA alignments (22,736 alignments)
- **Output**: 605 clusters

### Level 2: SUBCLUSTERS (667 subclusters for sample data)
- **Purpose**: Identify alternative splice variants within assembled transcripts
- **When**: After assembly
- **Algorithm**: Single-linkage clustering with ≥50% overlap requirement
- **Input**: PASA assemblies (853 assemblies from subcluster_link table)
- **Output**: 667 subclusters (1.09 subclusters per cluster)
- **Distribution**:
  - 559 clusters → 1 subcluster (91.5%)
  - 48 clusters → 2 subclusters
  - 4 clusters → 3 subclusters

## The Clustering Algorithm: Expanding Window

Perl's `Overlap_piler.pm` uses an **expanding window** approach:

```perl
# Pseudocode
sort alignments by lend (left coordinate)
window = first alignment coordinates
cluster = [first alignment]

for each subsequent alignment:
    if alignment overlaps window:
        # KEY: Expand window to include this alignment
        window.left = min(window.left, alignment.lend)
        window.right = max(window.right, alignment.rend)
        cluster.append(alignment)
    else:
        # No overlap with window - close this cluster
        save cluster
        window = alignment coordinates
        cluster = [alignment]

save final cluster
```

### Key Characteristic: Transitive Clustering

The expanding window creates **transitive clusters**:
- Alignment A overlaps B
- Alignment B overlaps C  
- Even if A doesn't directly overlap C, they're in the same cluster!

Example:
```
A: |-------|               (positions 100-200)
B:     |-------|           (positions 150-250)
C:         |-------|       (positions 200-300)
```

Expanding window process:
1. Start with A: window = 100-200
2. B overlaps window: expand to 100-250, cluster = [A, B]
3. C overlaps window: expand to 100-300, cluster = [A, B, C]
4. Result: ONE cluster containing A, B, C

## PyPASA's Current Algorithm: Union-Find

PyPASA uses **union-find with direct overlap**:

```python
for each pair of alignments:
    if they directly overlap:
        union(alignment1, alignment2)
```

This creates more restrictive clusters:
- A overlaps B → cluster [A, B]
- B overlaps C → cluster [B, C]
- If A doesn't directly overlap C, this still merges into one cluster via union-find

**Result**: Should be similar to Perl's expanding window!

## Why PyPASA Has 71% More Clusters (1,034 vs 605)

Possible reasons:

### 1. IntervalTree Overlap Detection
PyPASA uses IntervalTree which may have subtle differences in overlap semantics:
- Intervals are `[start, end)` half-open
- Perl uses `[lend, rend]` closed intervals

### 2. Scaffold Normalization Issue
PyPASA groups `gi|XXX` and `gi|revXXX` together, but may not be doing it correctly.

### 3. Validation Differences
- Perl: 22,736 valid (99.79%)
- PyPASA: 22,772 valid (99.95%)
- 36 extra alignments could create fragmentation

### 4. Algorithm Implementation Bug
Need to verify PyPASA's union-find is truly transitive.

## Workflow Comparison

### Perl PASA (Complete)
```
1. Import alignments (22,784)
2. Validate (22,736 valid)
3. Cluster by strand+overlap (605 clusters)          ← Expanding window
4. Assemble clusters (creates assemblies)
5. Subcluster assemblies (667 subclusters)           ← ≥50% overlap
```

### PyPASA (Current)
```
1. Import alignments (22,784)
2. Validate (22,772 valid)
3. Cluster by strand+overlap (1,034 clusters)        ← Union-find
```

## Recommendations

### Option 1: Implement Perl's Expanding Window Algorithm
Replace PyPASA's union-find with Perl's exact expanding window algorithm from `Overlap_piler.pm`.

**Pros:**
- Guaranteed to match Perl's cluster count
- Proven algorithm for genome annotation

**Cons:**
- O(n²) worst case vs O(n log n) for IntervalTree

### Option 2: Debug IntervalTree Implementation
Keep IntervalTree but fix the overlap logic to match Perl exactly.

**Pros:**
- Better performance (O(n log n))
- Modern data structure

**Cons:**
- May not perfectly match Perl due to subtle differences

### Option 3: Hybrid Approach
Use IntervalTree for initial overlap detection, then post-process with expanding window logic.

## Next Steps

1. **Pick one scaffold** and manually compare clustering output between Perl and PyPASA
2. **Implement expanding window** in PyPASA to match Perl exactly
3. **Verify cluster count** matches 605
4. **Document** that subclustering will be implemented after assembly step
