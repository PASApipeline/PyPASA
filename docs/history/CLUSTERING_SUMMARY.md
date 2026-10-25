# SUMMARY: Perl PASA Clustering Research Findings

## You Were Right!

Perl PASA does indeed use **two-pass clustering**:

### Pass 1: CLUSTERS (Strand-based grouping)
- **Algorithm**: Expanding window approach (`Overlap_piler.pm`)
- **Criteria**: Same strand + any overlap
- **Result**: 605 clusters from 22,736 valid alignments

### Pass 2: SUBCLUSTERS (Alternative splice variants)
- **Algorithm**: Single-linkage clustering with ≥50% overlap
- **Criteria**: Assemblies within a cluster that overlap ≥50%
- **When**: AFTER assembly (operates on assembled transcripts, not raw alignments)
- **Result**: 667 subclusters (most clusters have 1 subcluster, some split into 2-3)

## Key Insights

### The Expanding Window Algorithm

Perl's `Overlap_piler.pm` creates **transitive clusters** via an expanding window:

```
1. Sort alignments by left coordinate (lend)
2. Start window with first alignment
3. For each subsequent alignment:
   - If overlaps current window: EXPAND window, add to cluster
   - Else: Close cluster, start new window
4. Window expansion is KEY: once in a cluster, you stay even if 
   you don't overlap earlier members
```

This is different from a simple pairwise overlap approach!

### Why PyPASA Has More Clusters

**PyPASA: 1,034 clusters (71% more than Perl's 605)**

Potential causes:
1. **Algorithm difference**: Union-find may not perfectly replicate expanding window behavior
2. **IntervalTree semantics**: Half-open `[start, end)` vs closed `[lend, rend]` intervals
3. **Validation differences**: 36 more valid alignments (22,772 vs 22,736)

## Recommendation: Implement Expanding Window

To match Perl PASA's annotation behavior exactly, I recommend implementing Perl's expanding window algorithm in PyPASA. This is the most direct path to matching the 605 cluster count.

The algorithm is straightforward and well-defined in `Overlap_piler.pm` lines 84-150.

Would you like me to:
1. **Implement the expanding window algorithm** in PyPASA?
2. **First analyze a specific scaffold** to see exactly where clustering differs?
3. **Keep current approach** and document the clustering difference?

For genome annotation purposes where you want to stick closely to Perl's established behavior, option #1 (implement expanding window) is recommended.
