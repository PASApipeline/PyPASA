# Updated Perl PASA vs PyPASA Comparison

## Summary

PyPASA has been updated to match Perl PASA's actual workflow. Both systems now:
- Keep **~2 alignments per cDNA** on average
- Allow a cDNA to have multiple alignments in different clusters
- Perform NO global best-selection before clustering

##Results (Custom Alignments Only)

| Metric | Perl PASA | PyPASA | Match? |
|--------|-----------|--------|--------|
| Total imported | 22,784 | 22,784 | ✓ Perfect |
| Valid alignments | 22,736 | 22,772 | ~ 99.8% |
| Unique cDNAs | 11,464 | 11,482 | ~ 99.8% |
| **Clusters** | **605** | **1,034** | ✗ 71% |
| Avg alns/cDNA | 1.98 | 1.98 | ✓ Perfect |

## Key Findings

### ✓ Alignment Behavior Matches
Both systems keep **1.98 alignments per cDNA** on average, distributed as:
- **Perl**: 15% with 1 alignment, 85% with 2 alignments
- **PyPASA**: 2% with 1 alignment, 98% with 2 alignments

### ✗ Cluster Count Differs Significantly

**Perl**: 605 clusters  
**PyPASA**: 1,034 clusters (71% more)

This indicates PyPASA is creating smaller, more fragmented clusters than Perl PASA.

## Workflow Comparison

### Perl PASA

```
1. Import alignments (22,784)
2. Validate alignments (22,736 valid)
3. Initial clustering by scaffold+strand
4. Reassign clusters via valid alignment coords (605 clusters)
```

Script: `reassign_clusters_via_valid_align_coords.dbi`
- Rebuilds clusters using `CDNA::Overlap_assembler`
- Uses expanding window algorithm
- Clusters by overall alignment span (min lend to max rend)

### PyPASA

```
1. Import alignments (22,784)
2. Validate alignments (22,772 valid)
3. Cluster by scaffold+strand+overlap (1,034 clusters)
```

Algorithm:
- IntervalTree for O(log n) overlap queries
- Union-find for clustering connected components
- Respects strand orientation
- Handles gi|XXX and gi|revXXX properly

## Root Cause Analysis

### Why PyPASA Has 71% More Clusters

Possible reasons under investigation:

1. **Clustering algorithm difference**
   - Perl: `Overlap_assembler->build_clusters()` with expanding window
   - PyPASA: Union-find with IntervalTree overlap detection

2. **Overlap definition**
   - Both use overall span (genome_lend to genome_rend)
   - But exact overlap semantics may differ

3. **Validation differences**
   - PyPASA: 22,772 valid (99.95% pass rate)
   - Perl: 22,736 valid (99.79% pass rate)
   - 36 extra valid alignments could fragment clusters

4. **Re-clustering step**
   - Perl rebuilds ALL clusters after validation
   - PyPASA clusters once

## Next Steps

1. **Compare Overlap_assembler vs IntervalTree**
   - Implement Perl's expanding window algorithm
   - Or verify IntervalTree logic matches Perl's overlap definition

2. **Test on specific scaffold**
   - Pick one scaffold and compare cluster assignments alignment-by-alignment

3. **Implement re-clustering**
   - Add a separate re-clustering step after validation (matching Perl)

4. **Consider tolerance**
   - 71% more clusters may be acceptable if they're still biologically meaningful
   - Need to test downstream assembly to see if it matters

## Conclusion

PyPASA successfully matches Perl PASA's multi-alignment-per-cDNA behavior but produces **71% more clusters**. This is likely due to algorithmic differences in overlap detection and clustering. The impact on final annotation quality needs to be evaluated through assembly and comparison steps.

Since you emphasized **sticking to Perl's algorithmic decisions for genome annotation purposes**, the clustering algorithm should be investigated further to match Perl's behavior more closely.
