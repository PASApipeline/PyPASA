# Clustering Differences Between PyPASA and Perl PASA

## Summary

PyPASA creates **713 clusters** while Perl PASA creates **303 clusters** - a difference of **410 clusters (135% more in PyPASA)**.

## Root Cause: TBD - Algorithmic Difference in Clustering

Both PyPASA and Perl PASA:
- ✅ Treat each scaffold distinctly
- ✅ Group forward/reverse scaffold variants (`gi|68711` with `gi|rev68711`)
- ✅ Use strand-specific clustering
- ✅ Use transcribed (spliced) orientation when available

Yet PyPASA produces 2.4x more clusters than Perl.

## Detailed Analysis

### Cluster Size Distribution

**PyPASA** (713 clusters):
- Singletons: 106 (15%)
- 2-5 alignments: 248 (35%)
- 6-10 alignments: 105 (15%)
- 11-50 alignments: 196 (27%)
- >50 alignments: 58 (8%)
- Average: 17.5 alignments/cluster
- Largest: 424 alignments

**Perl PASA** (303 clusters):
- Singletons: 55 (18%)
- 2-5 alignments: 79 (26%)
- 6-10 alignments: 46 (15%)
- 11-50 alignments: 97 (32%)
- >50 alignments: 26 (9%)
- Average: 19.2 alignments/cluster (AFTER accounting for PyPASA's 214 extra valid alignments)
- Largest: 480 alignments

### Key Observations

1. **PyPASA has 2-3x more small clusters**: 248 clusters with 2-5 alignments vs Perl's 79
2. **Similar average cluster sizes**: 17.5 vs 19.2 alignments/cluster
3. **Both systems are strand-specific**: All clusters have consistent effective orientation
4. **PyPASA correctly uses transcribed orientation**: Prioritizes spliced_orient over aligned orient

The issue is NOT about:
- ❌ Coordinate normalization (both handle gi|XXX and gi|revXXX correctly)
- ❌ Strand-specificity (both cluster by effective orientation)
- ❌ Orientation detection (PyPASA correctly uses spliced orientation)

## Impact

**Validation (After aligned_length fix)**:
- Perl PASA: 12,251 valid alignments
- PyPASA: 12,465 valid alignments  
- Difference: 214 (1.7%) ✓ GOOD

**Clustering**:
- Perl PASA: 303 clusters (55 singletons, 248 multi-alignment)
- PyPASA: 713 clusters (106 singletons, 607 multi-alignment)
- Difference: 410 (135%) ❌ SIGNIFICANT

### Cluster Size Distribution

**Perl PASA**:
- Largest cluster: 480 alignments
- Average: 19.2 alignments/cluster
- More consolidated clusters

**PyPASA**:
- Largest cluster: 424 alignments  
- Average: 17.5 alignments/cluster
- More fragmented clusters

## Possible Causes

The clustering difference likely stems from one or more of these issues:

1. **Expanding window algorithm details**: PyPASA may be implementing the window expansion differently than Perl's Overlap_piler.pm
   
2. **Overlap detection sensitivity**: Small differences in how overlaps are detected (exact coordinate matching, off-by-one, etc.)

3. **Validation differences**: PyPASA has 214 more valid alignments (12,465 vs 12,251). These extra alignments may create fragmentation by falling between existing clusters.

4. **Cluster file interpretation**: The Perl cluster output file doesn't show scaffold names, making it difficult to verify exact cluster membership.

## Next Steps for Investigation

1. Compare Perl's `Overlap_piler.pm` implementation line-by-line with PyPASA's `_expanding_window_cluster()`
2. Check if the 214 extra valid alignments in PyPASA are causing cluster splits
3. Trace through a specific scaffold's clustering in both systems
4. Verify overlap detection logic handles edge cases identically

## References

- Perl PASA: `Overlap_piler.pm` - clustering algorithm
- PyPASA: `pypasa/alignment/clusterer.py` lines 195-268 - expanding window implementation
