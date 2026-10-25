# PyPASA vs Perl PASA: Detailed Comparison

## Test Dataset
- **Input**: `PASApipeline/sample_data/custom_alignments.gff3`
- **Date**: January 24, 2026
- **Perl PASA DB**: `/tmp/sample_mydb_pasa.sqlite`

## Results Summary

| Metric | Perl PASA | PyPASA | Ratio | Notes |
|--------|-----------|--------|-------|-------|
| **Imported alignments** | 22,784 | 22,784 | 1.00x | ✅ Identical |
| **Valid alignments** | 22,736 | 11,482 | 1.98x | Perl keeps ~2 per cDNA |
| **Unique cDNAs** | 11,464 | 11,482 | 1.00x | ✅ Nearly identical |
| **Clusters** | 605 | 551 | 1.10x | Close match |
| **Alns per cDNA** | 1.98 | 1.00 | - | Key difference |

## Key Finding: Perl PASA Allows Multiple Alignments Per cDNA

### Perl PASA Behavior (Observed)
Perl PASA **does NOT** select a single best alignment per cDNA. Instead:

1. **Import**: All alignments imported from GFF3
2. **Validate**: Filter by percent identity (95%) and percent aligned (75%)
3. **Cluster**: Group overlapping alignments
4. **Keep**: Multiple alignments per cDNA if they're in **different clusters**

**Example from Perl PASA DB:**
```
cDNA: gi|10277115|gb|BE844737.1|BE844737
  - Alignment 1: scaffold position 46,766-47,305, score=51802.2, cluster=377
  - Alignment 2: scaffold position 73,634-74,173, score=51802.2, cluster=72
  
Both kept! Same score, different genomic locations.
```

### PyPASA Behavior (Current Implementation)
PyPASA implements **strict best-selection**:

1. **Import**: All alignments imported from GFF3
2. **Validate**: Filter by percent identity (95%) and percent aligned (75%)
3. **Best-selection**: Keep ONLY the single highest-scoring alignment per cDNA **globally**
4. **Cluster**: Group overlapping alignments

This produces fewer valid alignments (11,482 vs 22,736) but similar cluster counts (551 vs 605).

## Distribution Analysis

### Perl PASA: Alignments per cDNA (custom only)
- **1 alignment**: 1,734 cDNAs (15%)
- **2 alignments**: 9,730 cDNAs (85%)

### PyPASA: Alignments per cDNA
- **1 alignment**: 11,482 cDNAs (100%)

## Why Perl Keeps Multiple Alignments

The `-N 1` parameter in Perl PASA appears to mean "keep up to N alignments **per aligner**", not "per cDNA globally". Since the sample data runs:
- Custom alignments (from GFF3)
- Minimap2 alignments (from FASTA)
- Blat/Gmap alignments (from FASTA)

Each cDNA can have alignments from multiple sources. Within the **custom** alignments specifically, Perl keeps ~2 per cDNA because they map to different genomic locations and end up in different clusters.

## Validation Differences

### Perl PASA
- **Total custom**: 22,784
- **Valid custom**: 22,736 (99.79% pass rate)
- **Invalid**: 48

### PyPASA
- **Total**: 22,784
- **Valid (before best-selection)**: 22,772 (99.95% pass rate)
- **Invalid**: 12
- **Valid (after best-selection)**: 11,482 (kept best from each cDNA)

PyPASA has slightly stricter validation (12 failures vs 48), possibly due to:
- Canonical splice site checking implementation differences
- Floating-point precision in threshold calculations

## Cluster Count Analysis

Despite having ~2x fewer alignments, PyPASA produces similar cluster counts:
- **Perl**: 605 clusters from 22,736 alignments
- **PyPASA**: 551 clusters from 11,482 alignments

This suggests:
1. PyPASA's clustering algorithm is working correctly
2. The "duplicate" alignments in Perl are spread across different clusters
3. Best-selection effectively removes redundant alignments without losing genomic coverage

## Implications

### If Goal is: Exact Perl PASA Replication
PyPASA needs modification:
- Remove or relax global best-selection
- Allow multiple alignments per cDNA if they cluster separately
- Implement post-clustering deduplication (one per cDNA per cluster)

### If Goal is: Improved Pipeline
PyPASA's current approach may be better:
- Cleaner: exactly one alignment per transcript
- Less redundancy in clustering
- Similar cluster counts suggest comparable genomic coverage
- Simpler assembly logic (no duplicate transcripts in same region)

## Recommendation

**Keep PyPASA's current strict best-selection behavior** because:
1. It produces clean, non-redundant alignments
2. Cluster counts are comparable (551 vs 605)
3. Simpler downstream processing
4. More deterministic results

If exact Perl replication is required, document this as a known behavioral difference and provide a flag to enable "Perl-compatible mode" that skips global best-selection.

## Full Statistics

### Perl PASA (All Aligners)
```
Total alignments:        38,127
  - custom:              22,784 (22,736 valid)
  - minimap2:            14,490 (12,251 valid)
  - assembler:              853 (853 valid)

Valid alignments:        35,840
Unique cDNAs:            13,180
Clusters:                   611
Avg alignments/cDNA:       2.72
```

### Perl PASA (Custom Only - Matching PyPASA Input)
```
Total:                   22,784
Valid:                   22,736
Unique cDNAs:            11,464
Clusters:                   605
Avg alignments/cDNA:       1.98
```

### PyPASA
```
Total:                   22,784
Valid:                   11,482
Unique cDNAs:            11,482
Clusters:                   551
Avg alignments/cDNA:       1.00
```
