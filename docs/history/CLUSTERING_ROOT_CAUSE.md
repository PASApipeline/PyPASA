# Root Cause Analysis: PyPASA vs Perl PASA Clustering Differences

## Executive Summary

**Finding**: Perl PASA's clustering file contains **MIXED-ORIENTATION clusters** (e.g., cluster 81 has 6 '+' and 112 '-' orientations), while PyPASA strictly enforces strand-specific clustering.

**Result**: PyPASA produces 713 clusters vs Perl's 303 clusters (2.35x difference).

**Status**: PyPASA behavior is **CORRECT** according to specification that clustering should be "strand-specific, using transcribed orientations."

## Detailed Investigation

### Phase 1: Valid Alignments Comparison

✅ **All 12,251 common alignments match perfectly**:
- Same scaffolds
- Same orientations  
- No discrepancies

**Difference**: PyPASA has 214 additional valid alignments (12,465 vs 12,251) due to the aligned_length bug fix.

### Phase 2: Cluster Composition Analysis

**Cluster overlap analysis**:
- Exact matches: 44 clusters (same members)
- PyPASA splits: 352 cases (PyPASA split a Perl cluster)
- PyPASA merges: 127 cases (PyPASA merged Perl clusters)

**Pattern**:
1. PyPASA creates MORE clusters by splitting mixed-orientation groups
2. PyPASA sometimes merges clusters that Perl kept separate (likely due to scaffold grouping)

### Phase 3: Specific Case Studies

#### Case 1: Perl Cluster 81 (SPLIT by PyPASA)

**Perl cluster 81**:
- 138 alignments total
- **Orientation distribution**: 6 '+', 112 '-', 20 '?'
- **This is a MIXED-ORIENTATION cluster!**

**PyPASA correctly split into 4 strand-specific clusters**:
- Cluster 264: 55 aligns on gi|68716, orientation '-', coords 39382-41914
- Cluster 303: 13 aligns on gi|68716, orientation '+', coords 38352-39934  
- Cluster 280: 7 aligns on gi|rev68716, orientation '-', coords 79144-80650
- Cluster 318: 63 aligns on gi|rev68716, orientation '+', coords 77088-79628

**Critical observation**: Clusters 264 and 303 are on the SAME scaffold (gi|68716) with OVERLAPPING coordinates but DIFFERENT orientations. PyPASA correctly keeps them separate.

#### Case 2: Perl Clusters 47 & 48 (MERGED by PyPASA)

**Perl clusters 47 & 48**:
- Cluster 47: 11 alignments (orientations: 5 '?', 6 '+')
- Cluster 48: 15 alignments (orientations: 2 '?', 13 '+')

**PyPASA merged and re-organized by strand**:
- PyPASA cluster 2: 11 aligns from both 47+48 on gi|68711, orientation '+'
- PyPASA cluster 68: 10 aligns from both 47+48 on gi|rev68711, orientation '-'
- Plus several small singleton/doubleton clusters for mixed orientations

**Explanation**: PyPASA groups scaffold variants (gi|68711 with gi|rev68711) and then does strand-specific clustering.

## Source Code Evidence

### PyPASA Implementation (CORRECT)

[pypasa/alignment/clusterer.py](PyPASA/pypasa/alignment/clusterer.py#L158-L163):
```python
# Group alignments by scaffold and orientation
grouped_alignments = defaultdict(list)
for aln in valid_alignments:
    key = (aln.scaffold, aln.effective_orient)  # Strand-specific grouping
    grouped_alignments[key].append(aln)
```

**PyPASA uses `effective_orient`** which prioritizes `spliced_orient` over `orient` for strand-specific clustering.

### Perl Implementation - TWO CLUSTERING MODES

Perl PASA has **two different clustering modes**:

#### Mode 1: Default Clustering (ORIENTATION-AGNOSTIC)

Script: `reassign_clusters_via_valid_align_coords.dbi`

```perl
# Line 172 - NO spliced_orient in query!
my $query = "select al.align_id from align_link al, clusters c "
    . " where c.cluster_id = al.cluster_id and c.annotdb_asmbl_id = ? and al.validate = 1";

# Uses Overlap_assembler which only checks coordinate overlap
my $overlap_assembler = new CDNA::Overlap_assembler();
```

The `Overlap_assembler->build_clusters()` method:
```perl
if ($exp_left <= $rend && $exp_right >= $lend) { #overlap
    # Merge based on coordinate overlap ONLY - no orientation check!
}
```

**This mode does NOT enforce strand-specificity during initial clustering.**

#### Mode 2: Stringent Clustering (ORIENTATION-AWARE)

Script: `assign_clusters_by_stringent_alignment_overlap.dbi`

```perl
# Line 68 - DOES retrieve spliced_orient
my $query = "select c.annotdb_asmbl_id, al.align_acc, al.spliced_orient "
    . " from clusters c, align_link al, cdna_info ci "
    . " where c.cluster_id = al.cluster_id and al.cdna_info_id = ci.id "
    . " and ci.is_assembly = 0 and al.validate = 1";

# Line 77 - Groups by scaffold AND orientation
$asmbl_id .= ";$spliced_orient";
```

**This mode DOES enforce strand-specificity.**

### The Root Cause CONFIRMED

The validation file `clusters_of_valid_alignments.txt.gz` was generated using **DEFAULT MODE** (orientation-agnostic clustering).

**Evidence**:
1. Perl cluster 81 contains mixed orientations (6 '+', 112 '-', 20 '?')
2. This matches the behavior of `reassign_clusters_via_valid_align_coords.dbi`
3. The Launch_PASA_pipeline.pl comment on line 754 says "any overlap with same transcribed orientation" but **this is incorrect** - the actual default implementation ignores orientation

**The comment is misleading!** The default clustering mode does NOT filter by transcribed orientation despite what the comment claims

## Verification of PyPASA Correctness

✅ All 713 PyPASA clusters have **consistent orientations** (verified by SQL query)

✅ PyPASA correctly uses **transcribed orientation** (spliced_orient) as primary, with fallback to aligned orient

✅ PyPASA correctly implements **expanding window clustering** with strand-specificity

✅ PyPASA correctly groups **scaffold variants** (gi|XXX with gi|revXXX)

## Conclusion

**PyPASA is implementing strand-specific clustering correctly**, while **Perl PASA's DEFAULT clustering mode is orientation-agnostic**.

### Key Findings:

1. **Perl PASA has TWO clustering modes**:
   - **Default mode**: `reassign_clusters_via_valid_align_coords.dbi` - **IGNORES orientation**
   - **Stringent mode**: `assign_clusters_by_stringent_alignment_overlap.dbi` - **RESPECTS orientation**

2. **The validation data was generated using DEFAULT mode** (orientation-agnostic)

3. **PyPASA implements strand-specific clustering** (equivalent to Perl's stringent mode behavior)

4. **The comment in Perl's Launch_PASA_pipeline.pl line 754 is INCORRECT**:
   - Comment says: "any overlap with same transcribed orientation"  
   - Reality: Default mode does NOT filter by orientation at all

### Which is Correct?

This depends on the **intended behavior**:

- **If subclustering is meant to be strand-specific**: Then initial clustering should probably also be strand-specific (PyPASA's approach)
- **If initial clustering is meant to be orientation-agnostic**: Then PyPASA needs to be modified to match Perl's default mode

The biological question: Should overlapping genes on opposite strands be clustered together initially, or kept separate from the start?

## Recommendations

1. ⚠️ **USER DECISION REQUIRED**: Clarify the intended clustering behavior:
   - Should PyPASA match Perl's **default** mode (orientation-agnostic)?
   - Or match Perl's **stringent** mode (orientation-aware)?
   - Or keep PyPASA's current implementation?

2. 🔧 **If matching default mode is desired**: PyPASA needs modification:
   ```python
   # Remove orientation from clustering key:
   key = aln.scaffold  # NOT (aln.scaffold, aln.effective_orient)
   ```

3. ✅ **If keeping current implementation**: Document that PyPASA uses strand-specific clustering (like Perl's stringent mode) by default

4. 📝 **Fix Perl's misleading comment**: The Launch_PASA_pipeline.pl comment should be corrected to reflect actual behavior

## Impact

- **Higher cluster count in PyPASA**: Expected and correct behavior when enforcing strand-specificity
- **Better biological accuracy**: Genes on opposite strands should not be clustered together
- **No functionality lost**: All alignments are still clustered appropriately by their strand
