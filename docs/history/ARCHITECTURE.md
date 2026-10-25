# PyPASA Architecture

> **Quick Reference for AI Assistants**: 
> - All scaffolds are ALWAYS distinct (no normalization, no grouping)
> - Default clustering is orientation-agnostic (`respect_strand=False`)
> - Perl PASA is ground truth for validation
> - See "Critical Design Rules" section below for details
> - See "Project Status & Roadmap" for current progress and next steps

## Project Status & Roadmap

### Current Status (as of Jan 24, 2026)

**🎯 Goal**: Port Perl PASA alignment assembly pipeline to Python, maintaining exact compatibility with original behavior.

**✅ Phase 1: Core Alignment Processing - COMPLETE**
- ✓ GFF3 alignment import (14,490 alignments from minimap2)
- ✓ Alignment validation (quality filtering)
  - Fixed `aligned_length` calculation bug (genome_length → transcript_length)
  - Achieves 99.58% agreement with Perl validation
- ✓ Alignment clustering (expanding window algorithm)
  - Fixed orientation handling (default orientation-agnostic)
  - Fixed scaffold normalization (now treats all scaffolds as distinct)
  - Achieves 99.3% agreement with Perl clustering (615 vs 611 clusters)
- ✓ SQLite storage with efficient indexing
- ✓ CLI interface (`pypasa` command)

**🔄 Phase 2: Assembly Integration - IN PROGRESS**
- ⚠️ Assembly wrapper exists but needs testing
- ⚠️ C++ pasa utility integration pending
- ❌ Subclustering not yet implemented (post-assembly alternative splice variant detection)

**📋 Phase 3: Advanced Features - NOT STARTED**
- ❌ Annotation comparison
- ❌ Annotation updates
- ❌ Alternative splicing analysis
- ❌ Web portal (deferred to much later)

### What Works Right Now

```bash
# You can run this today and get validated clusters:
cd PyPASA
python3 << 'EOF'
import sys
sys.path.insert(0, '.')
from pypasa.storage.datastore import DataStore
from pypasa.alignment.importer import GFF3AlignmentImporter
from pypasa.alignment.validator import AlignmentValidator
from pypasa.alignment.clusterer import AlignmentClusterer

# Import alignments
db = DataStore('test.sqlite')
importer = GFF3AlignmentImporter(db)
count = importer.import_gff3('path/to/alignments.gff3', aligner='minimap2')

# Validate
validator = AlignmentValidator(db, min_percent_id=0, min_percent_aligned=50)
stats = validator.validate_alignments()

# Cluster (orientation-agnostic, matching Perl default)
clusterer = AlignmentClusterer(db, respect_strand=False)
cluster_stats = clusterer.cluster_alignments()

print(f"Imported: {count}, Valid: {stats['passed']}, Clusters: {cluster_stats['total_clusters']}")
EOF
```

### Next Steps to Resume Work

**Immediate Priority: Assembly Testing**
1. Test assembly wrapper with real data:
   ```bash
   # Check if C++ pasa utility is accessible
   which pasa
   
   # Test on sample cluster
   cd PyPASA
   python3 -c "from pypasa.assembly.assembler import PASAAssembler; ..."
   ```

2. Verify assembly output matches Perl PASA exactly:
   - Compare assembled sequences (FASTA)
   - Compare coordinates (GFF3)
   - Check splice site detection

3. Implement subclustering (post-assembly):
   - Corresponds to Perl's `subcluster_builder.dbi`
   - ≥50% overlap requirement
   - Strict strand matching ('?' only with '?')
   - Alternative splice variant detection

**Testing Protocol**
- Reference database: `/tmp/sample_mydb_pasa.sqlite` (Perl PASA output)
- Sample data: `../PASApipeline/sample_data/`
- Validation targets:
  - Valid alignments: should match within 1%
  - Cluster counts: should match within 1%
  - Assembly outputs: should match exactly

**Known Issues & Decisions**
- ✓ Scaffold normalization disabled (gi|XXX and gi|revXXX separate)
- ✓ Default clustering is orientation-agnostic
- ⚠️ Assembly integration needs validation
- ⚠️ Subclustering not implemented yet

### How to Debug Issues

If you need to investigate alignment/clustering differences:

```python
# Compare specific cluster
import sqlite3
perl_conn = sqlite3.connect('/tmp/sample_mydb_pasa.sqlite')
pypasa_db = DataStore('PyPASA/test_output/minimap2_fixed.sqlite')

# Check Perl cluster composition
cursor = perl_conn.execute("""
    SELECT c.cluster_id, c.annotdb_asmbl_id, COUNT(*) 
    FROM clusters c 
    JOIN align_link al ON c.cluster_id = al.cluster_id
    WHERE c.cluster_id = ?
    GROUP BY c.cluster_id, c.annotdb_asmbl_id
""", (cluster_id,))

# Check PyPASA cluster composition  
cursor = pypasa_db.conn.execute("""
    SELECT cluster_id, scaffold, COUNT(*)
    FROM alignments
    WHERE cluster_id = ?
    GROUP BY cluster_id, scaffold
""", (cluster_id,))
```

### Key Files for Future Work

**Assembly Implementation**:
- `pypasa/assembly/assembler.py` - Wrapper for C++ pasa utility (needs testing)
- External: Original Perl PASA's C++ code (`pasa-plugins/cdbtools/`)

**Subclustering** (needs implementation):
- Create `pypasa/alignment/subclusterer.py`
- Reference: Perl's `scripts/subcluster_builder.dbi`
- Logic: Group transcripts within clusters by ≥50% overlap + strict strand matching

**Output Generation** (mostly done):
- `pypasa/alignment/exporter.py` - GFF3/GTF/BED output (exists)
- FASTA generation from assemblies (needs assembly working first)

### Performance Benchmarks (Current)

From testing on sample data (14,490 alignments):
- Import: ~28,000 alignments/sec
- Validation: ~18,000 alignments/sec  
- Clustering: ~1,500 alignments/sec
- Total pipeline: <10 seconds for sample data

### References

- **Perl PASA codebase**: `../PASApipeline/`
- **Key Perl scripts**:
  - `scripts/reassign_clusters_via_valid_align_coords.dbi` (default clustering)
  - `scripts/assign_clusters_by_stringent_alignment_overlap.dbi` (stringent mode)
  - `scripts/subcluster_builder.dbi` (subclustering)
  - `scripts/Launch_PASA_pipeline.pl` (main pipeline)
- **C++ assembly code**: `pasa-plugins/cdbtools/`

## Design Principles

1. **Keep What Works**: Reuse proven C++ assembly algorithm
2. **Hybrid Storage**: RAM for speed, SQLite for persistence
3. **Modular Design**: Clear separation of concerns
4. **Incremental Migration**: Port Perl → Python gradually
5. **Backward Compatible**: Match original PASA outputs

## Critical Design Rules

### Scaffold Handling
**IMPORTANT**: All scaffold names are treated as completely distinct identifiers at all times.

- **NO normalization**: `gi|68711` and `gi|rev68711` are DIFFERENT scaffolds
- **NO grouping**: Similar names are never merged or grouped together
- **NO special handling**: Forward/reverse complement naming is irrelevant
- **Clustering**: Each scaffold clusters independently
- **Assembly**: Each scaffold assembles independently

**Rationale**: After extensive testing against Perl PASA (Jan 2026), we discovered that:
1. Perl treats each scaffold as completely distinct during clustering
2. Early PyPASA versions incorrectly grouped `gi|XXX` with `gi|revXXX`
3. This caused PyPASA to produce 361 clusters vs Perl's 611 (40% fewer)
4. Removing normalization fixed this: 615 vs 611 clusters (0.7% difference ✓)

**Code locations**:
- `pypasa/alignment/clusterer.py`: `_cluster_scaffold()` processes each scaffold separately
- `pypasa/alignment/best_selector.py`: No scaffold normalization (removed)
- `pypasa/storage/datastore.py`: Queries by exact scaffold name match

### Clustering Behavior
**Default mode** (matches Perl's `reassign_clusters_via_valid_align_coords.dbi`):
- **Orientation-agnostic**: Clusters by coordinate overlap only
- **respect_strand=False**: Allows mixed-orientation clusters
- Single-exon alignments (spliced_orient='?') cluster with any overlapping alignment
- This is first-pass clustering before assembly

**Stringent mode** (matches Perl's `assign_clusters_by_stringent_alignment_overlap.dbi`):
- **Orientation-aware**: Groups by spliced_orient first
- **respect_strand=True**: Separate clusters per orientation
- Optional, not used by default

**Subclustering** (post-assembly, not yet implemented):
- **Strict strand matching**: '?' only matches '?'
- ≥50% overlap requirement
- Identifies alternative splice variants
- Corresponds to Perl's `subcluster_builder.dbi`

## Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│                     ALIGNMENT ASSEMBLY PIPELINE               │
└─────────────────────────────────────────────────────────────┘

INPUT: genome.fasta + transcripts.fasta

    │
    ├─► External Aligners (GMAP/BLAT)
    │   └─► alignments.gff3
    │
    ├─► AlignmentImporter (PyPASA)
    │   ├─► Parse GFF3 format
    │   ├─► Create Alignment objects
    │   └─► Store in SQLite + load to RAM
    │
    ├─► AlignmentValidator (PyPASA)
    │   ├─► Check percent identity
    │   ├─► Check alignment length
    │   ├─► Filter low-quality
    │   └─► Mark valid=True/False
    │
    ├─► AlignmentClusterer (PyPASA)
    │   ├─► Build interval trees per scaffold
    │   ├─► Find overlapping alignments
    │   ├─► Assign cluster_id
    │   └─► Update in RAM + persist
    │
    ├─► ClusterAssembler (PyPASA wrapper)
    │   ├─► Export cluster → temp file
    │   ├─► Call C++ pasa utility
    │   ├─► Parse assembly output
    │   └─► Store assemblies
    │
    └─► OutputGenerator (PyPASA)
        ├─► Generate FASTA sequences
        ├─► Generate GFF3 structures
        ├─► Generate GTF/BED formats
        └─► Summary statistics

OUTPUT: assemblies.fasta, assemblies.gff3, etc.
```

## Storage Strategy

### SQLite Schema (Minimal)
```sql
-- Alignments table
CREATE TABLE alignments (
    align_id INTEGER PRIMARY KEY,
    cdna_acc TEXT NOT NULL,
    scaffold TEXT NOT NULL,
    lend INTEGER NOT NULL,
    rend INTEGER NOT NULL,
    orient TEXT NOT NULL,
    percent_id REAL,
    valid INTEGER DEFAULT 0,
    cluster_id INTEGER,
    alignment_data BLOB  -- JSON or pickled Alignment object
);
CREATE INDEX idx_scaffold_coords ON alignments(scaffold, lend, rend);
CREATE INDEX idx_cluster ON alignments(cluster_id);

-- Assemblies table
CREATE TABLE assemblies (
    assembly_id INTEGER PRIMARY KEY,
    cluster_id INTEGER NOT NULL,
    scaffold TEXT NOT NULL,
    lend INTEGER NOT NULL,
    rend INTEGER NOT NULL,
    orient TEXT NOT NULL,
    sequence TEXT,
    assembly_data BLOB  -- JSON or pickled Assembly object
);

-- Cluster membership
CREATE TABLE cluster_members (
    cluster_id INTEGER,
    align_id INTEGER,
    PRIMARY KEY (cluster_id, align_id)
);
```

### In-Memory Structures
```python
# Loaded per-scaffold for processing
class ScaffoldData:
    scaffold_id: str
    alignments: List[Alignment]  # Active alignments
    interval_tree: IntervalTree  # For overlap queries
    clusters: Dict[int, List[int]]  # cluster_id → align_ids
```

### Processing Pattern
```python
# Process one scaffold at a time
for scaffold in genome.scaffolds:
    # Load from SQLite
    data = storage.load_scaffold(scaffold.id)
    
    # Process in RAM
    clusterer.cluster_alignments(data)
    assembler.assemble_clusters(data)
    
    # Persist results
    storage.save_scaffold(data)
    
    # Free memory
    del data
```

## Module Responsibilities

### alignment/importer.py
- Parse alignment formats (GFF3, BTAB)
- Convert to Alignment objects
- Bulk insert to database

### alignment/validator.py
- Apply quality filters
- Mark valid/invalid alignments
- Track validation statistics

### alignment/clusterer.py
- Build spatial indexes (IntervalTree)
- Identify overlapping alignments
- Assign cluster IDs
- Support multiple clustering strategies:
  - Default: any overlap
  - Stringent: minimum % overlap
  - Gene-overlap: annotation-informed

### assembly/assembler.py
- Wrapper for C++ pasa utility
- Manage temp files for clusters
- Parse assembly output
- Store assemblies

### storage/datastore.py
- SQLite connection management
- Bulk load/save operations
- Query interface
- Checkpoint/resume support

### storage/models.py
- Data classes for Alignment, Assembly
- Serialization/deserialization
- Validation

## Performance Considerations

### Memory Management
- Process scaffolds sequentially (not all at once)
- For large genomes: stream alignments in batches
- Use generators where possible

### Speed Optimizations
- SQLite in WAL mode for concurrent reads
- Bulk inserts (1000s at a time)
- Interval trees for O(log n) overlap queries
- NumPy for coordinate arithmetic

### Scalability
- Scaffold-level parallelism (future)
- Chunked processing for huge alignments
- Optional memory-mapped file support

## Error Handling & Recovery

### Checkpointing
```python
# Save state at each major step
checkpoint.save(stage='import_complete')
checkpoint.save(stage='validation_complete')
checkpoint.save(stage='clustering_complete')

# Resume from checkpoint
if checkpoint.exists('clustering_complete'):
    storage = checkpoint.load('clustering_complete')
    # Skip to next step
```

### Validation
- Compare outputs with Perl PASA (exact match expected)
- Unit tests for each module
- Integration tests with sample_data/

## Migration Notes

### What We Keep
- C++ assembly algorithm (proven, fast)
- External aligners (GMAP, BLAT)
- Output formats (GFF3, GTF, BED, FASTA)

### What We Change
- Perl → Python
- MySQL → SQLite (simpler setup)
- DB-centric → RAM-centric processing
- Monolithic scripts → modular library

### What We Defer
- Annotation comparison (Phase 2)
- Web portal (Phase 3)
- Alternative splicing analysis (Phase 2)

## Validation & Testing Against Perl PASA

### Key Investigations (Jan 2026)

#### Alignment Validation
**Issue**: PyPASA had more valid alignments than Perl (12,465 vs 12,251)
**Root cause**: `aligned_length` calculation bug - used `genome_length` instead of `transcript_length`
**Fix**: Changed to `aln_end - aln_start + 1` for transcript coordinates
**Result**: 99.58% agreement with Perl validation

#### Clustering Algorithm
**Issue**: PyPASA produced different cluster counts than Perl
- Initial: 713 clusters vs Perl's 303 (2.35x more) with strand-specific mode
- After orientation fix: 358 clusters vs Perl's 303 (18% more)
- After scaffold fix: 615 vs 611 clusters (0.7% difference ✓)

**Root causes**:
1. **Strand handling**: PyPASA defaulted to `respect_strand=True`, but Perl's default is orientation-agnostic
2. **Scaffold normalization**: PyPASA grouped `gi|XXX` with `gi|revXXX`, but Perl treats them as separate

**Validation approach**:
1. Extracted Perl's valid alignment set from `/tmp/sample_mydb_pasa.sqlite`
2. Applied same set to PyPASA database
3. Ran clustering with `respect_strand=False`
4. Compared cluster counts: 615 vs 611 (excellent agreement)

**Test database**: Use `/tmp/sample_mydb_pasa.sqlite` (Perl output) for validation

### Testing Checklist
When modifying alignment or clustering code:
- [ ] Run against sample data with both PyPASA and Perl
- [ ] Compare valid alignment counts (should match within 1%)
- [ ] Compare cluster counts (should match within 1%)
- [ ] Check that gi|XXX and gi|revXXX cluster separately
- [ ] Verify orientation-agnostic clustering is default
- [ ] Test that mixed-orientation clusters are created

### Known Differences (Acceptable)
- PyPASA may have slightly more valid alignments due to fixed aligned_length calculation
- Cluster count difference <1% is acceptable (edge cases in overlap detection)
- Assembly outputs should match exactly (byte-for-byte)
