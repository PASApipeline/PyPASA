# PyPASA Architecture

## Design Principles

1. **Keep What Works**: Reuse proven C++ assembly algorithm
2. **Hybrid Storage**: RAM for speed, SQLite for persistence
3. **Modular Design**: Clear separation of concerns
4. **Incremental Migration**: Port Perl → Python gradually
5. **Backward Compatible**: Match original PASA outputs

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
