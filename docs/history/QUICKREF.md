# PyPASA - Quick Reference

> **⚠️ Critical Notes for Developers/AI**:
> - **Scaffolds**: NEVER normalize, group, or merge scaffold names. `gi|68711` ≠ `gi|rev68711`
> - **Clustering**: Default is orientation-agnostic (`respect_strand=False`)
> - **Validation**: Always test against Perl PASA output in `/tmp/sample_mydb_pasa.sqlite`
> - See [ARCHITECTURE.md](ARCHITECTURE.md#critical-design-rules) for details

## Installation

```bash
cd /home/unix/bhaas/projects/PASA_port/PyPASA
pip install -e .
```

## Running PyPASA

### Full Pipeline (Command Line)
```bash
pypasa align-assemble \
    --alignments path/to/alignments.gff3 \
    --database output.sqlite \
    --genome path/to/genome.fasta \
    --min-percent-id 95 \
    --min-percent-aligned 75 \
    --stringent-clustering
```

### Check Statistics
```bash
pypasa stats --database output.sqlite
```

### Test on Sample Data
```bash
cd /home/unix/bhaas/projects/PASA_port/PyPASA
python scripts/test_sample_data.py
```

### Compare with Perl PASA
```bash
cd /home/unix/bhaas/projects/PASA_port/PyPASA
python scripts/compare_with_perl.py
```

## File Locations

### Source Code
- **Models**: `pypasa/storage/models.py`
- **Storage**: `pypasa/storage/datastore.py`
- **Import**: `pypasa/alignment/importer.py`
- **Validation**: `pypasa/alignment/validator.py`
- **Clustering**: `pypasa/alignment/clusterer.py`
- **Assembly**: `pypasa/assembly/assembler.py`
- **CLI**: `pypasa/cli.py`

### Test Data
- **Sample GFF3**: `../PASApipeline/sample_data/custom_alignments.gff3`
- **Genome**: `../PASApipeline/sample_data/genome_sample.fasta.gz`
- **Perl Output**: `../PASApipeline/sample_data/clusters_of_valid_alignments.txt.gz`

### Outputs
- **Test Database**: `test_output/pypasa_test.sqlite`
- **Test Results**: `test_output/`

### Documentation
- **README**: [README.md](README.md) - Getting started
- **Architecture**: [ARCHITECTURE.md](ARCHITECTURE.md) - Design details
- **Progress**: [PROGRESS.md](PROGRESS.md) - Development tracking
- **Comparison**: [COMPARISON_RESULTS.md](COMPARISON_RESULTS.md) - Perl vs Python
- **Status**: [STATUS.md](STATUS.md) - Current status report

## Current Status (Phase 1 Complete)

✅ **Working**:
- GFF3 import (28,000 alignments/sec)
- Quality validation (99.98% pass rate)
- Overlap clustering (1,500 alignments/sec)
- SQLite storage with WAL mode
- CLI interface
- Test infrastructure

⚠️ **Placeholder**:
- Assembly step (needs C++ pasa integration)

❌ **Not Yet Implemented**:
- GFF3/GTF/BED output writers
- FASTA sequence extraction
- Unit test suite
- FL-cDNA marking
- PolyA detection

## Key Findings

### Performance
- **PyPASA**: ~17 seconds on sample data (22,784 alignments)
- **Perl PASA**: ~30 seconds on same data
- **Result**: Python version is faster! ✅

### Clustering Differences
- **Perl PASA**: 303 clusters (5,811 alignments)
- **PyPASA**: 563 clusters (11,486 alignments)
- **Reason**: PyPASA keeps all alignments per transcript, Perl keeps only best

### Correctness
- ✅ GFF3 parsing correct
- ✅ Multi-exon alignment merging correct
- ✅ Validation filters working
- ✅ Overlap detection accurate
- ✅ No over-merging of clusters
- ⚠️ Different alignment filtering strategy than Perl

## Next Steps

1. **Decide alignment strategy**:
   - Match Perl (filter to best alignment)?
   - Keep current behavior (all alignments)?
   - Make it configurable? ← Recommended

2. **Integrate C++ pasa**:
   - Export cluster data in correct format
   - Call `bin/pasa` binary
   - Parse assembly results

3. **Add output writers**:
   - GFF3 for assemblies
   - GTF for assemblies
   - BED for visualization
   - FASTA sequences

## Configuration Options

### Validation Parameters
- `--min-percent-id`: Minimum alignment identity (default: 95%)
- `--min-percent-aligned`: Minimum aligned length (default: 75%)

### Clustering Parameters
- `--stringent-clustering`: Use stringent mode (50% overlap required)
- `--no-strand-specific`: Allow cross-strand clustering

### Performance Parameters
- Database cache: 64MB (hardcoded)
- WAL mode: Enabled by default
- Bulk operations: 1000 records per batch

## API Usage (Python)

```python
from pypasa.storage.datastore import DataStore
from pypasa.alignment.importer import import_gff3
from pypasa.alignment.validator import validate_alignments
from pypasa.alignment.clusterer import cluster_alignments

# Initialize database
datastore = DataStore("output.sqlite")

# Import alignments
import_gff3("alignments.gff3", datastore)

# Validate
validate_alignments(datastore, min_percent_id=95, min_percent_aligned=75)

# Cluster
cluster_alignments(datastore, stringent_alignment_overlap=False)

# Get statistics
stats = datastore.get_stats()
print(f"Total: {stats['total_alignments']}")
print(f"Valid: {stats['valid_alignments']}")
print(f"Clusters: {stats['num_clusters']}")

datastore.close()
```

## Common Tasks

### Query Database Directly

```python
import sqlite3

conn = sqlite3.connect("output.sqlite")
cursor = conn.cursor()

# Count alignments per scaffold
cursor.execute("""
    SELECT scaffold, COUNT(*) 
    FROM alignments 
    GROUP BY scaffold
""")
for scaffold, count in cursor.fetchall():
    print(f"{scaffold}: {count}")

conn.close()
```

### Export Clusters

```python
from pypasa.storage.datastore import DataStore

datastore = DataStore("output.sqlite")

# Get all clusters
for scaffold in datastore.get_scaffolds():
    alignments = datastore.get_alignments_by_scaffold(scaffold)
    clusters = {}
    for aln in alignments:
        if aln.cluster_id:
            clusters.setdefault(aln.cluster_id, []).append(aln)
    
    for cluster_id, alns in clusters.items():
        print(f"Cluster {cluster_id}: {len(alns)} alignments")
        for aln in alns:
            print(f"  - {aln.cdna_acc} at {aln.lend}-{aln.rend}")

datastore.close()
```

## Troubleshooting

### "No module named pypasa"
```bash
cd /home/unix/bhaas/projects/PASA_port/PyPASA
pip install -e .
```

### "Database is locked"
SQLite WAL mode should prevent this, but if it happens:
```python
# Close any open connections
datastore.close()

# Or check for other processes:
lsof output.sqlite
```

### Clustering too slow
For large datasets:
1. Process one scaffold at a time (already done)
2. Use stringent clustering (faster, fewer comparisons)
3. Consider parallel processing (future enhancement)

### Memory usage too high
The scaffold-based processing should keep memory manageable. If issues occur:
1. Check scaffold sizes (split large scaffolds)
2. Reduce database cache size
3. Process smaller batches

## Contact & Support

- **Project**: PyPASA
- **Location**: `/home/unix/bhaas/projects/PASA_port/PyPASA/`
- **Original**: PASA Pipeline (Perl version)
- **Status**: Phase 1 Complete, Phase 2 In Progress

## Quick Commands Cheat Sheet

```bash
# Install
pip install -e .

# Test
python scripts/test_sample_data.py

# Compare
python scripts/compare_with_perl.py

# Run pipeline
pypasa align-assemble --alignments INPUT.gff3 --database OUTPUT.sqlite --genome GENOME.fasta

# Stats
pypasa stats --database OUTPUT.sqlite

# Help
pypasa --help
pypasa align-assemble --help
```
