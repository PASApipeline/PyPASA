# PyPASA - Python Implementation of PASA Alignment Assembly Pipeline

## Overview

PyPASA is a Python port of the PASA (Program to Assemble Spliced Alignments) alignment assembly pipeline. This implementation focuses on the core alignment assembly functionality while maintaining compatibility with the original PASA workflow.

## Architecture

### Hybrid Processing Model
- **In-Memory Processing**: Fast operations on active data using NumPy/Pandas
- **Persistent Storage**: SQLite for state management and queryability
- **Native Integration**: Calls existing C++ assembly utility for proven algorithm

### Project Structure
```
PyPASA/
├── pypasa/
│   ├── alignment/        # Alignment import, validation, clustering
│   ├── assembly/         # Assembly pipeline and C++ wrapper
│   ├── storage/          # Hybrid storage manager (RAM + SQLite)
│   └── utils/            # Shared utilities
├── tests/                # Test suite
└── scripts/              # CLI tools
```

## Features

### Phase 1: Alignment Assembly (Current Focus)
- [x] Alignment import (GMAP/BLAT GFF3 format)
- [x] Alignment validation (quality filters)
- [x] Overlap-based clustering
- [x] Assembly via C++ utility
- [x] Output generation (FASTA, GFF3, GTF, BED)

### Future Phases
- [ ] Annotation comparison
- [ ] Annotation updates
- [ ] Alternative splicing analysis

## Installation

```bash
cd PyPASA
pip install -e .
```

## Usage

```bash
# Run alignment assembly pipeline
pypasa align-assemble \
  --genome genome.fasta \
  --transcripts transcripts.fasta \
  --output-dir pasa_out \
  --aligners gmap,blat \
  --cpu 8
```

## Performance

- **Memory**: Processes genome scaffolds independently
- **Storage**: SQLite for persistence, RAM for processing
- **Speed**: Comparable to original Perl implementation

## Compatibility

- Reuses existing C++ assembly utility
- Compatible with PASA sample data
- Output formats match original PASA

## Development Status

This is an active port from the original Perl implementation. See PROGRESS.md for migration status.
