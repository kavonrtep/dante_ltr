# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**dante_ltr** is a bioinformatics tool that identifies complete LTR (Long Terminal Repeat) retrotransposons from DANTE (Domain Based Annotation of Transposable Elements) output. It classifies elements into five ranks based on detected features: DLTP (complete), DLP, DLT, DL, D (domains only).

## Architecture

The project uses a **Python wrapper + R core** pattern:

- **Python** (`dante_ltr` script): Handles CLI, input validation, large genome chunking (splits FASTA/GFF3 into ≤100MB pieces), running the per-chunk R detector in a memory-gated process pool (single-threaded children, pool sized by `-c` and the memory budget from `--max_memory` / scheduler env / cgroup / `MemAvailable`), and coordinate remapping after R processing
- **R** (`utils/detect_putative_ltr.R` + `utils/ltr_utils.R`): Core genomic analysis — domain clustering, LTR detection, TSD/PBS identification, classification. Uses Bioconductor (GenomicRanges, Biostrings, rtracklayer)

### Entry Points

| Script | Language | Purpose |
|--------|----------|---------|
| `dante_ltr` | Python | Main tool — detects LTR retrotransposons from DANTE output |
| `dante_ltr_to_library` | Python | Creates non-redundant repeat libraries (calls extract_fasta.R → mmseq_clustering.R) |
| `dante_ltr_summary` | R | Generates HTML summary reports with embedded plots |
| `dante_reclassify` | R | Simplifies element classification |
| `clean_ltr.R` | R | Post-processing/validation of LTR annotations |

### Key Data Files

- `databases/lineage_domain_order.csv` — Defines per-lineage constraints (domain order, search offsets, LTR min length) for 37 lineages. This CSV is user-editable via `--te_constrains`
- `databases/feature_distances_model.RDS` — Pre-computed statistical model for validating domain distances
- `databases/tRNAscan-SE*.fasta` — tRNA BLAST databases for PBS detection

### Processing Pipeline

```
DANTE GFF3 + Reference FASTA
  → [dante_ltr] Python chunking/coordination
    → [detect_putative_ltr.R] domain clustering → LTR search → TSD/PBS detection → rank classification
      (chunks run concurrently in a memory-gated process pool; output is index-ordered so it stays identical to a serial run)
  → [coordinate recalculation if chunked]
  → Output: GFF3 + statistics CSV + FASTA sequences
```

## Commands

### Running Tests

```bash
# Run full test suite (defaults to 10 CPUs)
./tests.sh

# With specific CPU count
./tests.sh 4
```

Tests require conda environment `dante_ltr` to be active. Output goes to `tmp/` directory. There is no unit test framework — tests run the full pipeline on sample data in `test_data/`.

### Running the Tool

```bash
# Main LTR detection
./dante_ltr -g test_data/sample_DANTE.gff3 -s test_data/sample_genome.fasta -o output_prefix -c 10

# Create repeat library from results
./dante_ltr_to_library -g output_prefix.gff3 -s test_data/sample_genome.fasta -o library_dir

# Generate HTML summary
./dante_ltr_summary -g output_prefix.gff3 -o summary_prefix
```

### Environment Setup

```bash
conda create -n dante_ltr -c bioconda -c conda-forge -c petrnovak dante_ltr
conda activate dante_ltr
```

Dependencies are in `conda-deps.txt` (conda spec format, not pip — file is renamed from the conventional `requirements.txt` so GitHub's dependency-graph / Dependabot pip detector skips it). Key external tools: BLAST, bedtools, mmseqs2, mafft, tidehunter, parasail-python (used by `dante_ltr_refine`), pyfaidx.

## Code Conventions

- Version is centralized in `version.py` — imported by Python scripts
- Python scripts use only the standard library (argparse, subprocess, tempfile, etc.)
- R scripts use `optparse` for CLI and `suppressPackageStartupMessages()` for clean output
- GFF3 parsing in Python uses the `Gff3Feature` class (defined in `dante_ltr` script) — not a shared module
- All scripts are executable without extensions (except R utility scripts in `utils/`)
- No linting, formatting, or pre-commit hooks are configured
- No CI/CD pipeline — testing is manual via `tests.sh`


## Ignore files:

- directory `hermit` - this is not part of dante_ltr project
