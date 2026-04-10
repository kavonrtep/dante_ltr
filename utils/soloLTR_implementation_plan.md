# Solo LTR Detection — Implementation Plan

Reference design document: `soloLTR.md`

---

## Overview of new files

| File | Language | Role |
|------|----------|------|
| `dante_ltr_solo` | Python | Entry point: argument parsing, genome chunking, dispatch, output merging |
| `utils/build_ltr_library.R` | R | Step 1: extract LTRs, cluster, MAFFT consensus, 5'UTR tags, TSD map |
| `utils/detect_solo_ltr.R` | R | Steps 2–4: BLAST search, overlap filter, TSD check, 5'UTR junction check, GFF3 output |
| `utils/solo_ltr_utils.R` | R | Shared R helper functions sourced by both R scripts above |

Modified files: `requirements.txt`, `tests.sh`, `version.py`.

The architecture mirrors `dante_ltr` exactly: Python wrapper handles chunking and merging; R scripts do the biology per chunk.

---

## Phase 1 — Shared R utilities (`utils/solo_ltr_utils.R`)

Functions here are sourced by both `build_ltr_library.R` and `detect_solo_ltr.R`.

### 1.1 — TSD length map from GFF3

```
extract_tsd_length_map(gff) → named integer vector [lineage → modal_TSD_length]
```

- Filter input GFF to `transposable_element` features with `Rank` in `{DLT, DLTP}`
- For each, read the `TSD` attribute string (e.g. `TTGCA`); compute `nchar()`
- Group by `Final_Classification`; take modal length per lineage
- Return named vector; lineages with < 5 observations get `NA` (fallback: scan 4–6 bp)

### 1.2 — TSD check (ported from `evaluate_ltr()`)

```
check_tsd(hit_gr, genome, tsd_length) → list(TSD_sequence, TSD_length, TSD_positions, confirmed)
```

- Extract `tsd_length` bp immediately outside both ends of `hit_gr`
- Handle strand (reverse-complement for `−` strand hits), same as `evaluate_ltr()`
- Exact match → confirmed; ≤1 mismatch and length ≥ 5 bp → confirmed with mismatch
- Returns NA positions and `confirmed=FALSE` if no match

### 1.3 — Junction BLAST checks

Two functions covering both LTR boundaries:

```
check_utr5_junction(hit_gr, genome, utr5_db_path, strand) → logical
```
- Checks the **3' boundary** of an LTR hit (the 5' LTR → 5'UTR/PBS junction)
- Extract 30 bp immediately **after** the hit end (+ strand) / before the hit start (− strand)
- BLAST against `_5UTR_tags` database
- Return `TRUE` if any hit with `length ≥ 15` and `pident ≥ 80`
- A positive result means the candidate is the **5' LTR** of an unannotated complete element, with its 5'UTR still attached

```
check_ppt_junction(hit_gr, genome, ppt_db_path, strand) → logical
```
- Checks the **5' boundary** of an LTR hit (the PPT/3'UTR → 3' LTR junction)
- Extract 30 bp immediately **before** the hit start (+ strand) / after the hit end (− strand)
- BLAST against `_PPT_tags` database
- Return `TRUE` if any hit with `length ≥ 15` and `pident ≥ 80`
- A positive result means the candidate is the **3' LTR** of an unannotated complete element, with internal sequence still present upstream

Together these two checks discriminate both ends of any unannotated element:

```
5' boundary check         hit            3' boundary check
(PPT_tags)                               (5UTR_tags)
      ↓                                        ↓
[ PPT ][ TG... 3'LTR ...CA ]    vs    [ TG... 5'LTR ...CA ][ 5'UTR/PBS ]
 positive = 3'LTR of fragment          positive = 5'LTR of fragment
```

### 1.4 — GFF3 output helpers

```
make_solo_ltr_gff3(hits_df, tsd_results) → GRanges
```

- Assembles `solo_LTR` parent features and `target_site_duplication` child features
- Assigns sequential IDs `soloLTR_XXXXXXXX`
- Sets all required attributes: `Final_Classification`, `Identity`, `Coverage`, `TSD`, `Rank`, `UTR5_junction`, `PPT_junction`, `PBS_check`

### 1.5 — Statistics table

```
get_solo_ltr_statistics(gff_out, complete_elements_gff) → data.frame
```

- Counts SL and SL_noTSD per lineage
- Computes `Rsf = SL / n_complete_elements` per lineage (complete element count read from input GFF)

---

## Phase 2 — LTR library builder (`utils/build_ltr_library.R`)

Called once by the Python wrapper before chunking. Takes the full GFF3 and genome FASTA as input.

### CLI

```
utils/build_ltr_library.R -g <dante_ltr.gff3> -s <genome.fasta> -o <output_prefix> [-t <threads>]
```

Output files written to `<output_prefix>`:
- `_LTR_library.fasta` — consensus LTR sequences (BLAST queries)
- `_5UTR_tags.fasta` — 30 bp 5'UTR junction tags (BLAST database source); covers 3' boundary of 5' LTR hits
- `_PPT_tags.fasta` — 30 bp PPT/3'UTR junction tags (BLAST database source); covers 5' boundary of 3' LTR hits
- `_tsd_length_map.tsv` — two-column table `lineage \t tsd_length`

### 2.1 — Parse and rank-prioritise LTR features

- Import GFF3 with `rtracklayer::import()`
- Extract all features with `type == "long_terminal_repeat"`
- Join with parent `transposable_element` to get `Rank` and `Final_Classification`
- Assign priority: DLTP=4, DLT=3, DLP=2, DL=1

### 2.2 — Extract LTR sequences with flanking

- For each LTR feature, retrieve sequence ± 15 bp from genome using `getSeq()`
- Store pure LTR sequence (no flanking) separately for clustering
- Name sequences as `<seqname>_<start>_<end>_<lineage>`

### 2.3 — Per-lineage MMseqs2 clustering at 90 %

- Split LTR sequences by `Final_Classification`
- For each lineage with ≥ 2 sequences: call MMseqs2 (reuse `mmseq_clustering.R` logic) at 90 % identity
- For lineages with only 1 sequence: pass through directly
- Record cluster membership; retain highest-priority representative per cluster (by Rank, then by length)

### 2.4 — MAFFT MSA and consensus per cluster

For each cluster with ≥ 3 members:
- Write member sequences (with ± 15 bp flanking) to temp FASTA
- Run: `mafft --auto --thread <N> <input> > <output>`
- Parse MSA; call majority-vote consensus (gap threshold 50 %: column dropped if > 50 % gaps)
- Trim consensus to the region corresponding to pure LTR (remove flanking columns by locating TG/CA boundary in the consensus)

For clusters with < 3 members: use the representative sequence directly (no MSA).

Output: one consensus FASTA entry per cluster, named `>LTR_<XXXXXXXX> <Final_Classification>`.

### 2.5 — Build junction tag databases

Both databases are built from complete elements in the input GFF3 and are used in Step 4b to identify whether a no-TSD hit is actually a boundary LTR of an unannotated element.

**5'UTR tag database** (3' boundary of 5' LTR hits):
- For each element with a `5LTR` child feature:
  - Locate the `long_terminal_repeat` child with `LTR=5LTR`
  - Extract 30 bp immediately **3' of that LTR's end** (strand-aware)
  - This captures PBS and the start of the 5'UTR — conserved sequence that will be present if the hit is a 5' LTR still attached to internal sequence
  - Name: `>UTR5_<ID> <Final_Classification>`
- Write to `_5UTR_tags.fasta`; build BLAST DB with `makeblastdb -dbtype nucl`

**PPT/3'UTR tag database** (5' boundary of 3' LTR hits):
- For each element with a `3LTR` child feature:
  - Locate the `long_terminal_repeat` child with `LTR=3LTR`
  - Extract 30 bp immediately **5' of that LTR's start** (strand-aware)
  - This captures the PPT (polypurine tract) and flanking 3'UTR — conserved sequence that will be present if the hit is a 3' LTR still attached to internal sequence
  - Name: `>PPT_<ID> <Final_Classification>`
- Write to `_PPT_tags.fasta`; build BLAST DB with `makeblastdb -dbtype nucl`

### 2.6 — Build empirical TSD length map

- Call `extract_tsd_length_map()` from `solo_ltr_utils.R`
- Write result to `_tsd_length_map.tsv`

---

## Phase 3 — Core detection (`utils/detect_solo_ltr.R`)

Called once per genome chunk by the Python wrapper. Each invocation receives a FASTA chunk, the matching GFF3 slice, plus the pre-built library files from Phase 2 (paths passed as arguments).

### CLI

```
utils/detect_solo_ltr.R \
  -s <chunk.fasta> \
  -g <chunk_dante_ltr.gff3> \
  -l <LTR_library.fasta> \
  -u <5UTR_tags.fasta> \
  -p <PPT_tags.fasta> \
  -m <tsd_length_map.tsv> \
  -o <output_prefix> \
  -c <threads> \
  [-i <min_identity=80>] \
  [-C <min_coverage=0.8>]
```

### 3.1 — BLAST LTR library vs genome chunk

```bash
makeblastdb -in <chunk.fasta> -dbtype nucl -out <chunk_db>
blastn -task blastn -query LTR_library.fasta -db <chunk_db> \
  -num_threads <N> -dust no -perc_identity <min_identity> \
  -outfmt "6 qaccver saccver pident length qlen qstart qend sstart send sstrand evalue bitscore" \
  -out <blast_out.tsv>
```

Filter raw BLAST output:
- `length / qlen ≥ min_coverage`
- `evalue ≤ 1e-5`
- Parse `Final_Classification` from query name

Convert filtered hits to `GRanges` (handle `sstrand` for minus-strand hits: swap sstart/send, set strand `−`).

### 3.2 — Remove hits overlapping annotated complete elements

- Load input GFF3 chunk
- Extract all annotated features (all types, to catch LTR children too)
- Call `trim_gr(hits_gr, annotated_gr)` from `ltr_utils.R`; discard any hit overlapping by > 20 bp

### 3.3 — TSD check for each surviving hit

- Load `tsd_length_map.tsv`
- For each hit, call `check_tsd(hit_gr, genome_seq, tsd_length)` from `solo_ltr_utils.R`
- Partition hits into `tsd_confirmed` (SL) and `tsd_failed` (proceed to 3.4)

### 3.4 — Junction checks + PBS check (tsd_failed hits only)

Run in parallel with `mclapply`. Three independent checks per hit:

- `check_utr5_junction(hit_gr, genome_seq, utr5_db_path, strand)` → `UTR5_junction` attribute
  Checks 30 bp **after** the hit end (+ strand). Positive = hit is likely a **5' LTR** of an unannotated element.

- `check_ppt_junction(hit_gr, genome_seq, ppt_db_path, strand)` → `PPT_junction` attribute
  Checks 30 bp **before** the hit start (+ strand). Positive = hit is likely a **3' LTR** of an unannotated element.

- Reuse `add_pbs()` / `add_pbs_hemi()` from `ltr_utils.R` on the 30 bp after the hit end → `PBS_check` attribute
  Secondary confirmation for the 5'-LTR scenario.

All three results are recorded as attributes on `SL_noTSD` features; no hits are discarded based on them. A hit with both `UTR5_junction=negative` and `PPT_junction=negative` is more likely to be a genuine old insertion where the TSD has diverged, rather than a fragment of an unannotated element.

### 3.5 — Assemble and write output

- Call `make_solo_ltr_gff3(hits_df, tsd_results)` from `solo_ltr_utils.R`
- Export: `rtracklayer::export(gff_out, "<output_prefix>.gff3", format="gff3")`
- Call `get_solo_ltr_statistics(gff_out, input_gff)` → write `<output_prefix>_statistics.csv`

---

## Phase 4 — Python wrapper (`dante_ltr_solo`)

Mirrors `dante_ltr` in structure. The existing `Gff3Feature` class and all chunk-splitting functions can be imported directly from `dante_ltr` or copied into a shared module.

### 4.1 — Argument parsing

```
dante_ltr_solo [-h] -g GFF3 -s FASTA -o OUTPUT [-c CPU]
               [-i MIN_IDENTITY] [-C MIN_COVERAGE]
               [-S MAX_CHUNK_SIZE] [-v]
```

### 4.2 — Library build (once, before chunking)

```python
tool_path = os.path.dirname(os.path.realpath(__file__))
lib_prefix = os.path.join(tempdir, "ltr_lib")
subprocess.check_call([
    f"{tool_path}/utils/build_ltr_library.R",
    "-g", args.gff3, "-s", args.reference_sequence,
    "-o", lib_prefix, "-t", str(args.cpu)
])
```

Paths to `_LTR_library.fasta`, `_5UTR_tags.fasta`, `_tsd_length_map.tsv` are derived from `lib_prefix`.

### 4.3 — Genome chunking and per-chunk detection

Reuse `split_fasta_to_chunks()` and `recalculate_gff3_coordinates()` from `dante_ltr`. For each chunk pair (FASTA + GFF3 slice):

```python
subprocess.check_call([
    f"{tool_path}/utils/detect_solo_ltr.R",
    "-s", chunk_fasta, "-g", chunk_gff,
    "-l", ltr_library, "-u", utr5_tags, "-m", tsd_map,
    "-o", chunk_output_prefix,
    "-c", str(args.cpu),
    "-i", str(args.min_identity),
    "-C", str(args.min_coverage)
])
```

### 4.4 — Merge chunk outputs

- GFF3: recalculate coordinates back to original genome (`recalculate_gff3_back_to_original_coordinates()`), concatenate, remove duplicates at chunk boundaries using `get_unique_features()` logic (adapted for `solo_LTR` feature type)
- Statistics: sum counts per lineage across chunks (reuse `sum_up_stats_files()`); recompute Rsf from totals
- Prepend version comment to final GFF3 (`add_version_to_gff3()`)

---

## Phase 5 — Dependency and packaging

### 5.1 — Add MAFFT to `requirements.txt`

```
mafft>=7.490
```

Also add to conda recipe (`meta.yaml` if present).

### 5.2 — Version bump

Increment `version.py` when the feature is complete.

### 5.3 — `tests.sh` additions

```bash
echo "Running test 5: solo LTR detection"
./dante_ltr_solo \
  -g tmp/test_output1.gff3 \
  -s test_data/sample_genome.fasta \
  -o tmp/test_solo_output \
  -c $NCPU_TO_USE

cat tmp/test_solo_output_statistics.csv
```

---

## Implementation order and dependencies

```
Phase 1 (solo_ltr_utils.R)
    ↓
Phase 2 (build_ltr_library.R)    ←  depends on Phase 1
    ↓
Phase 3 (detect_solo_ltr.R)      ←  depends on Phases 1 & 2
    ↓
Phase 4 (dante_ltr_solo wrapper) ←  depends on Phases 2 & 3
    ↓
Phase 5 (packaging + tests)      ←  depends on Phase 4
```

Phases 2 and 3 can be developed and tested independently as standalone R scripts before the Python wrapper is written. Use `Rscript utils/build_ltr_library.R ...` and `Rscript utils/detect_solo_ltr.R ...` directly against `test_data/` during development.

---

## Validation approach

1. Run `dante_ltr` on `test_data/sample_genome.fasta` to get `test_output1.gff3` (already done by tests.sh test 1).
2. Run `dante_ltr_solo` on that output. Inspect:
   - At least some SL hits are found (the sample genome is known to contain LTR-RTs)
   - SL features have valid TSD child features
   - No SL feature overlaps a feature in `test_output1.gff3`
   - `_statistics.csv` contains plausible Rsf values (> 0)
3. Manual spot-check: for one SL call, verify the TSD sequence by extracting flanking sequence from the FASTA directly.
