# Implementation plan — 5'LTR-only flank-aware consensus for `dante_ltr_solo` library

**Scope:** `utils/build_ltr_library.R` only. All other pipeline stages (`dante_ltr`, `detect_solo_ltr.R`, `solo_ltr_utils.R`) consume the same file formats as today and do not change.

**Anchoring document:** `solo_ltr_library_consensus_design.md` (with user annotations).

---

## 1. How the user's comments narrow the design

| Comment location | User's note | Implication |
|---|---|---|
| §5.2 "clusters with only 5'LTRs" | "clustering with 5'LTR only is preferred approach" | Drop 3'LTRs from the library-build pipeline entirely. Cluster, align, and derive consensus from 5'LTRs only. 3'LTRs continue to feed the PPT-tag database (unchanged). |
| §5.3 solo LTRs in input | "this will not happen, normal pipeline returns only pairs" | Skip the solo-LTR filter — no input sanitation needed. |
| §5.4 truncated LTRs | "majority will outvote" | No special handling. |
| §5.6 insertion-site preference | "not important at this stage" | Single flat set of thresholds; no bimodality checks. |
| §5.7 unrelated element pairs | "given DANTE_LTR annotation this should not be issue" | No sibling-identity check. |

**Net effect:** the element-pair-preserving clustering (design §3.2), the anchored/joint MSA choice (design §3.3), and the 3'-subset conservation track (design §3.4) all collapse. The pipeline becomes: *cluster 5'LTRs → align 5'LTRs with flanks → detect 5'-boundary change-point → use median annotated 3' end → emit consensus between the two boundaries.*

---

## 2. Pipeline (final shape)

```
DANTE_LTR GFF3  +  reference FASTA
   │
   ▼
[A] Filter features:
     long_terminal_repeat with LTR=5LTR   →  used for library
     long_terminal_repeat with LTR=3LTR   →  used for PPT-tag DB only (unchanged)
   │
   ▼
[B] Per 5'LTR, extract body-only DNAStringSet (for clustering)
     and extended DNAStringSet with ±F bp flanks (for MSA/consensus).
     Record per-member offsets: bio_5flank_len, bio_body_len, bio_3flank_len.
   │
   ▼
[C] Partition by Final_Classification (lineage) — unchanged
   │
   ▼
[D] Per lineage: MMseqs2 easy-cluster on body-only sequences (identity 0.9)
   │
   ▼
[E] Per cluster with ≥ min_cluster_size members:
       • MAFFT on extended sequences (flanks included)
       • per-column conservation  c(i) = top_freq * (1 - gap_frac)
       • sliding-window mean (conservation_window = 5 bp)
       • map each member's annotated 5'/3' body ends to MSA columns
       • 5' boundary:  change-point scan on the 5' flank
       • 3' boundary:  median of mapped annotated 3' ends
       • sanity checks (shift ≤ F, length change ≤ 50 %)
       • consensus = column-wise majority of non-gap residues between
                     the corrected 5' boundary and the 3' boundary.
                     NO 50 %-gap filter inside the body.
   │
   ▼
[F] Small cluster (< min_cluster_size):
       unchanged — use highest-rank single LTR with its annotated coords.
   │
   ▼
[G] Write library + ID map + (new) per-cluster boundary-QC TSV.
```

---

## 3. Key algorithms

### 3.1 Flank-aware extraction

For each 5'LTR feature on strand `s`:

1. Build a stranded GRanges covering `[max(1, start - F), min(seqlength, end + F)]`.
2. `getSeq()` returns the sequence in biological orientation (reverse-complemented for `-` strand).
3. Compute offsets **in the biological orientation**:
   - `bio_5flank_len` = length of flank on the biological 5' side
     - `+` strand: `start(orig) - start(ext)`
     - `-` strand: `end(ext) - end(orig)`
   - `bio_body_len` = `width(orig)` (unchanged by extension)
   - `bio_3flank_len` = total length − `bio_5flank_len` − `bio_body_len`
4. Annotated boundary positions in the raw extended sequence (1-based):
   - `ann_5_pos_raw = bio_5flank_len + 1`
   - `ann_3_pos_raw = bio_5flank_len + bio_body_len`

### 3.2 Clustering

MMseqs2 `easy-cluster --min-seq-id 0.9 -c 0.8 --cov-mode 0`, driven by body-only sequences (same as today, just with 3'LTRs filtered out). Cluster memberships are applied to the extended sequences in [E].

### 3.3 MAFFT on extended sequences

Same `mafft --auto --thread N --quiet` as today. Input = extended sequences; the body anchors the alignment, the flanks fall into columns around it.

### 3.4 Mapping raw positions → MSA columns

For each aligned row, walk left-to-right counting non-gap characters until the target raw position is reached. Return the 1-based MSA column. Implemented in base R (no external deps).

### 3.5 Per-column conservation

```
gap_frac(i) = mean(col == "-")
non_gap     = col excluding "-" and "N"
top_freq(i) = max(table(non_gap)) / length(non_gap)   ; 0 if empty
c(i)        = top_freq(i) * (1 - gap_frac(i))
cw(i)       = mean of c over a window_size-wide window centred on i
```

Computed once per cluster over every MSA column.

### 3.6 5'-boundary change-point

`ann_5_col` = median of the mapped `ann_5_pos_raw` across all members.

Scan from `ann_5_col + F` (deep inside flank, high end of flank range) leftward:

1. Set `high_seen = FALSE`. Walk `i` from inside the body (`ann_5_col - 1`) outward to the flank.
2. Actually the design's directionality is: walk from **inside body** outward into flank. So go from `i = ann_5_col` leftward down to `i = ann_5_col - F`.
3. Track whether `cw(i)` has been above `T_high` (inside the body, expected).
4. Return the first column where, after the high-plateau, `cw(i) < T_low` — that column's right edge (`i + 1`) is the corrected 5' boundary.
5. Clip to `ann_5_col ± F`; on no-transition, fall back to `ann_5_col`.

### 3.7 3'-boundary (no flank-based detection)

`ann_3_col = median(mapped ann_3_pos_raw across members)`. That's the 3' boundary.

(Design note: for 5'LTRs the 3' flank is the conserved start of the internal region — a conservation drop is not expected at the LTR/internal junction. The median of annotated 3' ends still reduces noise compared to trusting any single member, which is the improvement over today.)

### 3.8 Sanity checks

- |corrected − annotated| ≤ `F` on the 5' side (guaranteed by the scan range, so mostly defensive).
- final body length ≥ 50 % of the cluster median annotated body length; otherwise fall back to annotated 5' boundary.
- resulting consensus ≥ 50 bp; otherwise fall back to highest-rank input LTR (current behaviour).

### 3.9 Consensus emission

For each MSA column `i` with `corrected_5_col ≤ i ≤ ann_3_col`:

- drop `-` and `N`
- if non-empty → majority base; ties broken by first-encountered order (R's `which.max(table(...))`)
- if empty → skip column

No 50 %-gap filter (design §3.5). Each retained column contributes one base to the consensus.

### 3.10 New QC output

`{output}_LTR_library_boundary_qc.tsv` with columns:

| col | meaning |
|---|---|
| `ltr_id` | library ID |
| `Final_Classification` | lineage |
| `n_members` | cluster size |
| `annotated_5_col` | median annotated 5' boundary (MSA col) |
| `corrected_5_col` | change-point-corrected 5' boundary (MSA col) |
| `shift_5` | corrected − annotated (bp, signed) |
| `annotated_3_col` | median annotated 3' boundary (MSA col) |
| `consensus_length` | final consensus width (bp) |
| `median_annot_body_len` | median body length across members (bp) |

One row per cluster that underwent MSA. Singletons and size-below-min clusters contribute a row with `NA` in the `corrected_*` / `shift_*` columns.

---

## 4. Parameter defaults

| Option flag | Default | Status |
|---|---|---|
| `--flank / -f` | **50** (was 15, unused) | Activate + increase. |
| `--min_cluster_size / -d` | 3 | Unchanged. |
| *(internal)* `conservation_window` | 5 | Design §4. |
| *(internal)* `T_high` | 0.7 | Design §4. |
| *(internal)* `T_low` | 0.4 | Design §4. |

Thresholds are hard-coded; they become CLI flags only if tuning demand arises.

---

## 5. Concrete code changes in `utils/build_ltr_library.R`

1. **Options block:** change `--flank` default `15L → 50L`; update help text ("flanking bp used for MAFFT and change-point boundary detection").
2. **Feature filter:** after computing `ltr_ltr`, restrict the clustering/consensus subset to `ltr_ltr == "5LTR"`. The 3'LTR subset stays in scope only for the PPT-tag section.
3. **Replace `mafft_consensus()`** with a new function `mafft_boundary_consensus(extended_seqs, body_starts_raw, body_ends_raw, flank, ...)` that:
   - writes input, runs MAFFT, reads alignment,
   - computes conservation profile,
   - maps per-member raw boundaries to MSA columns,
   - runs the 5' change-point scan,
   - derives 3' boundary from median,
   - emits majority consensus between boundaries,
   - returns `list(consensus = DNAStringSet, qc = data.frame_row)`.
4. **Extract both body-only and extended DNAStringSets** before the lineage loop. Keep names keyed on coordinates (as today) plus record the per-member offsets in parallel vectors.
5. **Lineage loop:**
   - Cluster using body sequences (unchanged call).
   - For each cluster ≥ `min_cluster_size`: call `mafft_boundary_consensus()` on the extended sequences.
   - For smaller clusters: current fallback path (highest-rank single annotated LTR).
6. **Accumulate** a `qc_rows` `data.frame` alongside `all_consensus`; write it to `{output}_LTR_library_boundary_qc.tsv` after the loop.
7. **Preserve all other outputs** (`_LTR_library.fasta`, `_LTR_library_map.tsv`, `_5UTR_tags.fasta`, `_PPT_tags.fasta`, `_tsd_length_map.tsv`) with identical schema.

### Helper functions (new, private to the script)

- `column_from_raw_pos(aligned_char_vec, raw_pos) -> int`: MSA-column of a given raw-sequence position in one aligned row.
- `conservation_profile(aln_matrix) -> numeric vector`: per-column `c(i)`.
- `sliding_mean(x, w) -> numeric vector`: length-preserving running mean (edges reuse a shorter window).
- `detect_5p_change_point(cw, ann_col, flank, T_high, T_low) -> int | NA`: returns corrected column or NA (NA → fall back).

---

## 6. Validation after implementation

1. `./tests.sh 4` must pass unchanged.
2. On `tmp/test_output1.gff3`, rerun `dante_ltr_to_library` and compare:
   - library FASTA record count (should match lineage × cluster count),
   - length distribution vs the prior run (expect tighter distribution, mean close to lineage-typical LTR length),
   - inspect `_LTR_library_boundary_qc.tsv` — most 5'-shifts should be ≤ a few bp.
3. Inspect one MSA in `alignments_dir` by eye for a sanity check that the flanks bracket a clear conservation jump on the 5' side.

Any regression vs the current output (empty library, dramatically different record count, consensuses shorter than 50 bp) is a fail.

---

## 7. Explicitly out of scope

- Option F from the design (per-LTR boundary back-projection into GFF3). Deferred.
- Re-clustering across lineages, aligner choice, or changes to the 5'UTR / PPT tag databases.
- CLI knobs for the internal thresholds.
