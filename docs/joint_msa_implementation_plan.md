# Implementation plan — joint 5'+3' LTR MSA for boundary detection

## 1. Motivation

The previous implementation (committed in `0.4.0.6`) clusters *and* MSAs 5'LTRs
only. That gives the change-point-corrected 5' boundary (from the random
genomic 5' flank), but leaves the 3' boundary stuck at the median annotated
position, because the 3' flank of a 5'LTR is the conserved start of the
internal region — no conservation drop to detect.

The original design
(`solo_ltr_library_consensus_design.md`, §2) handles this with asymmetric
flank conservation:

| LTR copy | 5' flank          | Body      | 3' flank          |
|----------|-------------------|-----------|-------------------|
| 5'LTR    | random (genomic)  | conserved | conserved (PBS)   |
| 3'LTR    | conserved (PPT)   | conserved | random (genomic)  |

So the 5' boundary is detectable from the 5'LTR subset's 5' flank, and the
3' boundary is detectable from the 3'LTR subset's 3' flank. This plan
restores that signal by including sibling 3'LTRs in the MSA.

## 2. What changes

**Keep unchanged**
- MMseqs2 clustering operates on 5'LTR bodies, grouped by lineage
  (`--min-seq-id 0.9 -c 0.8`).
- 5'UTR / PPT tag databases (built from the full LTR set as today).
- Output file set: `*_LTR_library.fasta`, `*_LTR_library_map.tsv`,
  `*_LTR_library_boundary_qc.tsv`, `*_5UTR_tags.fasta`, `*_PPT_tags.fasta`,
  `*_tsd_length_map.tsv`.
- Library ID scheme (`LTR_%06d`).

**Change**
- For each cluster, pull the sibling 3'LTR of every 5'LTR member (via the
  shared `Parent` TE id). Pass a **joint** DNAStringSet (5'LTRs + 3'LTRs)
  with `±flank` bp to MAFFT. Each row carries a role tag (`"5LTR"` or
  `"3LTR"`).
- Compute per-column base counts **twice** — once per subset. Derive two
  sliding-window conservation profiles.
- 5' boundary: change-point scan on the 5'LTR-subset profile, walking
  leftward from the annotated 5' column (unchanged algorithm).
- 3' boundary: change-point scan on the 3'LTR-subset profile, walking
  **rightward** from the annotated 3' column (new mirror of the 5p scan).
- Consensus is still a column-wise majority between the two corrected
  boundaries, but computed from the **full** base-count matrix (both
  subsets contribute, reinforcing the body consensus).
- QC TSV grows: now has both `corrected_5_col` / `shift_5` *and*
  `corrected_3_col` / `shift_3`, plus per-subset member counts.

## 3. Concrete code changes

All changes are in `utils/build_ltr_library.R`.

### 3.1 New helper — `detect_3p_change_point`

Mirror of `detect_5p_change_point`, walking rightward:

```r
detect_3p_change_point <- function(cw, ann_col, flank, T_high, T_low) {
  n <- length(cw)
  if (is.na(ann_col) || ann_col < 1L || ann_col > n) return(NA_integer_)
  hi <- min(n, ann_col + flank)
  high_seen <- FALSE
  for (i in seq.int(ann_col, hi)) {
    if (!high_seen && cw[i] >= T_high) high_seen <- TRUE
    if (high_seen && cw[i] < T_low) return(as.integer(i - 1L))
  }
  NA_integer_
}
```

The 5p variant is untouched.

### 3.2 Update extraction (main section)

Replace the `ltr_feat_5`-centric extraction with "keep everything, mark
which is which":

- `ltr_feat_all <- ltr_feat` (already includes both 5' and 3').
- `ltr_ltr_all <- ltr_ltr` (role string per feature).
- `ltr_parent_all <- ltr_parent` (TE id per feature).
- Compute body and extended sequences for **all** LTRs (same width+strand
  logic; just run over everything, not just 5'LTRs).
- Record `ann_5_pos_raw`, `ann_3_pos_raw`, `bio_body_len`,
  `bio_5flank_len`, `bio_3flank_len` for every LTR.
- For 3'LTRs: `ann_5_pos_raw` still marks the *biological* 5' end of the
  LTR body (the column immediately after the biological 5' flank). In the
  joint MSA, the 5'LTR's body 5' end and the 3'LTR's body 5' end should
  align (same family sequence, insertion-time identical).

### 3.3 Parent-indexed sibling lookup

Build a per-parent lookup once, outside the loop:

```r
by_parent_5 <- split(which(ltr_ltr_all == "5LTR"), ltr_parent_all[ltr_ltr_all == "5LTR"])
by_parent_3 <- split(which(ltr_ltr_all == "3LTR"), ltr_parent_all[ltr_ltr_all == "3LTR"])
```

For any 5'LTR at index `i`, its parent is `ltr_parent_all[i]` and the
sibling 3'LTR index is `by_parent_3[[parent]][1]` (or `NA` if missing).

### 3.4 Cluster loop changes

The clustering step still operates on 5'LTR bodies only — identical to
today:

```r
membership <- mmseqs_cluster(body_5_lin, identity = 0.9, threads = opt$threads)
```

For each cluster, assemble the MAFFT input:

```r
members_5_idx <- <indices of 5'LTRs in this cluster, in all-LTR frame>
# Pull sibling 3'LTR for each member
parents <- ltr_parent_all[members_5_idx]
members_3_idx <- unname(unlist(lapply(parents,
                       function(p) {
                         s <- by_parent_3[[p]]
                         if (is.null(s) || length(s) == 0L) NA_integer_ else s[1L]
                       })))

# Drop 5' members whose 3' sibling is missing? No — keep the 5' member,
# just omit the missing 3'. Result: MSA has all 5'LTRs + available 3'LTRs.
joint_idx   <- c(members_5_idx, members_3_idx[!is.na(members_3_idx)])
joint_roles <- c(rep("5LTR", length(members_5_idx)),
                 rep("3LTR", sum(!is.na(members_3_idx))))

ext_input   <- ltr_ext_seqs_all[joint_idx]
b5_input    <- ann_5_pos_raw_all[joint_idx]
b3_input    <- ann_3_pos_raw_all[joint_idx]
```

Call the boundary-aware consensus:

```r
res <- mafft_boundary_consensus(
  extended_seqs    = ext_input,
  body_starts_raw  = b5_input,
  body_ends_raw    = b3_input,
  roles            = joint_roles,      # NEW
  flank            = F_flank,
  threads          = opt$threads,
  save_aln_path    = aln_path,
  ltr_id           = ltr_id,
  lineage          = lin,
  timing_env       = timing_env
)
```

Small-cluster fallback (n < min_cluster_size) stays on the 5'LTR body
only — unchanged.

### 3.5 `mafft_boundary_consensus` — accept `roles`

Signature gains `roles = NULL` (backward-compat: if NULL, treat every row
as 5'LTR and degrade to the previous median-annotated 3' end).

Inside, after the MSA matrix is built:

```r
bc_all <- base_counts(mat)

if (!is.null(roles)) {
  is_5 <- roles == "5LTR"
  is_3 <- roles == "3LTR"
  bc_5 <- base_counts(mat[is_5, , drop = FALSE])
  bc_3 <- if (any(is_3)) base_counts(mat[is_3, , drop = FALSE]) else NULL
  cprof_5 <- conservation_profile(bc_5)
  cw_5    <- sliding_mean(cprof_5, conservation_window)
  if (!is.null(bc_3)) {
    cprof_3 <- conservation_profile(bc_3)
    cw_3    <- sliding_mean(cprof_3, conservation_window)
  } else {
    cw_3 <- NULL
  }
} else {
  bc_5 <- bc_all; cw_5 <- sliding_mean(conservation_profile(bc_5),
                                       conservation_window); cw_3 <- NULL
}
```

Annotated column medians (taken over the rows of the subset they refer
to):

```r
ann_5_col <- as.integer(median(start_cols[is_5]))     # 5'LTR rows
ann_3_col <- if (any(is_3)) as.integer(median(end_cols[is_3]))  # 3'LTR rows
             else as.integer(median(end_cols))
```

(When both subsets are present, the body columns of 5'LTRs and 3'LTRs
should align, so median choice is a minor bias; taking it from the
subset whose flank drives the scan keeps them consistent.)

Boundary detection:

```r
corrected_5 <- detect_5p_change_point(cw_5, ann_5_col, flank, T_high, T_low)
if (is.na(corrected_5)) corrected_5 <- ann_5_col

corrected_3 <- if (!is.null(cw_3)) {
  detect_3p_change_point(cw_3, ann_3_col, flank, T_high, T_low)
} else NA_integer_
if (is.na(corrected_3)) corrected_3 <- ann_3_col
```

Sanity on *total* body length:

```r
new_len <- corrected_3 - corrected_5 + 1L
if (is.na(new_len) || new_len < max_length_shrink * median_body_len) {
  corrected_5 <- ann_5_col
  corrected_3 <- ann_3_col
}
```

Consensus uses the **full** base-count matrix (both subsets vote on each
body column):

```r
consensus <- majority_from_counts(bc_all, corrected_5, corrected_3)
```

### 3.6 QC TSV — new columns

| column | meaning |
|---|---|
| `n_5ltr` | rows in the 5'LTR subset |
| `n_3ltr` | rows in the 3'LTR subset |
| `annotated_5_col` | median annotated 5' boundary (MSA column) |
| `corrected_5_col` | change-point-corrected 5' boundary |
| `shift_5` | corrected - annotated (bp, signed) |
| `annotated_3_col` | median annotated 3' boundary |
| `corrected_3_col` | change-point-corrected 3' boundary |
| `shift_3` | corrected - annotated (bp, signed) |
| `consensus_length` | final consensus width (bp) |
| `median_annot_body_len` | median body length across *5'LTR* members |

(`n_members` is retained for backwards readability: equals
`n_5ltr + n_3ltr`.)

### 3.7 No changes required

- `utils/detect_solo_ltr.R` — consumes `_LTR_library.fasta` and the map;
  format unchanged.
- `utils/select_solo_representatives.R` — untouched.
- `solo_ltr_utils.R` — untouched.
- Python wrapper — untouched.

## 4. Testing strategy

Compare current (committed `0.4.0.6`) vs new implementations on
`test_data/g1_dante_ltr.gff3` + `test_data/g1.fasta`.

### 4.1 Library build diffs

1. **Library FASTA**:
   - Record count should be identical (same cluster count).
   - Length distribution: expect most representatives to change slightly;
     median absolute delta under ~20 bp is reasonable.
2. **Boundary QC TSV**:
   - For multi-member clusters, compare `shift_3` distribution (was
     always 0 before; now should be mostly non-zero).
   - For clusters where both shifts are large and opposite-signed, inspect
     the alignment by eye.
3. **Alignment files** (`alignments/`): inspect 1-2 clusters for a clear
   conservation drop on both flanks.

### 4.2 Downstream solo-LTR diffs

Run `dante_ltr_solo` end-to-end with each library:
- Representative `solo_ltr.gff3` counts per lineage.
- SL vs SL_noTSD rate (ideally tiny improvement as library boundaries
  are tighter; BLAST now matches the LTR precisely).
- `Coverage` distribution in representatives: tighter boundaries should
  push mean coverage closer to 1.0.

### 4.3 Regression gate

Must pass:
- `./tests.sh 4` exits 0.
- No MAFFT / R error in the log.
- QC TSV has one row per cluster.

Acceptance:
- At least 50 % of multi-member clusters produce a non-zero `shift_3`
  (proving the new scan is firing).
- Library record count equals today's.
- Downstream `dante_ltr_solo` completes without error and the
  representative output sanity-check file counts are within ±10 % of the
  current implementation.

## 5. Out of scope

- Using the new `shift_3` to back-project corrected boundaries onto
  original GFF3 records (design §3.6 — still deferred).
- Any change to TSD handling, junction checks, or representative
  selection.
- Changing the asymmetric interpretation (i.e. we still do NOT detect the
  3' boundary from 5'LTR rows or vice versa — those rows' flanks are
  conserved-internal, not random).
