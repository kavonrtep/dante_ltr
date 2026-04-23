# Boundary validation via TSD and TG/CA motifs — implementation plan

## Goal

Add per-cluster validation signals for the MSA-driven boundary correction
in `utils/build_ltr_library.R`: does shifting the boundary from the
annotated column to the corrected column **preserve / gain / lose**
(a) the terminal TG/CA dinucleotide pair, and (b) the flanking Target Site
Duplication (TSD).

**Scope of this change is report-only.** Two new signal-pair measurements
are written into the boundary QC TSV, summarised at end-of-run, and
rendered in a standalone HTML report. No correction is rejected or
altered based on these signals in v1. Once we have data from Alyr plus
the larger **g2** genome we can re-evaluate whether either signal should
be promoted to a soft or hard gate.

## Context (why this is worth doing at all)

The change-point scan lands a boundary column using only the
conservation profile. That's noisy in clusters with short LTRs,
heterogeneous flanks, or degraded bodies. TG/CA and TSD are two
**independent** end-markers an LTR should have; they disagree with
conservation only when one is wrong.

Nuances worth keeping in mind:

- DANTE_LTR already biases per-element boundary search toward TG/CA, so
  annotated boundaries are *enriched* for TG/CA. The cluster MSA adds
  new information (evidence from many copies) that a single-element
  search cannot see.
- TSD in the annotation table is derived from flanking genomic
  sequence at detection time, per element. Reassessing it against the
  corrected cluster-level boundary, across multiple cluster members,
  can distinguish a bona-fide shift from an artifact.
- Some lineages genuinely lack TG/CA boundaries (e.g. certain lineages
  with different integrase specificity). We must not require TG/CA.
- Longer-term, the aggregated per-cluster signal could feed back into
  DANTE_LTR's original per-element annotation.

## Signals

### Signal A — TG/CA at the LTR terminal pair

LTR retrotransposons typically start with `TG` and end with `CA`
(5' → 3' on the forward biological strand of the element). The cluster
MSA gives us a consensus for this independent of the individual
genomic insertion contexts.

We compute at **two** locations:

- **At the annotated column pair** `(annotated_5_col, annotated_3_col)`
- **At the corrected column pair** `(corrected_5_col, corrected_3_col)`

And at **two granularities**:

- **Consensus**: majority base at each of the four positions involved
  (`5`, `5+1`, `3-1`, `3`) → is it exactly `TG…CA`? → boolean.
- **Per-member**: for each cluster row, what is the base (non-gap,
  majority vote if ambiguous) at those four positions on that row? →
  fraction of rows where the row has `TG…CA`.

Per-member fraction is the most informative for report-only. Consensus
is a useful secondary check (tells us whether the majority base across
the cluster gives TG/CA — a stricter signal than per-member fraction).

**Strictness (v1 default): exact `TG` and `CA`, 0 mismatches.** This is
the conservative choice and matches DANTE_LTR's own search semantics.
We can relax to 1 mismatch later if the data suggests it's too strict.

### Signal B — TSD around the element termini

A TSD is a short direct repeat (usually 4–6 bp) immediately flanking
the inserted element. For each cluster member (both 5'LTR-bearing and
3'LTR-bearing rows) we compare the `L` bp immediately left of the
corrected-5' boundary against the `L` bp immediately right of the
corrected-3' boundary on the **same genomic insertion**.

Because the MSA is built from joint 5'LTR + 3'LTR sequences, TSD
evaluation requires knowing which 5'LTR row is paired with which
3'LTR row (same parent TE) — we already have this mapping:
`parent_lin[midx] → sib_3[k]`.

Per cluster:

- For every `(5'LTR, 3'LTR)` sibling pair where both rows are present
  in the MSA, look up the raw sequence positions of the corrected-5'
  boundary on the 5'LTR member and the corrected-3' boundary on the
  3'LTR member using the existing `cum_mat` inverse mapping.
- Reuse the **same TSD semantics** as `check_tsd()` in
  `utils/solo_ltr_utils.R` (two-pass: exact match first, then
  1-mismatch for L ≥ 4; lineage-modal length from `tsd_map` ± 1,
  clipped to `[3, 8]`; fallback `4:6`).
- Tally the fraction of sibling pairs with a confirmed TSD.

Same dual computation as TG/CA: once at the annotated boundary pair,
once at the corrected boundary pair.

Some cluster members may have only a 5'LTR (no sibling 3'LTR). These
are skipped for TSD — denominator is `n_sibling_pairs`, not
`n_members`.

## New QC columns (four)

Added to the row written by `mafft_boundary_consensus()`:

| column | type | meaning |
|---|---|---|
| `tgca_frac_annotated` | numeric [0,1] or NA | fraction of cluster rows whose bases at `annotated_5_col, +1, annotated_3_col-1, annotated_3_col` form exactly `TG…CA` |
| `tgca_frac_corrected` | numeric [0,1] or NA | same, at `corrected_5_col, +1, corrected_3_col-1, corrected_3_col` |
| `tsd_frac_annotated` | numeric [0,1] or NA | fraction of sibling pairs in the cluster with a confirmed TSD at the annotated boundary positions (maps MSA col → raw pos via `cum_mat`) |
| `tsd_frac_corrected` | numeric [0,1] or NA | same, at the corrected boundary positions |

`NA` is used when the boundary column is out of range, the cluster has
no sibling pairs (TSD only), or the grow-cap reverted the correction
(in which case `annotated == corrected` and both columns carry the
same value).

**Consensus-level booleans are *not* added as columns** — they can be
derived by any consumer from the alignment + the two col indices, and
the per-member fraction already carries that information (a consensus
TG/CA is implied when `tgca_frac ≥ 0.5`). Keeping the QC schema tight.

## Where to plug in

Inside `mafft_boundary_consensus()` in `utils/build_ltr_library.R`,
**after** the grow-cap check (so we use the final
`corrected_5 / corrected_3` values that the consensus is actually
built from), and **before** the `qc <- data.frame(...)` block.

Signature additions:

```r
mafft_boundary_consensus(
  ...existing args...,
  tsd_map       = NULL,      # lineage -> modal TSD length
  sibling_pairs = NULL       # integer vector: for each 5'LTR row index,
                             # the row index of its paired 3'LTR row,
                             # or NA if no sibling in the MSA
)
```

Both are threaded down from the main loop. `tsd_map` is already
computed once per run at `build_ltr_library.R:466`. `sibling_pairs` is
cheap to construct alongside `joint_idx` in the main loop.

## Algorithm sketch (pseudocode)

```r
# After grow-cap, with ann_5_col, ann_3_col, corrected_5, corrected_3 all known.

tgca_frac <- function(mat, c5, c3) {
  if (c5 < 1 || c3 > ncol(mat) || c3 - 1 < c5 + 1) return(NA_real_)
  rows_ok <- mat[, c5]     == "T" & mat[, c5 + 1] == "G" &
             mat[, c3 - 1] == "C" & mat[, c3]     == "A"
  mean(rows_ok)
}

tgca_frac_ann  <- tgca_frac(mat, ann_5_col, ann_3_col)
tgca_frac_corr <- tgca_frac(mat, corrected_5, corrected_3)
```

For TSD, we need raw-position lookup:

```r
# cum_mat[r, c] = # non-gap bases in row r up to col c, inclusive.
# Raw position of the base at (r, c) in that row's underlying sequence
# is cum_mat[r, c] if the cell is non-gap; for a gap cell we use
# cum_mat[r, c] + 0 which gives "bases consumed so far" — for boundary
# lookup we require the specific column to be non-gap; otherwise skip
# that row.

raw_pos_at_col <- function(r, c) {
  if (!non_gap_mat[r, c]) return(NA_integer_)
  cum_mat[r, c]
}

tsd_frac <- function(c5, c3) {
  if (is.null(sibling_pairs) || c5 < 1 || c3 > ncol(mat)) return(NA_real_)
  ok <- 0L; denom <- 0L
  for (r5 in which(is_5)) {
    r3 <- sibling_pairs[r5]
    if (is.na(r3)) next
    p5 <- raw_pos_at_col(r5, c5)        # 1-based pos of 5' boundary in r5's raw seq
    p3 <- raw_pos_at_col(r3, c3)        # 1-based pos of 3' boundary in r3's raw seq
    if (is.na(p5) || is.na(p3)) next
    denom <- denom + 1L
    # Reuse check_tsd by constructing a synthetic DNAStringSet + GRanges
    # for the two sibling raw sequences. See "TSD adapter" below.
    if (check_tsd_pair(extended_seqs[r5], p5,
                       extended_seqs[r3], p3,
                       tsd_length = tsd_map[lineage])) {
      ok <- ok + 1L
    }
  }
  if (denom == 0L) NA_real_ else ok / denom
}
```

### TSD adapter

`check_tsd()` in `solo_ltr_utils.R` works on a single BLAST hit
against a genome. For per-cluster validation we aren't running BLAST;
we already have the extended LTR sequences `extended_seqs[r]` and raw
positions `p5` / `p3`. The cleanest approach is a small internal helper
`check_tsd_pair(seq5, p5, seq3, p3, tsd_length)` that reimplements the
core two-pass comparison (exact → 1-mismatch, longest-first within the
lineage-modal ± 1 range) without the GRanges wrapping. Concretely:

```r
check_tsd_pair <- function(seq5, p5, seq3, p3, tsd_length = NULL) {
  scan_lengths <- if (!is.null(tsd_length) && !is.na(tsd_length))
                    seq.int(max(3L, tsd_length - 1L),
                            min(8L, tsd_length + 1L))
                  else 4L:6L
  scan_lengths <- sort(unique(scan_lengths), decreasing = TRUE)
  max_len <- max(scan_lengths)

  l_end   <- p5 - 1L
  l_start <- max(1L, l_end - max_len + 1L)
  r_start <- p3 + 1L
  r_end   <- min(length(seq5), r_start + max_len - 1L)  # per-row length
  if (l_end < l_start || r_end < r_start) return(FALSE)

  l_full <- as.character(subseq(seq5, l_start, l_end))
  r_full <- as.character(subseq(seq3, r_start, r_end))

  for (len in scan_lengths) {
    if (nchar(l_full) < len || nchar(r_full) < len) next
    l <- substr(l_full, nchar(l_full) - len + 1L, nchar(l_full))
    r <- substr(r_full, 1L, len)
    if (l == r) return(TRUE)
  }
  for (len in scan_lengths) {
    if (len < 4L) next
    if (nchar(l_full) < len || nchar(r_full) < len) next
    l <- substr(l_full, nchar(l_full) - len + 1L, nchar(l_full))
    r <- substr(r_full, 1L, len)
    if (sum(strsplit(l, "")[[1L]] != strsplit(r, "")[[1L]]) <= 1L) return(TRUE)
  }
  FALSE
}
```

Inlined (not moved into `solo_ltr_utils.R`) to avoid coupling the
library builder to the solo-LTR pipeline. The logic is duplicated ~30
lines; that's fine.

**Strand note**: the extended sequences are already in biological
orientation (MAFFT input uses `getSeq` with strand). So the 5'LTR row's
left flank in its own extended sequence is the true upstream TSD, and
the 3'LTR row's right flank in its own extended sequence is the true
downstream TSD. No further orientation fixup is needed.

## `sibling_pairs` construction (caller-side)

In the main loop, where `joint_idx` / `joint_roles` are built:

```r
members_5 <- all_idx_lin[midx]              # indices in all-LTR frame
sib_3     <- unlist(lapply(parent_lin[midx], function(p) {
  s <- by_parent_3[[p]]; if (is.null(s) || length(s) == 0L) NA_integer_ else s[1L]
}))
members_3   <- sib_3[!is.na(sib_3)]
joint_idx   <- c(members_5, members_3)

# Row index within the MSA for each 5'LTR: positions 1..length(members_5).
# Its 3'LTR sibling's row index (if any): offset length(members_5) + rank in members_3.
sibling_pairs <- rep(NA_integer_, length(members_5))
pos3 <- which(!is.na(sib_3))
sibling_pairs[pos3] <- length(members_5) + seq_along(pos3)
# sibling_pairs is length == n_5 (not n_members). Pad with NA for 3'LTR rows
# before passing to mafft_boundary_consensus so it can be indexed uniformly:
sibling_pairs <- c(sibling_pairs, rep(NA_integer_, length(members_3)))
```

The wide-flank retry path (`extract_extended`) reuses the same
`joint_idx` (same members in the same order), so `sibling_pairs` is
valid unchanged for the retry.

## Summary log additions

End-of-run `[summary]` block (already in `build_ltr_library.R` around
line 719) gets two new rollups:

```
[summary]  TG/CA at corrected vs annotated boundary (n=<N> clusters with ≥2 members):
  corrected TG/CA frac    — min/med/max : 0.000 / 0.833 / 1.000
  annotated TG/CA frac    — min/med/max : 0.000 / 0.833 / 1.000
  clusters LOST  TG/CA after shift (corr < ann − 0.1) :  <k> / <N>
  clusters GAINED TG/CA after shift (corr > ann + 0.1) :  <k> / <N>
  clusters UNCHANGED (|corr-ann| ≤ 0.1)               :  <k> / <N>

[summary]  TSD at corrected vs annotated boundary (n=<N'> clusters with ≥1 sibling pair):
  corrected TSD frac      — min/med/max : 0.000 / 0.500 / 1.000
  annotated TSD frac      — min/med/max : 0.000 / 0.500 / 1.000
  clusters LOST  TSD (corr < ann − 0.1)  :  <k> / <N'>
  clusters GAINED TSD (corr > ann + 0.1) :  <k> / <N'>
  clusters UNCHANGED                    :  <k> / <N'>
```

`0.1` is a nominal jitter threshold to avoid counting single-row noise
as a shift. Tune after looking at first run.

Clusters with grow-cap-reverted corrections will show
`corrected == annotated` and fall into UNCHANGED for both signals.

## Inspector extension

`utils/inspect_ltr_alignment.R` (already prints QC + ASCII
conservation) gets a small appended block after the conservation plot:

```
Boundary validation:
  TG/CA (annotated) : 0.833   TG/CA (corrected) : 1.000   delta: +0.167
  TSD  (annotated) : 0.500   TSD  (corrected) : 0.667   delta: +0.167
```

Values are pulled directly from the QC TSV — no recomputation from the
alignment needed.

## HTML visual report

A self-contained HTML file
`<output>_LTR_library_boundary_report.html` is produced at the end of
the library-builder run. It complements the TSV and the end-of-run log
with visual summaries so patterns (systematic shifts, lineage
differences, grow-cap hits, TG/CA loss clusters) can be eyeballed in
one place.

### Generation

New utility script `utils/boundary_report.R`. It takes the QC TSV and
optionally the alignments directory, renders plots to PNG via base R
graphics (no external reporting framework — keeps the dependency
footprint tight), and inlines them as base64 `<img>` tags into a
single HTML file. No knitr / rmarkdown / pandoc dependency.

Invoked from `utils/build_ltr_library.R` at the very end of the run,
after the QC TSV is written and after the `[summary]` log block:

```r
if (!is.null(opt$alignments_dir)) {
  report_path <- paste0(opt$output, "_LTR_library_boundary_report.html")
  system2("Rscript", c(
    file.path(script_dir, "boundary_report.R"),
    "--qc",         qc_path,
    "--alignments", opt$alignments_dir,
    "--out",        report_path
  ))
}
```

Rendering a report only makes sense when alignments are kept, so the
step is conditional on `--alignments_dir` being set. If the script
fails we log a warning and continue — the report is an observability
aid, not a correctness dependency.

### Panels

All plots use base R graphics (or `lattice` if stratification is
cleaner there — but `base` suffices). Each panel is a PNG and embedded
into the HTML.

1. **Cluster size distribution**
   - Histogram of `n_members`, with a vertical line at
     `min_cluster_size`.
   - Side panel: count of clusters using MAFFT vs the single-member
     fallback.

2. **Shift distribution**
   - Two histograms side by side: `shift_5` and `shift_3` in bp.
     Overlay a line at 0. Annotate median and IQR.
   - Faceted (or separate panel) by `flank`: shows how often
     wide-flank retry was triggered and what shifts it found.

3. **TG/CA before vs after**
   - Scatter: x = `tgca_frac_annotated`, y = `tgca_frac_corrected`,
     one point per cluster, `alpha` scaled by `n_members`. Diagonal
     reference line `y = x`. Points above diagonal = GAINED, below =
     LOST.
   - Marginal density on both axes.
   - Colour by lineage family (e.g. Ty1/copia vs Ty3/gypsy), using a
     small palette. Legend placed in top-left corner.

4. **TSD before vs after**
   - Same format as panel 3 for `tsd_frac_annotated` vs
     `tsd_frac_corrected`.

5. **Per-lineage boxplots**
   - For each lineage (or lineage family if too many), side-by-side
     boxplots of `tgca_frac_corrected` and `tsd_frac_corrected`.
     Annotate `n_clusters` per box.
   - A small second row shows median `|shift_5|` and `|shift_3|` per
     lineage.

6. **Grow-cap activity**
   - Bar chart: count of clusters where the grow-cap reverted the
     correction (annotated == corrected on both sides despite a
     detection firing). Broken down by lineage family.

7. **Cluster-level "top cases" table**
   - Sortable plain HTML `<table>` of the 30 most interesting
     clusters, selected by a composite score (large `|shift|` AND
     `n_members ≥ 6` AND change in `tgca_frac ≥ 0.25` OR change in
     `tsd_frac ≥ 0.25`). Columns: ltr_id, lineage, n_members,
     shift_5, shift_3, flank, tgca_annotated, tgca_corrected,
     tsd_annotated, tsd_corrected, link-anchored to per-cluster panel
     (see next).

8. **Per-cluster mini panel** (optional, for the top ~20 only)
   - For each top-table cluster: a mini ASCII/matplotlib-style
     conservation-profile PNG (reusing the logic from
     `inspect_ltr_alignment.R` but in real plot form, not ASCII).
     Includes vertical lines for annotated vs corrected boundaries on
     both the 5' and 3' subsets.
   - Embedded via named anchors so the top-table rows link into them.

### HTML skeleton

```
<!doctype html>
<html><head>
  <meta charset="utf-8">
  <title>LTR library boundary report — &lt;output prefix&gt;</title>
  <style> /* simple, no framework: monospace headers, 900px max-width,
            table.striped, img.panel { display:block; margin:1em auto } */ </style>
</head><body>
  <h1>LTR library boundary report</h1>
  <p>Input: &lt;qc_path&gt; — &lt;N&gt; clusters — generated &lt;timestamp&gt;</p>

  <h2>1. Cluster sizes</h2>   <img class="panel" src="data:image/png;base64,...">
  <h2>2. Boundary shifts</h2> <img class="panel" src="data:image/png;base64,...">
  <h2>3. TG/CA: annotated vs corrected</h2> <img ...>
  <h2>4. TSD: annotated vs corrected</h2>    <img ...>
  <h2>5. Per-lineage breakdown</h2>          <img ...>
  <h2>6. Grow-cap activity</h2>              <img ...>
  <h2>7. Top 30 clusters of interest</h2>
    <table class="striped">...</table>
  <h2>8. Per-cluster profiles</h2>
    <section id="LTR_000145"><h3>LTR_000145</h3><img ...></section>
    ...
</body></html>
```

Plain CSS, no JS. Easy to diff by opening two reports from two runs
side by side.

### Why R + base64

- Already in the project's R/Bioconductor ecosystem, no new
  dependencies (no `rmarkdown`, no `knitr`, no `htmltools`).
- A single `.html` file works over scp / attachment / GitHub-rendered
  preview. No sibling `_files/` directory to keep in sync.
- Writing the HTML by hand from `paste0()` + `sprintf()` is trivial
  for a report with ~8 panels.

### Performance expectations

For Alyr (~200 clusters): < 2 s. For g2 (tens of thousands of
clusters expected): the bulk of the work is reading the QC TSV and
rendering ~6 aggregate PNGs + 20 per-cluster PNGs; well under 30 s
including alignment reads for the top-N panels.

### Optional CLI exposure

Add a standalone `Rscript utils/boundary_report.R --qc ... --out ...`
so the report can be regenerated against an existing QC TSV without
rerunning the builder. Useful for iterating on panel design.

## Out of scope (v1)

- **No gating.** We do not reject or alter any correction based on
  TG/CA or TSD in v1. The values only appear in the QC TSV and log.
- **No 1-mismatch TG/CA.** Exact match only for v1. If report numbers
  show many `TG/CA→TR/CA` near-misses we can add a second column
  later.
- **No consensus-level TG/CA columns.** Per-member fraction + raw
  alignment is sufficient; consumers can derive consensus from the
  stored alignment if needed.
- **No change to `utils/solo_ltr_utils.R`.** TSD logic is reused in
  spirit via a local helper, not by importing/refactoring `check_tsd`.
- **No retry or flank expansion based on TSD/TG/CA.** Wide-flank retry
  remains gated only on `detected_5 / detected_3` as today.

## Testing & validation

Correctness:

- Unit-level sanity: `tgca_frac` on a hand-built 4-row matrix (2 rows
  with TG…CA, 2 without) → 0.5.
- `check_tsd_pair` with known positive (`AAAAA…AAAAA`) → TRUE.
- `check_tsd_pair` with 1-mismatch at len=5 → TRUE; at len=3 → FALSE.

Integration:

- Run the full library builder on **Alyr** (`test_data/g1.fasta` +
  `g1_dante_ltr.gff3`; alignments already cached in `tmp/alyr_lib/`).
  Spot-check 5–10 clusters in the inspector:
  - A cluster with strong correction (`|shift| ≥ 20`) and high
    `n_members`. Expect `tgca_frac_corrected ≥ tgca_frac_annotated` if
    the scan found a real boundary.
  - The LTR_000145 grow-cap revert case. Expect annotated == corrected
    on both columns.
  - A wide-flank-retry case (from the run summary). Expect the
    corrected TG/CA frac at the retried boundary to be informative.
- Run on **g2** (`test_data/g2/genome.fasta` +
  `test_data/g2/DANTE_LTR.gff3`, ~1 GB genome, ~85 MB GFF3). This is
  our larger, more representative dataset. Expect many more clusters,
  broader lineage coverage, and a wider dynamic range of shifts —
  good for verifying that the TG/CA and TSD measurements scale
  sensibly and that report rendering stays fast.
- Drapa is intentionally **not** used for this evaluation — it's a
  distant species where the fallback classification path dominates,
  and it's not a good model for judging boundary-correction quality.

Performance:

- Per-cluster work is `O(n_members × 4)` for TG/CA and
  `O(n_sibling_pairs × max_tsd_len)` for TSD. Bounded by alignment
  size; sub-millisecond per cluster. Measure total added wall-clock on
  Alyr — should be < 1 s across the whole run.

## Rollout

1. Write this plan document. Commit on its own (`docs/`).
2. Implement the four QC columns + summary rollup + inspector block
   (one commit).
3. Implement `utils/boundary_report.R` + wire it into
   `build_ltr_library.R` (one commit — separate from step 2 so the QC
   logic and the visualisation land cleanly reviewable).
4. Run the new build on **Alyr** (cached alignments in
   `tmp/alyr_lib/`) and on **g2** (`test_data/g2/`). Examine the TSV,
   the `[summary]` log, and the HTML report. Capture a short findings
   note in `tmp/` (not committed).
5. **Move the `0.4.0.13` tag forward** to the implementation commit
   (matching the pattern from the earlier library-builder changes).
6. Re-evaluate based on the first-run numbers whether to promote
   TG/CA or TSD to a soft gate. That's a separate plan if we go
   there.

## Open questions to revisit after first run

- Is exact TG/CA too strict? (Baseline observation: DANTE already
  biases to TG/CA, so the annotated-column fraction should be near 1.
  If it isn't, that's signal.)
- Is per-member fraction noisier than consensus at the cluster level?
  (Compare the two on Alyr.)
- For lineages without TG/CA, should the QC TSV blank those columns
  rather than record low fractions? (Per-lineage masking would require
  a lineage list; simpler to report and let the consumer filter.)
