#!/usr/bin/env Rscript
# refine_mafft_fallback.R
# Per-cluster MAFFT change-point boundary refinement called as a
# subprocess by utils/refine_boundaries.py for clusters whose parasail
# validation rate falls below --mafft_fallback_threshold.
#
# Inputs:
#   --input_fasta   FASTA of cluster members in biological orientation,
#                   each row already extended by --flank bp on each side.
#   --members_tsv   per-row metadata (chrom, start, end, strand, role,
#                   ext_start, ext_end, body_start_raw, body_end_raw).
#                   Header line required.
#   --output_tsv    per-row output: name, corrected_g, motif_ok.
#   --flank         flank size used for the input extension (informational).
#   --boundary_motif  "TG/CA" or "none".
#   --threads       MAFFT threads.
#
# Logic mirrors utils/build_ltr_library.R::mafft_boundary_consensus, but
# the focus is per-member corrections — we report the corrected genomic
# position of each row's relevant boundary (5' for 5LTRs, 3' for 3LTRs)
# rather than building a consensus.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

option_list <- list(
  make_option("--input_fasta",   type = "character", default = NULL),
  make_option("--members_tsv",   type = "character", default = NULL),
  make_option("--output_tsv",    type = "character", default = NULL),
  make_option("--flank",         type = "integer",   default = 50L),
  make_option("--boundary_motif", type = "character", default = "TG/CA"),
  make_option("--threads",       type = "integer",   default = 1L)
)
opt <- parse_args(OptionParser(option_list = option_list))

for (k in c("input_fasta", "members_tsv", "output_tsv")) {
  if (is.null(opt[[k]]) || !nzchar(opt[[k]])) {
    stop(sprintf("Missing required arg --%s", k), call. = FALSE)
  }
}

# ---------- helpers (vendored from build_ltr_library.R, trimmed) ----------

base_counts <- function(mat) {
  list(
    A   = colSums(mat == "A"),
    C   = colSums(mat == "C"),
    G   = colSums(mat == "G"),
    T   = colSums(mat == "T"),
    gap = colSums(mat == "-"),
    n   = nrow(mat)
  )
}

conservation_profile <- function(bc) {
  non_gap_n <- bc$A + bc$C + bc$G + bc$T
  top_count <- pmax(bc$A, bc$C, bc$G, bc$T)
  top_freq  <- ifelse(non_gap_n > 0L, top_count / non_gap_n, 0)
  gap_frac  <- bc$gap / bc$n
  top_freq * (1 - gap_frac)
}

sliding_mean <- function(x, w) {
  n    <- length(x)
  half <- w %/% 2L
  cs   <- c(0, cumsum(x))
  i    <- seq_len(n)
  lo   <- pmax(1L, i - half)
  hi   <- pmin(n,  i + half)
  (cs[hi + 1L] - cs[lo]) / (hi - lo + 1L)
}

detect_5p_change_point <- function(cw, ann_col, flank, T_high, T_low,
                                   high_sustain = 2L) {
  n <- length(cw)
  if (is.na(ann_col) || ann_col < 1L || ann_col > n) return(NA_integer_)
  lo <- max(1L, ann_col - flank)
  high_streak <- 0L; high_seen <- FALSE
  for (i in seq.int(ann_col, lo, by = -1L)) {
    if (cw[i] >= T_high) {
      high_streak <- high_streak + 1L
      if (high_streak >= high_sustain) high_seen <- TRUE
    } else {
      high_streak <- 0L
    }
    if (high_seen && cw[i] < T_low) return(as.integer(i + 1L))
  }
  NA_integer_
}

detect_3p_change_point <- function(cw, ann_col, flank, T_high, T_low,
                                   high_sustain = 2L) {
  n <- length(cw)
  if (is.na(ann_col) || ann_col < 1L || ann_col > n) return(NA_integer_)
  hi <- min(n, ann_col + flank)
  high_streak <- 0L; high_seen <- FALSE
  for (i in seq.int(ann_col, hi)) {
    if (cw[i] >= T_high) {
      high_streak <- high_streak + 1L
      if (high_streak >= high_sustain) high_seen <- TRUE
    } else {
      high_streak <- 0L
    }
    if (high_seen && cw[i] < T_low) return(as.integer(i - 1L))
  }
  NA_integer_
}

revcomp_str <- function(s) {
  as.character(reverseComplement(DNAString(s)))
}

# ---------- main ----------

meta <- read.table(opt$members_tsv, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "", check.names = FALSE)

# Run MAFFT
tmp_aln <- tempfile(fileext = ".aln")
on.exit(unlink(tmp_aln), add = TRUE)
ret <- system2("mafft",
               args = c("--auto", "--thread", opt$threads, "--quiet",
                        opt$input_fasta),
               stdout = tmp_aln, stderr = NULL)
if (ret != 0L || !file.exists(tmp_aln) || file.size(tmp_aln) == 0L) {
  # Emit empty result with header so caller sees no rows (no fallback applied)
  writeLines("name\tcorrected_g\tmotif_ok", opt$output_tsv)
  quit(status = 0)
}

aln <- readDNAStringSet(tmp_aln)
if (length(aln) == 0L) {
  writeLines("name\tcorrected_g\tmotif_ok", opt$output_tsv)
  quit(status = 0)
}

# Re-order MAFFT rows to match input metadata order (MAFFT preserves
# order on --auto, but be defensive about whitespace in headers).
aln_names <- gsub(" .*", "", names(aln))
perm <- match(meta$name, aln_names)
if (any(is.na(perm))) {
  writeLines("name\tcorrected_g\tmotif_ok", opt$output_tsv)
  quit(status = 0)
}
aln <- aln[perm]
mat <- toupper(as.matrix(aln))
n_members <- nrow(mat); n_cols <- ncol(mat)

# Map per-row body boundaries from raw 1-based positions to MSA columns
non_gap_mat <- mat != "-"
cum_mat     <- t(apply(non_gap_mat, 1L, cumsum))

# body_start_raw / body_end_raw are in the per-row extended sequence,
# 1-based, in biological orientation.
body_start_raw <- as.integer(meta$body_start_raw)
body_end_raw   <- as.integer(meta$body_end_raw)
start_target <- matrix(body_start_raw, n_members, n_cols)
end_target   <- matrix(body_end_raw,   n_members, n_cols)
start_cond <- cum_mat >= start_target
end_cond   <- cum_mat >= end_target
start_cols <- ifelse(rowSums(start_cond) > 0L,
                     max.col(start_cond, ties.method = "first"), n_cols)
end_cols   <- ifelse(rowSums(end_cond) > 0L,
                     max.col(end_cond,   ties.method = "first"), n_cols)

is_5 <- meta$role == "5LTR"
is_3 <- meta$role == "3LTR"

# Need at least one 5LTR (for ann_5_col) — otherwise no scan possible.
ann_5_col <- if (any(is_5)) as.integer(median(start_cols[is_5])) else NA_integer_
ann_3_col <- if (any(is_3)) as.integer(median(end_cols[is_3])) else NA_integer_

# Conservation per-subset
cw_5 <- if (any(is_5)) {
  bc5 <- base_counts(mat[is_5, , drop = FALSE])
  sliding_mean(conservation_profile(bc5), 7L)
} else NULL
cw_3 <- if (any(is_3)) {
  bc3 <- base_counts(mat[is_3, , drop = FALSE])
  sliding_mean(conservation_profile(bc3), 7L)
} else NULL

corr_5 <- if (!is.null(cw_5))
  detect_5p_change_point(cw_5, ann_5_col, opt$flank, 0.7, 0.4) else NA_integer_
corr_3 <- if (!is.null(cw_3))
  detect_3p_change_point(cw_3, ann_3_col, opt$flank, 0.7, 0.4) else NA_integer_

if (is.na(corr_5)) corr_5 <- ann_5_col
if (is.na(corr_3)) corr_3 <- ann_3_col

# For each row, map the relevant corrected MSA column to a per-row raw
# position, and from raw position to a genomic coordinate.
# Side per row: 5LTR -> use corr_5; 3LTR -> use corr_3.
get_raw_pos <- function(row, col) {
  if (is.na(col) || col < 1L || col > n_cols) return(NA_integer_)
  # cum_mat row is integer; we want the raw position of the FIRST
  # non-gap base at or after `col`.  If the row has a gap at col,
  # take the next non-gap (mirrors mafft_boundary_consensus's
  # implicit cum_mat semantics).
  ng <- non_gap_mat[row, ]
  if (ng[col]) {
    return(as.integer(cum_mat[row, col]))
  }
  # walk forward to next non-gap
  hits <- which(ng[col:n_cols])
  if (length(hits) == 0L) return(NA_integer_)
  fwd_col <- col + hits[1L] - 1L
  as.integer(cum_mat[row, fwd_col])
}

# 5'-side raw position semantics: position of the first LTR base
# in bio orientation -> motif at raw_pos / raw_pos+1 should be "TG".
# 3'-side raw position semantics: position of the last LTR base
# in bio orientation -> motif at raw_pos-1 / raw_pos should be "CA".
# Note: change_point_3p returns the last LTR column, hence we map
# directly; but cum_mat[row, c] is the count of non-gap bases up to
# column c — so raw_pos for the 3' end == cum_mat[row, corr_3].

# Read per-row biological-orientation raw sequences from the FASTA
ext_seqs <- readDNAStringSet(opt$input_fasta)
ext_names <- gsub(" .*", "", names(ext_seqs))
perm_in <- match(meta$name, ext_names)
if (any(is.na(perm_in))) {
  writeLines("name\tcorrected_g\tmotif_ok", opt$output_tsv)
  quit(status = 0)
}
ext_seqs <- ext_seqs[perm_in]
ext_chars <- as.character(ext_seqs)
ext_lens  <- nchar(ext_chars)

raw_to_genomic <- function(raw_pos, strand, ext_start, ext_end) {
  if (is.na(raw_pos)) return(NA_integer_)
  if (strand == "+") return(as.integer(ext_start + raw_pos - 1L))
  return(as.integer(ext_end - raw_pos + 1L))
}

motif_at_raw <- function(seq_chr, raw_pos, side) {
  L <- nchar(seq_chr)
  if (is.na(raw_pos)) return(NA_character_)
  if (side == "5") {
    if (raw_pos < 1L || raw_pos + 1L > L) return(NA_character_)
    return(substr(seq_chr, raw_pos, raw_pos + 1L))
  }
  if (raw_pos - 1L < 1L || raw_pos > L) return(NA_character_)
  substr(seq_chr, raw_pos - 1L, raw_pos)
}

# Empirically the MSA change-point lands 5–15 bp inside the conserved
# region for typical LTR clusters; ±20 bp covers both ends of that
# range without producing false motif matches outside the LTR.
# (Validated on at + g2 in tmp/msa_validation_*.tsv — see
# docs/refine_msa_rescue_plan.md §3.)
snap_window <- 20L
target_for_side <- function(side) if (side == "5") "TG" else "CA"

snap_raw <- function(seq_chr, raw_pos, side, motif) {
  if (motif == "none" || is.na(raw_pos)) return(list(raw = raw_pos, ok = NA))
  tgt <- target_for_side(side)
  m0 <- motif_at_raw(seq_chr, raw_pos, side)
  if (!is.na(m0) && identical(m0, tgt)) return(list(raw = raw_pos, ok = TRUE))
  for (d in seq_len(snap_window)) {
    for (sgn in c(-1L, 1L)) {
      cand <- raw_pos + sgn * d
      mm <- motif_at_raw(seq_chr, cand, side)
      if (!is.na(mm) && identical(mm, tgt)) {
        return(list(raw = cand, ok = TRUE))
      }
    }
  }
  list(raw = raw_pos, ok = FALSE)
}

motif_policy <- opt$boundary_motif
out_corr_g  <- rep(NA_integer_, n_members)
out_motif   <- rep(NA, n_members)

for (i in seq_len(n_members)) {
  side <- if (is_5[i]) "5" else "3"
  col  <- if (side == "5") corr_5 else corr_3
  raw  <- get_raw_pos(i, col)
  if (is.na(raw)) next
  snap <- snap_raw(ext_chars[i], raw, side, motif_policy)
  raw  <- as.integer(snap$raw)
  # Raw -> genomic (in biological orientation, raw 1 = bio start)
  g <- raw_to_genomic(raw, meta$strand[i], meta$ext_start[i], meta$ext_end[i])
  out_corr_g[i] <- g
  if (motif_policy == "none") {
    out_motif[i] <- TRUE
  } else {
    out_motif[i] <- isTRUE(snap$ok)
  }
}

write.table(
  data.frame(
    name        = meta$name,
    corrected_g = ifelse(is.na(out_corr_g), "NA", as.character(out_corr_g)),
    motif_ok    = ifelse(is.na(out_motif), "NA",
                         ifelse(out_motif, "TRUE", "FALSE")),
    stringsAsFactors = FALSE
  ),
  opt$output_tsv,
  sep = "\t", quote = FALSE, row.names = FALSE
)
