#!/usr/bin/env Rscript
# build_ltr_library.R
# Build LTR reference library and junction tag databases from DANTE_LTR GFF3 output.
# Outputs: LTR consensus FASTA, ID→lineage map, 5'UTR tags, PPT tags, TSD length map.
# Called once (on the full genome) by dante_ltr_solo before chunk-level detection.

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name   <- "--file="
script_name     <- normalizePath(
  sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)])
)
script_dir <- dirname(script_name)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-g", "--gff3"), type = "character", default = NULL,
              help = paste("Refined DANTE_LTR GFF3 (output of dante_ltr_refine).",
                           "Cluster_ID, Cluster_Size and Refinement_Status",
                           "attributes drive library member grouping and",
                           "filter-then-rank ordering.")),
  make_option(c("-s", "--reference_sequence"), type = "character", default = NULL,
              help = "Reference genome FASTA"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file prefix"),
  make_option(c("-t", "--threads"), type = "integer", default = 1L,
              help = "Threads for MAFFT [default %default]"),
  make_option(c("-f", "--flank"), type = "integer", default = 50L,
              help = paste("Flanking bp added each side for the not_evaluated",
                           "fallback consensus path (change-point boundary",
                           "detection on small clusters refinement could not",
                           "evaluate). [default %default]")),
  make_option(c("--wide_flank"), type = "integer", default = 1000L,
              help = paste("Maximum retry-flank size for the not_evaluated",
                           "fallback path when the default scan cannot locate",
                           "a boundary within +-flank bp.  Proportional to the",
                           "cluster's annotated LTR length but capped above",
                           "by --wide_flank.  Set to 0 to disable the retry.",
                           "[default %default]")),
  make_option(c("-d", "--min_cluster_size"), type = "integer", default = 6L,
              help = paste("Min cluster members required for MAFFT consensus.",
                           "Below this, the highest-rank single annotated LTR",
                           "is used instead. [default %default]")),
  make_option(c("--alignments_dir"), type = "character", default = NULL,
              help = "Optional directory in which to save per-cluster MAFFT alignments (for debugging)")
)

parser <- OptionParser(option_list = option_list,
                       usage = "usage: %prog -g REFINED_GFF3 -s FASTA -o OUTPUT [OPTIONS]")
opt <- parse_args(parser)

if (is.null(opt$gff3) || !nzchar(opt$gff3)) {
  print_help(parser)
  stop("Mandatory argument missing: -g / --gff3 (refined DANTE_LTR GFF3)",
       call. = FALSE)
}
for (arg in c("reference_sequence", "output")) {
  if (is.null(opt[[arg]])) {
    print_help(parser)
    stop("Mandatory argument missing: --", arg, call. = FALSE)
  }
}
GFF_INPUT <- opt$gff3

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
})

source(file.path(script_dir, "solo_ltr_utils.R"))

RANK_PRIORITY <- c(DLTP = 4L, DLT = 3L, DLP = 2L, DL = 1L, D = 0L)

# ---- Helpers for flank-aware boundary detection ----

# Per-column base counts in an MSA matrix. Returns a list with integer vectors
# A, C, G, T, gap (length = ncol) and the row count n.
# This matrix is computed once per cluster and reused by conservation_profile()
# and majority_from_counts().
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

# Per-column conservation c(i) = top_freq * (1 - gap_frac), fully vectorized.
conservation_profile <- function(bc) {
  non_gap_n <- bc$A + bc$C + bc$G + bc$T
  top_count <- pmax(bc$A, bc$C, bc$G, bc$T)
  top_freq  <- ifelse(non_gap_n > 0L, top_count / non_gap_n, 0)
  gap_frac  <- bc$gap / bc$n
  top_freq * (1 - gap_frac)
}

# Centered sliding mean via cumsum — O(n), no loop.
sliding_mean <- function(x, w) {
  n    <- length(x)
  half <- w %/% 2L
  cs   <- c(0, cumsum(x))
  i    <- seq_len(n)
  lo   <- pmax(1L, i - half)
  hi   <- pmin(n,  i + half)
  (cs[hi + 1L] - cs[lo]) / (hi - lo + 1L)
}

# 5' boundary change-point. Walk from inside body (ann_col) outward into the 5' flank.
# Return the corrected MSA column of the first LTR base, or NA if no clear transition.
# `high_sustain` requires N consecutive smoothed columns ≥ T_high before the
# scan is "armed"; a single spurious spike (common in small clusters and in
# nested-context flanks) therefore can't trigger the scan.
detect_5p_change_point <- function(cw, ann_col, flank, T_high, T_low,
                                   high_sustain = 2L) {
  n <- length(cw)
  if (is.na(ann_col) || ann_col < 1L || ann_col > n) return(NA_integer_)
  lo <- max(1L, ann_col - flank)
  high_streak <- 0L
  high_seen   <- FALSE
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

# 3' boundary change-point — symmetric to the 5p version.
detect_3p_change_point <- function(cw, ann_col, flank, T_high, T_low,
                                   high_sustain = 2L) {
  n <- length(cw)
  if (is.na(ann_col) || ann_col < 1L || ann_col > n) return(NA_integer_)
  hi <- min(n, ann_col + flank)
  high_streak <- 0L
  high_seen   <- FALSE
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

# Majority consensus from a precomputed base-count matrix, vectorized.
# Gaps/Ns excluded; columns with no non-gap residue are skipped.
majority_from_counts <- function(bc, col_from, col_to) {
  if (col_to < col_from) return("")
  span      <- seq.int(col_from, col_to)
  count_mat <- cbind(bc$A[span], bc$C[span], bc$G[span], bc$T[span])
  non_gap_n <- rowSums(count_mat)
  maj_idx   <- max.col(count_mat, ties.method = "first")
  chars     <- c("A", "C", "G", "T")[maj_idx]
  chars[non_gap_n == 0L] <- NA_character_
  paste(chars[!is.na(chars)], collapse = "")
}

# ---- Flank-aware boundary-detecting consensus (not_evaluated fallback) ----
# Used ONLY for clusters refinement could not evaluate (Refinement_Status =
# not_evaluated for all members — typically because the cluster was too
# small to qualify for refinement).  All other clusters use refinement's
# validated coords directly via mafft_simple_consensus().
#
# extended_seqs : DNAStringSet of LTR ± flank (biological orientation)
# body_starts_raw, body_ends_raw : per-member 1-based positions in the raw extended
#   sequence marking the annotated LTR body (inclusive).
# Returns list(consensus = DNAStringSet | NULL, qc = data.frame_row, status = "ok"|"fallback").
mafft_changepoint_consensus <- function(extended_seqs, body_starts_raw, body_ends_raw,
                                        flank, threads = 1L,
                                        save_aln_path = NULL,
                                        ltr_id = "LTR_unknown",
                                        lineage = NA_character_,
                                        conservation_window = 7L,
                                        T_high = 0.7, T_low = 0.4,
                                        high_sustain = 2L,
                                        min_consensus_len = 50L,
                                        max_length_shrink = 0.5,
                                        max_length_grow   = 2.0,
                                        roles = NULL,
                                        timing_env = NULL) {

  add_time <- function(slot, dt) {
    if (!is.null(timing_env)) {
      cur <- timing_env[[slot]]
      timing_env[[slot]] <- (if (is.null(cur)) 0 else cur) + dt
    }
  }

  tmp_in  <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".fasta")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  writeXStringSet(extended_seqs, tmp_in)
  t0 <- proc.time()[["elapsed"]]
  ret <- system(
    paste("mafft --auto --thread", threads, "--quiet", tmp_in, ">", tmp_out, "2>/dev/null")
  )
  add_time("mafft", proc.time()[["elapsed"]] - t0)

  n_5ltr_in <- if (!is.null(roles)) sum(roles == "5LTR") else length(extended_seqs)
  n_3ltr_in <- if (!is.null(roles)) sum(roles == "3LTR") else 0L
  fallback_qc <- function() data.frame(
    ltr_id = ltr_id, Final_Classification = lineage,
    n_members = length(extended_seqs),
    n_5ltr = n_5ltr_in, n_3ltr = n_3ltr_in,
    flank = as.integer(flank),
    annotated_5_col = NA_integer_, corrected_5_col = NA_integer_,
    shift_5 = NA_integer_, detected_5 = FALSE,
    annotated_3_col = NA_integer_, corrected_3_col = NA_integer_,
    shift_3 = NA_integer_, detected_3 = FALSE,
    consensus_length = NA_integer_,
    median_annot_body_len = as.integer(median(body_ends_raw - body_starts_raw + 1L)),
    LibraryConfidence = "unrefined",
    stringsAsFactors = FALSE
  )

  if (ret != 0L || !file.exists(tmp_out) || file.size(tmp_out) == 0L) {
    return(list(consensus = NULL, qc = fallback_qc(), status = "fallback"))
  }

  if (!is.null(save_aln_path)) {
    try(file.copy(tmp_out, save_aln_path, overwrite = TRUE), silent = TRUE)
  }

  t0 <- proc.time()[["elapsed"]]
  aln <- readDNAStringSet(tmp_out)
  if (length(aln) == 0L) {
    return(list(consensus = NULL, qc = fallback_qc(), status = "fallback"))
  }

  # Reorder aligned rows to match input order so body_starts_raw / body_ends_raw align.
  in_names  <- names(extended_seqs)
  aln_names <- gsub(" .*", "", names(aln))
  perm <- match(gsub(" .*", "", in_names), aln_names)
  if (any(is.na(perm))) {
    return(list(consensus = NULL, qc = fallback_qc(), status = "fallback"))
  }
  aln <- aln[perm]

  mat <- toupper(as.matrix(aln))
  add_time("parse_aln", proc.time()[["elapsed"]] - t0)
  n_members <- nrow(mat)
  n_cols    <- ncol(mat)

  # Vectorized position mapping: one rowCumsums pass on the non-gap indicator,
  # then look up both start and end via max.col on the comparison matrix.
  t0 <- proc.time()[["elapsed"]]
  non_gap_mat <- mat != "-"
  cum_mat     <- t(apply(non_gap_mat, 1L, cumsum))
  start_target <- matrix(body_starts_raw, n_members, n_cols)
  end_target   <- matrix(body_ends_raw,   n_members, n_cols)
  start_cond <- cum_mat >= start_target
  end_cond   <- cum_mat >= end_target
  start_cols <- ifelse(rowSums(start_cond) > 0L,
                       max.col(start_cond, ties.method = "first"), n_cols)
  end_cols   <- ifelse(rowSums(end_cond) > 0L,
                       max.col(end_cond,   ties.method = "first"), n_cols)
  add_time("position_map", proc.time()[["elapsed"]] - t0)

  # Per-subset row masks. When no roles are supplied, treat every row as 5'LTR
  # and leave the 3' boundary uncorrected (legacy behaviour).
  is_5 <- if (!is.null(roles)) roles == "5LTR" else rep(TRUE, n_members)
  is_3 <- if (!is.null(roles)) roles == "3LTR" else rep(FALSE, n_members)

  # Annotated column medians — take each from the subset whose flank drives
  # the corresponding change-point scan, so they are consistent.
  ann_5_col <- as.integer(median(start_cols[is_5]))
  ann_3_col <- if (any(is_3)) as.integer(median(end_cols[is_3]))
               else          as.integer(median(end_cols))

  # Conservation: full matrix for the consensus vote; per-subset matrices
  # for the change-point scans.
  t0 <- proc.time()[["elapsed"]]
  bc_all <- base_counts(mat)
  bc_5   <- base_counts(mat[is_5, , drop = FALSE])
  cw_5   <- sliding_mean(conservation_profile(bc_5), conservation_window)
  cw_3   <- if (any(is_3)) {
    sliding_mean(conservation_profile(
      base_counts(mat[is_3, , drop = FALSE])
    ), conservation_window)
  } else NULL
  add_time("conservation", proc.time()[["elapsed"]] - t0)

  t0 <- proc.time()[["elapsed"]]
  corrected_5_raw <- detect_5p_change_point(
    cw_5, ann_5_col, flank, T_high, T_low, high_sustain = high_sustain
  )
  corrected_5 <- if (is.na(corrected_5_raw)) ann_5_col else corrected_5_raw

  corrected_3_raw <- if (!is.null(cw_3)) {
    detect_3p_change_point(
      cw_3, ann_3_col, flank, T_high, T_low, high_sustain = high_sustain
    )
  } else NA_integer_
  corrected_3 <- if (is.na(corrected_3_raw)) ann_3_col else corrected_3_raw

  # Sanity: total body length must stay within [max_length_shrink,
  # max_length_grow] × the annotated median body length.  Out-of-range
  # corrections revert to the annotated boundary on both sides.
  median_body_len <- as.integer(median(body_ends_raw - body_starts_raw + 1L))
  new_len <- corrected_3 - corrected_5 + 1L
  if (is.na(new_len) ||
      new_len < max_length_shrink * median_body_len ||
      new_len > max_length_grow   * median_body_len) {
    corrected_5 <- ann_5_col
    corrected_3 <- ann_3_col
  }
  add_time("change_point", proc.time()[["elapsed"]] - t0)

  t0 <- proc.time()[["elapsed"]]
  consensus <- majority_from_counts(bc_all, corrected_5, corrected_3)
  add_time("consensus_build", proc.time()[["elapsed"]] - t0)

  qc <- data.frame(
    ltr_id = ltr_id, Final_Classification = lineage,
    n_members = n_members,
    n_5ltr = sum(is_5), n_3ltr = sum(is_3),
    flank = as.integer(flank),
    annotated_5_col = ann_5_col,
    corrected_5_col = as.integer(corrected_5),
    shift_5 = as.integer(corrected_5 - ann_5_col),
    detected_5 = !is.na(corrected_5_raw),
    annotated_3_col = ann_3_col,
    corrected_3_col = as.integer(corrected_3),
    shift_3 = as.integer(corrected_3 - ann_3_col),
    detected_3 = if (!is.null(cw_3)) !is.na(corrected_3_raw) else FALSE,
    consensus_length = as.integer(nchar(consensus)),
    median_annot_body_len = median_body_len,
    LibraryConfidence = "unrefined",
    stringsAsFactors = FALSE
  )

  if (nchar(consensus) < min_consensus_len) {
    return(list(consensus = NULL, qc = qc, status = "fallback"))
  }

  list(
    consensus = DNAStringSet(setNames(consensus, ltr_id)),
    qc = qc,
    status = "ok"
  )
}

# ---- Refined-mode helpers ----
# Used when --refined_gff3 is supplied.  The cluster members are already
# at correct boundaries; build the consensus by majority vote over a
# straightforward MAFFT alignment of the (validated) bodies.

mafft_simple_consensus <- function(body_seqs, threads = 1L,
                                   ltr_id = "LTR_unknown",
                                   save_aln_path = NULL,
                                   min_consensus_len = 50L) {
  tmp_in  <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".fasta")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  writeXStringSet(body_seqs, tmp_in)
  ret <- system(
    paste("mafft --auto --thread", threads, "--quiet",
          tmp_in, ">", tmp_out, "2>/dev/null")
  )
  if (ret != 0L || !file.exists(tmp_out) || file.size(tmp_out) == 0L) {
    return(list(consensus = NULL, status = "fallback"))
  }
  if (!is.null(save_aln_path)) {
    try(file.copy(tmp_out, save_aln_path, overwrite = TRUE), silent = TRUE)
  }
  aln <- readDNAStringSet(tmp_out)
  if (length(aln) == 0L) {
    return(list(consensus = NULL, status = "fallback"))
  }
  mat <- toupper(as.matrix(aln))
  bc  <- base_counts(mat)
  cons_str <- majority_from_counts(bc, 1L, ncol(mat))
  if (nchar(cons_str) < min_consensus_len) {
    return(list(consensus = NULL, status = "fallback"))
  }
  list(consensus = DNAStringSet(setNames(cons_str, ltr_id)),
       status    = "ok")
}

# ============================================================
# MAIN
# ============================================================

if (!is.null(opt$alignments_dir)) {
  dir.create(opt$alignments_dir, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("Reading refined GFF3: %s\n", GFF_INPUT))
gff <- import(GFF_INPUT, format = "gff3")

cat("Reading genome ...\n")
s        <- readDNAStringSet(opt$reference_sequence)
names(s) <- gsub("\\s.+", "", names(s))
sl_vec   <- setNames(as.integer(nchar(s)), names(s))

# ---- Extract LTR features ----
ltr_feat <- gff[gff$type == "long_terminal_repeat"]
te_feat  <- gff[gff$type == "transposable_element"]

if (length(ltr_feat) == 0L) {
  stop("No long_terminal_repeat features in input GFF3. ",
       "Run dante_ltr first and use its output.", call. = FALSE)
}
cat(sprintf("Found %d long_terminal_repeat features\n", length(ltr_feat)))

# Parent metadata lookup
te_id    <- as.character(te_feat$ID)
te_rank  <- setNames(as.character(te_feat$Rank),               te_id)
te_class <- setNames(as.character(te_feat$Final_Classification), te_id)

ltr_parent <- as.character(unlist(ltr_feat$Parent))
ltr_rank   <- te_rank[ltr_parent]
ltr_class  <- te_class[ltr_parent]
ltr_ltr    <- as.character(ltr_feat$LTR)   # "5LTR" or "3LTR"
ltr_rank[is.na(ltr_rank)] <- "DL"

# Per-LTR refinement attributes — required.  The script assumes the input
# GFF3 came from dante_ltr_refine, so every LTR feature should carry
# Refinement_Status and Cluster_ID.  Members of small clusters that
# refinement could not evaluate carry Refinement_Status=not_evaluated and
# still have a Cluster_ID; they go to the change-point fallback path.
ltr_status <- if (!is.null(ltr_feat$Refinement_Status)) {
  as.character(ltr_feat$Refinement_Status)
} else {
  rep(NA_character_, length(ltr_feat))
}
ltr_cluster_id <- if (!is.null(ltr_feat$Cluster_ID)) {
  as.character(ltr_feat$Cluster_ID)
} else {
  rep(NA_character_, length(ltr_feat))
}

if (all(is.na(ltr_status) | ltr_status == "")) {
  stop("Input GFF3 has no Refinement_Status attributes — run ",
       "dante_ltr_refine first (or call dante_ltr_solo, which auto-",
       "refines if needed).", call. = FALSE)
}

# Refinement_Status values:
#   confirmed | refined     -> validated (Pool A)
#   unresolved              -> evaluated but not validated (Pool B inside
#                              evaluated cluster)
#   not_evaluated | <empty> -> cluster was too small for refinement
#                              (fallback path: change-point consensus)
ltr_validated <- ltr_status %in% c("confirmed", "refined")

# ---- Extract body + extended sequences for ALL LTRs (5' and 3') ----
# Refinement assigns the same Cluster_ID to a 5'LTR and its sibling
# 3'LTR (both LTRs of the same element), so grouping by Cluster_ID
# naturally pulls in both roles per cluster.
seqlengths(ltr_feat) <- sl_vec[seqlevels(ltr_feat)]

ltr5_mask <- ltr_ltr == "5LTR" & !is.na(ltr_ltr)
if (!any(ltr5_mask)) {
  stop("No 5'LTR features in input GFF3 — cannot build library.", call. = FALSE)
}

ltr_priority_all <- RANK_PRIORITY[ltr_rank]
ltr_priority_all[is.na(ltr_priority_all)] <- 0L

cat(sprintf("Found %d 5'LTRs and %d 3'LTRs; refinement assigned Cluster_ID to %d / %d LTRs\n",
            sum(ltr5_mask), sum(ltr_ltr == "3LTR" & !is.na(ltr_ltr)),
            sum(!is.na(ltr_cluster_id) & ltr_cluster_id != ""),
            length(ltr_feat)))

# Body-only sequences (used for MAFFT consensus building and the
# single-member fallback when MAFFT cannot run on a cluster)
ltr_body_seqs_all <- getSeq(s, ltr_feat)

# Extended ranges with ± flank, strand-aware via getSeq
F_flank  <- opt$flank
sn_v     <- as.character(seqnames(ltr_feat))
strand_v <- as.character(strand(ltr_feat))
sl_per   <- sl_vec[sn_v]
ext_start <- pmax(1L,     start(ltr_feat) - F_flank)
ext_end   <- pmin(sl_per, end(ltr_feat)   + F_flank)
ltr_feat_ext <- GRanges(sn_v, IRanges(ext_start, ext_end), strand = strand_v)
seqlengths(ltr_feat_ext) <- sl_vec[seqlevels(ltr_feat_ext)]
ltr_ext_seqs_all <- getSeq(s, ltr_feat_ext)

# Per-member offsets in biological orientation
gleft          <- start(ltr_feat) - ext_start    # genomic-left bases kept
gright         <- ext_end - end(ltr_feat)        # genomic-right bases kept
bio_5flank_len <- ifelse(strand_v == "-", gright, gleft)
bio_body_len   <- width(ltr_feat)

# Annotated boundary positions in the raw extended sequence (1-based)
ann_5_pos_raw <- as.integer(bio_5flank_len + 1L)
ann_3_pos_raw <- as.integer(bio_5flank_len + bio_body_len)

safe_coords <- paste0(seqnames(ltr_feat), "_", start(ltr_feat), "_", end(ltr_feat))
names(ltr_body_seqs_all) <- safe_coords
names(ltr_ext_seqs_all)  <- safe_coords

idx_5 <- which(ltr_ltr == "5LTR" & !is.na(ltr_ltr))
idx_3 <- which(ltr_ltr == "3LTR" & !is.na(ltr_ltr))

# ---- TSD length map ----
cat("Computing TSD length map ...\n")
tsd_map <- extract_tsd_length_map(gff)
tsd_map_path <- paste0(opt$output, "_tsd_length_map.tsv")
if (!is.null(tsd_map) && length(tsd_map) > 0L) {
  write.table(data.frame(lineage = names(tsd_map), tsd_length = as.integer(tsd_map),
                         stringsAsFactors = FALSE),
              tsd_map_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  %d lineages with TSD data\n", length(tsd_map)))
} else {
  cat("  No TSD data found; fallback scan (4-6 bp) will be used\n")
  write.table(data.frame(lineage = character(0L), tsd_length = integer(0L)),
              tsd_map_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Helper: re-extract a set of LTRs with a custom flank (used by the
# wide-flank retry below).  Returns a list(seqs = DNAStringSet,
# b5 = integer, b3 = integer) with the same structure the main
# extraction produces but for an arbitrary flank size.
extract_extended <- function(indices, flank_bp) {
  sn_x     <- sn_v[indices]
  strand_x <- strand_v[indices]
  sl_x     <- sl_per[indices]
  e_start <- pmax(1L, start(ltr_feat)[indices] - flank_bp)
  e_end   <- pmin(sl_x, end(ltr_feat)[indices] + flank_bp)
  ext_gr  <- GRanges(sn_x, IRanges(e_start, e_end), strand = strand_x)
  seqlengths(ext_gr) <- sl_vec[seqlevels(ext_gr)]
  seqs <- getSeq(s, ext_gr)
  gleft_x  <- start(ltr_feat)[indices] - e_start
  gright_x <- e_end - end(ltr_feat)[indices]
  bio5_x   <- ifelse(strand_x == "-", gright_x, gleft_x)
  body_x   <- width(ltr_feat)[indices]
  list(
    seqs = seqs,
    b5   = as.integer(bio5_x + 1L),
    b3   = as.integer(bio5_x + body_x)
  )
}

# ---- Build per-cluster consensus from refinement's clusters ----
# Members are grouped by Cluster_ID assigned by dante_ltr_refine.
# The path each cluster takes depends on the mix of Refinement_Status
# values among its members:
#
#   Pool A  : Refinement_Status in {confirmed, refined}  -> validated
#   Pool B  : Refinement_Status == "unresolved"          -> evaluated but
#                                                           not validated
#   None    : all Refinement_Status == "not_evaluated"   -> cluster was
#                                                           too small for
#                                                           refinement
#
#   Pool A non-empty               -> mafft_simple_consensus on Pool A
#                                     LibraryConfidence = "validated"
#   Pool A empty, Pool B non-empty -> mafft_simple_consensus on Pool B
#                                     LibraryConfidence = "mixed"
#   All not_evaluated              -> mafft_changepoint_consensus on all
#                                     members (with flank/wide_flank
#                                     boundary detection)
#                                     LibraryConfidence = "unrefined"
#
# Within Pool A or Pool B, members are kept in lineage order (as they
# appear in the GFF3); MAFFT then builds the consensus by majority vote.
# Single-member clusters fall back to the body sequence directly.

# Identify clusters and their member indices (into the all-LTR vectors)
cluster_keys <- ifelse(is.na(ltr_cluster_id) | ltr_cluster_id == "",
                       paste0("__nocluster__", seq_along(ltr_feat)),
                       ltr_cluster_id)
cluster_member_idx <- split(seq_along(ltr_feat), cluster_keys)

cat(sprintf("Building consensus for %d clusters ...\n",
            length(cluster_member_idx)))

all_consensus <- DNAStringSet()
id_to_lineage <- character(0L)
id_lib_conf   <- character(0L)    # validated | mixed | unrefined
id_n_total    <- integer(0L)
id_cluster    <- character(0L)    # original Cluster_ID (or empty for unclustered)
ltr_counter   <- 0L
qc_rows       <- list()

# Accumulate per-step wall-clock times across all clusters for profiling.
timing_env <- new.env(parent = emptyenv())
for (slot in c("mafft", "parse_aln", "position_map",
               "conservation", "change_point", "consensus_build")) {
  timing_env[[slot]] <- 0
}
t_main_start <- proc.time()[["elapsed"]]

# Per-cluster lineage (inherited from any member; all members share it
# unless the GFF3 is malformed).
for (cluster_key in names(cluster_member_idx)) {
  member_idx <- cluster_member_idx[[cluster_key]]
  if (length(member_idx) == 0L) next

  ltr_counter <- ltr_counter + 1L
  ltr_id <- sprintf("LTR_%06d", ltr_counter)

  # Lineage: take the modal lineage of the cluster's members.
  cluster_lineages <- ltr_class[member_idx]
  cluster_lineages <- cluster_lineages[!is.na(cluster_lineages)]
  lin <- if (length(cluster_lineages) > 0L) {
    names(sort(table(cluster_lineages), decreasing = TRUE))[1L]
  } else NA_character_
  if (is.na(lin)) next

  # Pool partition by Refinement_Status
  status_in_cluster <- ltr_status[member_idx]
  is_pool_a <- status_in_cluster %in% c("confirmed", "refined")
  is_pool_b <- status_in_cluster == "unresolved"
  is_unrefined <- !(is_pool_a | is_pool_b)  # not_evaluated / NA / ""

  cluster_lib_conf <- "unrefined"
  use_simple   <- FALSE      # refinement-validated path (mafft_simple_consensus)
  use_changept <- FALSE      # not_evaluated fallback (mafft_changepoint_consensus)
  pool_idx     <- integer(0L)

  if (any(is_pool_a)) {
    cluster_lib_conf <- "validated"
    pool_idx <- member_idx[is_pool_a]
    use_simple <- TRUE
  } else if (any(is_pool_b)) {
    cluster_lib_conf <- "mixed"
    pool_idx <- member_idx[is_pool_b]
    use_simple <- TRUE
  } else {
    cluster_lib_conf <- "unrefined"
    pool_idx <- member_idx
    use_changept <- TRUE
  }

  # Filter-then-rank ordering within the chosen pool: by Rank (DLTP > ...
  # > D), with NA priorities pushed to the back.
  pool_priority <- ltr_priority_all[pool_idx]
  pool_priority[is.na(pool_priority)] <- 0L
  pool_idx <- pool_idx[order(-pool_priority)]

  cons_done <- FALSE
  cons      <- NULL
  qc_row    <- NULL

  if (use_simple && length(pool_idx) >= opt$min_cluster_size) {
    aln_path <- if (!is.null(opt$alignments_dir)) {
      file.path(opt$alignments_dir, paste0(ltr_id, ".aln.fasta"))
    } else NULL
    body_seqs <- ltr_body_seqs_all[pool_idx]
    names(body_seqs) <- safe_coords[pool_idx]
    t0 <- proc.time()[["elapsed"]]
    res <- mafft_simple_consensus(
      body_seqs     = body_seqs,
      threads       = opt$threads,
      ltr_id        = ltr_id,
      save_aln_path = aln_path
    )
    timing_env$mafft <- timing_env$mafft + (proc.time()[["elapsed"]] - t0)
    if (res$status == "ok" && !is.null(res$consensus)) {
      cons      <- res$consensus
      cons_done <- TRUE
    }
    qc_row <- data.frame(
      ltr_id = ltr_id, Final_Classification = lin,
      n_members = length(pool_idx),
      n_5ltr = sum(ltr_ltr[pool_idx] == "5LTR", na.rm = TRUE),
      n_3ltr = sum(ltr_ltr[pool_idx] == "3LTR", na.rm = TRUE),
      flank = NA_integer_,
      annotated_5_col = NA_integer_, corrected_5_col = NA_integer_,
      shift_5 = NA_integer_, detected_5 = FALSE,
      annotated_3_col = NA_integer_, corrected_3_col = NA_integer_,
      shift_3 = NA_integer_, detected_3 = FALSE,
      consensus_length = if (cons_done) as.integer(width(cons)[1L]) else NA_integer_,
      median_annot_body_len = as.integer(median(width(ltr_body_seqs_all[pool_idx]))),
      LibraryConfidence = cluster_lib_conf,
      stringsAsFactors = FALSE
    )
  } else if (use_changept && length(pool_idx) >= opt$min_cluster_size) {
    aln_path <- if (!is.null(opt$alignments_dir)) {
      file.path(opt$alignments_dir, paste0(ltr_id, ".aln.fasta"))
    } else NULL
    joint_roles <- ltr_ltr[pool_idx]
    res <- mafft_changepoint_consensus(
      extended_seqs    = ltr_ext_seqs_all[pool_idx],
      body_starts_raw  = ann_5_pos_raw[pool_idx],
      body_ends_raw    = ann_3_pos_raw[pool_idx],
      flank            = F_flank,
      threads          = opt$threads,
      save_aln_path    = aln_path,
      ltr_id           = ltr_id,
      lineage          = lin,
      roles            = joint_roles,
      timing_env       = timing_env
    )
    if (res$status == "ok" &&
        opt$wide_flank > F_flank &&
        (!res$qc$detected_5 || !res$qc$detected_3)) {
      body_len <- as.integer(res$qc$median_annot_body_len)
      wide_flank_used <- min(opt$wide_flank, max(500L, body_len))
      ext <- extract_extended(pool_idx, wide_flank_used)
      names(ext$seqs) <- names(ltr_ext_seqs_all[pool_idx])
      res2 <- mafft_changepoint_consensus(
        extended_seqs   = ext$seqs,
        body_starts_raw = ext$b5,
        body_ends_raw   = ext$b3,
        flank           = wide_flank_used,
        threads         = opt$threads,
        save_aln_path   = aln_path,
        ltr_id          = ltr_id,
        lineage         = lin,
        roles           = joint_roles,
        timing_env      = timing_env
      )
      if (res2$status == "ok") {
        old_det <- sum(res$qc$detected_5, res$qc$detected_3, na.rm = TRUE)
        new_det <- sum(res2$qc$detected_5, res2$qc$detected_3, na.rm = TRUE)
        if (new_det > old_det) {
          cat(sprintf("    %s: wide-flank retry detected %d -> %d sides\n",
                      ltr_id, old_det, new_det))
          res <- res2
        }
      }
    }
    qc_row <- res$qc
    qc_row$LibraryConfidence <- cluster_lib_conf
    if (res$status == "ok" && !is.null(res$consensus)) {
      cons      <- res$consensus
      cons_done <- TRUE
    }
  }

  if (!cons_done) {
    # Single-member fallback: highest-rank annotated body in the pool.
    pool_ranks <- ltr_priority_all[pool_idx]
    pool_ranks[is.na(pool_ranks)] <- 0L
    best <- pool_idx[which.max(pool_ranks)]
    cons <- ltr_body_seqs_all[best]
    if (is.null(qc_row)) {
      qc_row <- data.frame(
        ltr_id = ltr_id, Final_Classification = lin,
        n_members = length(pool_idx),
        n_5ltr = sum(ltr_ltr[pool_idx] == "5LTR", na.rm = TRUE),
        n_3ltr = sum(ltr_ltr[pool_idx] == "3LTR", na.rm = TRUE),
        flank = NA_integer_,
        annotated_5_col = NA_integer_, corrected_5_col = NA_integer_,
        shift_5 = NA_integer_, detected_5 = FALSE,
        annotated_3_col = NA_integer_, corrected_3_col = NA_integer_,
        shift_3 = NA_integer_, detected_3 = FALSE,
        consensus_length = as.integer(width(cons)[1L]),
        median_annot_body_len = as.integer(median(width(ltr_body_seqs_all[pool_idx]))),
        LibraryConfidence = cluster_lib_conf,
        stringsAsFactors = FALSE
      )
    }
  }

  qc_rows[[length(qc_rows) + 1L]] <- qc_row
  names(cons)   <- ltr_id
  all_consensus <- c(all_consensus, cons)
  id_to_lineage[ltr_id] <- lin
  id_lib_conf <- c(id_lib_conf, cluster_lib_conf)
  id_n_total  <- c(id_n_total, length(member_idx))
  id_cluster  <- c(id_cluster,
                   if (startsWith(cluster_key, "__nocluster__")) "" else cluster_key)
}

t_main_total <- proc.time()[["elapsed"]] - t_main_start
cat("\n[timing] consensus/cluster loop breakdown (wall-clock s):\n")
cat(sprintf("  mafft (extern)  : %6.2f\n",  timing_env$mafft))
cat(sprintf("  parse alignment : %6.2f\n",  timing_env$parse_aln))
cat(sprintf("  position map    : %6.2f\n",  timing_env$position_map))
cat(sprintf("  conservation    : %6.2f\n",  timing_env$conservation))
cat(sprintf("  change-point    : %6.2f\n",  timing_env$change_point))
cat(sprintf("  consensus build : %6.2f\n",  timing_env$consensus_build))
cat(sprintf("  --- total loop  : %6.2f\n",  t_main_total))
cat(sprintf("LibraryConfidence summary : validated=%d  mixed=%d  unrefined=%d\n",
            sum(id_lib_conf == "validated"),
            sum(id_lib_conf == "mixed"),
            sum(id_lib_conf == "unrefined")))

if (length(all_consensus) == 0L) {
  stop("No LTR consensus sequences built. Check input GFF3.", call. = FALSE)
}

# Write library and ID map
ltr_lib_path <- paste0(opt$output, "_LTR_library.fasta")
writeXStringSet(all_consensus, ltr_lib_path)

ltr_map_path <- paste0(opt$output, "_LTR_library_map.tsv")
map_df <- data.frame(
  ltr_id               = names(id_to_lineage),
  Final_Classification = as.character(id_to_lineage),
  LibraryConfidence    = as.character(id_lib_conf),
  Cluster_ID           = as.character(id_cluster),
  Cluster_Size         = as.integer(id_n_total),
  stringsAsFactors     = FALSE
)
write.table(map_df, ltr_map_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

# Write boundary QC table
qc_path <- paste0(opt$output, "_LTR_library_boundary_qc.tsv")
if (length(qc_rows) > 0L) {
  qc_df <- do.call(rbind, qc_rows)
  write.table(qc_df, qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(data.frame(
    ltr_id = character(0L), Final_Classification = character(0L),
    n_members = integer(0L), n_5ltr = integer(0L), n_3ltr = integer(0L),
    flank = integer(0L),
    annotated_5_col = integer(0L), corrected_5_col = integer(0L),
    shift_5 = integer(0L), detected_5 = logical(0L),
    annotated_3_col = integer(0L), corrected_3_col = integer(0L),
    shift_3 = integer(0L), detected_3 = logical(0L),
    consensus_length = integer(0L), median_annot_body_len = integer(0L),
    LibraryConfidence = character(0L)
  ), qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("LTR library: %d consensus sequences written to %s\n",
            length(all_consensus), ltr_lib_path))

# ---- Cluster / boundary-adjustment summary -------------------------------
# Compact, human-readable report.  All numbers are also in qc_path.
if (length(qc_rows) > 0L) {
  qc_summary <- qc_df
  total_clusters <- nrow(qc_summary)
  size_bins <- cut(
    qc_summary$n_members,
    breaks = c(0, 1, 2, 5, 10, 20, Inf),
    labels = c("=1", "=2", "3-5", "6-10", "11-20", "21+"),
    right  = TRUE
  )
  cat("\n[summary] cluster / library report\n")
  cat(sprintf("  total clusters        : %d\n", total_clusters))
  cat("  size distribution     :\n")
  print(table(size_bins))
  if ("LibraryConfidence" %in% colnames(qc_summary)) {
    lc_tab <- table(qc_summary$LibraryConfidence)
    for (lvl in c("validated", "mixed", "unrefined")) {
      n <- if (lvl %in% names(lc_tab)) as.integer(lc_tab[[lvl]]) else 0L
      cat(sprintf("  LibraryConfidence %-9s : %d\n", lvl, n))
    }
  }
  # Change-point boundary detection stats (only meaningful for the
  # not_evaluated fallback path; LibraryConfidence == "unrefined")
  msa <- qc_summary[!is.na(qc_summary$corrected_5_col), ]
  if (nrow(msa) > 0L) {
    cat(sprintf("\n  not_evaluated fallback (change-point) clusters : %d\n",
                nrow(msa)))
    n_det5 <- sum(msa$detected_5)
    n_det3 <- sum(msa$detected_3)
    n_det_any  <- sum(msa$detected_5 | msa$detected_3)
    n_det_both <- sum(msa$detected_5 & msa$detected_3)
    n_blocked  <- sum(!msa$detected_5 | !msa$detected_3)
    cat(sprintf("  boundary detected     : 5'=%d, 3'=%d, both=%d, either=%d\n",
                n_det5, n_det3, n_det_both, n_det_any))
    cat(sprintf("  boundary blocked*     : %d  (* no transition found on at least one side)\n",
                n_blocked))
    if (n_det5 > 0L) {
      s5 <- msa$shift_5[msa$detected_5]
      cat(sprintf("  shift_5  min/med/max  : %d / %d / %d bp  (among %d detected)\n",
                  min(s5), as.integer(median(s5)), max(s5), length(s5)))
    }
    if (n_det3 > 0L) {
      s3 <- msa$shift_3[msa$detected_3]
      cat(sprintf("  shift_3  min/med/max  : %d / %d / %d bp  (among %d detected)\n",
                  min(s3), as.integer(median(s3)), max(s3), length(s3)))
    }
  }
  if (!is.null(opt$alignments_dir) && dir.exists(opt$alignments_dir)) {
    cat(sprintf("  per-cluster MSAs      : %s  (%d .aln.fasta files)\n",
                opt$alignments_dir,
                length(list.files(opt$alignments_dir,
                                  pattern = "\\.aln\\.fasta$"))))
  }
  cat(sprintf("  per-cluster QC TSV    : %s\n", qc_path))
}

# ---- HTML visual report --------------------------------------------------
# Only rendered when alignments were kept (report's per-cluster panels rely
# on them).  A failure here is non-fatal — the TSV is the source of truth.
if (!is.null(opt$alignments_dir) && dir.exists(opt$alignments_dir) &&
    length(qc_rows) > 0L) {
  report_script <- file.path(script_dir, "boundary_report.R")
  report_path   <- paste0(opt$output, "_LTR_library_boundary_report.html")
  rs_args <- c(report_script,
               "--qc",         qc_path,
               "--alignments", opt$alignments_dir,
               "--out",        report_path)
  cat("Rendering boundary HTML report ...\n")
  ret <- tryCatch(
    system2("Rscript", rs_args, stdout = "", stderr = ""),
    error = function(e) { cat("  warning: report render failed: ",
                              conditionMessage(e), "\n"); 1L })
  if (identical(ret, 0L) && file.exists(report_path)) {
    cat(sprintf("  boundary report      : %s\n", report_path))
  } else {
    cat("  (boundary report not produced)\n")
  }
}

# ---- 5'UTR tag database ----
cat("Building 5'UTR junction tag database ...\n")
ltr5_idx <- which(ltr_ltr == "5LTR" & !is.na(ltr_ltr))
utr5_tags <- DNAStringSet()

for (i in ltr5_idx) {
  sn_i <- as.character(seqnames(ltr_feat)[i])
  sl_i <- sl_vec[sn_i]
  s_i  <- as.character(strand(ltr_feat)[i])

  if (s_i == "+") {
    w_s <- end(ltr_feat)[i] + 1L
    w_e <- min(sl_i, end(ltr_feat)[i] + 30L)
  } else {
    w_s <- max(1L, start(ltr_feat)[i] - 30L)
    w_e <- start(ltr_feat)[i] - 1L
  }
  if (w_e < w_s) next

  tag <- getSeq(s, GRanges(sn_i, IRanges(w_s, w_e)))
  if (s_i == "-") tag <- reverseComplement(tag)
  names(tag) <- paste0("UTR5_", ltr_parent[i])
  utr5_tags <- c(utr5_tags, tag)
}

utr5_path <- paste0(opt$output, "_5UTR_tags.fasta")
writeXStringSet(utr5_tags, utr5_path)
if (length(utr5_tags) > 0L) {
  system(paste("makeblastdb -in", utr5_path, "-dbtype nucl -out", utr5_path,
               "-title UTR5_tags 2>/dev/null"))
  cat(sprintf("  %d 5'UTR tags\n", length(utr5_tags)))
} else {
  cat("  Warning: no 5'UTR tags extracted\n")
}

# ---- PPT/3'UTR tag database ----
cat("Building PPT junction tag database ...\n")
ltr3_idx  <- which(ltr_ltr == "3LTR" & !is.na(ltr_ltr))
ppt_tags  <- DNAStringSet()

for (i in ltr3_idx) {
  sn_i <- as.character(seqnames(ltr_feat)[i])
  sl_i <- sl_vec[sn_i]
  s_i  <- as.character(strand(ltr_feat)[i])

  if (s_i == "+") {
    w_s <- max(1L, start(ltr_feat)[i] - 30L)
    w_e <- start(ltr_feat)[i] - 1L
  } else {
    w_s <- end(ltr_feat)[i] + 1L
    w_e <- min(sl_i, end(ltr_feat)[i] + 30L)
  }
  if (w_e < w_s) next

  tag <- getSeq(s, GRanges(sn_i, IRanges(w_s, w_e)))
  if (s_i == "-") tag <- reverseComplement(tag)
  names(tag) <- paste0("PPT_", ltr_parent[i])
  ppt_tags <- c(ppt_tags, tag)
}

ppt_path <- paste0(opt$output, "_PPT_tags.fasta")
writeXStringSet(ppt_tags, ppt_path)
if (length(ppt_tags) > 0L) {
  system(paste("makeblastdb -in", ppt_path, "-dbtype nucl -out", ppt_path,
               "-title PPT_tags 2>/dev/null"))
  cat(sprintf("  %d PPT tags\n", length(ppt_tags)))
} else {
  cat("  Warning: no PPT tags extracted\n")
}

cat("\nLibrary build complete.\n")
cat(sprintf("  LTR library   : %s\n", ltr_lib_path))
cat(sprintf("  LTR ID map    : %s\n", ltr_map_path))
cat(sprintf("  Boundary QC   : %s\n", qc_path))
cat(sprintf("  5'UTR tags    : %s\n", utr5_path))
cat(sprintf("  PPT tags      : %s\n", ppt_path))
cat(sprintf("  TSD map       : %s\n", tsd_map_path))
