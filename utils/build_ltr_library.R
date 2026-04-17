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
              help = "DANTE_LTR GFF3 annotation file"),
  make_option(c("-s", "--reference_sequence"), type = "character", default = NULL,
              help = "Reference genome FASTA"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file prefix"),
  make_option(c("-t", "--threads"), type = "integer", default = 1L,
              help = "Threads for MMseqs2 and MAFFT [default %default]"),
  make_option(c("-f", "--flank"), type = "integer", default = 50L,
              help = "Flanking bp added each side for MSA and change-point boundary detection [default %default]"),
  make_option(c("-d", "--min_cluster_size"), type = "integer", default = 3L,
              help = "Min cluster members required for MAFFT consensus [default %default]"),
  make_option(c("--alignments_dir"), type = "character", default = NULL,
              help = "Optional directory in which to save per-cluster MAFFT alignments (for debugging)")
)

parser <- OptionParser(option_list = option_list,
                       usage = "usage: %prog -g GFF3 -s FASTA -o OUTPUT [OPTIONS]")
opt <- parse_args(parser)

for (arg in c("gff3", "reference_sequence", "output")) {
  if (is.null(opt[[arg]])) {
    print_help(parser)
    stop("Mandatory argument missing: --", arg, call. = FALSE)
  }
}

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
detect_5p_change_point <- function(cw, ann_col, flank, T_high, T_low) {
  n <- length(cw)
  if (is.na(ann_col) || ann_col < 1L || ann_col > n) return(NA_integer_)
  lo <- max(1L, ann_col - flank)
  high_seen <- FALSE
  for (i in seq.int(ann_col, lo, by = -1L)) {
    if (!high_seen && cw[i] >= T_high) high_seen <- TRUE
    if (high_seen && cw[i] < T_low) return(as.integer(i + 1L))
  }
  NA_integer_
}

# 3' boundary change-point. Walk from inside body (ann_col) outward into the 3' flank.
# Return the corrected MSA column of the last LTR base, or NA if no clear transition.
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

# ---- Flank-aware boundary-detecting consensus ----
# extended_seqs : DNAStringSet of LTR ± flank (biological orientation)
# body_starts_raw, body_ends_raw : per-member 1-based positions in the raw extended
#   sequence marking the annotated LTR body (inclusive).
# Returns list(consensus = DNAStringSet | NULL, qc = data.frame_row, status = "ok"|"fallback").
mafft_boundary_consensus <- function(extended_seqs, body_starts_raw, body_ends_raw,
                                     flank, threads = 1L,
                                     save_aln_path = NULL,
                                     ltr_id = "LTR_unknown",
                                     lineage = NA_character_,
                                     conservation_window = 5L,
                                     T_high = 0.7, T_low = 0.4,
                                     min_consensus_len = 50L,
                                     max_length_shrink = 0.5,
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
    annotated_5_col = NA_integer_, corrected_5_col = NA_integer_, shift_5 = NA_integer_,
    annotated_3_col = NA_integer_, corrected_3_col = NA_integer_, shift_3 = NA_integer_,
    consensus_length = NA_integer_,
    median_annot_body_len = as.integer(median(body_ends_raw - body_starts_raw + 1L)),
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
  corrected_5 <- detect_5p_change_point(cw_5, ann_5_col, flank, T_high, T_low)
  if (is.na(corrected_5)) corrected_5 <- ann_5_col

  corrected_3 <- if (!is.null(cw_3)) {
    detect_3p_change_point(cw_3, ann_3_col, flank, T_high, T_low)
  } else NA_integer_
  if (is.na(corrected_3)) corrected_3 <- ann_3_col

  # Sanity: total body length must stay within max_length_shrink of the annotated median
  median_body_len <- as.integer(median(body_ends_raw - body_starts_raw + 1L))
  new_len <- corrected_3 - corrected_5 + 1L
  if (is.na(new_len) || new_len < max_length_shrink * median_body_len) {
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
    annotated_5_col = ann_5_col,
    corrected_5_col = as.integer(corrected_5),
    shift_5 = as.integer(corrected_5 - ann_5_col),
    annotated_3_col = ann_3_col,
    corrected_3_col = as.integer(corrected_3),
    shift_3 = as.integer(corrected_3 - ann_3_col),
    consensus_length = as.integer(nchar(consensus)),
    median_annot_body_len = median_body_len,
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

# ---- MMseqs2 clustering (returns named vector: member -> representative) ----
mmseqs_cluster <- function(seqs, identity = 0.9, threads = 1L) {
  if (length(seqs) < 2L) return(setNames(names(seqs), names(seqs)))

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fa_in  <- file.path(tmp_dir, "input.fasta")
  cl_pfx <- file.path(tmp_dir, "cluster")
  tmp_mm <- file.path(tmp_dir, "tmp")
  dir.create(tmp_mm)

  writeXStringSet(seqs, fa_in)
  ret <- system(paste(
    "mmseqs easy-cluster", fa_in, cl_pfx, tmp_mm,
    "--min-seq-id", identity,
    "-c 0.8 --cov-mode 0",
    "--threads", threads,
    "-v 0 2>/dev/null"
  ))

  tsv <- paste0(cl_pfx, "_cluster.tsv")
  if (ret != 0L || !file.exists(tsv) || file.size(tsv) == 0L) {
    return(setNames(names(seqs), names(seqs)))  # fallback: each is its own cluster
  }

  cls <- read.table(tsv, as.is = TRUE, sep = "\t", col.names = c("rep", "member"))
  # MMseqs2 uses the first-space-truncated name; match back to full names
  name_map <- setNames(names(seqs), gsub(" .*", "", names(seqs)))
  rep_full    <- name_map[cls$rep]
  member_full <- name_map[cls$member]
  valid <- !is.na(rep_full) & !is.na(member_full)
  setNames(rep_full[valid], member_full[valid])
}

# ============================================================
# MAIN
# ============================================================

if (!is.null(opt$alignments_dir)) {
  dir.create(opt$alignments_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Reading GFF3 ...\n")
gff <- import(opt$gff3, format = "gff3")

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

# ---- Extract body + extended sequences for ALL LTRs (5' and 3') ----
# Clustering still uses 5'LTR bodies only, but the MAFFT step now also
# includes sibling 3'LTRs so the 3' boundary can be change-point-corrected.
seqlengths(ltr_feat) <- sl_vec[seqlevels(ltr_feat)]

ltr5_mask <- ltr_ltr == "5LTR" & !is.na(ltr_ltr)
if (!any(ltr5_mask)) {
  stop("No 5'LTR features in input GFF3 — cannot build library.", call. = FALSE)
}

ltr_priority_all <- RANK_PRIORITY[ltr_rank]
ltr_priority_all[is.na(ltr_priority_all)] <- 0L

cat(sprintf("Using %d 5'LTRs for clustering; %d 3'LTR siblings will join the MSA\n",
            sum(ltr5_mask), sum(ltr_ltr == "3LTR" & !is.na(ltr_ltr))))

# Body-only sequences (for MMseqs2 clustering and small-cluster fallback)
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

# Parent -> index lookup for sibling resolution
idx_5 <- which(ltr_ltr == "5LTR" & !is.na(ltr_ltr))
idx_3 <- which(ltr_ltr == "3LTR" & !is.na(ltr_ltr))
by_parent_3 <- setNames(as.list(idx_3), ltr_parent[idx_3])

# 5'LTR views used for clustering (clustering sees bodies only)
ltr_body_seqs   <- ltr_body_seqs_all[idx_5]
ltr_class_5     <- ltr_class[idx_5]
ltr_priority_5  <- ltr_priority_all[idx_5]
ltr_parent_5    <- ltr_parent[idx_5]
# Map cluster-local (5'LTR) index back into the all-LTR frame
all_idx_of_5    <- idx_5

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

# ---- Per-lineage clustering and boundary-aware consensus ----
cat("Clustering and building LTR consensus sequences ...\n")

lineages <- unique(ltr_class_5[!is.na(ltr_class_5)])
all_consensus <- DNAStringSet()
id_to_lineage <- character(0L)    # maps FASTA ID -> Final_Classification
ltr_counter   <- 0L
qc_rows       <- list()

# Accumulate per-step wall-clock times across all clusters for profiling.
timing_env <- new.env(parent = emptyenv())
for (slot in c("mmseqs", "mafft", "parse_aln", "position_map",
               "conservation", "change_point", "consensus_build")) {
  timing_env[[slot]] <- 0
}
t_main_start <- proc.time()[["elapsed"]]

for (lin in lineages) {
  idx <- which(ltr_class_5 == lin & !is.na(ltr_class_5))
  if (length(idx) == 0L) next

  body_lin     <- ltr_body_seqs[idx]
  parent_lin   <- ltr_parent_5[idx]
  all_idx_lin  <- all_idx_of_5[idx]           # index in the all-LTR frame
  priority_lin <- ltr_priority_5[idx]

  # Cluster on body-only sequences
  t0 <- proc.time()[["elapsed"]]
  if (length(body_lin) >= 2L) {
    membership <- mmseqs_cluster(body_lin, identity = 0.9, threads = opt$threads)
  } else {
    membership <- setNames(names(body_lin), names(body_lin))
  }
  timing_env$mmseqs <- timing_env$mmseqs + (proc.time()[["elapsed"]] - t0)

  reps <- unique(membership)
  cat(sprintf("  %-60s : %d 5'LTRs -> %d clusters\n",
              substr(lin, 1L, 60L), length(body_lin), length(reps)))

  for (rep_name in reps) {
    ltr_counter <- ltr_counter + 1L
    members     <- names(membership)[membership == rep_name]
    midx        <- match(members, names(body_lin))
    midx        <- midx[!is.na(midx)]
    if (length(midx) == 0L) next

    ltr_id <- sprintf("LTR_%06d", ltr_counter)

    use_mafft <- length(midx) >= opt$min_cluster_size
    cons_done <- FALSE

    if (use_mafft) {
      aln_path <- if (!is.null(opt$alignments_dir)) {
        file.path(opt$alignments_dir, paste0(ltr_id, ".aln.fasta"))
      } else NULL

      # Build joint 5'+3' input. For each 5'LTR in the cluster, pull its
      # sibling 3'LTR (same Parent TE id) if one exists.
      members_5 <- all_idx_lin[midx]
      sib_3     <- unlist(lapply(parent_lin[midx], function(p) {
        s <- by_parent_3[[p]]
        if (is.null(s) || length(s) == 0L) NA_integer_ else s[1L]
      }))
      members_3   <- sib_3[!is.na(sib_3)]
      joint_idx   <- c(members_5, members_3)
      joint_roles <- c(rep("5LTR", length(members_5)),
                       rep("3LTR", length(members_3)))

      res <- mafft_boundary_consensus(
        extended_seqs    = ltr_ext_seqs_all[joint_idx],
        body_starts_raw  = ann_5_pos_raw[joint_idx],
        body_ends_raw    = ann_3_pos_raw[joint_idx],
        flank            = F_flank,
        threads          = opt$threads,
        save_aln_path    = aln_path,
        ltr_id           = ltr_id,
        lineage          = lin,
        roles            = joint_roles,
        timing_env       = timing_env
      )
      qc_rows[[length(qc_rows) + 1L]] <- res$qc

      if (res$status == "ok" && !is.null(res$consensus)) {
        cons      <- res$consensus
        cons_done <- TRUE
      }
    }

    if (!cons_done) {
      # Fallback: use highest-rank single annotated LTR body
      best <- midx[which.max(priority_lin[midx])]
      cons <- body_lin[best]
      if (!use_mafft) {
        qc_rows[[length(qc_rows) + 1L]] <- data.frame(
          ltr_id = ltr_id, Final_Classification = lin,
          n_members = length(midx),
          n_5ltr = length(midx), n_3ltr = 0L,
          annotated_5_col = NA_integer_, corrected_5_col = NA_integer_,
          shift_5 = NA_integer_,
          annotated_3_col = NA_integer_, corrected_3_col = NA_integer_,
          shift_3 = NA_integer_,
          consensus_length = as.integer(width(cons)[1L]),
          median_annot_body_len = as.integer(median(width(body_lin[midx]))),
          stringsAsFactors = FALSE
        )
      }
    }

    names(cons)   <- ltr_id
    all_consensus <- c(all_consensus, cons)
    id_to_lineage[ltr_id] <- lin
  }
}

t_main_total <- proc.time()[["elapsed"]] - t_main_start
cat("\n[timing] consensus/cluster loop breakdown (wall-clock s):\n")
cat(sprintf("  mmseqs2         : %6.2f\n",  timing_env$mmseqs))
cat(sprintf("  mafft (extern)  : %6.2f\n",  timing_env$mafft))
cat(sprintf("  parse alignment : %6.2f\n",  timing_env$parse_aln))
cat(sprintf("  position map    : %6.2f\n",  timing_env$position_map))
cat(sprintf("  conservation    : %6.2f\n",  timing_env$conservation))
cat(sprintf("  change-point    : %6.2f\n",  timing_env$change_point))
cat(sprintf("  consensus build : %6.2f\n",  timing_env$consensus_build))
cat(sprintf("  --- total loop  : %6.2f\n",  t_main_total))

if (length(all_consensus) == 0L) {
  stop("No LTR consensus sequences built. Check input GFF3.", call. = FALSE)
}

# Write library and ID map
ltr_lib_path <- paste0(opt$output, "_LTR_library.fasta")
writeXStringSet(all_consensus, ltr_lib_path)

ltr_map_path <- paste0(opt$output, "_LTR_library_map.tsv")
write.table(data.frame(ltr_id = names(id_to_lineage),
                        Final_Classification = as.character(id_to_lineage),
                        stringsAsFactors = FALSE),
            ltr_map_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Write boundary QC table
qc_path <- paste0(opt$output, "_LTR_library_boundary_qc.tsv")
if (length(qc_rows) > 0L) {
  qc_df <- do.call(rbind, qc_rows)
  write.table(qc_df, qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(data.frame(
    ltr_id = character(0L), Final_Classification = character(0L),
    n_members = integer(0L), n_5ltr = integer(0L), n_3ltr = integer(0L),
    annotated_5_col = integer(0L), corrected_5_col = integer(0L), shift_5 = integer(0L),
    annotated_3_col = integer(0L), corrected_3_col = integer(0L), shift_3 = integer(0L),
    consensus_length = integer(0L), median_annot_body_len = integer(0L)
  ), qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("LTR library: %d consensus sequences written to %s\n",
            length(all_consensus), ltr_lib_path))

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
