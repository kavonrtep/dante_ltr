# solo_ltr_utils.R
# Shared utility functions for solo LTR detection.
# Sourced by build_ltr_library.R and detect_solo_ltr.R.

suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
})

# ============================================================
# TSD utilities
# ============================================================

# Build per-lineage modal TSD length from DLT/DLTP elements.
# Returns a named integer vector [lineage -> tsd_length], or NULL.
extract_tsd_length_map <- function(gff) {
  te_df <- as.data.frame(gff[gff$type == "transposable_element"])
  te_df <- te_df[!is.na(te_df$Rank) & te_df$Rank %in% c("DLT", "DLTP"), ]
  te_df <- te_df[!is.na(te_df$TSD) & te_df$TSD != "not_found" & nchar(te_df$TSD) > 0, ]
  if (nrow(te_df) == 0L) return(NULL)

  # TSD attribute can be "SEQ" (exact) or "SEQ1/SEQ2" (1-mismatch); take first part
  te_df$tsd_clean  <- gsub("/.+", "", te_df$TSD)
  te_df$tsd_length <- nchar(te_df$tsd_clean)
  te_df <- te_df[te_df$tsd_length >= 3L & te_df$tsd_length <= 8L, ]
  if (nrow(te_df) == 0L) return(NULL)

  lineages <- unique(te_df$Final_Classification)
  lineages <- lineages[!is.na(lineages)]

  tsd_map <- vapply(lineages, function(lin) {
    lens <- te_df$tsd_length[te_df$Final_Classification == lin]
    if (length(lens) < 3L) return(NA_integer_)
    as.integer(names(which.max(table(lens))))
  }, integer(1L))

  tsd_map[!is.na(tsd_map)]
}


# TSD check for a single solo LTR BLAST hit (GRanges of length 1).
# tsd_length: modal length from lineage map. When supplied, we scan (L-1):(L+1)
#   (clipped to [3, 8]); when NULL, we fall back to 4:6.
# Two-pass search: first try exact match (longest first), then 1-mismatch (for
#   len >= 4, longest first). Exact always wins over 1-mismatch.
# Returns: list(TSD_sequence, TSD_length, confirmed, TSD_L_position, TSD_R_position)
check_tsd <- function(hit_gr, s, tsd_length = NULL) {
  sn      <- as.character(seqnames(hit_gr))
  str_val <- as.character(strand(hit_gr))
  seq_len <- seqlengths(s)[sn]

  scan_lengths <- if (!is.null(tsd_length) && !is.na(tsd_length) && length(tsd_length) == 1L) {
    L <- as.integer(tsd_length)
    seq.int(max(3L, L - 1L), min(8L, L + 1L))
  } else {
    4L:6L
  }
  scan_lengths <- sort(unique(scan_lengths), decreasing = TRUE)
  max_len <- max(scan_lengths)

  l_start <- max(1L, start(hit_gr) - max_len)
  l_end   <- start(hit_gr) - 1L
  r_start <- end(hit_gr)   + 1L
  r_end   <- min(seq_len,   end(hit_gr) + max_len)

  not_found <- list(TSD_sequence = "not_found", TSD_length = 0L,
                    confirmed = FALSE, TSD_L_position = NULL, TSD_R_position = NULL)

  if (l_end < l_start || r_end < r_start) return(not_found)

  # Carry seqinfo so getSeq's internal merge doesn't warn about mismatched seqlevels.
  tsd_l_full <- as.character(getSeq(
    s, GRanges(sn, IRanges(l_start, l_end), seqinfo = seqinfo(s))
  ))
  tsd_r_full <- as.character(getSeq(
    s, GRanges(sn, IRanges(r_start, r_end), seqinfo = seqinfo(s))
  ))

  slice <- function(len) {
    if (nchar(tsd_l_full) < len || nchar(tsd_r_full) < len) return(NULL)
    l <- substr(tsd_l_full, nchar(tsd_l_full) - len + 1L, nchar(tsd_l_full))
    r <- substr(tsd_r_full, 1L, len)
    if (str_val == "-") {
      l <- as.character(reverseComplement(DNAString(l)))
      r <- as.character(reverseComplement(DNAString(r)))
    }
    list(l = l, r = r)
  }

  # Attach seqinfo to the TSD GRanges so later do.call(c, ...) across hits on
  # different contigs doesn't warn about disjoint seqlevels.
  si <- seqinfo(hit_gr)
  make_result <- function(seq_str, len) {
    list(
      TSD_sequence   = seq_str,
      TSD_length     = len,
      confirmed      = TRUE,
      TSD_L_position = GRanges(sn, IRanges(start(hit_gr) - len, start(hit_gr) - 1L),
                               seqinfo = si),
      TSD_R_position = GRanges(sn, IRanges(end(hit_gr) + 1L,    end(hit_gr) + len),
                               seqinfo = si)
    )
  }

  # Pass 1: exact match, longest first
  for (len in scan_lengths) {
    sl <- slice(len); if (is.null(sl)) next
    if (sl$l == sl$r) return(make_result(sl$l, len))
  }

  # Pass 2: 1-mismatch (len >= 4), longest first
  for (len in scan_lengths) {
    if (len < 4L) next
    sl <- slice(len); if (is.null(sl)) next
    n_mm <- sum(strsplit(sl$l, "")[[1L]] != strsplit(sl$r, "")[[1L]])
    if (n_mm <= 1L) return(make_result(paste0(sl$r, "/", sl$l), len))
  }

  not_found
}


# ============================================================
# Junction BLAST checks
# ============================================================

# Extract a genomic window for junction checking.
# direction: "downstream" (after LTR end, checks 5'UTR/PBS)
#            "upstream"   (before LTR start, checks PPT)
extract_junction_window <- function(hit_gr, s, direction, window = 30L) {
  sn      <- as.character(seqnames(hit_gr))
  str_val <- as.character(strand(hit_gr))
  seq_len <- seqlengths(s)[sn]

  if (direction == "downstream") {
    if (str_val == "+") {
      w_s <- end(hit_gr) + 1L
      w_e <- min(seq_len, end(hit_gr) + window)
      if (w_e < w_s) return(NULL)
      getSeq(s, GRanges(sn, IRanges(w_s, w_e)))
    } else {
      w_s <- max(1L, start(hit_gr) - window)
      w_e <- start(hit_gr) - 1L
      if (w_e < w_s) return(NULL)
      reverseComplement(getSeq(s, GRanges(sn, IRanges(w_s, w_e))))
    }
  } else {  # upstream
    if (str_val == "+") {
      w_s <- max(1L, start(hit_gr) - window)
      w_e <- start(hit_gr) - 1L
      if (w_e < w_s) return(NULL)
      getSeq(s, GRanges(sn, IRanges(w_s, w_e)))
    } else {
      w_s <- end(hit_gr) + 1L
      w_e <- min(seq_len, end(hit_gr) + window)
      if (w_e < w_s) return(NULL)
      reverseComplement(getSeq(s, GRanges(sn, IRanges(w_s, w_e))))
    }
  }
}


# BLAST a short window sequence against a nucleotide database.
# Returns TRUE if any hit passes min_length and min_pident thresholds.
blast_window_vs_db <- function(win_seq, db_path, min_length = 15L, min_pident = 80) {
  if (is.null(win_seq) || length(win_seq) == 0L || nchar(win_seq) < 10L) return(FALSE)
  tmp_q   <- tempfile()
  tmp_out <- tempfile()
  on.exit({ unlink(tmp_q); unlink(tmp_out) }, add = TRUE)

  names(win_seq) <- "query"
  writeXStringSet(DNAStringSet(as.character(win_seq)), tmp_q)

  cmd <- paste("blastn -task blastn -word_size 7 -dust no",
               "-perc_identity", min_pident,
               "-query", tmp_q, "-db", db_path, "-strand plus",
               '-outfmt "6 qaccver saccver pident length"',
               "-out", tmp_out, "2>/dev/null")
  system(cmd)

  if (!file.exists(tmp_out) || file.size(tmp_out) == 0L) return(FALSE)
  res <- tryCatch(
    read.table(tmp_out, as.is = TRUE, sep = "\t",
               col.names = c("qaccver", "saccver", "pident", "length")),
    error = function(e) data.frame()
  )
  nrow(res) > 0L && any(res$length >= min_length & res$pident >= min_pident)
}


# 5'UTR junction check: positive = hit is likely the 5' LTR of an unannotated element
check_utr5_junction <- function(hit_gr, s, utr5_db) {
  win <- extract_junction_window(hit_gr, s, "downstream")
  blast_window_vs_db(win, utr5_db)
}

# PPT junction check: positive = hit is likely the 3' LTR of an unannotated element
check_ppt_junction <- function(hit_gr, s, ppt_db) {
  win <- extract_junction_window(hit_gr, s, "upstream")
  blast_window_vs_db(win, ppt_db)
}


# ============================================================
# Batched junction + PBS checks
# ============================================================
# These replace the per-hit calls above for throughput.  One BLAST invocation
# per database instead of 3 per hit.  Query names encode the hit index (q<i>);
# any query with a qualifying hit in the BLAST output → that hit is flagged TRUE.

# Batched UTR5 / PPT junction check.
# direction: "downstream" (UTR5) or "upstream" (PPT).
# Returns a logical vector of length(hits_gr).
batch_check_junction <- function(hits_gr, s, direction, db_path,
                                 window = 30L, min_length = 15L,
                                 min_pident = 80, cpu = 1L) {
  n <- length(hits_gr)
  result <- logical(n)
  if (n == 0L) return(result)

  sn       <- as.character(seqnames(hits_gr))
  str_val  <- as.character(strand(hits_gr))
  seq_lens <- seqlengths(s)[sn]

  if (direction == "downstream") {
    # After the LTR (towards internal region): right flank on +, left flank on - (RC)
    w_s_p <- end(hits_gr) + 1L
    w_e_p <- pmin(seq_lens, end(hits_gr) + window)
    w_s_m <- pmax(1L, start(hits_gr) - window)
    w_e_m <- start(hits_gr) - 1L
    need_rc <- str_val == "-"
  } else {  # upstream
    w_s_p <- pmax(1L, start(hits_gr) - window)
    w_e_p <- start(hits_gr) - 1L
    w_s_m <- end(hits_gr) + 1L
    w_e_m <- pmin(seq_lens, end(hits_gr) + window)
    need_rc <- str_val == "-"
  }

  w_s   <- ifelse(str_val == "+", w_s_p, w_s_m)
  w_e   <- ifelse(str_val == "+", w_e_p, w_e_m)
  valid <- !is.na(seq_lens) & w_e >= w_s

  if (!any(valid)) return(result)

  idx    <- which(valid)
  win_gr <- GRanges(sn[idx], IRanges(w_s[idx], w_e[idx]), seqinfo = seqinfo(s))
  wins   <- getSeq(s, win_gr)
  if (any(need_rc[idx])) {
    rc_where <- which(need_rc[idx])
    wins[rc_where] <- reverseComplement(wins[rc_where])
  }

  # Drop windows shorter than 10 bp (BLAST word_size is 7; very short ones are noise)
  too_short <- width(wins) < 10L
  if (any(too_short)) {
    wins <- wins[!too_short]
    idx  <- idx[!too_short]
  }
  if (length(wins) == 0L) return(result)

  names(wins) <- sprintf("q%d", idx)

  tmp_q   <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".tsv")
  on.exit({ unlink(tmp_q); unlink(tmp_out) }, add = TRUE)
  writeXStringSet(wins, tmp_q)

  cmd <- paste("blastn -task blastn -word_size 7 -dust no",
               "-perc_identity", min_pident,
               "-query", tmp_q, "-db", db_path, "-strand plus",
               "-num_threads", cpu,
               '-outfmt "6 qaccver saccver pident length"',
               "-out", tmp_out, "2>/dev/null")
  system(cmd)

  if (!file.exists(tmp_out) || file.size(tmp_out) == 0L) return(result)
  res <- tryCatch(
    read.table(tmp_out, as.is = TRUE, sep = "\t",
               col.names = c("qaccver", "saccver", "pident", "length")),
    error = function(e) data.frame()
  )
  if (nrow(res) == 0L) return(result)

  pass <- res$length >= min_length & res$pident >= min_pident
  if (!any(pass)) return(result)
  hit_idx <- suppressWarnings(as.integer(sub("^q", "", res$qaccver[pass])))
  hit_idx <- hit_idx[!is.na(hit_idx) & hit_idx >= 1L & hit_idx <= n]
  result[hit_idx] <- TRUE
  result
}


# Batched PBS check across many hits.  For each + hit the window is the
# reverse-complement of the downstream flank; for each - hit it's the upstream
# flank as-is.
batch_check_pbs <- function(hits_gr, s, trna_db, cpu = 1L,
                            window = 31L, min_length = 12L,
                            min_pident = 90, min_qend = 20L) {
  n <- length(hits_gr)
  result <- logical(n)
  if (n == 0L) return(result)

  sn       <- as.character(seqnames(hits_gr))
  str_val  <- as.character(strand(hits_gr))
  seq_lens <- seqlengths(s)[sn]

  w_s_p <- end(hits_gr) + 1L
  w_e_p <- pmin(seq_lens, end(hits_gr) + window)
  w_s_m <- pmax(1L, start(hits_gr) - window)
  w_e_m <- start(hits_gr) - 1L
  w_s   <- ifelse(str_val == "+", w_s_p, w_s_m)
  w_e   <- ifelse(str_val == "+", w_e_p, w_e_m)
  valid <- !is.na(seq_lens) & w_e >= w_s
  need_rc <- str_val == "+"

  if (!any(valid)) return(result)

  idx    <- which(valid)
  win_gr <- GRanges(sn[idx], IRanges(w_s[idx], w_e[idx]), seqinfo = seqinfo(s))
  wins   <- getSeq(s, win_gr)
  if (any(need_rc[idx])) {
    rc_where <- which(need_rc[idx])
    wins[rc_where] <- reverseComplement(wins[rc_where])
  }

  too_short <- width(wins) < 10L
  if (any(too_short)) {
    wins <- wins[!too_short]
    idx  <- idx[!too_short]
  }
  if (length(wins) == 0L) return(result)

  names(wins) <- sprintf("q%d", idx)

  tmp_q   <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".tsv")
  on.exit({ unlink(tmp_q); unlink(tmp_out) }, add = TRUE)
  writeXStringSet(wins, tmp_q)

  cmd <- paste("blastn -task blastn -word_size 7 -dust no",
               "-perc_identity", min_pident,
               "-query", tmp_q, "-db", trna_db, "-strand plus",
               "-num_threads", cpu,
               '-outfmt "6 qaccver saccver pident length qstart qend sstart send"',
               "-out", tmp_out, "2>/dev/null")
  system(cmd)

  if (!file.exists(tmp_out) || file.size(tmp_out) == 0L) return(result)
  res <- tryCatch(
    read.table(tmp_out, as.is = TRUE, sep = "\t",
               col.names = c("qaccver", "saccver", "pident", "length",
                             "qstart", "qend", "sstart", "send")),
    error = function(e) data.frame()
  )
  if (nrow(res) == 0L) return(result)

  pass <- res$length >= min_length & res$pident >= min_pident & res$qend > min_qend
  if (!any(pass)) return(result)
  hit_idx <- suppressWarnings(as.integer(sub("^q", "", res$qaccver[pass])))
  hit_idx <- hit_idx[!is.na(hit_idx) & hit_idx >= 1L & hit_idx <= n]
  result[hit_idx] <- TRUE
  result
}

# Simplified PBS check using the tRNA database.
check_pbs_solo <- function(hit_gr, s, trna_db) {
  sn      <- as.character(seqnames(hit_gr))
  str_val <- as.character(strand(hit_gr))
  seq_len <- seqlengths(s)[sn]

  if (str_val == "+") {
    pbs_s <- end(hit_gr)   + 1L
    pbs_e <- min(seq_len,   end(hit_gr) + 31L)
    if (pbs_e < pbs_s) return(FALSE)
    pbs_seq <- reverseComplement(getSeq(s, GRanges(sn, IRanges(pbs_s, pbs_e))))
  } else {
    pbs_e <- start(hit_gr) - 1L
    pbs_s <- max(1L, start(hit_gr) - 31L)
    if (pbs_e < pbs_s) return(FALSE)
    pbs_seq <- getSeq(s, GRanges(sn, IRanges(pbs_s, pbs_e)))
  }

  tmp_q   <- tempfile()
  tmp_out <- tempfile()
  on.exit({ unlink(tmp_q); unlink(tmp_out) }, add = TRUE)
  names(pbs_seq) <- "pbs_region"
  writeXStringSet(DNAStringSet(as.character(pbs_seq)), tmp_q)

  cmd <- paste("blastn -task blastn -word_size 7 -dust no -perc_identity 90",
               "-query", tmp_q, "-db", trna_db, "-strand plus",
               '-outfmt "6 qaccver saccver pident length qstart qend sstart send"',
               "-out", tmp_out, "2>/dev/null")
  system(cmd)

  if (!file.exists(tmp_out) || file.size(tmp_out) == 0L) return(FALSE)
  res <- tryCatch(
    read.table(tmp_out, as.is = TRUE, sep = "\t",
               col.names = c("qaccver","saccver","pident","length","qstart","qend","sstart","send")),
    error = function(e) data.frame()
  )
  if (nrow(res) == 0L) return(FALSE)
  any(res$length >= 12L & res$pident >= 90 & res$qend > 20L)
}


# ============================================================
# GFF3 assembly
# ============================================================

# Assemble solo_LTR GRanges from BLAST hits + validation results.
# Vectorized: builds the solo GRanges and TSD children in whole-vector passes
# rather than per-hit concatenation (old O(N)-many GRanges appends was dominating
# the runtime on 3k+ hits).
# hits_gr       : GRanges with metadata: Final_Classification, Identity, Coverage, LibraryID
# tsd_list      : list of check_tsd() results, one per hit
# utr5_results  : logical vector (per hit)
# ppt_results   : logical vector (per hit)
# pbs_results   : logical vector (per hit)
# id_offset     : integer added to sequential IDs (for multi-chunk uniqueness)
make_solo_ltr_gff3 <- function(hits_gr, tsd_list,
                               utr5_results, ppt_results, pbs_results,
                               id_offset = 0L) {
  n <- length(hits_gr)
  if (n == 0L) return(GRanges())

  # ---- Extract per-hit TSD fields as vectors ----
  confirmed <- vapply(tsd_list, function(x) isTRUE(x$confirmed), logical(1L))
  tsd_seqs  <- vapply(tsd_list, function(x) as.character(x$TSD_sequence),
                      character(1L))
  rank      <- ifelse(confirmed, "SL", "SL_noTSD")
  ids       <- sprintf("soloLTR_%08d", seq_len(n) + id_offset)

  lib_id <- as.character(hits_gr$LibraryID)
  lib_id[is.na(lib_id)] <- ""

  utr5_ok <- as.logical(utr5_results); utr5_ok[is.na(utr5_ok)] <- FALSE
  ppt_ok  <- as.logical(ppt_results);  ppt_ok[is.na(ppt_ok)]   <- FALSE
  pbs_ok  <- as.logical(pbs_results);  pbs_ok[is.na(pbs_ok)]   <- FALSE

  attr_for <- function(is_no_tsd, ok) ifelse(
    is_no_tsd, ifelse(ok, "positive", "negative"), "not_applicable"
  )
  is_no_tsd <- rank == "SL_noTSD"

  # ---- Solo GRanges built in a single pass ----
  solo <- hits_gr
  mcols(solo) <- NULL
  solo$type               <- rep("solo_LTR",  n)
  solo$source             <- rep("dante_ltr", n)
  solo$ID                 <- ids
  solo$Final_Classification <- as.character(hits_gr$Final_Classification)
  solo$Identity           <- round(as.numeric(hits_gr$Identity), 2)
  solo$Coverage           <- round(as.numeric(hits_gr$Coverage), 3)
  solo$Rank               <- rank
  solo$LibraryID          <- lib_id
  solo$TSD                <- ifelse(confirmed, tsd_seqs, "not_found")
  solo$UTR5_junction      <- attr_for(is_no_tsd, utr5_ok)
  solo$PPT_junction       <- attr_for(is_no_tsd, ppt_ok)
  solo$PBS_check          <- attr_for(is_no_tsd, pbs_ok)

  # ---- TSD children: vectorized across confirmed hits ----
  tsd_idx <- which(confirmed)
  if (length(tsd_idx) > 0L) {
    tsd_l_list <- lapply(tsd_list[tsd_idx], `[[`, "TSD_L_position")
    tsd_r_list <- lapply(tsd_list[tsd_idx], `[[`, "TSD_R_position")
    # Drop any NULLs defensively
    keep_l <- !vapply(tsd_l_list, is.null, logical(1L))
    keep_r <- !vapply(tsd_r_list, is.null, logical(1L))
    keep   <- keep_l & keep_r
    if (any(keep)) {
      tsd_idx_ok <- tsd_idx[keep]
      tsd_l_gr   <- do.call(c, tsd_l_list[keep])
      tsd_r_gr   <- do.call(c, tsd_r_list[keep])
      tsd_all_gr <- c(tsd_l_gr, tsd_r_gr)
      # Align seqinfo with parent solo so the final c() has no seqlevel mismatch.
      seqlevels(tsd_all_gr) <- seqlevels(solo)
      seqinfo(tsd_all_gr)   <- seqinfo(solo)
      # Strand matches the parent solo_LTR
      str_vec              <- as.character(strand(hits_gr[tsd_idx_ok]))
      strand(tsd_all_gr)   <- c(str_vec, str_vec)
      tsd_all_gr$type      <- rep("target_site_duplication", length(tsd_all_gr))
      tsd_all_gr$source    <- rep("dante_ltr",               length(tsd_all_gr))
      tsd_all_gr$Parent    <- c(ids[tsd_idx_ok], ids[tsd_idx_ok])
    } else {
      tsd_all_gr <- GRanges()
    }
  } else {
    tsd_all_gr <- GRanges()
  }

  # Combining solo + TSD children emits a harmless seqinfo-merge warning when
  # one side has no seqlevels (e.g. empty chunk).  Suppress; output correctness
  # is checked by MD5 against the per-hit-loop implementation.
  suppressWarnings(c(solo, tsd_all_gr))
}


# ============================================================
# Statistics
# ============================================================

# Count SL / SL_noTSD per lineage, compute Rsf vs complete elements in input GFF.
get_solo_ltr_statistics <- function(gff_out, complete_gff) {
  solo     <- gff_out[!is.na(gff_out$type) & gff_out$type == "solo_LTR"]
  sl       <- solo[!is.na(solo$Rank) & solo$Rank == "SL"]
  sl_no    <- solo[!is.na(solo$Rank) & solo$Rank == "SL_noTSD"]
  complete <- if (length(complete_gff) > 0L) {
    complete_gff[!is.na(complete_gff$type) & complete_gff$type == "transposable_element"]
  } else GRanges()

  all_lin <- sort(unique(c(
    as.character(solo$Final_Classification),
    as.character(complete$Final_Classification)
  )))
  all_lin <- all_lin[!is.na(all_lin) & nchar(all_lin) > 0L]

  tbl <- do.call(rbind, lapply(all_lin, function(lin) {
    n_sl   <- sum(as.character(sl$Final_Classification)    == lin, na.rm = TRUE)
    n_no   <- sum(as.character(sl_no$Final_Classification) == lin, na.rm = TRUE)
    n_comp <- sum(as.character(complete$Final_Classification) == lin, na.rm = TRUE)
    rsf    <- if (n_comp > 0L) round(n_sl / n_comp, 3) else NA_real_
    data.frame(Classification = lin, SL = n_sl, SL_noTSD = n_no,
               Complete_elements = n_comp, Rsf = rsf,
               stringsAsFactors = FALSE)
  }))

  if (is.null(tbl) || nrow(tbl) == 0L) {
    return(data.frame(Classification = character(0), SL = integer(0),
                      SL_noTSD = integer(0), Complete_elements = integer(0),
                      Rsf = numeric(0)))
  }
  tbl
}
