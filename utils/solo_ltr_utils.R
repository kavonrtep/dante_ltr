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
# tsd_length: expected length from lineage map; NULL → scan 4:6 bp.
# Returns: list(TSD_sequence, TSD_length, confirmed, TSD_L_position, TSD_R_position)
check_tsd <- function(hit_gr, s, tsd_length = NULL) {
  sn      <- as.character(seqnames(hit_gr))
  str_val <- as.character(strand(hit_gr))
  seq_len <- seqlengths(s)[sn]

  scan_lengths <- if (!is.null(tsd_length) && !is.na(tsd_length) && length(tsd_length) == 1L) {
    as.integer(tsd_length)
  } else {
    4L:6L
  }
  max_len <- max(scan_lengths)

  l_start <- max(1L, start(hit_gr) - max_len)
  l_end   <- start(hit_gr) - 1L
  r_start <- end(hit_gr)   + 1L
  r_end   <- min(seq_len,   end(hit_gr) + max_len)

  not_found <- list(TSD_sequence = "not_found", TSD_length = 0L,
                    confirmed = FALSE, TSD_L_position = NULL, TSD_R_position = NULL)

  if (l_end < l_start || r_end < r_start) return(not_found)

  tsd_l_full <- as.character(getSeq(s, GRanges(sn, IRanges(l_start, l_end))))
  tsd_r_full <- as.character(getSeq(s, GRanges(sn, IRanges(r_start, r_end))))

  for (len in rev(sort(scan_lengths))) {
    if (nchar(tsd_l_full) < len || nchar(tsd_r_full) < len) next

    # innermost `len` bases on left; first `len` bases on right
    tsd_l <- substr(tsd_l_full, nchar(tsd_l_full) - len + 1L, nchar(tsd_l_full))
    tsd_r <- substr(tsd_r_full, 1L, len)

    if (str_val == "-") {
      tsd_l <- as.character(reverseComplement(DNAString(tsd_l)))
      tsd_r <- as.character(reverseComplement(DNAString(tsd_r)))
    }

    make_result <- function(seq_str) {
      list(
        TSD_sequence   = seq_str,
        TSD_length     = len,
        confirmed      = TRUE,
        TSD_L_position = GRanges(sn, IRanges(start(hit_gr) - len, start(hit_gr) - 1L)),
        TSD_R_position = GRanges(sn, IRanges(end(hit_gr) + 1L,    end(hit_gr) + len))
      )
    }

    if (tsd_l == tsd_r) return(make_result(tsd_l))

    if (len >= 5L) {
      n_mm <- sum(strsplit(tsd_l, "")[[1L]] != strsplit(tsd_r, "")[[1L]])
      if (n_mm <= 1L) return(make_result(paste0(tsd_r, "/", tsd_l)))
    }
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
# hits_gr       : GRanges with metadata: Final_Classification, Identity, Coverage
# tsd_list      : list of check_tsd() results, one per hit
# utr5_results  : logical vector
# ppt_results   : logical vector
# pbs_results   : logical vector
# id_offset     : integer added to sequential IDs (for multi-chunk uniqueness)
make_solo_ltr_gff3 <- function(hits_gr, tsd_list,
                               utr5_results, ppt_results, pbs_results,
                               id_offset = 0L) {
  out <- list()

  for (i in seq_along(hits_gr)) {
    hit  <- hits_gr[i]
    tsd  <- tsd_list[[i]]
    rank <- if (tsd$confirmed) "SL" else "SL_noTSD"
    id   <- sprintf("soloLTR_%08d", i + id_offset)

    solo                    <- hit
    mcols(solo)             <- NULL
    solo$type               <- "solo_LTR"
    solo$source             <- "dante_ltr"
    solo$ID                 <- id
    solo$Final_Classification <- as.character(hit$Final_Classification)
    solo$Identity           <- round(as.numeric(hit$Identity), 2)
    solo$Coverage           <- round(as.numeric(hit$Coverage), 3)
    solo$Rank               <- rank
    solo$TSD                <- if (tsd$confirmed) tsd$TSD_sequence else "not_found"
    solo$UTR5_junction      <- if (rank == "SL_noTSD") {
      if (utr5_results[[i]]) "positive" else "negative"
    } else "not_applicable"
    solo$PPT_junction       <- if (rank == "SL_noTSD") {
      if (ppt_results[[i]])  "positive" else "negative"
    } else "not_applicable"
    solo$PBS_check          <- if (rank == "SL_noTSD") {
      if (pbs_results[[i]])  "positive" else "negative"
    } else "not_applicable"

    out <- c(out, list(solo))

    if (tsd$confirmed && !is.null(tsd$TSD_L_position)) {
      for (tsd_pos in list(tsd$TSD_L_position, tsd$TSD_R_position)) {
        tsd_feat         <- tsd_pos
        tsd_feat$type    <- "target_site_duplication"
        tsd_feat$source  <- "dante_ltr"
        tsd_feat$Parent  <- id
        strand(tsd_feat) <- strand(hit)
        out <- c(out, list(tsd_feat))
      }
    }
  }

  if (length(out) == 0L) return(GRanges())
  do.call(c, out)
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
