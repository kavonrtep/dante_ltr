#!/usr/bin/env Rscript
# inspect_ltr_alignment.R --- quick-look for one per-cluster MSA produced
# by build_ltr_library.R.  Prints the cluster's QC numbers and an
# ASCII conservation profile with the annotated (A) and corrected
# (5/3) boundary columns marked, so you can see the flank -> body
# switchpoint at a glance.
#
# Usage:
#   Rscript utils/inspect_ltr_alignment.R <alignments_dir> <ltr_id>
#   Rscript utils/inspect_ltr_alignment.R <alignments_dir> <ltr_id> <qc.tsv>
#
# Defaults:
#   qc_tsv is <alignments_dir>/../*_LTR_library_boundary_qc.tsv
#   (first match)
#
# Example:
#   Rscript utils/inspect_ltr_alignment.R tmp/alyr_lib/alignments LTR_000001

suppressPackageStartupMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  cat("usage: Rscript utils/inspect_ltr_alignment.R ",
      "<alignments_dir> <LTR_id> [qc_tsv]\n", sep = "")
  quit(status = 2)
}
aln_dir <- args[[1L]]
ltr_id  <- args[[2L]]
qc_path <- if (length(args) >= 3L) args[[3L]] else {
  cand <- list.files(dirname(normalizePath(aln_dir)),
                     pattern = "_LTR_library_boundary_qc\\.tsv$",
                     full.names = TRUE)
  if (!length(cand)) stop("no QC tsv found; pass one explicitly as 3rd arg")
  cand[[1L]]
}

aln_path <- file.path(aln_dir, paste0(ltr_id, ".aln.fasta"))
if (!file.exists(aln_path)) stop("no alignment at ", aln_path)
if (!file.exists(qc_path))  stop("no QC tsv at ", qc_path)

qc  <- read.table(qc_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
row <- qc[qc$ltr_id == ltr_id, ]
if (nrow(row) != 1L) stop("no QC row for ", ltr_id)

ann5  <- as.integer(row$annotated_5_col)
corr5 <- as.integer(row$corrected_5_col)
ann3  <- as.integer(row$annotated_3_col)
corr3 <- as.integer(row$corrected_3_col)

aln <- readDNAStringSet(aln_path)
mat <- toupper(as.matrix(aln))
n   <- nrow(mat)
W   <- ncol(mat)

cat(sprintf("%s  (%d rows, %d cols, lineage=%s)\n",
            ltr_id, n, W, row$Final_Classification))
cat(sprintf("  subsets          : 5'LTR=%d, 3'LTR=%d\n",
            row$n_5ltr, row$n_3ltr))
cat(sprintf("  5'-boundary      : annotated col %d, corrected col %d  (shift %+d bp)\n",
            ann5, corr5, corr5 - ann5))
cat(sprintf("  3'-boundary      : annotated col %d, corrected col %d  (shift %+d bp)\n",
            ann3, corr3, corr3 - ann3))
cat(sprintf("  consensus length : %d bp  (median annotated body length %d bp)\n",
            row$consensus_length, row$median_annot_body_len))

# Per-column top-base frequency
top_freq <- apply(mat, 2L, function(col) {
  x <- col[!col %in% c("-", "N")]
  if (!length(x)) return(0)
  max(table(x)) / length(x)
})

# Bin to keep output narrow
BIN <- 5L
binned <- vapply(seq(1L, W, by = BIN), function(i) {
  mean(top_freq[i:min(i + BIN - 1L, W)])
}, numeric(1L))

marker <- vapply(binned, function(v) {
  if (is.na(v))    " " else
  if (v >= 0.9)    "#" else
  if (v >= 0.7)    "=" else
  if (v >= 0.5)    "-" else
  if (v >= 0.3)    "." else
                   " "
}, character(1L))

# Overlay boundary markers
mark_col <- function(col) floor((col - 1L) / BIN) + 1L
marker[mark_col(ann5)]  <- "A"
marker[mark_col(corr5)] <- if (corr5 == ann5) "A" else "5"
marker[mark_col(ann3)]  <- "A"
marker[mark_col(corr3)] <- if (corr3 == ann3) "A" else "3"

cat("\nconservation profile (1 char per 5 aln cols):\n")
cat(paste(marker, collapse = ""), "\n")
cat("legend:  # >=0.9   = >=0.7   - >=0.5   . >=0.3   ' ' < 0.3\n")
cat("         A=annotated boundary   5/3=corrected boundary (if shifted)\n")
