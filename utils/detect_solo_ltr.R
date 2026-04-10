#!/usr/bin/env Rscript
# detect_solo_ltr.R
# Detect solo LTRs in a single genome chunk.
# Steps: BLAST LTR library vs chunk, overlap filter, TSD check,
#        5'UTR/PPT junction checks for no-TSD hits, GFF3 output.
# Called per-chunk by the dante_ltr_solo Python wrapper.

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name   <- "--file="
script_name     <- normalizePath(
  sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)])
)
script_dir <- dirname(script_name)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-s", "--reference_sequence"), type = "character", default = NULL,
              help = "Genome chunk FASTA"),
  make_option(c("-g", "--gff3"), type = "character", default = NULL,
              help = "DANTE_LTR GFF3 for this chunk (used for overlap masking)"),
  make_option(c("-l", "--ltr_library"), type = "character", default = NULL,
              help = "LTR library FASTA (from build_ltr_library.R)"),
  make_option(c("-M", "--ltr_map"), type = "character", default = NULL,
              help = "LTR ID -> lineage map TSV (from build_ltr_library.R)"),
  make_option(c("-u", "--utr5_db"), type = "character", default = NULL,
              help = "5'UTR junction BLAST database prefix"),
  make_option(c("-p", "--ppt_db"), type = "character", default = NULL,
              help = "PPT junction BLAST database prefix"),
  make_option(c("-m", "--tsd_map"), type = "character", default = NULL,
              help = "TSD length map TSV (from build_ltr_library.R)"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file prefix"),
  make_option(c("-c", "--cpu"), type = "integer", default = 1L,
              help = "Number of CPUs [default %default]"),
  make_option(c("-i", "--min_identity"), type = "numeric", default = 80,
              help = "Minimum BLAST identity %% [default %default]"),
  make_option(c("-C", "--min_coverage"), type = "numeric", default = 0.8,
              help = "Minimum alignment coverage of query LTR [default %default]"),
  make_option(c("-T", "--trna_db"), type = "character", default = NULL,
              help = "tRNA BLAST database for PBS check (optional)"),
  make_option(c("--keep_blast"), type = "character", default = NULL,
              help = "Path at which to save the raw BLAST tabular output (for debugging)")
)

parser <- OptionParser(option_list = option_list,
                       usage = "usage: %prog [OPTIONS]")
opt <- parse_args(parser)

required_args <- c("reference_sequence", "gff3", "ltr_library", "ltr_map",
                   "utr5_db", "ppt_db", "tsd_map", "output")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
if (length(missing_args) > 0L) {
  print_help(parser)
  stop("Missing mandatory arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
  library(parallel)
})

source(file.path(script_dir, "solo_ltr_utils.R"))

# ---- Helpers ----

write_empty_outputs <- function(output) {
  writeLines("##gff-version 3", paste0(output, ".gff3"))
  write.table(
    data.frame(Classification = character(0L), SL = integer(0L),
               SL_noTSD = integer(0L), Complete_elements = integer(0L),
               Rsf = numeric(0L)),
    file = paste0(output, "_statistics.csv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

db_exists <- function(prefix) {
  any(file.exists(paste0(prefix, c(".nin", ".nhr", ".nsq", ".nsi", ".nsd"))))
}

# ============================================================
# MAIN
# ============================================================

cat("Reading genome chunk ...\n")
s        <- readDNAStringSet(opt$reference_sequence)
names(s) <- gsub("\\s.+", "", names(s))
sl_vec   <- setNames(as.integer(nchar(s)), names(s))

cat("Reading GFF3 ...\n")
gff <- if (file.exists(opt$gff3) && file.size(opt$gff3) > 0L) {
  tryCatch(import(opt$gff3, format = "gff3"), error = function(e) GRanges())
} else {
  GRanges()
}

# ---- Load helper tables ----
ltr_map_df <- tryCatch(
  read.table(opt$ltr_map, header = TRUE, sep = "\t", as.is = TRUE),
  error = function(e) data.frame(ltr_id = character(0L), Final_Classification = character(0L))
)
ltr_id_to_class <- setNames(ltr_map_df$Final_Classification, ltr_map_df$ltr_id)

tsd_map_df <- tryCatch(
  read.table(opt$tsd_map, header = TRUE, sep = "\t", as.is = TRUE),
  error = function(e) data.frame(lineage = character(0L), tsd_length = integer(0L))
)
tsd_map <- if (nrow(tsd_map_df) > 0L) {
  setNames(as.integer(tsd_map_df$tsd_length), tsd_map_df$lineage)
} else NULL

# ---- BLAST LTR library vs genome ----
cat("Building BLAST database from genome chunk ...\n")
blast_db <- tempfile()
ret_mkdb <- system(paste("makeblastdb -in", opt$reference_sequence,
                          "-dbtype nucl -out", blast_db,
                          "-title genome_chunk 2>/dev/null"))
if (ret_mkdb != 0L) stop("makeblastdb failed", call. = FALSE)

cat("Running BLAST ...\n")
blast_out <- tempfile()
ret_blast <- system(paste(
  "blastn -task blastn",
  "-query",       opt$ltr_library,
  "-db",          blast_db,
  "-num_threads", opt$cpu,
  "-dust no",
  "-perc_identity", opt$min_identity,
  '-outfmt "6 qaccver saccver pident length qlen qstart qend sstart send sstrand evalue bitscore"',
  "-out", blast_out, "2>/dev/null"
))

# Clean up BLAST DB temp files
invisible(lapply(
  list.files(dirname(blast_db),
             pattern = paste0("^", basename(blast_db), "\\."),
             full.names = TRUE),
  unlink
))

# Optionally keep raw BLAST output for debugging
if (!is.null(opt$keep_blast)) {
  if (file.exists(blast_out) && file.size(blast_out) > 0L) {
    try(file.copy(blast_out, opt$keep_blast, overwrite = TRUE), silent = TRUE)
  } else {
    try(file.create(opt$keep_blast), silent = TRUE)
  }
}

if (ret_blast != 0L || !file.exists(blast_out) || file.size(blast_out) == 0L) {
  cat("No BLAST hits found.\n")
  write_empty_outputs(opt$output)
  quit(save = "no", status = 0L)
}

hits_raw <- tryCatch(
  read.table(blast_out, as.is = TRUE, sep = "\t",
             col.names = c("qaccver", "saccver", "pident", "length", "qlen",
                           "qstart", "qend", "sstart", "send", "sstrand",
                           "evalue", "bitscore")),
  error = function(e) data.frame()
)
unlink(blast_out)

if (nrow(hits_raw) == 0L) {
  cat("No BLAST hits after parsing.\n")
  write_empty_outputs(opt$output)
  quit(save = "no", status = 0L)
}

# ---- Filter hits ----
hits_raw$coverage <- hits_raw$length / hits_raw$qlen
hits_raw <- hits_raw[hits_raw$coverage >= opt$min_coverage &
                     hits_raw$evalue   <= 1e-5, , drop = FALSE]

if (nrow(hits_raw) == 0L) {
  cat("No hits pass identity/coverage filters.\n")
  write_empty_outputs(opt$output)
  quit(save = "no", status = 0L)
}
cat(sprintf("%d hits after filtering\n", nrow(hits_raw)))

# Look up lineage from the ID map
hits_raw$Final_Classification <- ltr_id_to_class[hits_raw$qaccver]
hits_raw$Final_Classification[is.na(hits_raw$Final_Classification)] <- "Unknown"

# ---- Convert to GRanges ----
hits_raw$h_start  <- pmin(hits_raw$sstart, hits_raw$send)
hits_raw$h_end    <- pmax(hits_raw$sstart, hits_raw$send)
hits_raw$h_strand <- ifelse(hits_raw$sstrand == "plus", "+", "-")

hits_gr <- GRanges(
  seqnames = hits_raw$saccver,
  IRanges(start = hits_raw$h_start, end = hits_raw$h_end),
  strand               = hits_raw$h_strand,
  Final_Classification = hits_raw$Final_Classification,
  Identity             = round(hits_raw$pident, 2),
  Coverage             = round(hits_raw$coverage, 3)
)
seqlengths(hits_gr) <- sl_vec[seqlevels(hits_gr)]

# ---- Remove hits overlapping annotated features ----
cat("Removing hits overlapping annotated elements ...\n")
if (length(gff) > 0L && length(hits_gr) > 0L) {
  ovl <- findOverlaps(hits_gr, gff, ignore.strand = TRUE)
  if (length(ovl) > 0L) {
    hits_ext  <- hits_gr[from(ovl)]
    gff_ext   <- gff[to(ovl)]
    ovl_width <- pmin(end(hits_ext), end(gff_ext)) -
                 pmax(start(hits_ext), start(gff_ext)) + 1L
    remove_idx <- unique(from(ovl)[ovl_width > 20L])
    if (length(remove_idx) > 0L) hits_gr <- hits_gr[-remove_idx]
  }
}

if (length(hits_gr) == 0L) {
  cat("No hits remain after overlap filtering.\n")
  write_empty_outputs(opt$output)
  quit(save = "no", status = 0L)
}
cat(sprintf("%d hits after overlap filtering\n", length(hits_gr)))

# ---- TSD check ----
cat("Checking TSDs ...\n")
tsd_results <- mclapply(seq_along(hits_gr), function(i) {
  hit  <- hits_gr[i]
  lin  <- as.character(hit$Final_Classification)
  tlen <- if (!is.null(tsd_map)) tsd_map[lin] else NULL
  if (!is.null(tlen) && (length(tlen) == 0L || is.na(tlen))) tlen <- NULL
  check_tsd(hit, s, tlen)
}, mc.cores = opt$cpu, mc.preschedule = FALSE)

n_confirmed <- sum(sapply(tsd_results, `[[`, "confirmed"))
n_failed    <- length(tsd_results) - n_confirmed
cat(sprintf("TSD confirmed: %d | no TSD: %d\n", n_confirmed, n_failed))

# ---- Junction + PBS checks for no-TSD hits ----
failed_idx   <- which(!sapply(tsd_results, `[[`, "confirmed"))
utr5_results <- rep(FALSE, length(hits_gr))
ppt_results  <- rep(FALSE, length(hits_gr))
pbs_results  <- rep(FALSE, length(hits_gr))

utr5_ok <- db_exists(opt$utr5_db)
ppt_ok  <- db_exists(opt$ppt_db)
trna_ok <- !is.null(opt$trna_db) && db_exists(opt$trna_db)

if (length(failed_idx) > 0L) {
  cat(sprintf("Running junction/PBS checks on %d no-TSD hits ...\n", length(failed_idx)))

  junc <- mclapply(failed_idx, function(i) {
    hit <- hits_gr[i]
    list(
      utr5 = if (utr5_ok) check_utr5_junction(hit, s, opt$utr5_db) else FALSE,
      ppt  = if (ppt_ok)  check_ppt_junction(hit, s, opt$ppt_db)   else FALSE,
      pbs  = if (trna_ok) check_pbs_solo(hit, s, opt$trna_db)      else FALSE
    )
  }, mc.cores = opt$cpu, mc.preschedule = FALSE)

  for (k in seq_along(failed_idx)) {
    i <- failed_idx[k]
    utr5_results[i] <- junc[[k]]$utr5
    ppt_results[i]  <- junc[[k]]$ppt
    pbs_results[i]  <- junc[[k]]$pbs
  }
}

# ---- Assemble GFF3 ----
cat("Assembling output GFF3 ...\n")
gff_out <- make_solo_ltr_gff3(hits_gr, tsd_results,
                               utr5_results, ppt_results, pbs_results)

if (length(gff_out) == 0L) {
  cat("No features to write.\n")
  write_empty_outputs(opt$output)
  quit(save = "no", status = 0L)
}

# Append seqname to IDs to ensure uniqueness across chunks (mirrors dante_ltr)
solo_idx <- which(gff_out$type == "solo_LTR")
gff_out$ID[solo_idx] <- paste0(
  gff_out$ID[solo_idx], "_", as.character(seqnames(gff_out))[solo_idx]
)
tsd_idx <- which(gff_out$type == "target_site_duplication")
if (length(tsd_idx) > 0L) {
  gff_out$Parent[tsd_idx] <- paste0(
    gff_out$Parent[tsd_idx], "_", as.character(seqnames(gff_out))[tsd_idx]
  )
}

export(gff_out, con = paste0(opt$output, ".gff3"), format = "gff3")

# ---- Statistics ----
stats_tbl <- get_solo_ltr_statistics(gff_out, gff)
write.table(stats_tbl,
            file  = paste0(opt$output, "_statistics.csv"),
            sep   = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Done. %d solo_LTR features (%d SL, %d SL_noTSD) written to %s.gff3\n",
            length(solo_idx), n_confirmed, n_failed, opt$output))
