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
  make_option(c("-f", "--flank"), type = "integer", default = 15L,
              help = "Flanking bp added each side for MSA [default %default]"),
  make_option(c("-d", "--min_cluster_size"), type = "integer", default = 3L,
              help = "Min cluster members required for MAFFT consensus [default %default]")
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

# ---- MAFFT consensus from a set of sequences ----
mafft_consensus <- function(seqs, threads = 1L, gap_threshold = 0.5) {
  tmp_in  <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".fasta")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  writeXStringSet(seqs, tmp_in)
  ret <- system(
    paste("mafft --auto --thread", threads, "--quiet", tmp_in, ">", tmp_out, "2>/dev/null")
  )

  if (ret != 0L || !file.exists(tmp_out) || file.size(tmp_out) == 0L) {
    return(seqs[which.max(width(seqs))])  # fallback: longest member
  }

  aln <- readDNAStringSet(tmp_out)
  if (length(aln) == 0L) return(seqs[which.max(width(seqs))])

  mat <- as.matrix(aln)
  cons_chars <- apply(mat, 2L, function(col) {
    if (mean(col == "-") > gap_threshold) return(NA_character_)
    non_gap <- col[col != "-" & col != "N"]
    if (length(non_gap) == 0L) return(NA_character_)
    names(which.max(table(non_gap)))
  })
  consensus <- paste(cons_chars[!is.na(cons_chars)], collapse = "")
  if (nchar(consensus) < 50L) return(seqs[which.max(width(seqs))])
  DNAStringSet(setNames(consensus, paste0("consensus_", names(seqs)[1L])))
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

# ---- Extract sequences ----
seqlengths(ltr_feat) <- sl_vec[seqlevels(ltr_feat)]
ltr_seqs <- getSeq(s, ltr_feat)

# Safe names: no spaces (BLAST truncates at space), use # as separator for lineage
safe_coords <- paste0(seqnames(ltr_feat), "_", start(ltr_feat), "_", end(ltr_feat))
names(ltr_seqs) <- safe_coords   # will be enriched after clustering

ltr_priority <- RANK_PRIORITY[ltr_rank]
ltr_priority[is.na(ltr_priority)] <- 0L

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

# ---- Per-lineage clustering and consensus ----
cat("Clustering and building LTR consensus sequences ...\n")

lineages <- unique(ltr_class[!is.na(ltr_class)])
all_consensus  <- DNAStringSet()
id_to_lineage  <- character(0L)    # maps FASTA ID -> Full_Classification
ltr_counter    <- 0L

for (lin in lineages) {
  idx <- which(ltr_class == lin & !is.na(ltr_class))
  if (length(idx) == 0L) next

  seqs_lin     <- ltr_seqs[idx]
  priority_lin <- ltr_priority[idx]

  # Cluster within lineage
  if (length(seqs_lin) >= 2L) {
    membership <- mmseqs_cluster(seqs_lin, identity = 0.9, threads = opt$threads)
  } else {
    membership <- setNames(names(seqs_lin), names(seqs_lin))
  }

  reps <- unique(membership)
  cat(sprintf("  %-60s : %d LTRs -> %d clusters\n",
              substr(lin, 1L, 60L), length(seqs_lin), length(reps)))

  for (rep_name in reps) {
    ltr_counter <- ltr_counter + 1L
    members     <- names(membership)[membership == rep_name]
    midx        <- match(members, names(seqs_lin))
    midx        <- midx[!is.na(midx)]
    if (length(midx) == 0L) next

    if (length(midx) >= opt$min_cluster_size) {
      cons <- mafft_consensus(seqs_lin[midx], threads = opt$threads)
    } else {
      best <- midx[which.max(priority_lin[midx])]
      cons <- seqs_lin[best]
    }

    ltr_id        <- sprintf("LTR_%06d", ltr_counter)
    names(cons)   <- ltr_id
    all_consensus <- c(all_consensus, cons)
    id_to_lineage[ltr_id] <- lin
  }
}

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
cat(sprintf("  LTR library  : %s\n", ltr_lib_path))
cat(sprintf("  LTR ID map   : %s\n", ltr_map_path))
cat(sprintf("  5'UTR tags   : %s\n", utr5_path))
cat(sprintf("  PPT tags     : %s\n", ppt_path))
cat(sprintf("  TSD map      : %s\n", tsd_map_path))
