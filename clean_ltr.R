#!/usr/bin/env Rscript
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- normalizePath(sub(file_arg_name, "", initial_options[grep
                                                                    (file_arg_name,
                                                                     initial_options)]))
script_dir <- dirname(script_name)
library(optparse)

parser <- OptionParser()
option_list <- list(
  make_option(c("-g", "--gff3"), action = "store", type = "character",
              help = "gff3  with LTR Transposable elements", default = NULL),
  make_option(c("-s", "--reference_sequence"), action = "store", type = "character",
              help = "reference sequence as fasta",
              default = NULL),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output file prefix", default = NULL),
  make_option(c("-c", "--cpu"), type =
    "integer", default = 5, help = "Number of cpu to use [default %default]",
              metavar = "number")

)
description <- paste(strwrap(""))

epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description =
  description, usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))


# load packages
suppressPackageStartupMessages({ library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
  library(parallel)
   })

# CONFIGURATION
# load configuration files and functions:
lineage_file <- paste0(script_dir, "/databases/lineage_domain_order.csv")
ltr_utils_r <- paste0(script_dir, "/R/ltr_utils.R")
if (file.exists(lineage_file)) {
  lineage_info <- read.table(lineage_file, sep = "\t",
                             header = TRUE,
                             as.is = TRUE)
  source(ltr_utils_r)
}else {
  lineage_file <- paste0(script_dir, "/.
  ./share/dante_ltr/databases/lineage_domain_order.csv")
  ltr_utils_r <- paste0(script_dir, "/.
  ./share/dante_ltr/R/ltr_utils.R")
  if (file.exists(lineage_file)) {
    lineage_info <- read.table(lineage_file, sep = "\t",
                               header = TRUE,
                               as.is = TRUE)
    source(ltr_utils_r)
  }else {
    (stop("configuration files not found"))
  }
}


ncpus <- opt$cpu


# load data
cat("reading fasta...")
s <- readDNAStringSet(opt$reference_sequence)  # genome assembly
cat("done\n")
outfile <- opt$output
# clean sequence names:
names(s) <- gsub(" .+", "", names(s))
cat("reading gff...")
g <- rtracklayer::import(opt$gff3, format = 'gff3')  # DANTE gff3
cat("done\n")
# testing
if (FALSE) {
  g <- rtracklayer::import("./test_data/sample_ltr_annotation.gff3")
  s <- readDNAStringSet("./test_data/sample_genome.fasta")

  g <- rtracklayer::import("./test_data/DANTE_LTR_Vfaba_chr5.gff3")
  s <- readDNAStringSet("./test_data/211010_Vfaba_chr5.fasta")
  names(s) <- gsub(" .+", "", names(s))
  ncpus <- 10
  lineage_info <- read.table("databases/lineage_domain_order.csv", sep = "\t", header =
    TRUE, as.is = TRUE)
  source("./R/ltr_utils.R")
}


# get te sequence based on rank

# evaluate best TE -  DLTP grou
s_te <- get_te_sequences(g, s)  # split by 'element quality'
# best quality - split by lineage
word_size <- 15
# best TE
TE_DLTP_info <- analyze_TE(s_te$DLTP, word_size = word_size, ncpus = ncpus)

# TE rank 2:
TE_DLT_plus_DLP_info <- analyze_TE(c(s_te$DLP, s_te$DLT), word_size = word_size, ncpus
  = ncpus)
TE_DLT_plus_DLP_info_DLTP_verified <- compare_TE_datasets(c(s_te$DLT, s_te$DLP), ncpus
  = ncpus,
                                                          TE_DLTP_info$seqs_representative,
                                                          word_size = word_size
)

TE_DLT_plus_DLP_info_multiplicity <- verify_based_on_multiplicity(TE_DLT_plus_DLP_info)

# create additional library from rank 2 verified by multiplicity
id_for_additional_library <- setdiff(
  TE_DLT_plus_DLP_info_multiplicity$id_ok_mp_verified,
  TE_DLT_plus_DLP_info_DLTP_verified$id_ok_verified)

if (length(id_for_additional_library) > 1) {
  seqs_for_additional_library <- c(s_te$DLP, s_te$DLT)[names(c(s_te$DLP, s_te$DLT)) %in%
                                                         id_for_additional_library]
  seqs_additional_info <- analyze_TE(seqs_for_additional_library, word_size =
    word_size, ncpus = ncpus)
  seq_representative <- c(
    TE_DLTP_info$seqs_representative,
    seqs_additional_info$seqs_representative
  )
}else {
  if (length(id_for_additional_library) == 1) {
    seq_representative <- c(
      TE_DLTP_info$seqs_representative,
      c(s_te$DLP, s_te$DLT)[names(c(s_te$DLP, s_te$DLT)) %in% id_for_additional_library]
    )
  }else {
    seq_representative <- TE_DLTP_info$seqs_representative
  }
}


# TE  rank 3
TE_DL_info_DLTP_verified <- compare_TE_datasets(
  s_te$DL,
  TE_DLTP_info$seqs_representative, min_coverage = 0.98,
  ncpus = ncpus
)


R <- seq_diversity(seq_representative)$richness
SI <- seq_diversity(seq_representative)$shannon_index

# final RM library:
seq_representative_no_ssr <- seq_representative[R > 20]

ID <- g$ID[g$type == "transposable_element"]
names(ID) <- paste0(seqnames(g), "_",
                    start(g), "_",
                    end(g), "#",
                    g$Final_Classification
)[g$type %in% "transposable_element"]

# create clean gff3
id_of_good_te <- unique(c(
  TE_DLTP_info$te_conflict_info$ok,
  TE_DLT_plus_DLP_info_DLTP_verified$id_ok_verified,
  TE_DLT_plus_DLP_info_multiplicity$id_ok_mp_verified,
  TE_DL_info_DLTP_verified$id_ok_verified)
)

c1 <- g$ID %in% ID[id_of_good_te]
c2 <- sapply(g$Parent, function(x)ifelse(length(x) == 0, "", x)) %in% ID[id_of_good_te]

gff_out <- g[c1 | c2]


writeXStringSet(seq_representative, paste0(opt$output, "_RM_lib.fasta"))
export(gff_out, paste0(opt$output, "_clean.gff3"), format = "gff3")

