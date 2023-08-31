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
ltr_utils_r <- paste0(script_dir, "/utils/ltr_utils.R")
if (file.exists(lineage_file)) {
  lineage_info <- read.table(lineage_file, sep = "\t",
                             header = TRUE,
                             as.is = TRUE)
  source(ltr_utils_r)
}else {
  lineage_file <- paste0(script_dir, "/.
  ./share/dante_ltr/databases/lineage_domain_order.csv")
  ltr_utils_r <- paste0(script_dir, "/.
  ./share/dante_ltr/utils/ltr_utils.R")
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
  g <- rtracklayer::import("./test_data/big_test_data/dante_ltr_unfiltered_t.cacao.gff3")
  s <- readDNAStringSet("./test_data/big_test_data/T_cacao_chromosomes.fasta")

  g <- rtracklayer::import("./test_data/sample_ltr_annotation.gff3")
  g <- rtracklayer::import("./test_data/sample_DANTE_LTR_annotation.gff3")
  s <- readDNAStringSet("./test_data/sample_genome.fasta")

  g <- rtracklayer::import("./test_data/DANTE_LTR_Vfaba_chr5.gff3")
  s <- readDNAStringSet("./test_data/211010_Vfaba_chr5.fasta")
  names(s) <- gsub(" .+", "", names(s))
  ncpus <- 10
  lineage_info <- read.table("databases/lineage_domain_order.csv", sep = "\t", header =
    TRUE, as.is = TRUE)
  source("./R/ltr_utils.R")
}

## ID in g must be unique - this could be a problem if gff is concatenated from multiple files!
## id ID is renamed - rename parent to!
## add chromosom index to disctinguish same IDs
## do this only if IDs are not unique
if (any(duplicated(na.omit(g$ID)))){
  suffix <- as.numeric(seqnames(g))
  g$ID <- ifelse(is.na(g$ID), NA, paste0(g$ID,"_", suffix))
  g$Parent <- ifelse(is.na(g$Parent), NA, paste0(g$Parent,"_", suffix))
}

# get te sequence based on rank

# best quality - split by lineage
s_te <- get_te_sequences(g, s)  # split by 'element quality'
# evaluate best TE -  DLTP grou

# comparison parameters
word_size <- 11
task <- "blastn"
perc_identity <- 80

# best TE
TE_DLTP_info <- analyze_TE(s_te$DLTP,
                           word_size = word_size,
                           ncpus = ncpus,
                           perc_identity = perc_identity,
                           task = task)

# TE rank 2:
TE_DLT_plus_DLP_info <- analyze_TE(c(s_te$DLP, s_te$DLT),
                                   word_size = word_size,
                                   ncpus = ncpus,
                                   perc_identity = perc_identity,
                                   task = task)

TE_D_plus_DL_info <- analyze_TE(c(s_te$DL, s_te$D),
                                word_size = word_size,
                                ncpus = ncpus,
                                perc_identity = perc_identity,
                                task = task)

TE_DLT_plus_DLP_info_DLTP_verified <- compare_TE_datasets(
  c(s_te$DLT, s_te$DLP),
  ncpus = ncpus,
  TE_DLTP_info$seqs_representative,
  word_size = word_size,
  perc_identity = perc_identity,
  task = task)

TE_DLT_plus_DLP_info_multiplicity <- verify_based_on_multiplicity(TE_DLT_plus_DLP_info)
TE_D_plus_DL_info_multiplicity <- verify_based_on_multiplicity(TE_D_plus_DL_info)

# create additional library from rank 2 verified by multiplicity and DLTP
id_for_additional_library <- union(
  TE_DLT_plus_DLP_info_multiplicity$id_ok_mp_verified,
  TE_DLT_plus_DLP_info_DLTP_verified$id_ok_verified)

if (length(id_for_additional_library) > 1) {
  seqs_for_additional_library <- c(s_te$DLP, s_te$DLT)[names(c(s_te$DLP, s_te$DLT)) %in%
                                                         id_for_additional_library]
  seqs_additional_info <- analyze_TE(seqs_for_additional_library, word_size =
    word_size, ncpus = ncpus)
  seqs_representative <- c(
    TE_DLTP_info$seqs_representative,
    seqs_additional_info$seqs_representative
  )
}else {
  if (length(id_for_additional_library) == 1) {
    seqs_representative <- c(
      TE_DLTP_info$seqs_representative,
      c(s_te$DLP, s_te$DLT)[names(c(s_te$DLP, s_te$DLT)) %in% id_for_additional_library]
    )
  }else {
    seqs_representative <- TE_DLTP_info$seqs_representative
  }
}

# TE  rank 3 - verify agains good DLTP
TE_DL_info_DLTP_verified <- compare_TE_datasets(
  s_te$DL,
  seqs_representative, min_coverage = 0.90,
  ncpus = ncpus,
  word_size = word_size,
  task = task,
  perc_identity = perc_identity)

TE_D_info_DLTP_verified <- compare_TE_datasets(
  s_te$D,
  seqs_representative, min_coverage = 0.90,
  ncpus = ncpus,
  word_size = word_size,
  task = task,
  perc_identity = perc_identity)



R <- seq_diversity(seqs_representative)$richness
SI <- seq_diversity(seqs_representative)$shannon_index

# final RM library:
seqs_representative_no_ssr <- seqs_representative[R > 20]

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
  TE_DL_info_DLTP_verified$id_ok_verified,
  TE_D_info_DLTP_verified$id_ok_verified,
  TE_D_plus_DL_info_multiplicity$id_ok_mp_verified),
)

c1 <- g$ID %in% ID[id_of_good_te]
c2 <- sapply(g$Parent, function(x)ifelse(length(x) == 0, "", x)) %in% ID[id_of_good_te]

gff_out <- g[c1 | c2]

gff_te <- gff_out[gff_out$type %in% "transposable_element"]
# remove partial elements
gff_te_with_ltr <- gff_out[gff_out$type %in% "transposable_element" & gff_out$Rank != "D"]


gff_5ltr <- gff_out[gff_out$LTR %in% "5LTR"]
gff_3ltr <- gff_out[gff_out$LTR %in% "3LTR"]

full_te <- getSeqNamed(s, gff_te)
names(full_te) <- paste0(gff_te$ID,":",names(full_te))
ltr5 <-  getSeqNamed(s, gff_5ltr)
names(ltr5) <-  paste0(gff_5ltr$Parent,":",names(ltr5))
ltr3 <-  getSeqNamed(s, gff_3ltr)
names(ltr3) <- paste0(gff_3ltr$Parent,":",names(ltr3))
inc <- gff_te_with_ltr$Rank != "DL"

writeXStringSet(seqs_representative, paste0(opt$output, "_RM_lib_non_redundant.fasta"))
writeXStringSet(full_te, paste0(opt$output, "_RM_lib_full_TE.fasta"))
writeXStringSet(ltr5, paste0(opt$output, "_RM_lib_5LTR.fasta"))
writeXStringSet(ltr3, paste0(opt$output, "_RM_lib_3LTR.fasta"))

export(gff_out, paste0(opt$output, "_clean.gff3"), format = "gff3")

lv <- sort(unique(gff_te_with_ltr$Final_Classification))
te_count <- table(factor(gff_te_with_ltr$Final_Classification, levels=lv))

pdf(paste0(opt$output, "_summary.pdf"), width = 13, height=8, pointsize = 10)
par(mfrow=c(1,2), mar=c(5,7,2,0), las=1)
boxplot(width(gff_te_with_ltr) ~ factor(gff_te_with_ltr$Final_Classification, levels=lv),
        horizontal = TRUE, xlab="length[bp]", ylab="",
        names = paste0(gsub("^.+[|]", "", lv), " (", te_count, ")"),
        main = 'Full TE', at = seq_along(lv)*4
)
boxplot(width(gff_te_with_ltr[inc]) ~ factor(gff_te_with_ltr$Final_Classification[inc], levels=lv),
        horizontal = TRUE, xlab="length[bp]", ylab="",
        names = rep("", length(lv)),
        main = 'Full TE', at = seq_along(lv)*4-1, add=TRUE, col=2
)
par(mar=c(5,0,2,7))
boxplot(width(gff_5ltr) ~ factor(gff_5ltr$Final_Classification, levels=lv),
        horizontal = TRUE, xlab="length[bp]", ylab="",
        names = rep("", length(lv)),
        main = "5'LTR", at = seq_along(lv)*4
)
boxplot(width(gff_5ltr[inc]) ~ factor(gff_5ltr$Final_Classification[inc], levels=lv),
        horizontal = TRUE, xlab="length[bp]", ylab="",
        names = rep("", length(lv)),
        main = "5'LTR", at = seq_along(lv)*4-1, add=TRUE, col=2
)
legend('bottomright', col=c("grey","2"), legend=c("All TE", "TE with PBS/TSD"), pch=15)
dev.off()
