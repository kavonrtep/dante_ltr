#!/usr/bin/env Rscript
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- normalizePath(sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)]))
script_dir <- dirname(script_name)
library(optparse)

parser <- OptionParser()
option_list <- list(
  make_option(c("-g", "--gff3"), action = "store", type = "character",
              help = "gff3 with dante results", default = NULL),
  make_option(c("-s", "--reference_sequence"), action = "store", type = "character",
              help = "reference sequence as fasta", default = NULL),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output file path and prefix", default = NULL),
  make_option(c("-c", "--cpu"), type = "integer", default = 5,
              help = "Number of cpu to use [default %default]", metavar = "number"),
  make_option(c("-M", "--max_missing_domains"), type = "integer", default = 0,
              help = "Maximum number of missing domains is retrotransposon [default %default]",
              metavar = "number")


)
description <- paste(strwrap(""))

epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))


# load packages
suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
  library(parallel)
})


# CONFIGURATION
OFFSET <- 15000
# load configuration files and functions:
lineage_file <- paste0(script_dir, "/databases/lineage_domain_order.csv")
trna_db <- paste0(script_dir, "/databases/tRNAscan-SE_ALL_spliced-no_plus-old-tRNAs_UC_unique-3ends.fasta")
ltr_utils_r <- paste0(script_dir, "/R/ltr_utils.R")
if (file.exists(lineage_file) & file.exists(trna_db)) {
  lineage_info <- read.table(lineage_file, sep = "\t", header = TRUE, as.is = TRUE)
  source(ltr_utils_r)
}else {
  lineage_file <- paste0(script_dir, "/../share/dante_ltr/databases/lineage_domain_order.csv")
  trna_db <- paste0(script_dir, "/../share/dante_ltr/databases/tRNAscan-SE_ALL_spliced-no_plus-old-tRNAs_UC_unique-3ends.fasta")
  ltr_utils_r <- paste0(script_dir, "/../share/dante_ltr/R/ltr_utils.R")
  if (file.exists(lineage_file) & file.exists((trna_db))) {
    lineage_info <- read.table(lineage_file, sep = "\t", header = TRUE, as.is = TRUE)
    source(ltr_utils_r)
  }else(
    stop("configuration files not found")
  )
}


# for testing
if (FALSE) {
  g <- rtracklayer::import("/mnt/raid/454_data/cuscuta/Ceuropea_assembly_v4/0_final_asm_hifiasm+longstitch/repeat_annotation/DANTE_on_CEUR_filtered_short_names.gff3")
  s <- readDNAStringSet("/mnt/raid/454_data/cuscuta/Ceuropea_assembly_v4/0_final_asm_hifiasm+longstitch/asm.bp.p.ctg_scaffolds.short_names.fa")
  lineage_info <- read.table("/mnt/raid/users/petr/workspace/ltr_finder_test/lineage_domain_order.csv", sep = "\t", header = TRUE, as.is = TRUE)

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/ltr_finder_test/test_data/DANTE_filtered_part.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/ltr_finder_test/test_data/Rbp_part.fa")

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data
  /DANTE_Vfaba_chr5.gff3")
  s <- readDNAStringSet("/mnt/ceph/454_data/Vicia_faba_assembly/assembly/ver_210910
  /fasta_parts/211010_Vfaba_chr5.fasta")

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/big_test_data//Cocoa_theobroma_DANTE_filtered.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/dante_ltr/test_data/big_test_data/Cocoa_theobroma_chr1.fasta.gz")

  source("R/ltr_utils.R")

  g <- rtracklayer::import("./test_data/sample_DANTE.gff3")
  s <- readDNAStringSet("./test_data/sample_genome.fasta")
  outfile <- "/mnt/raid/users/petr/workspace/ltr_finder_test/te_with_domains_2.gff3"
  lineage_info <- read.table("databases/lineage_domain_order.csv", sep = "\t", header =
    TRUE, as.is = TRUE)
  trna_db <- "./databases/tRNAscan-SE_ALL_spliced-no_plus-old-tRNAs_UC_unique-3ends.fasta"

}

# MAIN #############################################################

# load data:

cat("reading gff...")
g <- rtracklayer::import(opt$gff3, format = 'gff3')  # DANTE gff3
cat("done\n")
cat("reading fasta...")
s <- readDNAStringSet(opt$reference_sequence)  # genome assembly
cat("done\n")
outfile <- opt$output
# clean sequence names:
names(s) <- gsub(" .+", "", names(s))
lineage_domain <- lineage_info$Domains.order
lineage_offset5prime <- lineage_info$offset5prime
lineage_offset3prime <- lineage_info$offset3prime
ln <- gsub("ss/I", "ss_I", gsub("_", "/", gsub("/", "|", lineage_info$Lineage)))
names(lineage_offset3prime) <-  ln
names(lineage_offset5prime) <-  ln
names(lineage_domain) <- ln


seqlengths(g) <- seqlengths(s)[names(seqlengths(g))]
g <- add_coordinates_of_closest_neighbor(g)

gcl <- get_domain_clusters(g)
gcl_clean <- clean_domain_clusters(gcl)
# glc annotation
gcl_clean_lineage <- sapply(gcl_clean, function(x)  x$Final_Classification[1])
gcl_clean_domains <- sapply(gcl_clean, function(x) ifelse(x$strand[1] == "-",
                                                paste(rev(x$Name), collapse = " "),
                                                paste(x$Name, collapse = " "))
)

dd <- mapply(domain_distance,
             d_query = gcl_clean_domains,
             d_reference = lineage_domain[gcl_clean_lineage])

# get lineages which has correct number and order of domains
# gcl_clean2 <- gcl_clean[gcl_clean_domains == lineage_domain[gcl_clean_lineage]]
gcl_clean2 <- gcl_clean[dd <= opt$max_missing_domains]

gcl_clean_with_domains <- gcl_clean2[check_ranges(gcl_clean2, s)]
gr <- get_ranges(gcl_clean_with_domains)


cat('Number of analyzed regions:\n')
cat('Total number of domain clusters             : ', length(gcl), '\n')
cat('Number of clean clusters                    : ', length(gcl_clean), '\n')
cat('Number of clusters with complete domain set : ', length(gcl_clean_with_domains), '\n')


te_strand <- sapply(gcl_clean_with_domains, function(x)x$strand[1])
te_lineage <- sapply(gcl_clean_with_domains, function(x)x$Final_Classification[1])

max_left_offset <- ifelse(te_strand == "+", lineage_offset5prime[te_lineage], lineage_offset3prime[te_lineage])
max_right_offset <- ifelse(te_strand == "-", lineage_offset5prime[te_lineage], lineage_offset3prime[te_lineage])

grL <- get_ranges_left(gcl_clean_with_domains, max_left_offset)
grR <- get_ranges_right(gcl_clean_with_domains, max_right_offset)

s_left <- getSeq(s, grL)
s_right <- getSeq(s, grR)

# for statistics
RT <- g[g$Name == "RT" & substring(g$Final_Classification, 1, 11) == "Class_I|LTR"]
# cleanup
rm(g)
rm(gcl)
rm(gcl_clean)
rm(gcl_clean2)

names(te_strand) <- paste(seqnames(gr), start(gr), end(gr), sep = "_")
names(s_left) <- paste(seqnames(grL), start(grL), end(grL), sep = "_")
names(s_right) <- paste(seqnames(grR), start(grR), end(grR), sep = "_")
cat('Identification of LTRs...')
TE <- mclapply(seq_along(gr), function(x)get_TE(s_left[x],
                                                s_right[x],
                                                gcl_clean_with_domains[[x]],
                                                gr[x],
                                                grL[x],
                                                grR[x]),
               mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE
)

cat('done.\n')

good_TE <- TE[!sapply(TE, is.null)]
cat('Number of putative TE with identified LTR   :', length(good_TE), '\n')


ID <- paste0('TE_', sprintf("%08d", seq(good_TE)))
gff3_list <- mcmapply(get_te_gff3, g = good_TE, ID = ID, mc.cores = opt$cpu)

cat('Identification of PBS ...')
gff3_list2 <- mclapply(gff3_list, FUN = add_pbs, s = s, trna_db = trna_db, mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE)
cat('done\n')
gff3_out <- do.call(c, gff3_list2)

# define new source
src <- as.character(gff3_out$source)
src[is.na(src)] <- "dante_ltr"
gff3_out$source <- src
gff3_out$Rank <- get_te_rank(gff3_out)
# TODO add attributte specifying individual groups DL, DLT, DLP DLPT gff3

export(gff3_out, con = paste0(outfile, ".gff3"), format = 'gff3')

all_tbl <- get_te_statistics(gff3_out, RT)
write.table(all_tbl, file = paste0(outfile, "_statistics.csv"), sep = "\t", quote = FALSE, row.names = TRUE)
# export fasta files:
s_te <- get_te_sequences(gff3_out, s)

for (i in seq_along(s_te)) {
  outname <- paste0(outfile, "_", names(s_te)[i], ".fasta")
  writeXStringSet(s_te[[i]], filepath = outname)
}

