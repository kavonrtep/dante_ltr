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
              metavar = "number"),
  make_option(c("-L", "--min_relative_length"), type = "numeric", default = 0.6,
              help = "Minimum relative length of protein domain to be considered for retrostransposon detection [default %default]",
              metavar = "number"),
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Debug mode [default %default]")

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
FDM_file <- paste0(script_dir, "/databases/feature_distances_model.RDS")
trna_db <- paste0(script_dir, "/databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_unique-3ends.fasta")
trna_db_hemi <- paste0(script_dir, "/databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_numbered_unique-half-tRNA-20nt.fasta")
ltr_utils_r <- paste0(script_dir, "/R/ltr_utils.R")
if (file.exists(lineage_file) & file.exists(trna_db)) {
  lineage_info <- read.table(lineage_file, sep = "\t", header = TRUE, as.is = TRUE)
  FDM <- readRDS(FDM_file)
  source(ltr_utils_r)
}else {
  # this destination work is installed using conda
  lineage_file <- paste0(script_dir, "/../share/dante_ltr/databases/lineage_domain_order.csv")
  FDM_file <- paste0(script_dir, "/../share/dante_ltr/databases/feature_distances_model.RDS")
  trna_db <- paste0(script_dir, "/../share/dante_ltr/databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_unique-3ends.fasta")
  trna_db_half <- paste0(script_dir, "/../share/dante_ltr/databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_numbered_unique-half-tRNA-20nt.fasta")
  ltr_utils_r <- paste0(script_dir, "/../share/dante_ltr/R/ltr_utils.R")
  if (file.exists(lineage_file) & file.exists((trna_db))) {
    lineage_info <- read.table(lineage_file, sep = "\t", header = TRUE, as.is = TRUE)
    source(ltr_utils_r)
    FDM <- readRDS(FDM_file)
  }else(
    stop("configuration files not found")
  )
}

if (opt$debug) {
  # create additional directory with extra information
  debug_dir <- paste0(opt$output, "_debug")
  dir.create(debug_dir, showWarnings = FALSE)
}

# for testing
if (FALSE) {
  g <- rtracklayer::import("/mnt/raid/454_data/cuscuta/Ceuropea_assembly_v4/0_final_asm_hifiasm+longstitch/repeat_annotation/DANTE_on_CEUR_filtered_short_names.gff3")
  s <- readDNAStringSet("/mnt/raid/454_data/cuscuta/Ceuropea_assembly_v4/0_final_asm_hifiasm+longstitch/asm.bp.p.ctg_scaffolds.short_names.fa")
  lineage_info <- read.table("/mnt/raid/users/petr/workspace/ltr_finder_test/lineage_domain_order.csv", sep = "\t", header = TRUE, as.is = TRUE)

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/sample_DANTE_unfiltered.gff3")
  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/ltr_finder_test/test_data/DANTE_filtered_part.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/ltr_finder_test/test_data/Rbp_part.fa")

  # oriza
  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/big_test_data/DANTE_full_oryza.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/dante_ltr/test_data/big_test_data/o_sativa_msu7.0.fasta")

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data
  /DANTE_Vfaba_chr5.gff3")
  s <- readDNAStringSet("/mnt/ceph/454_data/Vicia_faba_assembly/assembly/ver_210910
  /fasta_parts/211010_Vfaba_chr5.fasta")

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/big_test_data//Cocoa_theobroma_DANTE_filtered.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/dante_ltr/test_data/big_test_data/Cocoa_theobroma_chr1.fasta.gz")
  # test on bigger data:

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/tmp/DANTE_unfiltered/1.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/dante_ltr/test_data/tmp/fasta_parts/1.fasta")

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/ikeros_test_dante.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/dante_ltr/test_data/ikeros_test.fasta")

  g <- rtracklayer::import("/mnt/raid/454_data/dante/Intact_LTR_dataset/assembly_o_sativa/Oryza_sativa_msu7_dante.gff3")
  s <- readDNAStringSet("/mnt/raid/454_data/dante/Intact_LTR_dataset/assembly_o_sativa/Oryza_sativa_msu7.fasta")


  # test - TE insertion in satellite:
  g <- rtracklayer::import("/mnt/raid/454_data/dante/Intact_LTR_dataset/assembly_o_sativa/crm_test_dante.gff3")
  s <- readDNAStringSet("/mnt/raid/454_data/dante/Intact_LTR_dataset/assembly_o_sativa/crm_test.fasta")


  source("R/ltr_utils.R")
  ## feature distance model
  FDM <- readRDS("./databases/feature_distances_model.RDS")
  g <- rtracklayer::import("./test_data/sample_DANTE.gff3")
  s <- readDNAStringSet("./test_data/sample_genome.fasta")
  outfile <- "/mnt/raid/users/petr/workspace/ltr_finder_test/te_with_domains_2.gff3"
  lineage_info <- read.table("databases/lineage_domain_order.csv", sep = "\t", header =
    TRUE, as.is = TRUE)
  trna_db <- "./databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_unique-3ends.fasta"
  trna_db_hemi <- "./databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_unique-3ends_hemi.fasta"
  opt <- list(min_relative_length=0.6, cpu = 8, max_missing_domains = 0, debug = FALSE)

}

# MAIN #############################################################

# load data:
cat("reading gff...")
g <- rtracklayer::import(opt$gff3, format = 'gff3')  # DANTE gff3


g <- CHD_CHDCR_correction(g)
cat("done\n")
cat("reading fasta...")
s <- readDNAStringSet(opt$reference_sequence)  # genome assembly
cat("done\n")
outfile <- opt$output
# clean sequence names:
names(s) <- gsub(" .+", "", names(s))
lineage_domain <- lineage_info$Domains.order
lineage_domain_span <- lineage_info$domain_span
lineage_ltr_mean_length <- lineage_info$ltr_length
lineage_offset5prime <- lineage_info$offset5prime
lineage_offset3prime <- lineage_info$offset3prime
ln <- gsub("ss/I", "ss_I", gsub("_", "/", gsub("/", "|", lineage_info$Lineage)))
names(lineage_offset3prime) <-  ln
names(lineage_offset5prime) <-  ln
names(lineage_domain) <- ln
names(lineage_domain_span) <- ln
names(lineage_ltr_mean_length) <- ln
lineage_domains_sequence <- unlist(mapply(function(d,l) {
  paste(strsplit(d, " ")[[1]], ":", l, sep = "")
}, d = lineage_domain, l = names(lineage_domain)))

# this repeat block is run just once
# it can breaks in eny point if zero TE is found
repeat{
  # remove reduntaunt attributes:
  g$Region_Hits_Classifications <- NULL
  g$Region_Hits_Classifications_ <- NULL

  # prefilering analysis gff3
  # clusters of domain are identified, domains which are close to expected other domains
  # are filterer with less strict criteria
  g <- sort(g, by = ~ seqnames * start)
  # add info about domain order:
  g$domain_order <- as.numeric(factor(paste(g$Name, g$Final_Classification, sep = ":"), levels = lineage_domains_sequence))
  # set NA to 0 in  g$domain_order ( some domains are not fromm ClassI elements
  g$domain_order[is.na(g$domain_order)] <- 0

  cls_prefilter <- paste(get_domain_clusters_alt(g, FDM), get_domain_clusters(g))
  g$cls_prefilter <- cls_prefilter
  neighbors_count <- c(head(cls_prefilter, -1) == cls_prefilter[-1], 0)
  neighbors_count <- c(0, head(neighbors_count,-1)) + neighbors_count
  # neighbors_count is number of direct neighbors witch are in the same cluster
  # it is used to improve domain filtering. If there are 2 neigbors, filtering can
  # be less strict

  # filtering
  g$neighbors_count <- neighbors_count
  g <- dante_filtering(g, Relative_Length = opt$min_relative_length) # default
  seqlengths(g) <- seqlengths(s)[names(seqlengths(g))]
  g <- add_coordinates_of_closest_neighbor(g)

  # add info about domain order:
  g$domain_order <- as.numeric(factor(paste(g$Name, g$Final_Classification, sep = ":"), levels = lineage_domains_sequence))
  # set NA to 0 in  g$domain_order ( some domains are not fromm ClassI elements
  g$domain_order[is.na(g$domain_order)] <- 0

  # NOTE - some operation is faster of GrangesList but some on list of data.frames
  # this is primary clusteing
  cls <- get_domain_clusters(g)
  gcl <- split(as.data.frame(g), cls)
  # gcl_as_GRL <- split(g, cls)  # delete?
  cls_alt <- get_domain_clusters_alt(g, FDM)
  g$Cluster <- as.numeric(factor(cls_alt))
  gcl_alt <- split(as.data.frame(g), cls_alt)

  TE_partial <-  GRanges(seqnames =  sapply(gcl_alt, function(x) x$seqnames[1]),
                         Name = sapply(gcl_alt, function(x) x$Final_Classification[1]),
                         Final_Classification = sapply(gcl_alt, function(x) x$Final_Classification[1]),
                         ID = sapply(gcl_alt, function(x) paste0("TE_partial_", sprintf("%08d", x$Cluster[1]))),
                         strand = sapply(gcl_alt, function(x) x$strand[1]),
                         Ndomains = sapply(gcl_alt, function(x) nrow(x)),
                         type = "transposable_element",
                         source = "dante_ltr",
                         Rank="D",
                         IRanges(start = sapply(gcl_alt, function(x) min(x$start)),
                                 end = sapply(gcl_alt, function(x) max(x$end)))
  )
  if (opt$debug) {
    export(TE_partial, con = paste0(debug_dir,"/TE_partial_all.gff3"), format = "gff3")
  }
  g$Ndomains_in_cluster <- count_occurences_for_each_element(g$Cluster)
  g$Parent <- paste0("TE_partial_", sprintf("%08d", g$Cluster))
  g$Rank <- "D"
  # for statistics
  RT <- g[g$Name == "RT" & substring(g$Final_Classification, 1, 11) == "Class_I|LTR"]

  # keep only partial TE with more than one domain
  TE_partial_with_more_than_one_domain <- TE_partial[TE_partial$Ndomains > 1]
  g_with_more_than_one_domain <- g[as.vector(g$Ndomains_in_cluster > 1)]

  # first filtering  - remove TEs with low number of domains
  gcl_clean <- clean_domain_clusters(gcl, lineage_domain_span, min_domains = 5 - opt$max_missing_domains)

  # glc annotation
  gcl_clean_lineage <- sapply(gcl_clean, function(x)  x$Final_Classification[1])
  gcl_clean_domains <- sapply(gcl_clean, function(x) ifelse(x$strand[1] == "-",
                                                  paste(rev(x$Name), collapse = " "),
                                                  paste(x$Name, collapse = " "))
  )
  # compare detected domains with domains in lineages from REXdb database
  dd <- mapply(domain_distance,
               d_query = gcl_clean_domains,
               d_reference = lineage_domain[gcl_clean_lineage])

  # get lineages which has correct number and order of domains
  # gcl_clean2 <- gcl_clean[gcl_clean_domains == lineage_domain[gcl_clean_lineage]]


  gcl_clean2 <- gcl_clean[dd <= opt$max_missing_domains]
  if(length(gcl_clean2) == 0) {
    cat("No complete TE found\n")
    good_TE <- list()
    break
  }

  cat("Number of complete TE found:")
  cat(length(gcl_clean2), "\n")
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
  save.image(file = "debug.RData")
  grR <- get_ranges_right(gcl_clean_with_domains, offset=max_right_offset, SL = seqlengths(s))
  s_left <- getSeq(s, grL)
  s_right <- getSeq(s, grR)


  expected_ltr_length <- lineage_ltr_mean_length[sapply(gcl_clean_with_domains, function (x)x$Final_Classification[1])]

  # cleanup
  #rm(g)
  rm(gcl)
  rm(gcl_clean)
  rm(gcl_clean2)

  names(te_strand) <- paste(seqnames(gr), start(gr), end(gr), sep = "_")
  names(s_left) <- paste(seqnames(grL), start(grL), end(grL), sep = "_")
  names(s_right) <- paste(seqnames(grR), start(grR), end(grR), sep = "_")
  cat('Identification of LTRs...')
  if (opt$debug){
    save.image(file = paste0("debug_dante_ltr.RData"))
  }
  TE <- mclapply(seq_along(gr), function(x)get_TE(s_left[x],
                                                  s_right[x],
                                                  gcl_clean_with_domains[[x]],
                                                  gr[x],
                                                  grL[x],
                                                  grR[x],
                                                  expected_ltr_length[x]),
                 mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE
  )
  cat('done.\n')
  good_TE <- TE[!sapply(TE, is.null)]
  cat('Number of putative TE with identified LTR   :', length(good_TE), '\n')
  if (opt$debug) {
    saveRDS(good_TE, paste0(debug_dir, "/good_TE.rds"))
  }
  break
  }

if (length(good_TE)>0){   # handle empty list
  ID <- paste0('TE_', sprintf("%08d", seq(good_TE)))
  gff3_list <- mcmapply(get_te_gff3, g = good_TE, ID = ID, mc.cores = opt$cpu)
  cat('Identification of PBS ...')
  gff3_list2 <- mclapply(gff3_list, FUN = add_pbs, s = s, trna_db = trna_db, mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE)
  # gff3_list2 <- lapply(gff3_list, FUN = add_pbs, s = s, trna_db = trna_db)
  cat('done\n')


  gff3_list2_pbs_positive <- gff3_list2[sapply(gff3_list2, function(x) "primer_binding_site" %in% x$type)]

  gff3_list2_pbs_negative <- gff3_list[!sapply(gff3_list2, function(x) "primer_binding_site" %in% x$type)]

  cat("Identification of PBS - hemi")
  gff3_list3 <- mclapply(gff3_list2_pbs_negative, FUN = add_pbs_hemi, s = s,
                         trna_db = trna_db_hemi,
                         mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE)
  cat("done\n")

  gff3_out <- do.call(c, append(gff3_list3, gff3_list2_pbs_positive))

  # define new source
  src <- as.character(gff3_out$source)
  src[is.na(src)] <- "dante_ltr"
  gff3_out$source <- src
  gff3_out$Rank <- get_te_rank(gff3_out)
  # add partial TEs but first remove all ovelaping
  if (opt$debug) {
    save.image(file = paste0("debug_dante_ltr.RData"))
  }
  # use complete TE as mask for partial TE

  TE_partial_parent_part <- trim_gr(TE_partial_with_more_than_one_domain, gff3_out)

  if (!is.null(TE_partial_parent_part)) {
    TE_partial_domain_part <-  g[g$Parent %in% TE_partial_parent_part$ID]
    # and trim it
    TE_partial_domain_part <-  trim_gr(TE_partial_domain_part, gff3_out)
    gff3_out <- sort(merge_gr(gff3_out, TE_partial_parent_part, TE_partial_domain_part), by = ~ seqnames * start)
  }else{
    gff3_out <- sort(gff3_out, by = ~ seqnames * start)

  }


  # this is to convert Parent from characterList
  gff3_out$Parent <- as.character(gff3_out$Parent)

}else{
  # but this could be a problem if there are no TEs in the sequence
  if (length(TE_partial_with_more_than_one_domain)>0){
    TE_partial_parent_part <- TE_partial_with_more_than_one_domain
    TE_partial_domain_part <-  g[g$Parent %in% TE_partial_parent_part$ID]
    gff3_out <- sort(c(TE_partial_domain_part, TE_partial_parent_part), by = ~ seqnames * start)
  }else{
    gff3_out <- NULL
  }
}

if (is.null(gff3_out)){
  cat('No TEs found.\n')
}else{
# modify ID and Parent - add seqname - this will ensure it is unique is done on chunk level
  gff3_out$ID[!is.na(gff3_out$ID)] <- paste0(gff3_out$ID[!is.na(gff3_out$ID)], "_", seqnames(gff3_out)[!is.na(gff3_out$ID)])
  gff3_out$Parent[!is.na(gff3_out$Parent)] <- paste0(gff3_out$Parent[!is.na(gff3_out$Parent)], "_", seqnames(gff3_out)[!is.na(gff3_out$Parent)])
  gff3_out <- convert_gr_Lists_to_Vectors(gff3_out)
  gff3_out <- revert_CHDCR_correction(gff3_out)
  export(gff3_out, con = paste0(outfile, ".gff3"), format = 'gff3')
  all_tbl <- get_te_statistics(gff3_out, RT)
  all_tbl <- cbind(Classification = rownames(all_tbl), all_tbl)
  write.table(all_tbl, file = paste0(outfile, "_statistics.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # export fasta files:
  s_te <- get_te_sequences(gff3_out, s)
  for (i in seq_along(s_te)) {
    outname <- paste0(outfile, "_", names(s_te)[i], ".fasta")
    writeXStringSet(s_te[[i]], filepath = outname)
  }
}

