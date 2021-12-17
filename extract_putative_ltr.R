#!/usr/bin/env Rscript

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- normalizePath(sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)]))
script_dir <- dirname(script_name)
library(optparse)


parser = OptionParser()
option_list = list(
  make_option(c("-g", "--gff3"), action = "store", type = "character",
              help = "gff3 with dante results", default = NULL),
  make_option(c("-s", "--reference_sequence"), action = "store", type = "character",
              help = "reference sequence as fasta", default = NULL),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output gff3 file name", default = NULL),
  make_option(c("-c", "--cpu"), type = "integer", default = 5,
              help = "Number of cpu to use [default %default]", metavar = "number")

)

## params
OFFSET = 15000

description = paste(strwrap(""))

epilogue = ""
parser = OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                      usage = "usage: %prog COMMAND [OPTIONS]")
opt = parse_args(parser, args = commandArgs(TRUE))


# load packages
suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
  library(parallel)
})

cat("reading gff...")
g <- rtracklayer::import(opt$gff3)  # DANTE gff3
cat("done\n")

cat("reading fasta...")
s <- readDNAStringSet(opt$reference_sequence)  # genome assembly
cat("done\n")
lineage_info <- read.table(paste0(script_dir, "/lineage_domain_order.csv"), sep = "\t", header = TRUE, as.is = TRUE)
outfile <- opt$output


# for testing
if (FALSE) {
  g <- rtracklayer::import("/mnt/raid/454_data/cuscuta/Ceuropea_assembly_v4/0_final_asm_hifiasm+longstitch/repeat_annotation/DANTE_on_CEUR_filtered_short_names.gff3")
  s <- readDNAStringSet("/mnt/raid/454_data/cuscuta/Ceuropea_assembly_v4/0_final_asm_hifiasm+longstitch/asm.bp.p.ctg_scaffolds.short_names.fa")
  lineage_info <- read.table("/mnt/raid/users/petr/workspace/ltr_finder_test/lineage_domain_order.csv", sep = "\t", header = TRUE, as.is = TRUE)

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/ltr_finder_test/test_data/DANTE_filtered_part.gff3")
  s <- readDNAStringSet("/mnt/raid/users/petr/workspace/ltr_finder_test/test_data/Rbp_part.fa")

  g <- rtracklayer::import("/mnt/raid/users/petr/workspace/dante_ltr/test_data/DANTE_Vfaba_chr5.gff3")
  s <- readDNAStringSet("/mnt/ceph/454_data/Vicia_faba_assembly/assembly/ver_210910/fasta_parts/211010_Vfaba_chr5.fasta")
  lineage_info <- read.table("/mnt/raid/users/petr/workspace/ltr_finder_test/lineage_domain_order.csv", sep = "\t", header = TRUE, as.is = TRUE)
  outfile <- "/mnt/raid/users/petr/workspace/ltr_finder_test/te_with_domains_2.gff3"


}
if (file.exists(outfile)) {
  stop("output file already exists")
}

# clean sequence names:

names(s) <- gsub(" .+", "", names(s))
lineage_domain <- lineage_info$Domains.order
names(lineage_domain) <- gsub("ss/I", "ss_I", gsub("_", "/", gsub("/", "|", lineage_info$Lineage)))

# functions
get_coordinates_of_closest_neighbor <- function(gff) {
  gff <- gff[order(seqnames(gff), start(gff))]
  # split to chromosomes:
  gff_parts <- split(gff, seqnames(gff))
  upstreams <- c(sapply(gff_parts, function(x) c(1, head(end(x), -1))))
  downstreams <- c(sapply(gff_parts, function(x) c(start(x)[-1], seqlengths(x)[runValue(seqnames(x[1]))])))
  gff_updated <- unlist(gff_parts)
  gff_updated$upstream_domain <- unlist(upstreams)
  gff_updated$downstream_domain <- unlist(downstreams)
  names(gff_updated) <- NULL
  return(gff_updated)
}

get_domain_clusters <- function(gff) {
  # get consecutive domains from same linage
  # must be sorted!
  # TODO - sorting at this point is not probably necessary - remove and test
  # gff <- gff[order(seqnames(gff), start(gff))]
  # split before GAG (+) or after (- strand)
  gag_plus <- as.numeric(cumsum(gff$Name == "GAG" & strand(gff) == '+'))
  gag_minus <- rev(as.numeric(cumsum(rev(gff$Name == "GAG" & strand(gff) == '-'))))
  # split based on classification - must be same and consecutive
  x <- rle(gff$Final_Classification)
  # split on strand change
  s <- rep(seq_along(runLength(strand(gff))), runLength(strand(gff)))
  domain_cluster <- paste0(rep(seq_along(x$lengths), x$lengths), "_", seqnames(gff),
                           "_", gag_plus, "_", gag_minus, "_", s)
  gff_clusters <- split(as.data.frame(gff), factor(domain_cluster, levels = unique(domain_cluster)))
  gff_clusters
}

clean_domain_clusters <- function(gcl) {
  ## remove clusters wich does not have enough domains or domains
  ## are on different strand
  N_domains <- sapply(gcl, nrow)
  N_unique_domains <- sapply(gcl, function(x)length(unique(x$Name)))
  S <- sapply(gcl, function(x)paste(sort(unique(x$strand)), collapse = " "))
  S_OK <- S %in% c("+", "-")
  min_domains <- 5
  maxlength <- 15000 # max span between domains
  span <- sapply(gcl, function(x)max(x$end) - min(x$start))
  cond1 <- S_OK &
    N_unique_domains == N_domains &
    N_domains >= min_domains &
    span <= maxlength
  return(gcl[cond1])
}


check_ranges <- function(gx, s, offset = OFFSET) {
  START <- sapply(gx, function(x)min(x$start)) - offset
  END <- sapply(gx, function(x)max(x$end)) + offset
  MAX <- seqlengths(s)[sapply(gx, function(x)as.character(x$seqnames[1]))]
  good_ranges <- (START > 0) & (END <= MAX)
  return(good_ranges)
}

get_ranges <- function(gx, offset = OFFSET) {
  S <- sapply(gx, function(x)min(x$start))
  E <- sapply(gx, function(x)max(x$end))
  gr <- GRanges(seqnames = sapply(gx, function(x)x$seqnames[1]), IRanges(start = S - offset, end = E + offset))
}

get_ranges_left <- function(gx, offset = OFFSET, offset2 = 300) {
  S <- sapply(gx, function(x)min(x$start))
  max_offset <- S - sapply(gx, function(x)min(x$upstream_domain))
  offset_adjusted <- ifelse(max_offset < offset, max_offset, offset)
  gr <- GRanges(seqnames = sapply(gx, function(x)x$seqnames[1]), IRanges(start = S - offset_adjusted, end = S + offset2))
  return(gr)
}

get_ranges_right <- function(gx, offset = OFFSET, offset2 = 300) {
  E <- sapply(gx, function(x)max(x$end))
  max_offset <- sapply(gx, function(x)max(x$downstream_domain)) - E
  offset_adjusted <- ifelse(max_offset < offset, max_offset, offset)
  gr <- GRanges(seqnames = sapply(gx, function(x)x$seqnames[1]), IRanges(start = E - offset2, end = E + offset_adjusted))
  return(gr)
}

firstTG <- function(ss) {
  x <- matchPattern("TG", ss)
  if (length(x) == 0) {
    return(0)
  }else {
    return(min(start(x)))
  }
}

lastCA <- function(ss) {
  x <- matchPattern("CA", ss)
  if (length(x) == 0) {
    return(0)
  }else {
    return(max(start(x)))
  }
}

trim2TGAC <- function(bl) {
  for (i in 1:nrow(bl)) {
    tg_L <- firstTG(bl$qseq[i])
    tg_R <- firstTG(bl$sseq[i])
    ca_L <- lastCA(bl$qseq[i])
    ca_R <- lastCA(bl$sseq[i])
    e_dist <- bl$length[i] - ca_R
    no_match <- any(tg_L == 0, tg_R == 0, ca_L == 0, ca_R == 0)
    if (!no_match &
      tg_L == tg_R &
      ca_L == ca_R &
      tg_L < 8 &
      e_dist < 8) {
      ## trim alignment
      bl[i,] <- trim_blast_table(bl[i,], tg_L, e_dist - 1)
    }
  }
  return(bl)
}

trim_blast_table <- function(b, T1, T2) {
  b$qstart <- b$qstart + T1
  b$sstart <- b$sstart + T1
  b$qend <- b$qend - T2
  b$send <- b$send - T2
  b$sseq <- substring(b$sseq, T1, b$length - T2)
  b$qseq <- substring(b$qseq, T1, b$length - T2)
  b$length <- nchar(b$sseq)
  return(b)
}

blast <- function(s1, s2) {
  tmp1 <- tempfile()
  tmp2 <- tempfile()
  tmp_out <- tempfile()
  writeXStringSet(DNAStringSet(s1), tmp1)
  writeXStringSet(DNAStringSet(s2), tmp2)
  # alternative blast:
  cmd <- paste("blastn -task blastn -word_size 7 -dust no -gapextend 1 -gapopen 2 -reward 1 -penalty -1",
               " -query ", tmp1, ' -subject ', tmp2, ' -strand plus ',
               '-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"',
               '  -out', tmp_out)

  system(cmd)
  out_raw <- read.table(tmp_out, as.is = TRUE, sep = "\t",
                        col.names = strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq", split = ' ')[[1]])

  if (nrow(out_raw) == 0) {
    return(out_raw)
  }
  out <- trim2TGAC(out_raw)
  # remove alingment shorted that
  out <- out[out$length > 100,]
  if (nrow(out) == 0) {
    return(out)
  }
  ## filter for TGCA
  TG_L <- substring(out$qseq, 1, 2) == "TG"
  TG_R <- substring(out$sseq, 1, 2) == "TG"
  CA_L <- substring(out$qseq, out$length - 1, out$length) == "CA"
  CA_R <- substring(out$sseq, out$length - 1, out$length) == "CA"
  cond <- TG_L & TG_R & CA_L & CA_R
  out <- out[cond, , drop = FALSE]
  # TODO - detele all temp files!
  return(out)
}

get_correct_coordinates <- function(b) {
  do.call(rbind, strsplit(b$qaccver, split = "_"))
}

evaluate_ltr <- function(GR, GRL, GRR, blast_line, Lseq, Rseq) {
  LTR_L <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRL) + blast_line$qstart - 1,
                                               end = start(GRL) + blast_line$qend - 1))
  LTR_R <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRR) + blast_line$sstart - 1,
                                               end = start(GRR) + blast_line$send - 1))

  TSD_L <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRL) + blast_line$qstart - 3:10,
                                               end = start(GRL) + blast_line$qstart - 3))
  TSD_R <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRR) + blast_line$send,
                                               end = start(GRR) + blast_line$send + 0:7))

  TSD_L_seq <- DNAStringSet(substring(Lseq, blast_line$qstart - 2:9, blast_line$qstart - 2))
  TSD_R_seq <- DNAStringSet(substring(Rseq, blast_line$send + 1, blast_line$send + 1:8))

  matching_tsd <- TSD_R_seq == TSD_L_seq
  matching_tsd[1:3] <- FALSE # remove short tsd
  p <- which(matching_tsd)
  if (length(p) > 0) {
    TSD_Length <- max(p)
    TSD_sequence <- TSD_L_seq[TSD_Length]
    TSD_position <- append(TSD_R[TSD_Length], TSD_L[TSD_Length])
  }else {
    TSD_Length <- 0
    TSD_sequence <- ""
    TSD_position <- NA
  }

  TE_Length <- end(LTR_R) - start(LTR_L)
  LTR_Identity <- blast_line$pident
  out <- list(TSD_position = TSD_position, TSD_sequence = TSD_sequence, TSD_Length = TSD_Length,
              LTR_R_position = LTR_R, LTR_L_position = LTR_L, TE_Length = TE_Length, LTR_Identity = LTR_Identity)
  return(out)
}

get_best_ltr <- function(x) {
  tsd_ok <- sapply(x, function(k)k$TSD_Length > 3)
  te_length_ok <- sapply(x, function(k)k$TE_Length < 30000)
  ltr_length_ok <- sapply(x, function(k)width(k$LTR_R_position) >= 100 & width(k$LTR_L_position) >= 100)
  if (sum(tsd_ok & te_length_ok & ltr_length_ok) >= 1) {
    # return the first one (best bitscore)
    return(x[tsd_ok & te_length_ok][1])
  }
  if (any(te_length_ok & ltr_length_ok)) {
    return(x[te_length_ok & ltr_length_ok][1])
  }else {
    return(NULL)
  }
}

get_te_gff3 <- function(g, ID) {
  D <- makeGRangesFromDataFrame(g$domain, keep.extra.columns = TRUE)
  sn <- seqnames(D)[1]
  S <- strand(D)[1]
  TE <- GRanges(seqnames = sn,
                IRanges(start = start(g$ltr_info[[1]]$LTR_L_position),
                        end = end(g$ltr_info[[1]]$LTR_R_position)), strand = S)
  TE$type <- "transposable_element"
  TE$ID <- ID

  if (as.character(S) == "+") {
    LTR_5 <- g$ltr_info[[1]]$LTR_L
    LTR_3 <- g$ltr_info[[1]]$LTR_R
  }else {
    LTR_3 <- g$ltr_info[[1]]$LTR_L
    LTR_5 <- g$ltr_info[[1]]$LTR_R
  }

  LTR_5$type <- "long_terminal_repeat"
  LTR_3$type <- "long_terminal_repeat"
  strand(LTR_3) <- S
  strand(LTR_5) <- S
  LTR_3$Parent <- ID
  LTR_5$Parent <- ID
  LTR_3$Final_Classification <- D$Final_Classification[1]
  LTR_5$Final_Classification <- D$Final_Classification[1]
  LTR_5$LTR_Identity <- g$ltr_info[[1]]$LTR_Identity
  LTR_3$LTR_Identity <- g$ltr_info[[1]]$LTR_Identity

  TE$LTR_Identity <- g$ltr_info[[1]]$LTR_Identity
  TE$LTR5_length <- width(LTR_5)
  TE$LTR3_length <- width(LTR_3)

  if (is.na(g$ltr_info[[1]]$TSD_position)[1]) {
    # no TSD found
    TSD <- NULL
    TE$TSD <- 'not_found'
  }else {
    TSD <- g$ltr_info[[1]]$TSD_position
    TSD$type <- "target_site_duplication"
    TSD$Parent <- ID
    TE$TSD <- as.character(g$ltr_info[[1]]$TSD_sequence)
  }


  TE$Final_Classification <- D$Final_Classification[1]

  D$Parent <- ID
  out <- c(TE, LTR_3, LTR_5, D, TSD)
  return(out)
}

get_TE <- function(Lseq, Rseq, domains_gff, GR, GRL, GRR) {
  xx <- blast(Lseq, Rseq)
  if (nrow(xx) == 0) {
    return(NULL)
  }else {
    ltr_tmp <- list()
    for (j in 1:nrow(xx)) {
      ltr_tmp[[j]] <- evaluate_ltr(GR, GRL, GRR, xx[j, , drop = FALSE], Lseq, Rseq)
    }
    ltr <- get_best_ltr(ltr_tmp)
    if (length(ltr) == 0) {
      return(NULL)
      ## add good ltr
    }else {
      return(list(domain = domains_gff, ltr_info = ltr, blast_out = xx))
    }
  }
}


# MAIN #############################################################
# sort and add upstrean and downstream

seqlengths(g) <- seqlengths(s)[names(seqlengths(g))]
g <- get_coordinates_of_closest_neighbor(g)
gcl <- get_domain_clusters(g)

gcl_clean <- clean_domain_clusters(gcl)
# glc annotation
lineage <- sapply(gcl_clean, function(x)  x$Final_Classification[1])
domains <- sapply(gcl_clean, function(x) ifelse(x$strand[1] == "-",
                                                paste(rev(x$Name), collapse = " "),
                                                paste(x$Name, collapse = " "))
)
# get lineages which has correct number and order of domains
gcl_clean2 <- gcl_clean[domains == lineage_domain[lineage]]
gcl_clean_with_domains <- gcl_clean2[check_ranges(gcl_clean2, s)]
gr <- get_ranges(gcl_clean_with_domains)


print('dataset size')
print(length(gcl))
print(length(gcl_clean))
print(length(gcl_clean_with_domains))


te_strand <- sapply(gcl_clean_with_domains, function(x)x$strand[1])
grL <- get_ranges_left(gcl_clean_with_domains)
grR <- get_ranges_right(gcl_clean_with_domains)

s_left <- getSeq(s, grL)
s_right <- getSeq(s, grR)

# for statistics
RT <- g[g$Name == "RT" & substring(g$Final_Classification, 1, 11) == "Class_I|LTR"]
input_TE_RT_domain_tbl <- sort(table(RT$Final_Classification), decreasing = TRUE)

# cleanup
#gc()
#print(sort( sapply(ls(),function(x){object.size(get(x))})))
rm(s)
rm(g)
rm(gcl)
rm(gcl_clean)
rm(gcl_clean2)
rm(RT)
gc()
#print(sort( sapply(ls(),function(x){object.size(get(x))})))

names(te_strand) <- paste(seqnames(gr), start(gr), end(gr), sep = "_")
names(s_left) <- paste(seqnames(grL), start(grL), end(grL), sep = "_")
names(s_right) <- paste(seqnames(grR), start(grR), end(grR), sep = "_")


print(opt$cpu)
TE <- mclapply(seq_along(gr), function(x)get_TE(s_left[x],
                                                s_right[x],
                                                gcl_clean_with_domains[[x]],
                                                gr[x],
                                                grL[x],
                                                grR[x]),
               mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE
)
good_TE <- TE[!sapply(TE, is.null)]
print(length(good_TE))

# testing
#x=6
# get_TE(s_left[x],s_right[x], gcl_clean_with_domains[[x]],gr[x],grL[x],grR[x])
# get_te_gff3(get_TE(s_left[x],s_right[x], gcl_clean_with_domains[[x]],gr[x],grL[x],grR[x]), 'testid')


ID <- paste0('TE_', sprintf("%08d", seq(good_TE)))
gff3_list <- mcmapply(get_te_gff3, g = good_TE, ID = ID, mc.cores = opt$cpu)
gff3_out <- do.call(c, gff3_list)

print(table(gff3_out$type))

export(gff3_out, con = outfile, format = 'gff3')


# TODO remove bad TE containing extra domain

# summary statistics
input_TE_complete_domain_tbl <- table(sapply(gcl_clean_with_domains, function(x)x$Final_Classification[1]))
good_TE_tbl <- table(sapply(good_TE, function(x)x$domain$Final_Classification[1]))

all_tbl <- data.frame(
  lineage = names(input_TE_RT_domain_tbl),
  RT = as.integer(input_TE_RT_domain_tbl),
  complete_TE_domains = as.integer(input_TE_complete_domain_tbl[names(input_TE_RT_domain_tbl)]),
  TE_with_domains_ltr_tsd = as.integer(good_TE_tbl[names(input_TE_RT_domain_tbl)])
)

write.table(all_tbl, file = paste0(outfile, "_statistics.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
#  <- read_tsv(".txt")
