#!/usr/bin/env Rscript
# parse arguments
library(optparse)
opt_list <- list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, help="Input fasta file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="output directory", metavar="character"),
  make_option(c("-m", "--min_coverage"), type="numeric", default=3, help="minimal number of sequences in a cluster", metavar="numeric"),
  make_option(c("-t", "--threads"), type="numeric", default=1, help="number of threads", metavar="numeric"),
  make_option(c("-p", "--proportion_min"), type="numeric", default=0.95, help="minimal proportion of the main class in a cluster", metavar="numeric"),
  make_option(c("-a", "--annotation_conflict"), type="character", default="strict",
              help=paste("how to treat a cluster whose members carry different",
                         "classifications: 'strict' (default, historical) drops",
                         "anything whose label is not its majority; 'nested'",
                         "keeps ancestor/descendant chains and may promote them",
                         "to the deepest supported label [default %default]"),
              metavar="character"),
  make_option(c("--lineage_promotion_min_elements"), type="numeric", default=2,
              help="under 'nested', minimum distinct elements carrying the deepest label [default %default]",
              metavar="numeric"),
  make_option(c("--lineage_promotion_min_share"), type="numeric", default=0.25,
              help="under 'nested', minimum share of the cluster's elements carrying it [default %default]",
              metavar="numeric")
)

calculate_segments <- function(LEN, seqname,  overlap = 100, approx_window_size = 1000) {
  N <- round(LEN/approx_window_size)
  if (N == 0) N <- 1
  # Calculate the segment length
  L <- round((LEN + (N - 1) * overlap) / N)
  # Calculate the start positions of each segment
  starts <- 1 + (0:(N - 1)) * (L - overlap)
  # Calculate the end positions of each segment
  ends <- starts + L - 1
  # Make sure the last segment doesn't go over the total length
  ends[N] <- LEN
  # Make sure all segments have the same length (may reduce the length of the last segment)
  diff <- max(diff(c(starts, LEN + 1))) - 1
  ends <- starts + diff - 1
  ends[N] <- min(ends[N], LEN)
  # Return the start and end positions
  return(data.frame(seqname = seqname, start = starts, end = ends))
}


slide_and_cut <- function(s){
  L <- nchar(s)
  N <- names(s)
  dfs <- mapply(calculate_segments, LEN = L, seqname = N, SIMPLIFY = FALSE)
  gr <- makeGRangesFromDataFrame(do.call(rbind, dfs))
  s_parts <- subseqs <- getSeq(s, gr)
  names(s_parts) <- paste(seqnames(gr), "_sliding:", start(gr), "-", end(gr), sep = "")
  return(s_parts)
}


resolve_name <- function(x){
  if (length(x)==1){
    # no conflict
    return(x)
  } else{
    y <- sapply(x, strsplit, split="|", fixed = TRUE)
    ny <- table(unlist(sapply(y, function(x)paste(seq_along(x), x))))
    if (max(ny)<length(x)){
      return("Unknown")
    }else{
      k <- which(ny==length(x))
      r <- max(as.numeric((gsub(" .+", "", names(k)))))
      out <- paste(y[[1]][1:r], collapse="|")
      return(out)
    }
  }
}

initial_options <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("--file=", "",
  initial_options[grep("--file=", initial_options)])))
source(file.path(script_dir, "library_policy.R"))

opt_parser <- OptionParser(option_list=opt_list)
opt <- parse_args(opt_parser)

# check mandatory arguments
if (is.null(opt$fasta) | is.null(opt$output_dir)){
  message("Missing arguments")
  print_help(opt_parser)
  q(status=0)
}
if (!opt$annotation_conflict %in% c("strict", "nested")) {
  stop("--annotation_conflict must be 'strict' or 'nested', got '",
       opt$annotation_conflict, "'")
}

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(rtracklayer))

dir.create(opt$output_dir, showWarnings = FALSE)
message("Reading fasta file")
s <- readDNAStringSet(opt$fasta)

# Canonical order for reproducible clustering. `mmseqs easy-cluster` is
# order-sensitive: the same set of sequences in a different order elects a
# different set of representatives and a different cluster count. TE_all.fasta
# arrives in DANTE_LTR.gff3 row order, which is not canonical on the multi-chunk
# path -- chunks are concatenated in chunk-index order and the chunk count
# depends on genome size / open-file limit / machine -- so the library, and the
# downstream RepeatMasker annotation, would otherwise vary run-to-run and
# across machines. Sort by sequence content (then name, to break ties) to make
# the library a deterministic function of the input SET, independent of record
# order and robust to upstream coordinate jitter. method = "radix" sorts in the
# C locale, so the order does not depend on the machine's LC_COLLATE.
s <- s[order(as.character(s), names(s), method = "radix")]

size_total <- sum(nchar(s))

message("Partitioning sequences")
s_parts <- slide_and_cut(s)
rm(s)
fasta_parts <- paste(opt$output_dir, "partitioned_s900_w1000.fasta", sep="/")
writeXStringSet(s_parts, fasta_parts)
rm(s_parts)
# run mmseqs2 clustering
# example command:
# mmseqs easy-cluster fasta/TE_all_partitioned_s900_w1000.fasta TE_all_partitioned_s900_w1000_clustered /ssd.scratch/tmp --threads 20 --spaced-kmer-mode 0
#
# --spaced-kmer-mode 0 is load-bearing, not a tuning knob. The canonical sort
# above fixes input-ORDER sensitivity; this fixes run-to-run variation on
# byte-identical input. `easy-cluster` runs cascaded clustering whose first step
# is linclust/kmermatcher, and linclust requests spaced k-mers. With no pattern
# supplied mmseqs generates one AT RANDOM per process, so kmermatcher emits a
# different prefilter every run and the whole cascade inherits it: one fixed
# input here gave 5146 / 5162 / 5167 / 5169 / 5173 clusters over repeated runs.
# It is not thread-related (it varies at --threads 1) and not memory-related (it
# varies with --split-memory-limit pinned); kmermatcher is deterministic as soon
# as this flag is set. Consecutive k-mers cost ~1% more representatives (5221 vs
# 5167, +1.1% bp on that input) and in exchange the library becomes a pure
# function of the input set -- identical across repeated runs and across
# different --threads values.

tempdir <- tempdir()
message("Running mmseqs2 clustering")
cmd <- paste("mmseqs easy-cluster", fasta_parts, paste(opt$output_dir, "mmseqs", sep="/"),
             tempdir, "--threads", opt$threads, "--spaced-kmer-mode 0",
             "-v 1 2>&1" , sep=" ")
out <- system(cmd, intern=TRUE)

cls <- read.table(paste0(opt$output_dir,"/mmseqs_cluster.tsv"), as.is=TRUE, comment.char = "")

# detect conflicting annotations
annot1 <- gsub(".+#","",
               gsub("_sliding.+", "", cls$V1))
annot2 <- gsub(".+#","",
               gsub("_sliding.+", "", cls$V2))

# calculate how many sequences are in each cluster
cls_count <- table(cls$V1)
cls_remove <- names(cls_count)[cls_count < opt$min_coverage]

cls_clusters <- split(cls$V2, cls$V1)
annot_in_clusters <- split(annot2, cls$V1)
size_of_clusters <- sapply(cls_clusters, length)
# remove clusters with less than min_coverage sequences
cls_clusters <- cls_clusters[size_of_clusters >= opt$min_coverage]
annot_in_clusters <- annot_in_clusters[size_of_clusters >= opt$min_coverage]
size_of_clusters <- size_of_clusters[size_of_clusters >= opt$min_coverage]
# Source element of every member, for the element-level counting the nested
# policy uses.  Member names are <element_id>#<classification>_sliding:<s>-<e>,
# and element ids contain neither "#" nor "_sliding" -- assert it rather than
# silently mis-grouping if that ever changes.
element_of <- function(x) sub("#.*", "", x)
elements_in_clusters <- lapply(cls_clusters, element_of)
stopifnot(!any(grepl("_sliding", unlist(elements_in_clusters), fixed = TRUE)))

# One decision per cluster.  lapply, not sapply: sapply simplifies to a matrix
# when every cluster happens to carry the same number of distinct labels, after
# which the per-cluster vectors silently go out of step.
decisions <- lapply(names(annot_in_clusters), function(k)
  classify_cluster(annot_in_clusters[[k]], elements_in_clusters[[k]],
                   proportion_min = opt$proportion_min,
                   policy = opt$annotation_conflict,
                   promote_min_elements = opt$lineage_promotion_min_elements,
                   promote_min_share = opt$lineage_promotion_min_share,
                   resolve_name = resolve_name))
names(decisions) <- names(annot_in_clusters)

kept <- vapply(decisions, function(d) d$keep, logical(1))
final_name <- vapply(decisions[kept], function(d) d$label, character(1))

outcomes <- vapply(decisions, function(d) d$outcome, character(1))
message(sprintf(
  "annotation policy %s: %d clusters | kept %d (majority %d, lca %d, recovered %d, promoted %d) | dropped %d",
  opt$annotation_conflict, length(decisions), sum(kept),
  sum(outcomes == "majority"), sum(outcomes == "lca"),
  sum(outcomes == "recovered"), sum(outcomes == "promoted"),
  sum(outcomes == "dropped")))

final_names_rm_compatible <- gsub("|", "/", gsub("/","_", final_name, fixed=TRUE), fixed=TRUE)
uniq_id <- paste0(gsub("#.+", "", names(final_name)),"_", gsub( ".+sliding:","", names(final_name)))



## first column is representative sequence - how many times there is a conflicts:
rep_conflict <- cls$V1[annot1 != annot2]
rep_count  <- length(unique(cls$V1))
conflict_count <-  length(unique(rep_conflict))

message('Proportion of conflicting annotations : ',
        round(100*conflict_count/rep_count, 2), "%")


rep_seq <- readDNAStringSet(paste0(opt$output_dir,"/mmseqs_rep_seq.fasta"))
names(rep_seq) <- gsub(" .*", "", names(rep_seq))

rep_seq_clean <- rep_seq[match( names(final_name), names(rep_seq))]


names(rep_seq_clean) <- paste(uniq_id, final_name, sep = "#")
rep_seq_clean2 <- rep_seq_clean; names(rep_seq_clean2) <- paste(uniq_id, final_names_rm_compatible, sep = "#")

message("Writing representative sequences")
writeXStringSet(rep_seq_clean, paste0(opt$output_dir,"/mmseqs_representative_seq_clean.fasta"))
writeXStringSet(rep_seq_clean2, paste0(opt$output_dir,"/mmseqs_representative_seq_clean_rm_compatible.fasta"))

size_reduced <- sum(nchar(rep_seq_clean))
message("Input library size          : ", size_total)
message("Representative library size : ", size_reduced)
message("Size reduction              : ", round(100*(1-size_reduced/size_total),2), "%")

