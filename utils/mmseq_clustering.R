#!/usr/bin/env Rscript
# parse arguments
library(optparse)
opt_list <- list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, help="Input fasta file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="output directory", metavar="character"),
  make_option(c("-m", "--min_coverage"), type="numeric", default=3, help="minimal number of sequences in a cluster", metavar="numeric"),
  make_option(c("-t", "--threads"), type="numeric", default=1, help="number of threads", metavar="numeric"),
  make_option(c("-p", "--proportion_min"), type="numeric", default=0.95, help="minimal proportion of the main class in a cluster", metavar="numeric")
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

opt_parser <- OptionParser(option_list=opt_list)
opt <- parse_args(opt_parser)

# check mandatory arguments
if (is.null(opt$fasta) | is.null(opt$output_dir)){
  message("Missing arguments")
  print_help(opt_parser)
  q(status=0)
}

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(rtracklayer))

dir.create(opt$output_dir, showWarnings = FALSE)
message("Reading fasta file")
s <- readDNAStringSet(opt$fasta)
size_total <- sum(nchar(s))

message("Partitioning sequences")
s_parts <- slide_and_cut(s)
rm(s)
fasta_parts <- paste(opt$output_dir, "partitioned_s900_w1000.fasta", sep="/")
writeXStringSet(s_parts, fasta_parts)
rm(s_parts)
# run mmseqs2 clustering
# example command:
# mmseqs easy-cluster fasta/TE_all_partitioned_s900_w1000.fasta TE_all_partitioned_s900_w1000_clustered /ssd.scratch/tmp --threads 20

tempdir <- tempdir()
message("Running mmseqs2 clustering")
cmd <- paste("mmseqs easy-cluster", fasta_parts, paste(opt$output_dir, "mmseqs", sep="/"),
             tempdir, "--threads", opt$threads, "-v 1 2>&1" , sep=" ")
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
main_class_proportion <- sapply(annot_in_clusters, function(x) max(table(x))/length(x))
main_class_name <- sapply(annot_in_clusters, function(x) names(which.max(table(x))))
all_names <- sapply(annot_in_clusters, function(x) unique(x))
consensus_names <- sapply(all_names, resolve_name)


final_name <- ifelse(main_class_proportion > opt$proportion_min,
                     main_class_name, consensus_names)
resoved_names_l <- final_name == main_class_name

final_name <- final_name[resoved_names_l]

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

