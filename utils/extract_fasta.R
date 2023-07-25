#!/usr/bin/env Rscript
# parse command line arguments
library(optparse)
opt_list <- list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, help="File with fasta sequences annotated by DANTE_LTR", metavar="character"),
  make_option(c("-g", "--gff3"), type="character", default=NULL, help="GFF3 with DANTE_LTR annotation", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=opt_list)
opt <- parse_args(opt_parser)

# all arguments are mandatory
if (is.null(opt$fasta) | is.null(opt$gff3) | is.null(opt$output)){
  message("Missing arguments")
  print_help(opt_parser)
  q(status=0)
}

suppressPackageStartupMessages({
    library(Biostrings)
    library(rtracklayer)
    library(BSgenome)
})


getSeqNamed <- function(s, gr, name = NULL) {
  spart <- getSeq(s, gr)
  if (is.null(name)){
    id1 <- paste0(seqnames(gr), '_', start(gr), "_", end(gr))
  }else{
    id1 <- mcols(gr)[,name]
  }
  id2 <- gr$Final_Classification
  names(spart) <- paste0(id1, "#", id2)
  spart
}

s <- readDNAStringSet(opt$fasta)
# rename - keep only part of name before whitespace - space or tab
names(s) <- gsub("([^\t ]+).*", "\\1", names(s))



g <- import(opt$gff3)

gparts <- list(
  TE_all= g[g$type == "transposable_element"],
  TE_DLplus = g[g$type == "transposable_element" & g$Rank!="D"],
  TE_DLTP = g[g$type == "transposable_element" & g$Rank == "DLTP"],
  TE_DL = g[g$type == "transposable_element" & g$Rank == "DL"],
  TE_DLP = g[g$type == "transposable_element" & g$Rank == "DLP"],
  TE_DLT = g[g$type == "transposable_element" & g$Rank == "DLT"]
)

# remove empty elements
gparts <- gparts[sapply(gparts, length) > 0]
s_te <- lapply(gparts, function(x) sp = getSeqNamed(s, x))
dir.create(opt$output, showWarnings = FALSE)
x <- lapply(names(s_te), function(x)writeXStringSet(s_te[[x]], filepath = paste0(opt$output,"/",x,".fasta")))
