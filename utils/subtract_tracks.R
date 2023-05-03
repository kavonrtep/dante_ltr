#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer);
  library(optparse)
})

# parse command line arguments
option_list <- list(
  make_option(c("-A", "--A_input"), type="character", default=NULL, help="Input  GFF3/BED file A", metavar="character"),
  make_option(c("-B", "--B_input"), type="character", default=NULL, help="Input GFF3/BED file B", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output (BED)", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


A <- import(opt$A_input)

B <- import(opt$B_input)

# merge all overlaping features
A <- reduce(A)
B <- reduce(B)

A_total_size <- sum(width(A))
B_total_size <- sum(width(B))

# use bedtools to subtract B from A
# write to temp files
tmpA <- tempfile(fileext = ".bed")
tmpB <- tempfile(fileext = ".bed")
export(A, con = tmpA, format = "bed")
export(B, con = tmpB, format = "bed")
tmpOut <- tempfile(fileext = ".bed")

suppressMessages({
  system(paste("bedtools subtract -a", tmpA, "-b", tmpB, ">", tmpOut))
})
out <- import(tmpOut, format = "bed")
A_minus_B_size <- sum(width(out))
output <- paste0(opt$output, "_A_minus_B.bed")
export(out, con = output, format = "bed")

suppressMessages({
  system(paste("bedtools subtract -b", tmpA, "-a", tmpB, ">", tmpOut))
})
out <- import(tmpOut, format = "bed")
B_minus_A_size <- sum(width(out))
output <- paste0(opt$output, "_B_minus_A.bed")
export(out, con = output, format = "bed")



# intersect A and B
suppressMessages({
  system(paste("bedtools intersect -b", tmpA, "-a", tmpB, ">", tmpOut))
})
out <- import(tmpOut, format = "bed")
Ovelap_size <- sum(width(out))
output <- paste0(opt$output, "_AB_intersect.bed")
export(out, con = output, format = "bed")

# remove temp files
suppressMessages({
  s <- file.remove(tmpB)
  s <- file.remove(tmpA)
  s <- file.remove(tmpOut)
})

# export statistics
output <- paste0(opt$output, "_stats.txt")
con <- file(output, open = "w")
cat(paste0("A file: ", opt$A_input, "\n"), file = con)
cat(paste0("B file: ", opt$B_input, "\n"), file = con)
cat(paste0("----------------------------------------------------------------------\n"), file = con)
cat(paste0("A_total_size:\t", A_total_size,"\n"), file = con)
cat(paste0("B_total_size:\t", B_total_size,"\n"), file = con)
cat(paste0("A_minus_B_size:\t", A_minus_B_size,"\n"), file = con)
cat(paste0("B_minus_A_size:\t", B_minus_A_size,"\n"), file = con)
cat(paste0("Ovelap_size:\t", Ovelap_size,"\n"), file = con)
close(con)

