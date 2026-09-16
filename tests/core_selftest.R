#!/usr/bin/env Rscript
# tests/core_selftest.R
#
# Unit tests for utils/core_ltr_utils.R.  Builds synthetic GRanges
# directly, so it needs no genome and no BLAST and runs in well under a
# second -- which is what makes the core-mode R code iterable.
#
# Run directly, or via tests/core.sh.

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", initial_options[grep(file_arg, initial_options)])
script_dir <- dirname(normalizePath(script_path))
root <- dirname(script_dir)
source(file.path(root, "utils", "core_ltr_utils.R"))

CONSTRAINTS <- read_core_constraints(
  file.path(root, "databases", "core_domain_order.csv"))

FAILURES <- 0L
CHECKS <- 0L

ok <- function(label, condition, detail = "") {
  CHECKS <<- CHECKS + 1L
  if (isTRUE(condition)) {
    cat("  ok   ", label, "\n", sep = "")
  } else {
    FAILURES <<- FAILURES + 1L
    cat("  FAIL ", label, if (nzchar(detail)) paste0("  [", detail, "]") else "",
        "\n", sep = "")
  }
}

section <- function(title) cat("\n== ", title, " ==\n", sep = "")

#' Build a synthetic DANTE-like GRanges.
#' `spec` is a data.frame with start, end, name and optional strand,
#' classification, db, similarity, relat_length.
gr_of <- function(spec, seqname = "chr1") {
  n <- nrow(spec)
  fill <- function(col, default) {
    if (col %in% names(spec)) as.character(spec[[col]]) else rep(default, n)
  }
  g <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = spec$start, end = spec$end),
    strand = fill("strand", "+"))
  mcols(g)$Name <- as.character(spec$name)
  mcols(g)$Final_Classification <- fill("classification",
                                        "Class_I|LTR|Ty3/gypsy")
  mcols(g)$Best_Hit_DB_Pos <- fill("db", "1:100of100")
  mcols(g)$Similarity <- fill("similarity", "0.8")
  mcols(g)$Relat_Length <- fill("relat_length", "0.9")
  mcols(g)$Region_Hits_Classifications <- CharacterList(
    as.list(fill("region_hits", "RT|Class_I|LTR|Ty3/gypsy[100bp]")))
  g
}

seeds_of <- function(g) find_core_seeds(core_candidates(g), CONSTRAINTS)

# gypsy core on the plus strand at a fixed offset
gypsy_spec <- function(base = 10000, gap = 200) {
  data.frame(
    start = c(base, base + 500 + gap, base + 1000 + 2 * gap),
    end   = c(base + 499, base + 999 + gap, base + 1499 + 2 * gap),
    name  = c("RT", "RH", "INT"), stringsAsFactors = FALSE)
}


# --- seeding ---------------------------------------------------------

section("core seeding (design 6.2)")

s <- seeds_of(gr_of(gypsy_spec()))
ok("single gypsy core -> one seed", nrow(s) == 1L)
ok("superfamily from order alone",
   identical(s$superfamily[1], "Class_I/LTR/Ty3_gypsy"))

copia <- gypsy_spec()
copia$name <- c("INT", "RT", "RH")
s <- seeds_of(gr_of(copia))
ok("copia order -> copia superfamily",
   nrow(s) == 1L && identical(s$superfamily[1], "Class_I/LTR/Ty1_copia"))

# tandem array: two complete cores in one run
tandem <- rbind(gypsy_spec(base = 10000), gypsy_spec(base = 13000))
s <- seeds_of(gr_of(tandem))
ok("tandem array RT RH INT RT RH INT -> two seeds", nrow(s) == 2L,
   paste("got", nrow(s)))

# minus strand: same coordinates, reversed element orientation
minus <- gypsy_spec()
minus$name <- c("INT", "RH", "RT")   # reads RT,RH,INT right-to-left
minus$strand <- "-"
s <- seeds_of(gr_of(minus))
ok("minus strand gypsy -> one gypsy seed",
   nrow(s) == 1L && identical(s$superfamily[1], "Class_I/LTR/Ty3_gypsy"))

minus_bad <- gypsy_spec()
minus_bad$name <- c("RT", "RH", "INT")  # reads INT,RH,RT in element order
minus_bad$strand <- "-"
ok("minus strand with plus-strand order -> no seed",
   nrow(seeds_of(gr_of(minus_bad))) == 0L)

two <- gypsy_spec()[1:2, ]
ok("RT + RH only -> no seed", nrow(seeds_of(gr_of(two))) == 0L)

wrong <- gypsy_spec()
wrong$name <- c("RH", "RT", "INT")
ok("wrong core order -> no seed", nrow(seeds_of(gr_of(wrong))) == 0L)

far <- gypsy_spec(gap = 5000)   # exceeds core_max_gap of 2000
ok("gap beyond core_max_gap -> no seed", nrow(seeds_of(gr_of(far))) == 0L)

non_ltr <- gypsy_spec()
s <- seeds_of(gr_of(cbind(non_ltr,
                          classification = "Class_II|Subclass_1|TIR",
                          region_hits = "RT|Class_II|Subclass_1|TIR[100bp]")))
ok("non-LTR classification -> no seed", nrow(s) == 0L)

secondary <- gypsy_spec()
s <- seeds_of(gr_of(cbind(
  secondary,
  classification = "Class_I|pararetrovirus",
  region_hits = "RT|Class_I|LTR|Ty3/gypsy|chromovirus[100bp]")))
ok("non-LTR best hit but LTR in Region_Hits -> seeds (design 6.1)",
   nrow(s) == 1L)

# a duplicated RT that is NOT a split: the seed forms, the spare remains
dup <- gypsy_spec()
dup <- rbind(data.frame(start = 8000, end = 8499, name = "RT",
                        stringsAsFactors = FALSE), dup)
s <- seeds_of(gr_of(dup))
ok("duplicated RT -> exactly one seed", nrow(s) == 1L, paste("got", nrow(s)))

# Determinism: the same input in a different row order must give the
# same seeds.  i1/i2/i3 index the candidate object, so they legitimately
# differ; compare the positional result.
g <- gr_of(rbind(gypsy_spec(base = 10000), gypsy_spec(base = 13000)))
loc <- function(x) {
  x <- x[, c("seqnames", "strand", "superfamily", "start", "end")]
  rownames(x) <- NULL
  x
}
ok("seeding is order-independent",
   identical(loc(seeds_of(g)),
             loc(seeds_of(g[order(end(g), decreasing = TRUE)]))))
ok("seeding is repeatable", identical(seeds_of(g), seeds_of(g)))


# --- fragment merging ------------------------------------------------

section("fragment merging (design 6.0.1)")

split_spec <- data.frame(
  start = c(10000, 10260), end = c(10250, 10600),
  name = c("RT", "RT"),
  db = c("1:60of120", "63:120of120"),
  relat_length = c("0.35", "0.30"),
  stringsAsFactors = FALSE)
m <- merge_split_domains(gr_of(split_spec))
ok("complementary tiling within max_gap -> merged", length(m) == 1L,
   paste("got", length(m)))
ok("merged feature spans both fragments",
   length(m) == 1L && start(m)[1] == 10000 && end(m)[1] == 10600)
ok("merged Relat_Length is the sum",
   length(m) == 1L && abs(as.numeric(m$Relat_Length[1]) - 0.65) < 1e-9)
ok("merged records Nfragments", length(m) == 1L && m$Nfragments[1] == 2L)

overlap_spec <- split_spec
overlap_spec$db <- c("1:120of120", "1:120of120")
ok("two genuine copies (overlapping DB hits) -> not merged",
   length(merge_split_domains(gr_of(overlap_spec))) == 2L)

far_spec <- split_spec
far_spec$start <- c(10000, 12000); far_spec$end <- c(10250, 12340)
ok("split beyond max_gap -> not merged",
   length(merge_split_domains(gr_of(far_spec))) == 2L)

diff_name <- split_spec
diff_name$name <- c("RT", "RH")
ok("different domain names -> not merged",
   length(merge_split_domains(gr_of(diff_name))) == 2L)

diff_strand <- split_spec
diff_strand$strand <- c("+", "-")
ok("different strands -> not merged",
   length(merge_split_domains(gr_of(diff_strand))) == 2L)

minus_split <- split_spec
minus_split$strand <- c("-", "-")
minus_split$db <- c("63:120of120", "1:60of120")   # reference order flips
ok("minus-strand split -> merged",
   length(merge_split_domains(gr_of(minus_split))) == 1L)

minus_wrong <- split_spec
minus_wrong$strand <- c("-", "-")   # plus-strand DB order on minus strand
ok("minus-strand split in plus-strand DB order -> not merged",
   length(merge_split_domains(gr_of(minus_wrong))) == 2L)

# merge-then-filter: two sub-threshold halves become one passing domain.
# Filtering first would delete both -- this ordering is the whole point.
source(file.path(root, "utils", "ltr_utils.R"))
frag <- gr_of(split_spec)
mcols(frag)$Identity <- "0.5"
mcols(frag)$Relat_Interruptions <- "0"
mcols(frag)$neighbors_count <- 0
filtered_first <- dante_filtering(frag, Relative_Length = 0.6)
merged_first <- dante_filtering(merge_split_domains(frag),
                                Relative_Length = 0.6)
ok("filter-then-merge loses both fragments", length(filtered_first) == 0L)
ok("merge-then-filter rescues the domain (design 5.3)",
   length(merged_first) == 1L)


# --- summary ---------------------------------------------------------

cat("\n", CHECKS - FAILURES, "/", CHECKS, " checks passed\n", sep = "")
if (FAILURES > 0L) {
  cat("FAILED:", FAILURES, "\n")
  quit(save = "no", status = 1)
}
quit(save = "no", status = 0)
