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


# --- search-space delimitation ---------------------------------------

section("blocking walk (design 6.3)")

# A gypsy core with room on both sides; accessory domains are added
# relative to it.  Element orientation is plus, so 5' is to the left.
walk_case <- function(extra = NULL) {
  spec <- gypsy_spec(base = 40000)
  if (!is.null(extra)) {
    spec <- rbind(spec[, c("start", "end", "name")],
                  extra[, c("start", "end", "name")])
    for (col in c("strand", "classification", "db", "similarity",
                  "relat_length")) {
      if (col %in% names(extra)) {
        spec[[col]] <- c(rep(NA, 3), as.character(extra[[col]]))
      }
    }
    for (col in intersect(names(spec), c("strand", "classification"))) {
      d <- if (col == "strand") "+" else "Class_I|LTR|Ty3/gypsy"
      spec[[col]][is.na(spec[[col]])] <- d
    }
  }
  g <- gr_of(spec)
  s <- find_core_seeds(core_candidates(g), CONSTRAINTS)
  seed_search_limits(s, g, CONSTRAINTS)
}

acc <- function(start, end, name, ...) {
  data.frame(start = start, end = end, name = name, ...,
             stringsAsFactors = FALSE)
}

# PROT is accessory_5 for gypsy: the walk passes through it
s <- walk_case(acc(39000, 39500, "PROT"))
ok("PROT 5' of a gypsy core -> traversed",
   nrow(s) == 1L && s$accessory_5[1] == "PROT" && is.na(s$left_limit[1]))

# PROT then GAG, inward->outward: both traversed
s <- walk_case(rbind(acc(39000, 39500, "PROT"), acc(37000, 37500, "GAG")))
ok("PROT then GAG going outward -> both traversed",
   nrow(s) == 1L && s$accessory_5[1] == "PROT GAG")

# GAG then PROT: out of canonical order, so the walk stops at PROT
s <- walk_case(rbind(acc(39000, 39500, "GAG"), acc(37000, 37500, "PROT")))
ok("GAG before PROT -> walk stops (rule 5)",
   nrow(s) == 1L && s$accessory_5[1] == "GAG" && s$left_limit[1] == 37500)

# a second GAG going outward belongs to the neighbouring element
s <- walk_case(rbind(acc(39000, 39500, "GAG"), acc(37000, 37500, "GAG")))
ok("second GAG going outward -> walk stops (rule 4)",
   nrow(s) == 1L && s$accessory_5[1] == "GAG" && s$left_limit[1] == 37500)

# GAG is accessory_5, so on the 3' side it is out of place
s <- walk_case(acc(43000, 43500, "GAG"))
ok("GAG 3' of a gypsy core -> blocks (rule 3)",
   nrow(s) == 1L && s$accessory_3[1] == "" && s$right_limit[1] == 43000)

# CHD is accessory_3 for gypsy
s <- walk_case(acc(43000, 43500, "CHD"))
ok("CHD 3' of a gypsy core -> traversed",
   nrow(s) == 1L && s$accessory_3[1] == "CHD" && is.na(s$right_limit[1]))

s <- walk_case(acc(39000, 39500, "CHD"))
ok("CHD 5' of a gypsy core -> blocks (rule 3)",
   nrow(s) == 1L && s$left_limit[1] == 39500)

s <- walk_case(acc(39000, 39500, "PROT", strand = "-"))
ok("opposite-strand domain -> blocks (rule 1)",
   nrow(s) == 1L && s$left_limit[1] == 39500)

s <- walk_case(acc(39000, 39500, "TPase",
                   classification = "Class_II|Subclass_1|TIR"))
ok("Class II transposase -> blocks (rule 2)",
   nrow(s) == 1L && s$left_limit[1] == 39500)

s <- walk_case(acc(39000, 39500, "RT"))
ok("another core domain -> blocks (rule 3)",
   nrow(s) == 1L && s$left_limit[1] == 39500)

# rule 6: a GAG confidently called copia next to a gypsy core is an
# element boundary
s <- walk_case(acc(39000, 39500, "PROT",
                   classification = "Class_I|LTR|Ty1/copia|Ivana"))
ok("copia-classified PROT 5' of a gypsy core -> blocks (rule 6)",
   nrow(s) == 1L && s$left_limit[1] == 39500 &&
     s$accessory_sf_mismatch[1] == 1L)

# ...but rule 6 must not fire on missing information
s <- walk_case(acc(39000, 39500, "PROT", classification = "Class_I|LTR"))
ok("PROT resolving only to Class_I|LTR -> traversed (rule 6 vacuous)",
   nrow(s) == 1L && s$accessory_5[1] == "PROT" &&
     s$accessory_sf_mismatch[1] == 0L)

# no accessory domains at all: the window runs to the offset cap
s <- walk_case()
ok("no accessory domains -> no limit on either side",
   nrow(s) == 1L && is.na(s$left_limit[1]) && is.na(s$right_limit[1]))

SL <- c(chr1 = 200000L)
grL <- core_ranges_left(s, CONSTRAINTS)
grR <- core_ranges_right(s, CONSTRAINTS, SL)
ok("unblocked window uses the table offset",
   start(grL)[1] == s$start[1] - 17000 && end(grR)[1] == s$end[1] + 14000)
ok("window overlaps the core by offset2",
   end(grL)[1] == s$start[1] + 300 && start(grR)[1] == s$end[1] - 300)

s_blocked <- walk_case(acc(39000, 39500, "RT"))
grL <- core_ranges_left(s_blocked, CONSTRAINTS)
ok("blocked window extends 100 bp into the blocker",
   start(grL)[1] == 39500 - 100)

# minus strand: the element's 5' end is at higher coordinates, so the
# offsets swap
minus_seed <- s
minus_seed$strand <- "-"
grL <- core_ranges_left(minus_seed, CONSTRAINTS)
grR <- core_ranges_right(minus_seed, CONSTRAINTS, SL)
ok("minus strand swaps the 5'/3' offsets",
   start(grL)[1] == minus_seed$start[1] - 14000 &&
     end(grR)[1] == minus_seed$end[1] + 17000)

near_start <- s
near_start$start <- 500L; near_start$end <- 2000L
ok("window clamped at the sequence start",
   start(core_ranges_left(near_start, CONSTRAINTS))[1] == 1L)
near_end <- s
near_end$start <- 190000L; near_end$end <- 199000L
ok("window clamped at the sequence end",
   end(core_ranges_right(near_end, CONSTRAINTS, SL))[1] == 200000L)


# --- summary ---------------------------------------------------------

cat("\n", CHECKS - FAILURES, "/", CHECKS, " checks passed\n", sep = "")
if (FAILURES > 0L) {
  cat("FAILED:", FAILURES, "\n")
  quit(save = "no", status = 1)
}
quit(save = "no", status = 0)
