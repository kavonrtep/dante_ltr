#!/usr/bin/env Rscript
# select_solo_representatives.R
# Collapse overlapping solo_LTR BLAST hits into one representative per locus.
# Overlap rule: reciprocal overlap >= --overlap_frac of the shorter member.
# Representative preference: SL (with TSD) over SL_noTSD; then longest;
#   ties by Identity, then LibraryID lexicographic.
# Flags (attributes):
#   boundary_uncertain=true  - SL representative strictly contained (>= nest_frac)
#                              inside a longer SL_noTSD member of the same cluster
#   class_conflict=true      - >=1 within-cluster pair has reciprocal overlap
#                              above threshold but different Final_Classification
# Added attributes on every representative: LibraryID, ClusterSize, SupportingHits

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input raw solo GFF3 (with LibraryID on each solo_LTR)"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output representative GFF3"),
  make_option(c("-a", "--annotation_gff3"), type = "character", default = NULL,
              help = "Optional DANTE_LTR GFF3 with complete elements; when supplied, representative-based statistics CSV is written"),
  make_option(c("--stats_out"), type = "character", default = NULL,
              help = "Path for the statistics CSV (default: <output>_statistics.csv)"),
  make_option(c("--overlap_frac"), type = "numeric", default = 0.5,
              help = "Reciprocal overlap threshold of the shorter member [default %default]"),
  make_option(c("--nest_frac"), type = "numeric", default = 0.8,
              help = "SL-inside-SL_noTSD coverage threshold for boundary_uncertain [default %default]")
)

parser <- OptionParser(option_list = option_list,
                       usage = "usage: %prog -i RAW.gff3 -o REP.gff3 [OPTIONS]")
opt <- parse_args(parser)

if (is.null(opt$input) || is.null(opt$output)) {
  print_help(parser)
  stop("Both --input and --output are required.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

# ---- Union-find with path compression, in place via environment scope ----
# Avoids copying `parent` on every edge (copy-on-modify kills performance here).
make_uf <- function(n) {
  env <- new.env(parent = emptyenv())
  env$parent <- seq_len(n)
  env$find <- function(x) {
    p <- env$parent
    while (p[x] != x) {
      p[x] <- p[p[x]]
      x <- p[x]
    }
    env$parent <- p
    x
  }
  env$union <- function(a, b) {
    ra <- env$find(a); rb <- env$find(b)
    if (ra != rb) env$parent[ra] <- rb
  }
  env
}

# ---- Load GFF3 ----
if (!file.exists(opt$input) || file.size(opt$input) == 0L) {
  writeLines("##gff-version 3", opt$output)
  cat("Empty input; wrote empty output.\n")
  quit(save = "no", status = 0L)
}

t0 <- proc.time()[["elapsed"]]
gff <- tryCatch(import(opt$input, format = "gff3"), error = function(e) GRanges())
solo <- gff[!is.na(gff$type) & gff$type == "solo_LTR"]
tsd  <- gff[!is.na(gff$type) & gff$type == "target_site_duplication"]
t_read <- proc.time()[["elapsed"]] - t0

if (length(solo) == 0L) {
  writeLines("##gff-version 3", opt$output)
  cat("No solo_LTR features in input; wrote empty output.\n")
  quit(save = "no", status = 0L)
}

n <- length(solo)
cat(sprintf("Read %d solo_LTR, %d target_site_duplication features\n", n, length(tsd)))

# ---- Build reciprocal-overlap graph (vectorized) ----
t0 <- proc.time()[["elapsed"]]
hits_ol <- findOverlaps(solo, ignore.strand = TRUE,
                        drop.self = TRUE, drop.redundant = TRUE)
t_overlap <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
uf <- make_uf(n)
edges_kept <- 0L

if (length(hits_ol) > 0L) {
  qi <- from(hits_ol); si <- to(hits_ol)
  a_start <- start(solo)[qi]; a_end <- end(solo)[qi]; a_w <- width(solo)[qi]
  b_start <- start(solo)[si]; b_end <- end(solo)[si]; b_w <- width(solo)[si]
  ov_w   <- pmin(a_end, b_end) - pmax(a_start, b_start) + 1L
  ov_w[ov_w < 0L] <- 0L
  frac   <- ov_w / pmin(a_w, b_w)
  keep   <- which(frac >= opt$overlap_frac)
  for (k in keep) uf$union(qi[k], si[k])
  edges_kept <- length(keep)
}

# Component labels
roots      <- vapply(seq_len(n), uf$find, integer(1L))
cluster_id <- as.integer(factor(roots))
n_clusters <- max(cluster_id)
t_graph <- proc.time()[["elapsed"]] - t0

cat(sprintf("Reciprocal-overlap edges kept: %d\n", edges_kept))
cat(sprintf("Clusters: %d (singletons: %d)\n",
            n_clusters, sum(table(cluster_id) == 1L)))

# ---- Pre-compute metadata vectors ----
rank_v   <- as.character(solo$Rank)
class_v  <- as.character(solo$Final_Classification)
id_v     <- as.character(solo$ID)
libid_v  <- as.character(solo$LibraryID)
libid_v[is.na(libid_v)] <- ""
ident_v  <- as.numeric(solo$Identity)
width_v  <- width(solo)

# ---- Pick representative per cluster + flag demotions ----
rep_idx       <- integer(n_clusters)
flag_bound    <- logical(n_clusters)
flag_conflict <- logical(n_clusters)
cluster_sizes <- integer(n_clusters)
supp_lists    <- vector("list", n_clusters)

# Group member indices by cluster once
members_by_cluster <- split(seq_len(n), cluster_id)

# Pre-index edges by cluster for conflict and nesting checks
edge_table <- if (edges_kept > 0L) {
  data.frame(a = qi[keep], b = si[keep], frac = frac[keep])
} else data.frame(a = integer(0), b = integer(0), frac = numeric(0))

# Cluster id for each edge endpoint (same for both, picking a)
edge_cluster <- if (nrow(edge_table) > 0L) cluster_id[edge_table$a] else integer(0)
edges_by_cluster <- split(seq_len(nrow(edge_table)), edge_cluster)

t0 <- proc.time()[["elapsed"]]
n_singleton <- 0L; n_multi <- 0L
for (c in seq_len(n_clusters)) {
  mems <- members_by_cluster[[as.character(c)]]
  cluster_sizes[c] <- length(mems)

  # Fast-path singletons: no rep competition, no flags, no supporting hits.
  if (length(mems) == 1L) {
    rep_idx[c]      <- mems
    supp_lists[[c]] <- character(0L)
    n_singleton     <- n_singleton + 1L
    next
  }
  n_multi <- n_multi + 1L

  pool_sl <- mems[rank_v[mems] == "SL"]
  pool    <- if (length(pool_sl) > 0L) pool_sl else mems
  ord <- order(-width_v[pool], -ident_v[pool], libid_v[pool])
  r <- pool[ord[1L]]
  rep_idx[c] <- r

  # SupportingHits = unique LibraryIDs of the OTHER members
  others <- mems[mems != r]
  sup <- unique(libid_v[others])
  sup <- sup[nzchar(sup) & sup != libid_v[r]]
  supp_lists[[c]] <- sup

  # class_conflict: any intra-cluster edge with differing Final_Classification
  e_idx <- edges_by_cluster[[as.character(c)]]
  if (!is.null(e_idx) && length(e_idx) > 0L) {
    aa <- edge_table$a[e_idx]; bb <- edge_table$b[e_idx]
    if (any(class_v[aa] != class_v[bb])) flag_conflict[c] <- TRUE
  }

  # boundary_uncertain: rep is SL, some longer SL_noTSD in cluster covers
  # >= nest_frac of rep's width
  if (rank_v[r] == "SL") {
    noTSD_longer <- mems[rank_v[mems] == "SL_noTSD" & width_v[mems] > width_v[r]]
    if (length(noTSD_longer) > 0L) {
      r_start <- start(solo)[r]; r_end <- end(solo)[r]; r_w <- width_v[r]
      ov <- pmin(r_end, end(solo)[noTSD_longer]) -
            pmax(r_start, start(solo)[noTSD_longer]) + 1L
      ov[ov < 0L] <- 0L
      if (max(ov) / r_w >= opt$nest_frac) flag_bound[c] <- TRUE
    }
  }
}
t_repsel <- proc.time()[["elapsed"]] - t0

cat(sprintf("Representatives: %d  (singleton clusters: %d, multi-member: %d)\n",
            n_clusters, n_singleton, n_multi))
cat(sprintf("  boundary_uncertain flagged: %d\n", sum(flag_bound)))
cat(sprintf("  class_conflict     flagged: %d\n", sum(flag_conflict)))

# ---- Build representative GRanges ----
rep_gr <- solo[rep_idx]

# Preserve existing attributes; add/overwrite the new ones.
rep_gr$ClusterSize    <- cluster_sizes
sup_strs <- vapply(supp_lists,
                   function(x) if (length(x) == 0L) NA_character_ else paste(x, collapse = ","),
                   character(1L))
rep_gr$SupportingHits <- sup_strs
# Flags: NA when false so rtracklayer omits the attribute
rep_gr$boundary_uncertain <- ifelse(flag_bound,    "true", NA_character_)
rep_gr$class_conflict     <- ifelse(flag_conflict, "true", NA_character_)

# ---- Carry TSD children of representatives ----
kept_ids <- as.character(rep_gr$ID)
if (length(tsd) > 0L) {
  parent_raw <- tsd$Parent
  parent_str <- if (is.list(parent_raw) || methods::is(parent_raw, "CompressedList")) {
    vapply(parent_raw, function(x) if (length(x) == 0L) NA_character_ else as.character(x[[1L]]),
           character(1L))
  } else {
    as.character(parent_raw)
  }
  tsd_keep <- tsd[parent_str %in% kept_ids]
} else {
  tsd_keep <- GRanges()
}

# ---- Assemble + write ----
out <- suppressWarnings(c(rep_gr, tsd_keep))
# Stable sort by seqname then start
o <- order(as.character(seqnames(out)), start(out))
out <- out[o]

t0 <- proc.time()[["elapsed"]]
export(out, con = opt$output, format = "gff3")
t_export <- proc.time()[["elapsed"]] - t0

cat(sprintf("Wrote %d solo_LTR + %d TSD to %s\n",
            length(rep_gr), length(tsd_keep), opt$output))

cat("\n[timing] select_solo_representatives breakdown (wall-clock s):\n")
cat(sprintf("  read GFF3          : %6.2f\n", t_read))
cat(sprintf("  overlap detection  : %6.2f\n", t_overlap))
cat(sprintf("  graph + union-find : %6.2f\n", t_graph))
cat(sprintf("  rep selection loop : %6.2f\n", t_repsel))
cat(sprintf("  export GFF3        : %6.2f\n", t_export))

# ---- Statistics (representative-based) ----
if (!is.null(opt$annotation_gff3)) {
  initial_options <- commandArgs(trailingOnly = FALSE)
  file_arg_name   <- "--file="
  script_name     <- normalizePath(
    sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)])
  )
  source(file.path(dirname(script_name), "solo_ltr_utils.R"))

  ann <- tryCatch(import(opt$annotation_gff3, format = "gff3"),
                  error = function(e) GRanges())
  stats <- get_solo_ltr_statistics(out, ann)

  stats_path <- if (!is.null(opt$stats_out)) opt$stats_out
                else sub("\\.gff3?$", "_statistics.csv", opt$output, ignore.case = TRUE)
  write.table(stats, file = stats_path, sep = "\t",
              quote = FALSE, row.names = FALSE)
  cat(sprintf("Stats (representative-based): %s\n", stats_path))
}
