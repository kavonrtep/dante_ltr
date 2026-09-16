#!/usr/bin/env Rscript
# detect_core_ltr.R -- per-chunk detector for core-domain mode.
#
# Mirrors detect_putative_ltr.R, but seeds on the ordered RT/RH/INT core
# instead of a lineage-keyed full domain complement, and classifies only
# after the element boundaries are fixed.
# See docs/core_domain_mode_design.md.
#
# Everything from the BLAST search onward -- LTR pair choice, TSD, PBS,
# ranking, statistics, FASTA export -- is reused verbatim from
# ltr_utils.R; only candidate generation and classification differ.

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- normalizePath(sub(file_arg_name, "",
                                 initial_options[grep(file_arg_name, initial_options)]))
script_dir <- dirname(script_name)
library(optparse)

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
              help = paste("Maximum number of domains an element may carry that a",
                           "lineage does not have, when testing lineage",
                           "compatibility [default %default]"),
              metavar = "number"),
  make_option(c("-L", "--min_relative_length"), type = "numeric", default = 0.6,
              help = "Minimum relative length of a protein domain [default %default]",
              metavar = "number"),
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Debug mode [default %default]"),
  make_option(c("-t", "--te_constrains"), type = "character", default = NULL,
              help = "core constraints table [default %default]"),
  make_option(c("--min_similarity"), type = "numeric", default = 0.4,
              help = "Minimum domain Similarity [default %default]"),
  make_option(c("--min_relative_length_core"), type = "numeric", default = 0.3,
              help = paste("Relaxed minimum relative length for RT/RH/INT when",
                           "seeding.  The ordered triplet is strong joint",
                           "evidence, so individually weak core domains are",
                           "admitted [default %default]")),
  make_option(c("--core_max_gap"), type = "integer", default = NA_integer_,
              help = "Override core_max_gap from the constraints table"),
  make_option(c("--core_max_span"), type = "integer", default = NA_integer_,
              help = "Override core_max_span from the constraints table"),
  make_option(c("--min_ltr_length"), type = "integer", default = 100,
              help = "Minimum LTR length [default %default]"),
  make_option(c("--max_te_length"), type = "integer", default = 35000,
              help = "Maximum element length [default %default]"),
  make_option(c("--core_require_tsd"), action = "store_true", default = FALSE,
              help = "Reject elements without a TSD [default %default]"),
  make_option(c("--split_max_gap"), type = "integer", default = 50,
              help = paste("Max bp between two annotations of one frameshifted",
                           "domain [default %default]")),
  make_option(c("-n", "--no_ambiguous_domains"), action = "store_true",
              default = FALSE, help = "Remove ambiguous domains from analysis")
)

parser <- OptionParser(option_list = option_list,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
  library(parallel)
})

core_constraints_file <- if (is.null(opt$te_constrains)) {
  paste0(script_dir, "/../databases/core_domain_order.csv")
} else {
  opt$te_constrains
}
lineage_file <- paste0(script_dir, "/../databases/lineage_domain_order.csv")
FDM_file <- paste0(script_dir, "/../databases/feature_distances_model.RDS")
trna_db <- paste0(script_dir, "/../databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_unique-3ends.fasta")
trna_db_hemi <- paste0(script_dir, "/../databases/tRNAscan-SE_ALL_spliced-yes_2022-12-14_plus-old-tRNAs_UC_numbered_unique-half-tRNA-20nt.fasta")

if (!all(file.exists(core_constraints_file, lineage_file, trna_db))) {
  stop("configuration files not found")
}
source(paste0(script_dir, "/../utils/ltr_utils.R"))
source(paste0(script_dir, "/../utils/core_ltr_utils.R"))

constraints <- read_core_constraints(core_constraints_file)
if (!is.na(opt$core_max_gap)) {
  constraints$core_max_gap <- opt$core_max_gap
}
if (!is.na(opt$core_max_span)) {
  constraints$core_max_span <- opt$core_max_span
}
lineage_info <- read.table(lineage_file, sep = "\t", header = TRUE, as.is = TRUE)
lineage_domain <- lineage_domain_map(lineage_info)
FDM <- readRDS(FDM_file)

outfile <- opt$output


#' Write the empty output set and exit cleanly.
exit_empty <- function(msg) {
  cat(msg, "\n", sep = "")
  cat("##gff-version 3\n", file = paste0(outfile, ".gff3"))
  file.create(paste0(outfile, "_statistics.csv"))
  quit(save = "no", status = 0, runLast = FALSE)
}


# MAIN #############################################################

cat("reading gff...")
if (file.size(opt$gff3) == 0) {
  exit_empty("No TE domains on input.")
}
g <- rtracklayer::import(opt$gff3, format = "gff3")
if (length(g) < 3) {
  exit_empty("Less than 3 domains found in input GFF3 file, exiting")
}

if (opt$no_ambiguous_domains) {
  ambiguous_names <- c(
    "Class_I|LTR|Ty1/copia",
    "Class_I|LTR|Ty3/gypsy",
    "Class_I|LTR|Ty3/gypsy|chromovirus",
    "Class_I|LTR|Ty3/gypsy|non-chromovirus",
    "Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA",
    "Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Tat")
  g <- g[!g$Final_Classification %in% ambiguous_names]
}

ori_seqlevels <- seqlevels(g)
decode <- !all(URLencode(seqlevels(g), reserved = TRUE) == seqlevels(g))
g <- CHD_CHDCR_correction(g)
cat("done\n")

g <- gff_cleanup_overlaps(g)
if (length(g) < 3) {
  exit_empty("Less than 3 domains found in input GFF3 file, exiting")
}

cat("reading fasta...")
s <- readDNAStringSet(opt$reference_sequence)
names(s) <- gsub("\\s.+", "", names(s))
cat("done\n")
if (!all(seqlevels(g) %in% names(s))) {
  stop("\nSequence names in input GFF3 do not match sequence names in FASTA file\n\n")
}

# Region_Hits_Classifications carries the per-domain candidate list and
# is needed by core_candidates() and by classification, so unlike
# lineage mode we keep it rather than dropping it here.
region_hits_all <- flatten_region_hits(g)
mcols(g)$Region_Hits_Flat <- region_hits_all
# ...but the CharacterList form cannot be serialised by export(), which
# is why lineage mode drops it at the top of its main block.  Keep only
# the flattened copy.
mcols(g)$Region_Hits_Classifications <- NULL
mcols(g)$Region_Hits_Classifications_ <- NULL

# Fragment merging must precede filtering: each half of a frameshifted
# domain has a reduced Relat_Length by construction (design 6.0.1).
n_before <- length(g)
g <- merge_split_domains(g, max_gap = opt$split_max_gap)
cat("fragment merging: ", n_before, " -> ", length(g), " domains\n", sep = "")

# dante_filtering()'s second arm rescues shorter domains that have a
# same-cluster neighbour.  Lineage mode derives neighbours from its
# lineage-keyed clustering, which core mode does not have; approximate it
# by spatial adjacency on the same strand within the largest core gap.
g <- g[order(as.character(seqnames(g)), start(g))]
adj_gap <- max(constraints$core_max_gap)
same_run <- as.character(seqnames(g))[-1] == head(as.character(seqnames(g)), -1) &
  as.character(strand(g))[-1] == head(as.character(strand(g)), -1) &
  (start(g)[-1] - head(end(g), -1)) <= adj_gap
mcols(g)$neighbors_count <- c(0, as.integer(same_run)) + c(as.integer(same_run), 0)

# Two thresholds, not one (design 6.0.2): seeds may be admitted on weak
# individual evidence because the ordered triplet carries the weight, but
# only solidly-real domains are allowed to truncate a search window.
g_block <- dante_filtering(g, min_similarity = opt$min_similarity,
                           Relative_Length = opt$min_relative_length)
g_seed_pool <- dante_filtering(g, min_similarity = opt$min_similarity,
                               Relative_Length = opt$min_relative_length_core)
g_seed_pool <- g_seed_pool[as.character(mcols(g_seed_pool)$Name) %in% CORE_DOMAINS]

if (length(g_block) == 0) {
  exit_empty("No domains passed filtering, exiting")
}
seqlengths(g_block) <- seqlengths(s)[seqlevels(g_block)]
SL <- seqlengths(s)

candidates <- core_candidates(g_seed_pool)
n_secondary <- sum(mcols(candidates)$Domain_LTR_Support == "secondary")
cat("core candidates (RT/RH/INT)        : ", length(candidates),
    "   (", n_secondary, " via secondary LTR support)\n", sep = "")

seeds <- find_core_seeds(candidates, constraints)
n_gypsy <- sum(seeds$superfamily == "Class_I/LTR/Ty3_gypsy")
n_copia <- sum(seeds$superfamily == "Class_I/LTR/Ty1_copia")
cat("ordered core seeds                 : ", nrow(seeds),
    "   (gypsy ", n_gypsy, " / copia ", n_copia, ")\n", sep = "")

good_TE <- list()
seed_meta <- list()

if (nrow(seeds) > 0) {
  seeds <- seed_search_limits(seeds, g_block, constraints)
  n_trunc <- sum(!is.na(seeds$left_limit) | !is.na(seeds$right_limit))
  cat("windows truncated by blocking rule : ", n_trunc, " (of ", nrow(seeds),
      " seeds)\n", sep = "")

  grL <- core_ranges_left(seeds, constraints)
  grR <- core_ranges_right(seeds, constraints, SL)
  gr <- GRanges(seqnames = seeds$seqnames,
                ranges = IRanges(start = seeds$start, end = seeds$end),
                strand = seeds$strand)

  # drop any seed whose windows fell outside the sequence
  usable <- start(grL) >= 1 & end(grR) <= SL[seeds$seqnames] &
    width(grL) > 0 & width(grR) > 0
  seeds <- seeds[usable, , drop = FALSE]
  grL <- grL[usable]; grR <- grR[usable]; gr <- gr[usable]

  if (nrow(seeds) > 0) {
    s_left <- getSeq(s, grL)
    s_right <- getSeq(s, grR)
    names(s_left) <- paste(seqnames(grL), start(grL), end(grL), sep = "_")
    names(s_right) <- paste(seqnames(grR), start(grR), end(grR), sep = "_")

    expected_ltr <- pmax(
      opt$min_ltr_length,
      constraints$ltr_length[match(seeds$superfamily, constraints$Superfamily)])

    # the seed's own core domains -- get_TE only carries them through to
    # its result; the element's full domain set is collected afterwards,
    # once the boundaries are known (design 7)
    seed_domains <- lapply(seq_len(nrow(seeds)), function(i) {
      as.data.frame(candidates[c(seeds$i1[i], seeds$i2[i], seeds$i3[i])])
    })

    cat("Identification of LTRs...")
    TE <- mclapply(seq_len(nrow(seeds)), function(x) {
      get_TE(s_left[x], s_right[x], seed_domains[[x]], gr[x], grL[x], grR[x],
             expected_ltr[x])
    }, mc.set.seed = TRUE, mc.cores = opt$cpu, mc.preschedule = FALSE)
    cat("done.\n")

    keep <- !sapply(TE, is.null)
    good_TE <- TE[keep]
    seed_meta <- split(seeds[keep, , drop = FALSE], seq_len(sum(keep)))
    cat("seeds with LTR pair                : ", length(good_TE), "\n", sep = "")
  }
}


# --- structural gates (design 6.4) -----------------------------------

gate_counts <- c(G2 = 0L, G3 = 0L, G4 = 0L, G5 = 0L, G6 = 0L)
if (length(good_TE) > 0) {
  block_strand <- as.character(strand(g_block))
  block_sf_ok <- is_ltr_classification(mcols(g_block)$Final_Classification) |
    (!is.na(mcols(g_block)$Region_Hits_Flat) &
       grepl("\\|Class_I\\|LTR", mcols(g_block)$Region_Hits_Flat))
  block_index <- build_domain_index(g_block)

  passes <- logical(length(good_TE))
  for (i in seq_along(good_TE)) {
    info <- good_TE[[i]]$ltr_info[[1]]
    meta <- seed_meta[[i]]
    te_start <- start(info$LTR_L_position)
    te_end <- end(info$LTR_R_position)
    te_len <- te_end - te_start + 1
    core_span <- meta$end - meta$start + 1

    if (te_len < core_span + 2 * opt$min_ltr_length) {
      gate_counts["G2"] <- gate_counts["G2"] + 1L
      next
    }
    if (te_len > opt$max_te_length) {
      gate_counts["G3"] <- gate_counts["G3"] + 1L
      next
    }
    inside <- domains_within(block_index, meta$seqnames, te_start, te_end)
    # G4: an element that swallows an opposite-strand or non-LTR domain
    # is not one element
    bad <- inside[block_strand[inside] != meta$strand | !block_sf_ok[inside]]
    if (length(bad) > 0) {
      gate_counts["G4"] <- gate_counts["G4"] + 1L
      next
    }
    # G5: two seeds inside one LTR pair means the outer repeat is not
    # this element's LTR (design 6.4, O4 resolved: reject)
    n_seeds_inside <- sum(seeds$seqnames == meta$seqnames &
                            seeds$start >= te_start & seeds$end <= te_end)
    if (n_seeds_inside > 1) {
      gate_counts["G5"] <- gate_counts["G5"] + 1L
      next
    }
    if (opt$core_require_tsd && info$TSD_Length <= 3) {
      gate_counts["G6"] <- gate_counts["G6"] + 1L
      next
    }
    passes[i] <- TRUE
  }
  good_TE <- good_TE[passes]
  seed_meta <- seed_meta[passes]
  cat("rejected by structural gates       : ", sum(!passes), "   (",
      paste(sprintf("%s %d", names(gate_counts), gate_counts), collapse = ", "),
      ")\n", sep = "")
}


# --- element assembly and classification -----------------------------

gff3_out <- NULL
n_lineage_call <- 0L
n_demoted <- 0L
n_conflict <- 0L

if (length(good_TE) > 0) {
  # Re-collect the element's full domain set now that the boundaries are
  # known, and hand *that* to get_te_gff3() (design 7).
  block_index <- build_domain_index(g_block)
  for (i in seq_along(good_TE)) {
    info <- good_TE[[i]]$ltr_info[[1]]
    meta <- seed_meta[[i]]
    te_start <- start(info$LTR_L_position)
    te_end <- end(info$LTR_R_position)
    inside <- domains_within(block_index, meta$seqnames, te_start, te_end,
                             strand = meta$strand)
    dom <- g_block[inside]
    dom <- dom[order(start(dom))]
    if (meta$strand == "-") {
      dom <- rev(dom)
    }
    good_TE[[i]]$domain <- as.data.frame(dom)
    good_TE[[i]]$core_meta <- meta
  }

  ID <- paste0("TE_", sprintf("%08d", seq_along(good_TE)))
  gff3_list <- mcmapply(get_te_gff3, g = good_TE, ID = ID, mc.cores = opt$cpu)

  # get_te_gff3() takes Name/Final_Classification from the first domain;
  # core mode overrides both with the computed classification.
  for (i in seq_along(gff3_list)) {
    x <- gff3_list[[i]]
    meta <- good_TE[[i]]$core_meta
    dom <- good_TE[[i]]$domain
    cls <- classify_core_element(
      domain_names = as.character(dom$Name),
      domain_classifications = as.character(dom$Final_Classification),
      region_hits = as.character(dom$Region_Hits_Flat),
      superfamily = meta$superfamily,
      lineage_domain = lineage_domain,
      max_missing_domains = opt$max_missing_domains)

    is_te <- x$type == "transposable_element"
    is_ltr <- x$type == "long_terminal_repeat"
    x$Final_Classification[is_te | is_ltr] <- cls$Final_Classification
    x$Name[is_te] <- cls$Final_Classification

    n <- length(x)
    blank <- rep(NA_character_, n)
    put <- function(v, value) {
      v[is_te] <- value
      v
    }
    x$Superfamily_Evidence <- put(blank, cls$Superfamily_Evidence)
    x$Core_Domains <- put(blank, cls$Core_Domains)
    x$Accessory_Domains <- put(blank, paste(
      c(meta$accessory_5, meta$accessory_3)[nzchar(c(meta$accessory_5,
                                                     meta$accessory_3))],
      collapse = " "))
    x$Accessory_Superfamily_Mismatch <- put(
      blank, as.character(meta$accessory_sf_mismatch))
    x$Lineage_Call <- put(blank, cls$Lineage_Call)
    x$Lineage_Candidates <- put(blank, cls$Lineage_Candidates)
    x$Lineage_Support <- put(blank, cls$Lineage_Support)
    x$Classification_Demoted <- put(
      blank, tolower(as.character(cls$Classification_Demoted)))
    x$Classification_Conflict <- put(
      blank, tolower(as.character(cls$Classification_Conflict)))

    if (!is.na(cls$Lineage_Call)) n_lineage_call <- n_lineage_call + 1L
    if (isTRUE(cls$Classification_Demoted)) n_demoted <- n_demoted + 1L
    if (isTRUE(cls$Classification_Conflict)) n_conflict <- n_conflict + 1L
    gff3_list[[i]] <- x
  }

  cat("Identification of PBS ...")
  gff3_list2 <- mclapply(gff3_list, FUN = add_pbs, s = s, trna_db = trna_db,
                         mc.set.seed = TRUE, mc.cores = opt$cpu,
                         mc.preschedule = FALSE)
  cat("done\n")
  pbs_pos <- gff3_list2[sapply(gff3_list2,
                               function(x) "primer_binding_site" %in% x$type)]
  pbs_neg <- gff3_list[!sapply(gff3_list2,
                               function(x) "primer_binding_site" %in% x$type)]
  cat("Identification of PBS - half-molecule tRNA ...")
  gff3_list3 <- mclapply(pbs_neg, FUN = add_pbs_hemi, s = s,
                         trna_db = trna_db_hemi, mc.set.seed = TRUE,
                         mc.cores = opt$cpu, mc.preschedule = FALSE)
  cat(" done\n")

  gff3_out <- do.call(c, append(gff3_list3, pbs_pos))
  src <- as.character(gff3_out$source)
  src[is.na(src) | src == "dante_ltr"] <- "dante_ltr_core"
  gff3_out$source <- src
  gff3_out$Rank <- get_te_rank(gff3_out)
}


# --- rank D track (design 6.5) ---------------------------------------

# Everything that passed the filter and is not inside a called element is
# reported at rank D, reusing the lineage-mode machinery unchanged.
g_block$domain_order <- 0
lineage_domains_sequence <- unlist(mapply(function(d, l) {
  paste(strsplit(d, " ")[[1]], ":", l, sep = "")
}, d = lineage_domain, l = names(lineage_domain)))
dom_ord <- as.numeric(factor(
  paste(g_block$Name, g_block$Final_Classification, sep = ":"),
  levels = lineage_domains_sequence))
dom_ord[is.na(dom_ord)] <- 0
g_block$domain_order <- dom_ord

cls_alt <- get_domain_clusters_alt(g_block, FDM)
g_block$Cluster <- as.numeric(factor(cls_alt))
gcl_alt <- split(as.data.frame(g_block), cls_alt)
TE_partial <- GRanges(
  seqnames = sapply(gcl_alt, function(x) x$seqnames[1]),
  Name = sapply(gcl_alt, function(x) x$Final_Classification[1]),
  Final_Classification = sapply(gcl_alt, function(x) x$Final_Classification[1]),
  ID = sapply(gcl_alt, function(x) paste0("TE_partial_",
                                          sprintf("%08d", x$Cluster[1]))),
  strand = sapply(gcl_alt, function(x) x$strand[1]),
  Ndomains = sapply(gcl_alt, function(x) nrow(x)),
  type = "transposable_element",
  source = "dante_ltr_core",
  Rank = "D",
  IRanges(start = sapply(gcl_alt, function(x) min(x$start)),
          end = sapply(gcl_alt, function(x) max(x$end))))
g_block$Ndomains_in_cluster <- count_occurences_for_each_element(g_block$Cluster)
g_block$Parent <- paste0("TE_partial_", sprintf("%08d", g_block$Cluster))
g_block$Rank <- "D"
RT <- g_block[g_block$Name == "RT" &
                substring(g_block$Final_Classification, 1, 11) == "Class_I|LTR"]
TE_partial_multi <- TE_partial[TE_partial$Ndomains > 1]

# export() fails on a zero-length GRanges ("differing number of rows:
# 0, 1"), and trim_gr() exports its input, so guard the empty case.
# Lineage mode has the same latent issue but rarely reaches it.
if (!is.null(gff3_out) && length(TE_partial_multi) > 0) {
  TE_partial_parent_part <- trim_gr(TE_partial_multi, gff3_out)
  if (!is.null(TE_partial_parent_part)) {
    TE_partial_domain_part <- g_block[g_block$Parent %in% TE_partial_parent_part$ID]
    TE_partial_domain_part <- trim_gr(TE_partial_domain_part, gff3_out)
    gff3_out <- sort(merge_gr(gff3_out, TE_partial_parent_part,
                              TE_partial_domain_part), by = ~ seqnames * start)
  } else {
    gff3_out <- sort(gff3_out, by = ~ seqnames * start)
  }
  gff3_out$Parent <- as.character(gff3_out$Parent)
} else if (!is.null(gff3_out)) {
  gff3_out <- sort(gff3_out, by = ~ seqnames * start)
  gff3_out$Parent <- as.character(gff3_out$Parent)
} else if (length(TE_partial_multi) > 0) {
  TE_partial_domain_part <- g_block[g_block$Parent %in% TE_partial_multi$ID]
  gff3_out <- sort(c(TE_partial_domain_part, TE_partial_multi),
                   by = ~ seqnames * start)
}

cat("elements with lineage-level call   : ", n_lineage_call, "\n", sep = "")
cat("elements demoted to superfamily    : ", n_demoted, "\n", sep = "")
cat("classification conflicts           : ", n_conflict, "\n", sep = "")


# --- export ----------------------------------------------------------

if (is.null(gff3_out)) {
  exit_empty("No TEs found.")
}

gff3_out$Region_Hits_Flat <- NULL
gff3_out$ID[!is.na(gff3_out$ID)] <-
  paste0(gff3_out$ID[!is.na(gff3_out$ID)], "_",
         seqnames(gff3_out)[!is.na(gff3_out$ID)])
gff3_out$Parent[!is.na(gff3_out$Parent)] <-
  paste0(gff3_out$Parent[!is.na(gff3_out$Parent)], "_",
         seqnames(gff3_out)[!is.na(gff3_out$Parent)])
gff3_out <- convert_gr_Lists_to_Vectors(gff3_out)
gff3_out <- revert_CHDCR_correction(gff3_out)
gff3_out <- add_info_about_flanking_sequences(gff3_out, s, 10)

seqlevels_changed <- !all(seqlevels(gff3_out) %in% ori_seqlevels)
if (decode & seqlevels_changed) {
  export_gff3_without_url_encoding(gff3_out, con = paste0(outfile, ".gff3"))
} else {
  export(gff3_out, con = paste0(outfile, ".gff3"), format = "gff3")
}

all_tbl <- get_core_te_statistics(gff3_out, RT)
all_tbl <- cbind(Classification = rownames(all_tbl), all_tbl)
write.table(all_tbl, file = paste0(outfile, "_statistics.csv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

s_te <- get_te_sequences(gff3_out, s)
for (i in seq_along(s_te)) {
  writeXStringSet(s_te[[i]],
                  filepath = paste0(outfile, "_", names(s_te)[i], ".fasta"))
}
