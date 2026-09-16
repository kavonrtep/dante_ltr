# core_ltr_utils.R
#
# Support functions for core-domain mode (`dante_ltr --mode core`).
# See docs/core_domain_mode_design.md.
#
# Core mode seeds on the ordered RT/RH/INT core -- whose relative order
# is the structural difference between the two LTR superfamilies -- and
# classifies afterwards, instead of requiring a lineage-keyed full domain
# complement up front.  This file holds the parts with no lineage-mode
# analogue:
#
#   1. fragment merging          (design 6.0.1)
#   2. core candidate selection  (design 6.1)
#   3. ordered core seeding      (design 6.2)
#   4. search-space delimitation (design 6.3)
#   5. reannotation and LCA classification (design 7)
#
# Everything downstream of the search window -- BLAST, LTR pair choice,
# TSD, PBS, ranking, statistics -- is reused verbatim from ltr_utils.R.

CORE_DOMAINS <- c("RT", "RH", "INT")

# Ordered core triplet, read in element orientation, per superfamily.
# This mapping is the premise core mode rests on; it held on 51875 of
# 51875 validated elements across three genomes (design 5.2).
CORE_ORDER <- list(
  "Class_I/LTR/Ty3_gypsy" = c("RT", "RH", "INT"),
  "Class_I/LTR/Ty1_copia" = c("INT", "RT", "RH")
)

# Two annotations overlapping on the genome by more than this are rival
# hits to one locus, not fragments of one domain.
MAX_GENOME_OVERLAP <- 30


# --- helpers ---------------------------------------------------------

#' Numeric metadata column, coerced and NA-safe.
num_mcol <- function(g, name, default = 0) {
  if (!name %in% names(mcols(g))) {
    return(rep(default, length(g)))
  }
  v <- suppressWarnings(as.numeric(as.character(mcols(g)[[name]])))
  v[is.na(v)] <- default
  v
}


#' Parse DANTE's Best_Hit_DB_Pos ("1:120of121") into a 3-column matrix.
#' Rows that do not match give NA.
parse_db_pos <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexec("^([0-9]+):([0-9]+)of([0-9]+)$", x))
  out <- matrix(NA_real_, nrow = length(x), ncol = 3,
                dimnames = list(NULL, c("start", "end", "length")))
  ok <- lengths(m) == 4
  if (any(ok)) {
    parts <- do.call(rbind, m[ok])
    out[ok, ] <- as.numeric(parts[, 2:4])
  }
  out
}


#' Does this classification string denote an LTR retrotransposon?
is_ltr_classification <- function(cls) {
  cls <- as.character(cls)
  !is.na(cls) & grepl("^Class_I[|/]LTR", cls)
}


#' Flatten Region_Hits_Classifications, which rtracklayer splits into a
#' CharacterList because DANTE separates the hits with commas.
flatten_region_hits <- function(g) {
  if (!"Region_Hits_Classifications" %in% names(mcols(g))) {
    return(rep(NA_character_, length(g)))
  }
  v <- mcols(g)$Region_Hits_Classifications
  if (is(v, "List")) {
    return(as.character(unstrsplit(v, sep = ",")))
  }
  as.character(v)
}


#' Superfamily of a REXdb classification, or NA when it is not an LTR-RT.
superfamily_of <- function(cls) {
  cls <- gsub("\\|", "/", as.character(cls))
  out <- rep(NA_character_, length(cls))
  ltr <- !is.na(cls) & grepl("^Class_I/LTR/", cls)
  third <- sapply(strsplit(cls, "/"), function(p) if (length(p) >= 3) p[3] else NA)
  out[ltr & grepl("^Ty1", third)] <- "Class_I/LTR/Ty1_copia"
  out[ltr & grepl("^Ty3", third)] <- "Class_I/LTR/Ty3_gypsy"
  out
}


# --- 1. fragment merging (design 6.0.1) ------------------------------

#' Collapse annotations that are one protein domain split in two.
#'
#' A frameshift or in-frame stop can break a domain into two adjacent
#' DANTE hits, each with a reduced Relat_Length.  Lineage mode discards
#' any cluster with a repeated domain name outright
#' (clean_domain_clusters(), ltr_utils.R:206), so such elements never
#' reach its output at all.  Core mode rejoins them instead, because a
#' leftover half is a core domain outside the seed, which the blocking
#' walk of design 6.3 would treat as another element's core and stop at.
#'
#' The signature that separates a split from a genuine duplication is
#' that the two hits tile *complementary* parts of the same reference
#' domain, in element order.  Two real copies hit overlapping parts.
#'
#' Must run BEFORE dante_filtering: each fragment's Relat_Length is
#' reduced by construction, so filtering first destroys exactly what we
#' are trying to rejoin (design 5.3).
#'
#' @param g GRanges of DANTE domains.
#' @param max_gap max bp between fragments.  The frameshift excess is
#'   concentrated below ~20 bp and complementary tiling stops being
#'   specific beyond a few hundred bp (design 5.3, 5.4), hence 50.
#' @param max_db_overlap max overlap of the two reference intervals,
#'   as a fraction of the shorter.
#' @return GRanges, sorted, with an Nfragments column.
merge_split_domains <- function(g, max_gap = 50, max_db_overlap = 0.3) {
  if (length(g) == 0) {
    return(g)
  }
  g <- g[order(as.character(seqnames(g)), start(g), end(g))]
  mcols(g)$Nfragments <- 1L
  if (length(g) < 2) {
    return(g)
  }

  n <- length(g)
  i <- seq_len(n - 1)
  j <- i + 1

  same <- as.character(seqnames(g))[i] == as.character(seqnames(g))[j] &
    as.character(mcols(g)$Name)[i] == as.character(mcols(g)$Name)[j] &
    as.character(strand(g))[i] == as.character(strand(g))[j]

  gap <- start(g)[j] - end(g)[i] - 1L
  close_enough <- gap >= -MAX_GENOME_OVERLAP & gap <= max_gap

  db <- parse_db_pos(mcols(g)$Best_Hit_DB_Pos)
  a_start <- db[i, "start"]; a_end <- db[i, "end"]
  b_start <- db[j, "start"]; b_end <- db[j, "end"]
  parsed <- !is.na(a_start) & !is.na(b_start)

  # fragments of one domain appear in reference order along the element;
  # on the minus strand the element runs right-to-left, so the test flips
  plus <- as.character(strand(g))[i] == "+"
  ordered_ok <- ifelse(plus, b_start >= a_start, a_start >= b_start)

  overlap <- pmin(a_end, b_end) - pmax(a_start, b_start) + 1
  shorter <- pmin(a_end - a_start + 1, b_end - b_start + 1)
  frac <- ifelse(shorter > 0, pmax(0, overlap) / shorter, 1)

  mergeable <- same & close_enough & parsed & ordered_ok & frac <= max_db_overlap
  mergeable[is.na(mergeable)] <- FALSE

  if (!any(mergeable)) {
    return(g)
  }

  grp <- cumsum(c(TRUE, !mergeable))
  keep_idx <- which(!duplicated(grp))
  out <- g[keep_idx]

  sizes <- as.integer(table(grp)[as.character(grp[keep_idx])])
  multi <- which(sizes > 1)
  if (length(multi) > 0) {
    rel <- num_mcol(g, "Relat_Length")
    sim <- num_mcol(g, "Similarity")
    sp <- split(seq_len(n), grp)
    for (k in multi) {
      members <- sp[[as.character(grp[keep_idx[k]])]]
      # keep the dominant fragment's annotation, widen to span both
      dominant <- members[which.max(rel[members])]
      mcols(out)[k, ] <- mcols(g)[dominant, ]
      start(out)[k] <- min(start(g)[members])
      end(out)[k] <- max(end(g)[members])
      mcols(out)$Relat_Length[k] <- as.character(min(1, sum(rel[members])))
      mcols(out)$Similarity[k] <- as.character(max(sim[members]))
      mcols(out)$Nfragments[k] <- length(members)
    }
  }
  out
}


# --- 2. core candidate selection (design 6.1) ------------------------

#' Domains eligible to seed an element.
#'
#' RT/RH/INT with LTR-retrotransposon support, where support is either
#' the domain's own Final_Classification or -- and this matters on
#' distant genomes, where a genuine LTR RT can lose its best hit to a
#' pararetrovirus or a LINE -- any entry of Region_Hits_Classifications.
#' The latter are marked Domain_LTR_Support = "secondary".
core_candidates <- function(g) {
  if (length(g) == 0) {
    return(g)
  }
  is_core <- as.character(mcols(g)$Name) %in% CORE_DOMAINS
  primary <- is_ltr_classification(mcols(g)$Final_Classification)
  rhc <- flatten_region_hits(g)
  secondary <- !is.na(rhc) & grepl("\\|Class_I\\|LTR", rhc)

  keep <- is_core & (primary | secondary)
  out <- g[keep]
  if (length(out) > 0) {
    mcols(out)$Domain_LTR_Support <-
      ifelse(primary[keep], "primary", "secondary")
  }
  out
}


# --- 3. ordered core seeding (design 6.2) ----------------------------

#' Enumerate ordered RT/RH/INT seeds.
#'
#' Per sequence and strand: build runs of core candidates separated by at
#' most core_max_gap, enumerate index triples spelling a valid core order
#' in element orientation within the gap and span caps, then greedily
#' accept non-overlapping triples, most compact first.
#'
#' Ordering is fully deterministic -- (span, -sum(Similarity), start) --
#' so repeated runs on the same input are byte-identical.  The positional
#' key is not cosmetic: without it ties resolve by whatever order order()
#' happens to produce.
#'
#' @param candidates GRanges from core_candidates().
#' @param constraints data.frame from read_core_constraints().
#' @return data.frame with one row per seed: seqnames, strand,
#'   superfamily, start, end, and i1/i2/i3 indexing `candidates`.
find_core_seeds <- function(candidates, constraints) {
  empty <- data.frame(seqnames = character(0), strand = character(0),
                      superfamily = character(0), start = integer(0),
                      end = integer(0), i1 = integer(0), i2 = integer(0),
                      i3 = integer(0), stringsAsFactors = FALSE)
  if (length(candidates) == 0) {
    return(empty)
  }

  nm <- as.character(mcols(candidates)$Name)
  sim <- num_mcol(candidates, "Similarity")
  sq <- as.character(seqnames(candidates))
  st <- as.character(strand(candidates))
  S <- start(candidates)
  E <- end(candidates)

  # caps are per superfamily; a run is examined against both, so take the
  # widest of the two to build runs and re-check per candidate triple
  max_gap_any <- max(constraints$core_max_gap)

  seeds <- list()
  groups <- split(seq_along(candidates), paste(sq, st, sep = "\r"))
  for (gi in groups) {
    gi <- gi[order(S[gi], E[gi])]
    if (length(gi) < 3) {
      next
    }
    # split into runs on gaps larger than the cap
    gaps <- S[gi][-1] - E[gi][-length(gi)] - 1L
    run_id <- cumsum(c(1L, as.integer(gaps > max_gap_any)))
    for (run in split(gi, run_id)) {
      if (length(run) < 3) {
        next
      }
      strand_here <- st[run[1]]
      ord <- if (strand_here == "+") run else rev(run)
      cand <- .enumerate_core_triples(ord, nm, S, E, sim, strand_here,
                                      constraints)
      if (nrow(cand) == 0) {
        next
      }
      # greedy: most compact first, then best summed Similarity, then
      # leftmost -- see the determinism note above
      cand <- cand[order(cand$span, -cand$qual, cand$start, cand$end), ,
                   drop = FALSE]
      used <- integer(0)
      for (r in seq_len(nrow(cand))) {
        members <- c(cand$i1[r], cand$i2[r], cand$i3[r])
        if (any(members %in% used)) {
          next
        }
        used <- c(used, members)
        seeds[[length(seeds) + 1]] <- data.frame(
          seqnames = sq[members[1]], strand = strand_here,
          superfamily = cand$superfamily[r],
          start = cand$start[r], end = cand$end[r],
          i1 = members[1], i2 = members[2], i3 = members[3],
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(seeds) == 0) {
    return(empty)
  }
  out <- do.call(rbind, seeds)
  out[order(out$seqnames, out$start, out$end), , drop = FALSE]
}


#' Candidate triples within one run, in element orientation.
#' `ord` is already ordered 5'->3' for the element.
.enumerate_core_triples <- function(ord, nm, S, E, sim, strand_here,
                                    constraints) {
  n <- length(ord)
  rows <- list()
  # element-orientation gap between two positions in `ord`
  eg <- function(a, b) {
    if (strand_here == "+") S[b] - E[a] - 1L else S[a] - E[b] - 1L
  }
  for (a in seq_len(n - 2)) {
    for (b in seq(a + 1, n - 1)) {
      gap_ab <- eg(ord[a], ord[b])
      if (gap_ab > max(constraints$core_max_gap)) {
        break
      }
      for (cc in seq(b + 1, n)) {
        gap_bc <- eg(ord[b], ord[cc])
        if (gap_bc > max(constraints$core_max_gap)) {
          break
        }
        members <- c(ord[a], ord[b], ord[cc])
        trio <- nm[members]
        sf <- NA_character_
        for (k in names(CORE_ORDER)) {
          if (identical(as.character(trio), CORE_ORDER[[k]])) {
            sf <- k
            break
          }
        }
        if (is.na(sf)) {
          next
        }
        row <- constraints[constraints$Superfamily == sf, , drop = FALSE]
        if (nrow(row) == 0) {
          next
        }
        if (gap_ab > row$core_max_gap[1] || gap_bc > row$core_max_gap[1]) {
          next
        }
        lo <- min(S[members]); hi <- max(E[members])
        if (hi - lo + 1 > row$core_max_span[1]) {
          next
        }
        rows[[length(rows) + 1]] <- data.frame(
          superfamily = sf, start = lo, end = hi, span = hi - lo + 1,
          qual = sum(sim[members]),
          i1 = members[1], i2 = members[2], i3 = members[3],
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(superfamily = character(0), start = integer(0),
                      end = integer(0), span = integer(0), qual = numeric(0),
                      i1 = integer(0), i2 = integer(0), i3 = integer(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}


#' Read databases/core_domain_order.csv (design 5).
read_core_constraints <- function(path) {
  t <- read.table(path, sep = "\t", header = TRUE, as.is = TRUE,
                  comment.char = "", quote = "")
  required <- c("Superfamily", "core_order", "offset5prime", "offset3prime",
                "core_max_gap", "core_max_span", "ltr_length",
                "accessory_5", "accessory_3")
  missing <- setdiff(required, names(t))
  if (length(missing) > 0) {
    stop("core constraints table is missing column(s): ",
         paste(missing, collapse = ", "))
  }
  unknown <- setdiff(t$Superfamily, names(CORE_ORDER))
  if (length(unknown) > 0) {
    stop("unknown superfamily in constraints table: ",
         paste(unknown, collapse = ", "))
  }
  t
}


#' Accessory domain names allowed on one side of the core, in element
#' orientation.  "-" means none are (a copia core's 3' side).
accessory_names <- function(constraints, superfamily, side) {
  col <- if (side == "5") "accessory_5" else "accessory_3"
  v <- constraints[constraints$Superfamily == superfamily, col][1]
  if (is.na(v) || v == "-" || v == "") {
    return(character(0))
  }
  strsplit(v, " +")[[1]]
}
