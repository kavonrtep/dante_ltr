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


# --- 4. search-space delimitation (design 6.3) -----------------------

# This is the only part with no lineage-mode analogue.  get_ranges_left()
# and get_ranges_right() (ltr_utils.R:241,254) clamp the BLAST window at
# the nearest neighbouring domain, which works there because the
# element's own GAG/PROT are inside the cluster.  In core mode they are
# outside it, so the naive clamp would stop the window dead at the
# element's own GAG.  Measurement says this is not a corner case: 49% of
# Pisum seeds and 99% of Drapa seeds have a domain within 500 bp of the
# core edge, almost all of them PROT/CHD/GAG (design 5.3, 5.4).
#
# The replacement is an outward walk governed by one principle: a domain
# that survived dante_filtering is real and blocks, unless it can be
# positively explained as part of this element.

#' Per-seed limits for the outward LTR search.
#'
#' @param seeds data.frame from find_core_seeds().
#' @param g_block GRanges of domains at the *standard* filter threshold.
#'   Seed candidates admitted only by the relaxed core threshold are
#'   deliberately absent: a domain too weak to be trusted should not
#'   truncate a window (design 6.0.2).
#' @param constraints core constraints table.
#' @return the seeds data.frame with added columns left_limit,
#'   right_limit (genomic coordinates of the nearest blocking feature
#'   edge, or NA when the walk ran to the offset cap),
#'   accessory_5, accessory_3, accessory_sf_mismatch.
seed_search_limits <- function(seeds, g_block, constraints) {
  seeds$left_limit <- NA_integer_
  seeds$right_limit <- NA_integer_
  seeds$accessory_5 <- ""
  seeds$accessory_3 <- ""
  seeds$accessory_sf_mismatch <- 0L
  if (nrow(seeds) == 0) {
    return(seeds)
  }

  ok_ltr <- is_ltr_classification(mcols(g_block)$Final_Classification)
  rhc <- flatten_region_hits(g_block)
  ok_ltr <- ok_ltr | (!is.na(rhc) & grepl("\\|Class_I\\|LTR", rhc))
  dom_sf <- superfamily_of(mcols(g_block)$Final_Classification)
  dom_name <- as.character(mcols(g_block)$Name)
  dom_strand <- as.character(strand(g_block))
  dom_seq <- as.character(seqnames(g_block))
  dom_start <- start(g_block)
  dom_end <- end(g_block)

  by_seq <- split(seq_along(g_block), dom_seq)
  by_seq <- lapply(by_seq, function(ix) ix[order(dom_start[ix], dom_end[ix])])

  for (r in seq_len(nrow(seeds))) {
    ix <- by_seq[[seeds$seqnames[r]]]
    if (is.null(ix)) {
      next
    }
    plus <- seeds$strand[r] == "+"
    core_lo <- seeds$start[r]
    core_hi <- seeds$end[r]
    sf <- seeds$superfamily[r]

    # walking left in coordinate space is walking 5' on the plus strand
    left <- .walk_side(
      ix = rev(ix[dom_end[ix] < core_lo]), direction = "left",
      side = if (plus) "5" else "3", sf = sf, seed_strand = seeds$strand[r],
      dom_name = dom_name, dom_strand = dom_strand, dom_sf = dom_sf,
      ok_ltr = ok_ltr, constraints = constraints)
    right <- .walk_side(
      ix = ix[dom_start[ix] > core_hi], direction = "right",
      side = if (plus) "3" else "5", sf = sf, seed_strand = seeds$strand[r],
      dom_name = dom_name, dom_strand = dom_strand, dom_sf = dom_sf,
      ok_ltr = ok_ltr, constraints = constraints)

    if (!is.na(left$blocker)) {
      seeds$left_limit[r] <- dom_end[left$blocker]
    }
    if (!is.na(right$blocker)) {
      seeds$right_limit[r] <- dom_start[right$blocker]
    }
    acc_left <- paste(left$traversed, collapse = " ")
    acc_right <- paste(right$traversed, collapse = " ")
    if (plus) {
      seeds$accessory_5[r] <- acc_left
      seeds$accessory_3[r] <- acc_right
    } else {
      seeds$accessory_5[r] <- acc_right
      seeds$accessory_3[r] <- acc_left
    }
    seeds$accessory_sf_mismatch[r] <- left$sf_mismatch + right$sf_mismatch
  }
  seeds
}


#' Walk outward on one side until something blocks.
#' `ix` is already ordered from the core outward.
.walk_side <- function(ix, direction, side, sf, seed_strand, dom_name,
                       dom_strand, dom_sf, ok_ltr, constraints) {
  allowed <- accessory_names(constraints, sf, side)
  traversed <- character(0)
  last_rank <- 0L
  sf_mismatch <- 0L

  for (k in ix) {
    # 1. same strand
    if (dom_strand[k] != seed_strand) {
      break
    }
    # 2. LTR-retrotransposon support
    if (!ok_ltr[k]) {
      break
    }
    # 3. an accessory type allowed on this side, in element orientation
    rank <- match(dom_name[k], allowed)
    if (is.na(rank)) {
      break
    }
    # 4. + 5. not already traversed, and in canonical inward->outward
    #    order -- a second GAG going outward belongs to the neighbouring
    #    element, and GAG before PROT means the walk has crossed into one
    if (rank <= last_rank) {
      break
    }
    # 6. superfamily must not contradict the seed's order-derived call.
    #    Fires only on positive disagreement: a domain resolving no
    #    deeper than Class_I|LTR passes vacuously, so this cannot bite on
    #    the shallow-classification case core mode exists for.
    if (!is.na(dom_sf[k]) && dom_sf[k] != sf) {
      sf_mismatch <- sf_mismatch + 1L
      break
    }
    traversed <- c(traversed, dom_name[k])
    last_rank <- rank
  }

  blocker <- NA_integer_
  n_pass <- length(traversed)
  if (length(ix) > n_pass) {
    blocker <- ix[n_pass + 1L]
  }
  list(blocker = blocker, traversed = traversed,
       traversed_idx = if (n_pass > 0) ix[seq_len(n_pass)] else integer(0),
       sf_mismatch = sf_mismatch)
}


#' Left-hand BLAST window, mirroring get_ranges_left() (ltr_utils.R:241)
#' including its +100 over-extension into the blocking feature and its
#' offset2 overlap into the core, with the per-seed limit substituted for
#' `upstream_domain`.
core_ranges_left <- function(seeds, constraints, offset2 = 300) {
  offs <- .side_offsets(seeds, constraints)
  # Anchored on the core, not on the outermost accessory domain the walk
  # reached.  Anchoring on the accessory domain was tried, on the theory
  # that it stops the window reaching into the element's own GAG/PROT
  # region where get_TE()'s innermost-pair preference can prefer a
  # spurious internal repeat.  Measured on Pisum it repaired 2 elements
  # and broke 6: when the traversed domain actually belongs to the
  # neighbouring element, the anchor jumps past the true LTR and the
  # element is lost entirely.  See design §10.
  S <- seeds$start
  limit <- ifelse(is.na(seeds$left_limit), 1L, seeds$left_limit)
  max_offset <- S - limit + 100
  adjusted <- pmin(max_offset, offs$left)
  starts <- pmax(1L, as.integer(S - adjusted))
  GRanges(seqnames = seeds$seqnames,
          ranges = IRanges(start = starts, end = as.integer(S + offset2)))
}


#' Right-hand BLAST window, mirroring get_ranges_right()
#' (ltr_utils.R:254).  `SL` is seqlengths of the reference.
core_ranges_right <- function(seeds, constraints, SL, offset2 = 300) {
  offs <- .side_offsets(seeds, constraints)
  E <- seeds$end
  seq_end <- as.integer(SL[seeds$seqnames])
  limit <- ifelse(is.na(seeds$right_limit), seq_end, seeds$right_limit)
  max_offset <- limit - E + 100
  adjusted <- pmin(max_offset, offs$right)
  ends <- pmin(seq_end, as.integer(E + adjusted))
  GRanges(seqnames = seeds$seqnames,
          ranges = IRanges(start = as.integer(E - offset2), end = ends))
}


#' Table offsets mapped onto coordinate sides.  On the minus strand the
#' element's 5' end is at higher coordinates, so the offsets swap.
.side_offsets <- function(seeds, constraints) {
  i <- match(seeds$superfamily, constraints$Superfamily)
  o5 <- constraints$offset5prime[i]
  o3 <- constraints$offset3prime[i]
  plus <- seeds$strand == "+"
  list(left = ifelse(plus, o5, o3), right = ifelse(plus, o3, o5))
}


# --- 5. reannotation and classification (design 7) -------------------

# Core mode calls the superfamily structurally, from domain order, and
# only then asks what the classifications say.  The order is the
# authority: it is structural evidence, while a classification is a best
# hit against a database that may be phylogenetically distant.

# The constraints table names superfamilies with "/" separators; DANTE
# and REXdb use "|" and spell the superfamily Ty1/copia, Ty3/gypsy.
SUPERFAMILY_LABEL <- c(
  "Class_I/LTR/Ty1_copia" = "Class_I|LTR|Ty1/copia",
  "Class_I/LTR/Ty3_gypsy" = "Class_I|LTR|Ty3/gypsy"
)


#' Lowest common ancestor of REXdb classification paths.
#' Returns "" when there is no shared root.
lca_classification <- function(labels) {
  labels <- unique(as.character(labels[!is.na(labels) & nzchar(labels)]))
  if (length(labels) == 0) {
    return(NA_character_)
  }
  parts <- strsplit(labels, "|", fixed = TRUE)
  common <- parts[[1]]
  for (p in parts[-1]) {
    n <- min(length(common), length(p))
    same <- which(common[seq_len(n)] != p[seq_len(n)])
    keep <- if (length(same) == 0) n else same[1] - 1L
    common <- common[seq_len(keep)]
    if (length(common) == 0) {
      return("")
    }
  }
  paste(common, collapse = "|")
}


#' Is `child` the same as, or below, `parent` in the REXdb tree?
is_at_or_below <- function(child, parent) {
  if (is.na(child) || is.na(parent) || !nzchar(child)) {
    return(FALSE)
  }
  child == parent || startsWith(child, paste0(parent, "|"))
}


#' Lineage labels named in a Region_Hits_Classifications string.
#' Entries look like "RT|Class_I|LTR|Ty3/gypsy|chromovirus|Tekay[524bp]".
region_hit_labels <- function(rhc) {
  if (is.na(rhc) || !nzchar(rhc)) {
    return(character(0))
  }
  items <- strsplit(rhc, ",", fixed = TRUE)[[1]]
  items <- sub("\\[[0-9]+bp\\]$", "", trimws(items))
  items <- sub("^[^|]+\\|", "", items)      # drop the leading domain name
  items[nzchar(items)]
}


#' Is an observed domain order compatible with a lineage's canonical one?
#'
#' Compatible means the observed domains that the lineage knows about
#' appear in the lineage's order (a subsequence), with at most
#' `max_extra` observed domains the lineage does not have at all.
#'
#' This deliberately does NOT use domain_distance() (ltr_utils.R:285),
#' which lineage mode uses for the same purpose.  That function compares
#' `d_query_p == d_reference_p[d_reference_p %in% d_query_p]` without
#' checking lengths, so whenever the query carries more domains than the
#' reference the comparison recycles -- R warns, and the returned
#' distance is meaningless.  In lineage mode that is rare, because a
#' cluster is already keyed to one lineage.  In core mode an incomplete
#' or unexpected domain complement is the normal case, which is the
#' whole point of the mode, so the unsound path would be the common one.
#'
#' Missing domains are not penalised: an element carrying only RT/RH/INT
#' is genuinely compatible with every gypsy lineage, and reporting that
#' breadth honestly is what Lineage_Candidates is for.
order_compatible_with <- function(observed, reference, max_extra = 0) {
  ref <- strsplit(as.character(reference), " +")[[1]]
  obs <- as.character(observed)
  extra <- sum(!(obs %in% ref))
  if (extra > max_extra) {
    return(FALSE)
  }
  obs <- obs[obs %in% ref]
  # greedy subsequence match
  at <- 0L
  for (d in obs) {
    nxt <- match(d, ref[seq.int(at + 1L, length(ref))])
    if (is.na(nxt)) {
      return(FALSE)
    }
    at <- at + nxt
  }
  TRUE
}


#' Classify one delimited element (design 7).
#'
#' @param domain_names domain names in element orientation.
#' @param domain_classifications their Final_Classification values.
#' @param region_hits their Region_Hits_Classifications, flattened.
#' @param superfamily order-derived superfamily, "/" form.
#' @param lineage_domain named character vector: DANTE-form lineage name
#'   -> its canonical domain order (as built in detect_putative_ltr.R).
#' @param max_missing_domains tolerance passed to domain_distance().
#' @return list of the attributes described in design 8.
classify_core_element <- function(domain_names, domain_classifications,
                                  region_hits, superfamily, lineage_domain,
                                  max_missing_domains = 0) {
  sf_label <- unname(SUPERFAMILY_LABEL[superfamily])

  # lineages of this superfamily whose canonical domain order is
  # compatible with what the element actually carries
  under_sf <- names(lineage_domain)[
    vapply(names(lineage_domain), is_at_or_below, logical(1), parent = sf_label)]
  order_compatible <- character(0)
  if (length(under_sf) > 0) {
    keep <- vapply(under_sf, function(ln) {
      order_compatible_with(domain_names, lineage_domain[[ln]],
                            max_missing_domains)
    }, logical(1))
    order_compatible <- under_sf[keep]
  }

  cls <- as.character(domain_classifications)
  cls <- cls[!is.na(cls) & nzchar(cls)]
  n_total <- length(domain_names)

  # Final_Classification: LCA of the per-domain calls, clipped at the
  # structural superfamily.  If the LCA is not at or below the
  # order-derived superfamily -- the domains disagree with each other
  # across superfamilies, or with the order -- the order wins.
  lca <- lca_classification(cls)
  conflict <- !is_at_or_below(lca, sf_label)
  final <- if (conflict) sf_label else lca

  # a lineage-depth call needs to be order-compatible *and* carried by a
  # majority of the element's domains
  lineage_call <- NA_character_
  if (!conflict && length(order_compatible) > 0) {
    tab <- table(cls[cls %in% order_compatible])
    if (length(tab) > 0) {
      top <- names(tab)[which.max(tab)]
      if (max(tab) > n_total / 2 && sum(tab == max(tab)) == 1) {
        lineage_call <- top
        if (is_at_or_below(top, final) && top != final) {
          # a clear majority justifies reporting the deeper label
          final <- top
        }
      }
    }
  }

  # candidates: order-compatible lineages named anywhere in the domains'
  # hit lists, most-supported first.  Advisory only -- these never set
  # Final_Classification.
  seen <- unlist(lapply(region_hits, region_hit_labels), use.names = FALSE)
  seen <- seen[seen %in% order_compatible]
  candidates <- character(0)
  if (length(seen) > 0) {
    tab <- sort(table(seen), decreasing = TRUE)
    candidates <- names(tab)
  }

  supporting <- sum(vapply(cls, is_at_or_below, logical(1), parent = final))
  deepest <- if (length(cls) == 0) 0L else max(lengths(strsplit(cls, "|", fixed = TRUE)))
  reported_depth <- if (is.na(final) || !nzchar(final)) 0L else
    length(strsplit(final, "|", fixed = TRUE)[[1]])

  list(
    Final_Classification = final,
    Lineage_Call = lineage_call,
    Lineage_Candidates = if (length(candidates) > 0)
      paste(candidates, collapse = ",") else NA_character_,
    Lineage_Support = paste0(supporting, "/", n_total),
    Classification_Demoted = reported_depth < deepest,
    Classification_Conflict = conflict,
    Superfamily_Evidence = paste0("domain_order:",
                                  paste(CORE_ORDER[[superfamily]],
                                        collapse = ",")),
    Core_Domains = paste(CORE_ORDER[[superfamily]], collapse = " ")
  )
}


#' Lineage name -> canonical domain order, in DANTE's label form.
#' Mirrors the conversion in detect_putative_ltr.R so both modes read the
#' same lineage_domain_order.csv identically.
lineage_domain_map <- function(lineage_info) {
  nm <- gsub("ss/I", "ss_I",
             gsub("_", "/", gsub("/", "|", lineage_info$Lineage)))
  setNames(lineage_info$Domains.order, nm)
}


#' Per-rank element counts, same shape as get_te_statistics().
#'
#' get_te_statistics() (ltr_utils.R:988) takes its row labels from the RT
#' domains' classifications.  In lineage mode an element's classification
#' is by construction one of those, so nothing is lost.  Core mode
#' computes the element's classification instead (design 7), and a
#' demoted or conflict-resolved label need not appear among the RT
#' domains at all -- such elements would silently vanish from the table
#' and from its Total row.
#'
#' Columns are identical (D, DL, DLT, DLP, DLTP, RT_domain) so the
#' Python sum_up_stats_files() merges chunk outputs unchanged.
get_core_te_statistics <- function(gr, RT) {
  ranks <- c("D", "DL", "DLT", "DLP", "DLTP")
  te <- gr[gr$type == "transposable_element"]
  rt_class <- as.character(RT$Final_Classification)
  te_class <- as.character(te$Final_Classification)

  rt_tab <- sort(table(rt_class), decreasing = TRUE)
  te_tab <- sort(table(te_class), decreasing = TRUE)
  all_class <- unique(c(names(rt_tab), names(te_tab)))

  rank_table <- lapply(ranks, function(r) {
    as.integer(table(factor(te_class[te$Rank == r], levels = all_class)))
  })
  names(rank_table) <- ranks
  out <- cbind(do.call(cbind, rank_table),
               RT_domain = as.integer(table(factor(rt_class,
                                                   levels = all_class))))
  out <- rbind(out, Total = colSums(out))
  rownames(out) <- c(all_class, "Total")
  out
}


#' Indices of domains fully contained in [lo, hi] on one sequence.
#'
#' Built for repeated queries: the caller precomputes the per-sequence
#' index once, and each lookup is a binary search rather than a scan.
#' The naive form -- re-deriving as.character(seqnames(g)) inside a loop
#' over elements -- is quadratic in (domains x elements), which is
#' tolerable on a test fixture and intractable on a real chunk carrying
#' 10^5 domains and 10^3 seeds.
#'
#' @param index list from build_domain_index().
#' @param seqname sequence to query.
#' @param lo,hi inclusive bounds.
#' @param strand optional; restrict to this strand.
domains_within <- function(index, seqname, lo, hi, strand = NULL) {
  ix <- index$by_seq[[seqname]]
  if (is.null(ix) || length(ix) == 0) {
    return(integer(0))
  }
  s <- index$start[ix]
  first <- findInterval(lo - 1L, s) + 1L
  last <- findInterval(hi, s)
  if (first > last || first > length(ix)) {
    return(integer(0))
  }
  cand <- ix[first:last]
  cand <- cand[index$end[cand] <= hi]
  if (!is.null(strand)) {
    cand <- cand[index$strand[cand] == strand]
  }
  cand
}


#' Precompute the coordinate index used by domains_within().
build_domain_index <- function(g) {
  st <- start(g)
  ord <- order(as.character(seqnames(g)), st)
  by_seq <- split(ord, as.character(seqnames(g))[ord])
  list(by_seq = by_seq, start = st, end = end(g),
       strand = as.character(strand(g)),
       seqnames = as.character(seqnames(g)))
}
