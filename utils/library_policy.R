# library_policy.R -- the keep/label decision for one cluster of the repeat
# library, factored out of utils/mmseq_clustering.R so it can be tested without
# running mmseqs over a real library.
#
# See docs/dante_ltr_core_library_policy_implementation_plan.md.
#
# Two policies:
#
#   strict   the historical rule, byte-identical.  Counts sliding windows.
#            A cluster survives only when the label it is given equals its
#            majority label -- so a cluster whose members are an ancestor and
#            its descendant (Ty3/gypsy and Ty3/gypsy|chromovirus) is discarded,
#            because the LCA is the ancestor while the majority may be the
#            descendant.  That is the behaviour core mode trips over constantly.
#
#   nested   counts distinct source elements rather than windows, keeps a
#            cluster whose labels form a single ancestor chain, and may promote
#            it to the deepest label when enough elements carry it.
#
# TIE-BREAKS.  `which.max(table(x))` picks the first label in *sorted* order,
# not the first seen.  That is load-bearing: the library must be a function of
# the input set, not of record order (the same reason mmseq_clustering.R sorts
# canonically before clustering).  Every majority here goes through
# `which.max(table(...))` for exactly that reason -- do not replace it with
# something that depends on encounter order.


#' Split a DANTE classification into its levels.
.levels_of <- function(x) strsplit(x, "|", fixed = TRUE)[[1]]


#' Is `a` the same as, or an ancestor of, `b`?  Compared at "|" boundaries, so
#' "Ty3/gypsy|chromo" is not treated as an ancestor of "Ty3/gypsy|chromovirus".
is_ancestor_or_equal <- function(a, b) {
  la <- .levels_of(a); lb <- .levels_of(b)
  length(la) <= length(lb) && identical(la, lb[seq_along(la)])
}


#' Do these distinct labels form a single ancestor chain -- is every pair
#' ancestor-or-descendant?  This is stricter than "the LCA is one of the
#' labels": it rejects {gypsy, chromovirus|Tekay, chromovirus|Reina}, whose LCA
#' is present but which would put two named lineages under one representative.
is_single_chain <- function(labels) {
  u <- unique(labels)
  if (length(u) <= 1) return(TRUE)
  u <- u[order(lengths(lapply(u, .levels_of)))]
  for (i in seq_len(length(u) - 1)) {
    if (!is_ancestor_or_equal(u[i], u[i + 1])) return(FALSE)
  }
  TRUE
}


#' Shallowest / deepest label of a set already known to form a chain.
.shallowest <- function(labels) {
  u <- unique(labels); u[which.min(lengths(lapply(u, .levels_of)))]
}
.deepest <- function(labels) {
  u <- unique(labels); u[which.max(lengths(lapply(u, .levels_of)))]
}


#' The historical rule, unchanged.
#'
#' `resolve_name` is passed in rather than duplicated, so the strict path keeps
#' calling the function that produced the shipped output.
.classify_strict <- function(labels, proportion_min, resolve_name) {
  tab  <- table(labels)
  prop <- max(tab) / length(labels)
  main <- names(which.max(tab))
  label <- if (prop > proportion_min) main else resolve_name(unique(labels))
  if (identical(label, main)) {
    list(keep = TRUE, label = label,
         outcome = if (prop > proportion_min) "majority" else "lca")
  } else {
    list(keep = FALSE, label = NA_character_, outcome = "dropped")
  }
}


#' Ancestor-chain recovery and lineage promotion.
.classify_nested <- function(labels, elements, proportion_min,
                             promote_min_elements, promote_min_share) {
  # one vote per source element, not per 1 kb window: a long element is cut
  # into ~8 windows and would otherwise outvote several short ones
  keep_first <- !duplicated(elements)
  el_labels  <- labels[keep_first]
  n <- length(el_labels)

  tab  <- table(el_labels)
  prop <- max(tab) / n
  main <- names(which.max(tab))
  chain <- is_single_chain(el_labels)

  if (!(prop > proportion_min || chain)) {
    return(list(keep = FALSE, label = NA_character_, outcome = "dropped"))
  }

  if (prop > proportion_min) {
    # Promotion is deliberately not evaluated here, and cannot apply: a
    # majority above the threshold leaves every other label below
    # 1 - proportion_min, far under promote_min_share, and a deepest label that
    # *is* the majority is already this cluster's label.  Measured on the
    # request's genome: 0 of the 69 promotions came from this branch.
    return(list(keep = TRUE, label = main, outcome = "majority"))
  }

  base <- .shallowest(el_labels)
  deep <- .deepest(el_labels)
  n_deep <- sum(el_labels == deep)

  if (deep != base && n_deep >= promote_min_elements &&
        n_deep / n >= promote_min_share) {
    # refinement along one path, never a reclassification
    stopifnot(is_ancestor_or_equal(base, deep))
    return(list(keep = TRUE, label = deep, outcome = "promoted"))
  }
  list(keep = TRUE, label = base,
       outcome = if (identical(base, main)) "lca" else "recovered")
}


#' Decide one cluster.
#'
#' @param labels    classification of every cluster member (one per window)
#' @param elements  source element id of every member, same length as `labels`
#' @param policy    "strict" or "nested"
#' @return list(keep, label, outcome); outcome is one of
#'   "majority", "lca", "recovered", "promoted", "dropped"
classify_cluster <- function(labels, elements, proportion_min,
                             policy = c("strict", "nested"),
                             promote_min_elements = 2,
                             promote_min_share = 0.25,
                             resolve_name = NULL) {
  policy <- match.arg(policy)
  if (policy == "strict") {
    stopifnot(is.function(resolve_name))
    .classify_strict(labels, proportion_min, resolve_name)
  } else {
    stopifnot(length(labels) == length(elements))
    .classify_nested(labels, elements, proportion_min,
                     promote_min_elements, promote_min_share)
  }
}
