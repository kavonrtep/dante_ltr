#!/usr/bin/env Rscript
# tests/library_policy_selftest.R -- unit tests for utils/library_policy.R.
#
# Drives classify_cluster() with synthetic label/element vectors, so no genome,
# no mmseqs and no library build are needed.  That is not only for speed: the
# repository's fixtures are too small to form a cluster of three windows, and
# the element-vs-window counting rule cannot be exercised end-to-end at all on
# the request's genome, where both conventions give the same answer.  These
# tests are its only coverage.

initial_options <- commandArgs(trailingOnly = FALSE)
root <- dirname(dirname(normalizePath(sub("--file=", "",
  initial_options[grep("--file=", initial_options)]))))
source(file.path(root, "utils", "library_policy.R"))

# resolve_name() as shipped, so the strict cases test the real thing
resolve_name <- function(x){
  if (length(x)==1){ return(x) } else{
    y <- sapply(x, strsplit, split="|", fixed = TRUE)
    ny <- table(unlist(sapply(y, function(x)paste(seq_along(x), x))))
    if (max(ny)<length(x)){ return("Unknown") } else{
      k <- which(ny==length(x)); r <- max(as.numeric((gsub(" .+", "", names(k)))))
      return(paste(y[[1]][1:r], collapse="|"))
    }
  }
}

FAIL <- 0L; N <- 0L
ok <- function(label, cond, detail = "") {
  N <<- N + 1L
  if (isTRUE(cond)) cat("  ok   ", label, "\n", sep = "")
  else { FAIL <<- FAIL + 1L
         cat("  FAIL ", label, if (nzchar(detail)) paste0("  [", detail, "]") else "",
             "\n", sep = "") }
}
section <- function(x) cat("\n== ", x, " ==\n", sep = "")

G  <- "Class_I|LTR|Ty3/gypsy"
C  <- "Class_I|LTR|Ty3/gypsy|chromovirus"
TK <- "Class_I|LTR|Ty3/gypsy|chromovirus|Tekay"
RE <- "Class_I|LTR|Ty3/gypsy|chromovirus|Reina"
CO <- "Class_I|LTR|Ty1/copia"

# one window per element unless a test says otherwise
strict <- function(lab) classify_cluster(lab, seq_along(lab), 0.95, "strict",
                                         resolve_name = resolve_name)
nested <- function(lab, el = seq_along(lab), ...)
  classify_cluster(lab, el, 0.95, "nested", ...)


section("label relations")
ok("ancestor recognised",            is_ancestor_or_equal(G, C))
ok("descendant is not an ancestor",  !is_ancestor_or_equal(C, G))
ok("equal counts as ancestor",       is_ancestor_or_equal(C, C))
ok("boundary respected: gypsy|chromo is not an ancestor of gypsy|chromovirus",
   !is_ancestor_or_equal("Class_I|LTR|Ty3/gypsy|chromo", C))
ok("chain of two",   is_single_chain(c(G, C)))
ok("chain of three", is_single_chain(c(G, C, TK)))
ok("single label is a chain", is_single_chain(c(C, C, C)))
ok("siblings are not a chain",        !is_single_chain(c(TK, RE)))
ok("siblings plus parent not a chain", !is_single_chain(c(G, TK, RE)))
ok("cross-superfamily is not a chain", !is_single_chain(c(G, CO)))


section("strict -- the historical rule, unchanged")
r <- strict(rep(C, 10))
ok("uniform cluster kept", r$keep && r$label == C && r$outcome == "majority")
r <- strict(c(rep(C, 39), G))
ok("majority above 95% kept, labelled with the majority",
   r$keep && r$label == C && r$outcome == "majority")
r <- strict(c(rep(C, 6), rep(G, 4)))
ok("ancestor+descendant, majority is the descendant -> DROPPED (the bug)",
   !r$keep)
r <- strict(c(rep(G, 6), rep(C, 4)))
ok("ancestor+descendant, majority is the ancestor -> kept as LCA",
   r$keep && r$label == G && r$outcome == "lca")
r <- strict(c(rep(TK, 6), rep(RE, 4)))
ok("siblings dropped", !r$keep)
r <- strict(c(rep(G, 6), rep(TK, 2), rep(RE, 2)))
ok("siblings plus parent, parent in majority -> kept by strict",
   r$keep && r$label == G)
r <- strict(c(rep(G, 6), rep(CO, 4)))
ok("cross-superfamily dropped", !r$keep)


section("nested -- R2, recover ancestor chains")
# Recovery and promotion usually coincide: if strict dropped the cluster then
# the majority was the *descendant*, which therefore already clears min_share.
r <- nested(c(rep(C, 6), rep(G, 4)))
ok("the strict bug case is recovered, and promotion carries it to the descendant",
   r$keep && r$label == C && r$outcome == "promoted")
# Recovered without promotion needs the deepest label to be too thin: majority
# is a descendant (so strict drops), but only one element carries the deepest.
r <- nested(c(rep(G, 3), rep(C, 5), TK))
ok("recovered and left at the LCA when the deepest label is too thin",
   r$keep && r$label == G && r$outcome == "recovered")
r <- nested(c(rep(TK, 6), rep(RE, 4)))
ok("siblings still dropped", !r$keep)
r <- nested(c(rep(G, 6), rep(CO, 4)))
ok("cross-superfamily still dropped", !r$keep)
# decided behaviour: nested follows R2 literally, so it is NOT a superset of
# strict.  This cluster is in the library today and is not under nested.
r <- nested(c(rep(G, 6), rep(TK, 2), rep(RE, 2)))
ok("siblings plus parent DROPPED under nested, though strict keeps it",
   !r$keep)


section("nested -- R4, promotion")
r <- nested(c(rep(C, 6), rep(TK, 4)))          # 40% carry the deepest
ok("promoted to the deepest label",
   r$keep && r$label == TK && r$outcome == "promoted")
r <- nested(c(rep(C, 9), TK))                  # 10% -- below min_share
ok("not promoted when the share is too small",
   r$keep && r$label == C && r$outcome == "lca")
r <- nested(c(rep(C, 2), TK), promote_min_share = 0.25)  # 33% but 1 element
ok("not promoted when only one element carries it",
   r$keep && r$label == C)
r <- nested(c(rep(C, 6), rep(TK, 4)), promote_min_elements = 5)
ok("min_elements is honoured", r$keep && r$label == C)
r <- nested(c(rep(C, 6), rep(TK, 4)), promote_min_share = 0.5)
ok("min_share is honoured", r$keep && r$label == C)
r <- nested(c(rep(TK, 39), C))
ok("majority above 95% is kept as the majority, never promoted",
   r$keep && r$label == TK && r$outcome == "majority")


section("nested -- R3, elements not windows")
# One long element cut into 10 windows, plus 2 short elements.  The deepest
# label is carried by 2 of 12 windows (17%, below min_share) but by 2 of 3
# elements (67%, above it) -- so the counting unit decides the label.
lab <- c(rep(G, 10), rep(C, 2))
el  <- c(rep("long", 10), "s1", "s2")
r <- nested(lab, el)
ok("by elements the cluster promotes to the descendant",
   r$keep && r$label == C && r$outcome == "promoted")
r2 <- nested(lab, seq_along(lab))
ok("counting windows instead leaves it at the LCA",
   r2$keep && r2$label == G && r2$label != r$label)


section("tie-breaks must stay sorted-order, not first-seen")
# "Class_I|LTR|Ty1/copia" sorts before "Class_I|LTR|Ty3/gypsy"; a 5/5 tie must
# resolve to copia regardless of which appears first in the vector.
a <- strict(c(rep(G, 5), rep(CO, 5)))
b <- strict(c(rep(CO, 5), rep(G, 5)))
ok("strict tie-break is independent of member order",
   identical(a$keep, b$keep) && identical(a$outcome, b$outcome))
tie <- c(rep(C, 5), rep(TK, 5))
x <- nested(tie); y <- nested(rev(tie))
ok("nested tie-break is independent of member order", identical(x$label, y$label))


cat("\n", N - FAIL, "/", N, " checks passed\n", sep = "")
if (FAIL > 0L) quit(save = "no", status = 1)
quit(save = "no", status = 0)
