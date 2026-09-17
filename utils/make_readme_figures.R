#!/usr/bin/env Rscript
# make_readme_figures.R -- generate the three detection-mode diagrams used in
# README.md.  Base R graphics only, so it needs nothing beyond the runtime
# dante_ltr already depends on.
#
#   utils/make_readme_figures.R [output_dir]     (default: repository root)
#
# Produces:
#   dante_ltr_workflow.png    the default lineage-keyed pipeline (Principle)
#   dante_ltr_fallback.png    what --fallback_mode changes
#   dante_ltr_core.png        what --mode core changes
#
# The three share one visual vocabulary -- same domain glyph, same colours,
# same stage layout -- so a reader moving between README sections can compare
# them directly.
#
# Sizing: these are read inside a GitHub README column, roughly 900 px wide,
# so they are drawn at 2x that and every text size is chosen to stay legible
# after the browser halves them.  The stage caption sits *above* its track
# rather than in a left-hand gutter, which hands the genome track the full
# width -- the gutter was costing a fifth of the figure and shrinking every
# glyph with it.
#
# Domain orders follow databases/lineage_domain_order.csv:
#   Ty1/copia          GAG PROT INT RT RH
#   Ty3/gypsy chromo   GAG PROT RT RH INT CHD

args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- if (length(args) >= 1) args[1] else "."

# --- palette --------------------------------------------------------------
# Okabe-Ito, which stays distinguishable under the common forms of colour
# blindness.  Classification colours and structural-feature colours are kept
# disjoint so that a colour never means two things.
COL <- list(
  lin1    = "#0072B2",  # blue            a REXdb lineage
  lin2    = "#56B4E9",  # sky blue        another lineage
  lin3    = "#D55E00",  # vermillion      another lineage
  lin4    = "#CC79A7",  # reddish purple  another lineage
  sfam    = "#009E73",  # bluish green    superfamily-level call
  shallow = "#8A8A8A",  # grey            not resolved to lineage depth
  ltr     = "#3A3A3A",  # structural: long terminal repeat
  tsd     = "#E69F00",  # structural: target site duplication
  pbs     = "#7D3C98",  # structural: primer binding site
  rule    = "#CFCFCF",  # genome line
  ink     = "#1A1A1A",
  muted   = "#6E6E6E",
  ok      = "#0B7A4B",
  no      = "#B3261E"
)

FADE <- 0.16   # opacity of context that is not being acted on

# --- geometry -------------------------------------------------------------
TRACK_X <- 2      # the track now spans essentially the whole figure
TRACK_W <- 96
DH      <- 3.6    # domain glyph height
LAB_DY  <- 6.6    # caption sits this far above its track
NOTE_DY <- 5.4    # annotations sit this far below it

CEX_TITLE <- 0.82
CEX_SUB   <- 0.68
CEX_DOM   <- 0.70
CEX_NOTE  <- 0.62

# --- primitives -----------------------------------------------------------

new_panel <- function() {
  par(mar = c(0.1, 0.1, 0.1, 0.1), xaxs = "i", yaxs = "i", family = "sans")
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), type = "n",
       axes = FALSE, xlab = "", ylab = "")
}

#' y centre of stage i (1 = top), for a panel of n stages.
row_y <- function(n, top, bottom) {
  if (n == 1) return((top + bottom) / 2)
  top - (seq_len(n) - 1) * (top - bottom) / (n - 1)
}

#' Caption line above a track: "n  Title  sub-label", with an optional
#' right-aligned outcome.
stage_label <- function(i, y, title, sub = NULL, verd = NULL, ok = TRUE) {
  yy <- y + LAB_DY
  num <- as.character(i)
  text(TRACK_X, yy, num, adj = c(0, 0.5), cex = CEX_TITLE, font = 2,
       col = COL$muted)
  x <- TRACK_X + strwidth(num, cex = CEX_TITLE, font = 2) + 1.4
  text(x, yy, title, adj = c(0, 0.5), cex = CEX_TITLE, font = 2,
       col = COL$ink)
  if (!is.null(sub)) {
    x <- x + strwidth(title, cex = CEX_TITLE, font = 2) + 1.6
    text(x, yy, sub, adj = c(0, 0.5), cex = CEX_SUB, col = COL$muted)
  }
  if (!is.null(verd)) {
    text(TRACK_X + TRACK_W, yy,
         paste0(if (ok) "✓  " else "✕  ", verd),
         adj = c(1, 0.5), cex = CEX_SUB, col = if (ok) COL$ok else COL$no,
         font = 2)
  }
}

genome_line <- function(y, x0 = TRACK_X, x1 = TRACK_X + TRACK_W) {
  segments(x0, y, x1, y, col = COL$rule, lwd = 1.4)
}

#' Readable label colour for a given fill.
ink_on <- function(col) {
  v <- col2rgb(col) / 255
  lum <- 0.2126 * v[1] + 0.7152 * v[2] + 0.0722 * v[3]
  if (lum > 0.6) "#12212B" else "#FFFFFF"
}

#' One protein domain, drawn as a chevron pointing in the strand direction.
domain <- function(x, w, y, label, col, dir = 1, alpha = 1, cex = CEX_DOM) {
  tip <- min(w * 0.30, 1.8)
  if (dir >= 0) {
    xs <- c(x, x + w - tip, x + w, x + w - tip, x)
  } else {
    xs <- c(x + w, x + tip, x, x + tip, x + w)
  }
  ys <- c(y - DH / 2, y - DH / 2, y, y + DH / 2, y + DH / 2)
  polygon(xs, ys, col = adjustcolor(col, alpha),
          border = adjustcolor(col, alpha), lwd = 0.8)
  text(x + w / 2 + (if (dir >= 0) -tip / 3 else tip / 3), y, label, cex = cex,
       col = adjustcolor(ink_on(col), alpha), font = 2)
}

#' Draw a run of domains from a data.frame(x, w, label, col, dir).
domains <- function(d, y, alpha = 1, ...) {
  for (i in seq_len(nrow(d))) {
    domain(d$x[i], d$w[i], y, d$label[i], d$col[i],
           dir = if ("dir" %in% names(d)) d$dir[i] else 1, alpha = alpha, ...)
  }
}

#' The window searched for an LTR.
#'
#' The constraints table allows the search to run several kb out from the
#' element -- `x_offset` -- but it stops at the first annotated domain in the
#' way.  Both extents are drawn: the dashed outline is what the table allows,
#' the filled part is what is actually searched, and a bar marks where a
#' domain cut it short.  Accessory domains the walk may pass lie *inside* the
#' filled part.
search_window <- function(x_inner, x_offset, y, x_block = NULL) {
  h <- DH * 1.15
  lo <- min(x_inner, x_offset); hi <- max(x_inner, x_offset)
  rect(lo, y - h, hi, y + h, col = NA, border = adjustcolor(COL$ltr, 0.40),
       lty = 2, lwd = 0.8)
  if (is.null(x_block)) {
    slo <- lo; shi <- hi
  } else if (x_offset < x_inner) {
    slo <- x_block; shi <- x_inner
  } else {
    slo <- x_inner; shi <- x_block
  }
  rect(slo, y - h, shi, y + h, col = adjustcolor(COL$ltr, 0.10), border = NA)
  if (!is.null(x_block)) {
    segments(x_block, y - h, x_block, y + h, col = COL$no, lwd = 2.0)
  }
}

#' A long terminal repeat.
ltr <- function(x, w, y, label = NULL, alpha = 1) {
  tip <- min(w * 0.24, 1.4)
  xs <- c(x, x + w - tip, x + w, x + w - tip, x)
  ys <- c(y - DH / 2, y - DH / 2, y, y + DH / 2, y + DH / 2)
  polygon(xs, ys, col = adjustcolor(COL$ltr, alpha),
          border = adjustcolor(COL$ltr, alpha), lwd = 0.8)
  if (!is.null(label)) {
    text(x + w / 2 - tip / 3, y, label, cex = CEX_DOM * 0.9,
         col = adjustcolor("white", alpha), font = 2)
  }
}

tsd_mark <- function(x, y, w = 1.5) {
  rect(x - w / 2, y - DH / 2, x + w / 2, y + DH / 2,
       col = COL$tsd, border = COL$tsd)
}

pbs_mark <- function(x, y, w = 1.9) {
  rect(x - w / 2, y - DH / 2 * 0.72, x + w / 2, y + DH / 2 * 0.72,
       col = COL$pbs, border = COL$pbs)
}

note <- function(x, y, text, cex = CEX_NOTE, col = COL$muted,
                 adj = c(0.5, 0.5), font = 1) {
  text(x, y, text, cex = cex, col = col, adj = adj, font = font)
}

#' Bracket under a span, used to mark the core triplet.
brace <- function(x0, x1, y, label, col = COL$ink) {
  yy <- y - DH / 2 - 1.4
  segments(x0, yy, x1, yy, col = col, lwd = 1.4)
  segments(c(x0, x1), yy, c(x0, x1), yy + 1.0, col = col, lwd = 1.4)
  text((x0 + x1) / 2, yy - 2.6, label, cex = CEX_NOTE, col = col, font = 2)
}

open_png <- function(file, width_px, height_px, res = 170) {
  png(file.path(OUTDIR, file), width = width_px, height = height_px,
      res = res, type = "cairo", bg = "white")
}

#' Lay a run of domains out left to right from x0.
run <- function(x0, labels, widths, col, gap = 1.2, dir = 1) {
  xs <- numeric(length(widths))
  xs[1] <- x0
  for (i in seq_along(widths)[-1]) xs[i] <- xs[i - 1] + widths[i - 1] + gap
  data.frame(x = xs, w = widths, label = labels, col = col, dir = dir,
             stringsAsFactors = FALSE)
}

W_GAG <- 5.9; W_PROT <- 4.7; W_INT <- 5.4; W_RT <- 5.1; W_RH <- 4.1
W_CHD <- 4.9; W_LTR <- 7.4


# =========================================================================
# Figure 1 -- the default, lineage-keyed pipeline
# =========================================================================
fig_workflow <- function() {
  open_png("dante_ltr_workflow.png", 1800, 900)
  new_panel()

  ys <- row_y(5, top = 84, bottom = 14)

  h <- run(37, c("GAG", "PROT", "INT", "RT", "RH"),
           c(W_GAG, W_PROT, W_INT, W_RT, W_RH), COL$lin1)
  hx0 <- min(h$x); hx1 <- max(h$x + h$w)              # 37 .. 66.8
  nb <- run(14.5, c("RT", "RH"), c(W_RT, W_RH), COL$lin4)
  nb_end <- max(nb$x + nb$w)                          # 24.9

  w5_off <- hx0 - 24; w5_in <- hx0 - 0.7; w5_block <- nb_end + 0.7
  w3_in  <- hx1 + 0.7; w3_off <- hx1 + 24

  x_ltr5 <- 27.0; x_ltr3 <- 70.5

  # 1 -- DANTE input
  stage_label(1, ys[1], "DANTE protein domains",
              "colour = the lineage DANTE assigned")
  genome_line(ys[1]); domains(nb, ys[1]); domains(h, ys[1])

  # 2 -- cluster
  stage_label(2, ys[2], "Cluster of domains",
              "one lineage, expected order and spacing",
              verd = "complement complete", ok = TRUE)
  genome_line(ys[2]); domains(nb, ys[2], alpha = FADE); domains(h, ys[2])

  # 3 -- search window
  stage_label(3, ys[3], "Search window for LTRs",
              "several kb out, from the constraints table")
  genome_line(ys[3])
  search_window(w5_in, w5_off, ys[3], x_block = w5_block)
  search_window(w3_in, w3_off, ys[3])
  domains(nb, ys[3], alpha = FADE); domains(h, ys[3])
  note((w5_off + w5_block) / 2, ys[3] - NOTE_DY,
       "stops at the neighbouring domain", col = COL$no)
  note((w3_in + w3_off) / 2, ys[3] - NOTE_DY,
       "nothing in the way: the full offset is searched")

  # 4 -- LTRs
  stage_label(4, ys[4], "LTR detection",
              "closest direct repeat pair, TG ... CA")
  genome_line(ys[4])
  search_window(w5_in, w5_off, ys[4], x_block = w5_block)
  search_window(w3_in, w3_off, ys[4])
  ltr(x_ltr5, W_LTR, ys[4], "5'LTR"); ltr(x_ltr3, W_LTR, ys[4], "3'LTR")
  domains(nb, ys[4], alpha = FADE); domains(h, ys[4])
  note(x_ltr5 + 0.6, ys[4] - NOTE_DY, "TG", col = COL$ink)
  note(x_ltr3 + W_LTR - 0.6, ys[4] - NOTE_DY, "CA", col = COL$ink)

  # 5 -- TSD / PBS
  stage_label(5, ys[5], "TSD and PBS",
              "presence of these sets the rank",
              verd = "rank DLTP", ok = TRUE)
  genome_line(ys[5]); domains(nb, ys[5], alpha = FADE)
  ltr(x_ltr5, W_LTR, ys[5], "5'LTR"); ltr(x_ltr3, W_LTR, ys[5], "3'LTR")
  domains(h, ys[5])
  tsd_mark(x_ltr5 - 1.3, ys[5]); tsd_mark(x_ltr3 + W_LTR + 1.3, ys[5])
  pbs_mark(hx0 - 1.4, ys[5])
  note(x_ltr5 - 1.3, ys[5] - NOTE_DY, "TSD", col = COL$tsd, font = 2)
  note(x_ltr3 + W_LTR + 1.3, ys[5] - NOTE_DY, "TSD", col = COL$tsd, font = 2)
  note(hx0 - 1.4, ys[5] - NOTE_DY, "PBS", col = COL$pbs, font = 2)

  dev.off()
}


# =========================================================================
# Figure 2 -- --fallback_mode
# =========================================================================
fig_fallback <- function() {
  open_png("dante_ltr_fallback.png", 1800, 610)
  new_panel()

  ys <- row_y(3, top = 80, bottom = 26)

  # one Ty3/gypsy chromovirus element: two closely related lineages plus
  # domains DANTE could not resolve below chromovirus
  el <- run(31, c("GAG", "PROT", "RT", "RH", "INT"),
            c(W_GAG, W_PROT, W_RT, W_RH, W_INT), COL$lin1)
  el$col <- c(COL$lin1, COL$shallow, COL$lin1, COL$lin4, COL$shallow)
  calls  <- c("Tekay", "chromovirus", "Tekay", "Reina", "chromovirus")
  x_mid  <- el$x + el$w / 2
  ex0 <- min(el$x); ex1 <- max(el$x + el$w)

  # 1 -- the problem
  stage_label(1, ys[1], "DANTE, genome far from REXdb",
              "one element, ambiguous classification",
              verd = "no cluster forms", ok = FALSE)
  genome_line(ys[1]); domains(el, ys[1])
  for (i in seq_len(nrow(el))) note(x_mid[i], ys[1] - NOTE_DY, calls[i])

  # 2 -- demotion
  stage_label(2, ys[2], "--fallback_mode",
              "all calls demoted to one depth")
  genome_line(ys[2])
  el2 <- el; el2$col <- COL$sfam
  domains(el2, ys[2])
  for (i in seq_len(nrow(el))) {
    note(x_mid[i], ys[2] - NOTE_DY, "Ty3/gypsy", col = COL$sfam)
  }

  # 3 -- cluster forms, pipeline continues
  stage_label(3, ys[3], "Cluster forms",
              "the rest of the pipeline is unchanged",
              verd = "element detected", ok = TRUE)
  genome_line(ys[3]); domains(el2, ys[3])
  ltr(ex0 - W_LTR - 2.6, W_LTR, ys[3], "5'LTR")
  ltr(ex1 + 2.6, W_LTR, ys[3], "3'LTR")

  note(TRACK_X, 6, paste("The domain complement must still be complete and",
                         "in the expected order — only its classification",
                         "is relaxed."), adj = c(0, 0.5))

  dev.off()
}


# =========================================================================
# Figure 3 -- --mode core
# =========================================================================
fig_core <- function() {
  open_png("dante_ltr_core.png", 1800, 790)
  new_panel()

  ys <- row_y(4, top = 82, bottom = 22)

  # a chromovirus with no detectable GAG and no call below superfamily
  el <- run(29, c("PROT", "RT", "RH", "INT", "CHD"),
            c(W_PROT, W_RT, W_RH, W_INT, W_CHD), COL$shallow)
  core_i <- 2:4
  cx0 <- min(el$x[core_i]); cx1 <- max(el$x[core_i] + el$w[core_i])
  blocker <- data.frame(x = 72.0, w = 6.2, label = "TPase", col = COL$lin3,
                        dir = -1, stringsAsFactors = FALSE)

  # the window runs out from the *core*, and is free to pass over the
  # accessory domains on the way
  w5_off <- cx0 - 22; w5_in <- cx0 - 0.7
  w3_in  <- cx1 + 0.7; w3_off <- cx1 + 22; w3_block <- blocker$x - 0.7

  el2 <- el; el2$col[core_i] <- COL$sfam

  # 1 -- the problem
  stage_label(1, ys[1], "DANTE, no close REXdb reference",
              "no GAG; calls stop at superfamily",
              verd = "complement incomplete", ok = FALSE)
  genome_line(ys[1]); domains(el, ys[1]); domains(blocker, ys[1])

  # 2 -- core seed
  stage_label(2, ys[2], "Core seed: RT → RH → INT",
              "the order alone gives the superfamily")
  genome_line(ys[2])
  domains(el, ys[2], alpha = FADE); domains(blocker, ys[2], alpha = FADE)
  domains(el2[core_i, ], ys[2])
  brace(cx0, cx1, ys[2], "RT RH INT  ⇒  Ty3/gypsy")

  # 3 -- outward walk
  stage_label(3, ys[3], "Search window, measured from the core",
              "accessory domains are passed over")
  genome_line(ys[3])
  search_window(w5_in, w5_off, ys[3])
  search_window(w3_in, w3_off, ys[3], x_block = w3_block)
  domains(el2, ys[3]); domains(blocker, ys[3])
  note(el$x[1] + el$w[1] / 2, ys[3] - NOTE_DY, "passed over")
  note(el$x[5] + el$w[5] / 2, ys[3] - NOTE_DY, "passed over")
  note(blocker$x + blocker$w / 2, ys[3] - NOTE_DY, "blocks", col = COL$no,
       font = 2)

  # 4 -- outcome
  stage_label(4, ys[4], "LTRs, then classification",
              "from the domains inside the element",
              verd = "rank DLTP, Ty3/gypsy", ok = TRUE)
  genome_line(ys[4]); domains(blocker, ys[4], alpha = FADE)
  x_ltr5 <- el$x[1] - W_LTR - 3.0
  x_ltr3 <- max(el$x + el$w) + 1.4
  ltr(x_ltr5, W_LTR, ys[4], "5'LTR"); ltr(x_ltr3, W_LTR, ys[4], "3'LTR")
  domains(el2, ys[4])
  tsd_mark(x_ltr5 - 1.3, ys[4]); tsd_mark(x_ltr3 + W_LTR + 1.3, ys[4])
  pbs_mark(el$x[1] - 1.4, ys[4])

  note(TRACK_X, 6, paste("Demoting the classification would not rescue this",
                         "element: the accessory domains are missing, not",
                         "merely named too deeply."), adj = c(0, 0.5))

  dev.off()
}

fig_workflow()
fig_fallback()
fig_core()
cat("figures written to ", normalizePath(OUTDIR), "\n", sep = "")
