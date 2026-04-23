#!/usr/bin/env Rscript
# boundary_report.R --- Render a self-contained HTML report of the
# boundary-adjustment / TG-CA / TSD signals emitted by build_ltr_library.R.
#
# Input : one *_LTR_library_boundary_qc.tsv (from the builder) and,
#         optionally, the alignments directory (for per-cluster
#         conservation-profile mini-panels).
# Output: a single .html with PNGs inlined as base64, no sibling files.
#
# Usage:
#   Rscript utils/boundary_report.R \
#     --qc   tmp/alyr_lib/lib_LTR_library_boundary_qc.tsv \
#     --out  tmp/alyr_lib/lib_LTR_library_boundary_report.html \
#    [--alignments tmp/alyr_lib/alignments] \
#    [--top_n 30]

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

# Pure-R base64 encoder — avoids adding base64enc / openssl / jsonlite
# as dependencies.  ~200 KB/s on my laptop; plenty for a handful of PNGs.
b64_table <- c(LETTERS, letters, as.character(0L:9L), "+", "/")
base64_encode_raw <- function(x) {
  n <- length(x)
  pad <- (3L - n %% 3L) %% 3L
  if (pad > 0L) x <- c(x, as.raw(rep(0L, pad)))
  ints <- as.integer(x)
  n3   <- length(ints) / 3L
  m    <- matrix(ints, nrow = 3L)
  b1 <- bitwShiftR(m[1L, ], 2L)
  b2 <- bitwOr(bitwShiftL(bitwAnd(m[1L, ], 0x03L), 4L),
               bitwShiftR(m[2L, ], 4L))
  b3 <- bitwOr(bitwShiftL(bitwAnd(m[2L, ], 0x0FL), 2L),
               bitwShiftR(m[3L, ], 6L))
  b4 <- bitwAnd(m[3L, ], 0x3FL)
  out <- character(4L * n3)
  out[seq.int(1L, by = 4L, length.out = n3)] <- b64_table[b1 + 1L]
  out[seq.int(2L, by = 4L, length.out = n3)] <- b64_table[b2 + 1L]
  out[seq.int(3L, by = 4L, length.out = n3)] <- b64_table[b3 + 1L]
  out[seq.int(4L, by = 4L, length.out = n3)] <- b64_table[b4 + 1L]
  s <- paste0(out, collapse = "")
  if (pad > 0L) substr(s, 1L, nchar(s) - pad) <- substr(s, 1L, nchar(s) - pad)  # no-op, keep
  if (pad > 0L) {
    s <- paste0(substr(s, 1L, nchar(s) - pad), strrep("=", pad))
  }
  s
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--qc",         type = "character", default = NULL),
  make_option("--out",        type = "character", default = NULL),
  make_option("--alignments", type = "character", default = NULL),
  make_option("--top_n",      type = "integer",   default = 30L)
)))

for (arg in c("qc", "out")) {
  if (is.null(opt[[arg]])) stop("missing --", arg, call. = FALSE)
}
if (!file.exists(opt$qc)) stop("no QC tsv at ", opt$qc, call. = FALSE)

qc <- read.table(opt$qc, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 quote = "", comment.char = "")

# Backward-compatibility: if run against an older QC TSV without the
# TG/CA / TSD columns, fill them with NA so the report renders the
# aggregate panels and shows "(no data)" for the validation scatters.
for (cn in c("tgca_frac_annotated", "tgca_frac_corrected",
             "tsd_frac_annotated",  "tsd_frac_corrected")) {
  if (!(cn %in% names(qc))) qc[[cn]] <- NA_real_
}

# Normalise: some columns are NA for the single-member fallback path;
# that's fine, the plots filter them out.
n_total <- nrow(qc)

# ----------------------------------------------------------------------
# Small PNG-to-base64 helper.  Render with base R graphics into a
# tempfile, read bytes, encode, wrap in <img>.
# ----------------------------------------------------------------------
render_png <- function(expr, width = 900, height = 450, res = 100) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  png(tmp, width = width, height = height, res = res)
  on_err <- function(e) { dev.off(); stop(e) }
  tryCatch(expr, error = on_err)
  dev.off()
  b <- readBin(tmp, what = "raw", n = file.info(tmp)$size)
  paste0('<img class="panel" src="data:image/png;base64,',
         base64_encode_raw(b), '">')
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;",  x, fixed = TRUE)
  x <- gsub("<", "&lt;",   x, fixed = TRUE)
  x <- gsub(">", "&gt;",   x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

# ----------------------------------------------------------------------
# Lineage family (top-level grouping) — best-effort split of
# Final_Classification strings of the form "A|B|C|Copia|...":
# we return the first segment that matches copia/gypsy/other.
# ----------------------------------------------------------------------
lineage_family <- function(cls) {
  cls_lc <- tolower(as.character(cls))
  ifelse(grepl("copia",   cls_lc), "Ty1/Copia",
  ifelse(grepl("gypsy|ty3", cls_lc), "Ty3/Gypsy",
  ifelse(is.na(cls) | cls == "",    "unknown",
                                    "other")))
}
qc$family <- lineage_family(qc$Final_Classification)

# Subset that actually went through the MAFFT/change-point path.
msa <- qc[!is.na(qc$corrected_5_col), , drop = FALSE]

# ----------------------------------------------------------------------
# Panel 1 — cluster size distribution
# ----------------------------------------------------------------------
panel_sizes <- render_png({
  old <- par(mfrow = c(1L, 2L), mar = c(4, 4, 3, 1))
  on.exit(par(old), add = TRUE)
  h <- hist(qc$n_members, breaks = "FD", plot = FALSE)
  plot(h, main = "cluster size (n_members)", xlab = "n_members",
       col = "#cccccc", border = "#555555")
  abline(v = 6L, col = "red", lty = 2L)
  mtext("dashed: MAFFT threshold", side = 3, line = 0.1, cex = 0.8)
  counts <- c(MAFFT = nrow(msa), fallback = n_total - nrow(msa))
  barplot(counts, col = c("#6baed6", "#cccccc"), main = "path taken",
          ylab = "# clusters")
}, width = 900, height = 380)

# ----------------------------------------------------------------------
# Panel 2 — shift distribution, faceted by flank (50 vs wide-flank retry)
# ----------------------------------------------------------------------
panel_shifts <- render_png({
  old <- par(mfrow = c(2L, 2L), mar = c(4, 4, 3, 1))
  on.exit(par(old), add = TRUE)
  msa$flank_label <- ifelse(msa$flank > 50L, "wide-flank retry",
                            "default flank (50)")
  for (flk in c("default flank (50)", "wide-flank retry")) {
    sub <- msa[msa$flank_label == flk, , drop = FALSE]
    if (!nrow(sub)) {
      plot.new(); title(paste0("shift_5  [", flk, ": n=0]")); next
    }
    s5 <- sub$shift_5[sub$detected_5]
    if (length(s5) > 0L) {
      hist(s5, breaks = 40L, col = "#6baed6", border = "#3079a8",
           main = sprintf("shift_5   %s  (n=%d)", flk, length(s5)),
           xlab = "shift_5 (bp)")
      abline(v = 0, col = "red", lty = 2L)
    } else plot.new()
  }
  for (flk in c("default flank (50)", "wide-flank retry")) {
    sub <- msa[msa$flank_label == flk, , drop = FALSE]
    if (!nrow(sub)) {
      plot.new(); title(paste0("shift_3  [", flk, ": n=0]")); next
    }
    s3 <- sub$shift_3[sub$detected_3]
    if (length(s3) > 0L) {
      hist(s3, breaks = 40L, col = "#fd8d3c", border = "#b85a10",
           main = sprintf("shift_3   %s  (n=%d)", flk, length(s3)),
           xlab = "shift_3 (bp)")
      abline(v = 0, col = "red", lty = 2L)
    } else plot.new()
  }
}, width = 900, height = 600)

# ----------------------------------------------------------------------
# Generic "annotated vs corrected" scatter panel (used for TG/CA + TSD).
# ----------------------------------------------------------------------
scatter_ann_vs_corr <- function(df, col_ann, col_corr, title,
                                palette = c("Ty1/Copia" = "#1f77b4",
                                            "Ty3/Gypsy" = "#ff7f0e",
                                            "other"     = "#7f7f7f",
                                            "unknown"   = "#bbbbbb")) {
  render_png({
    old <- par(mar = c(4, 4, 3, 1))
    on.exit(par(old), add = TRUE)
    a <- df[[col_ann]]; c <- df[[col_corr]]
    keep <- !is.na(a) & !is.na(c)
    if (!any(keep)) {
      plot.new(); title(paste0(title, " [no data]"))
    } else {
      a <- a[keep]; cc <- c[keep]
      fam <- df$family[keep]; col_v <- palette[fam]; col_v[is.na(col_v)] <- "#7f7f7f"
      nm <- df$n_members[keep]
      cex_v <- pmax(0.4, pmin(2.0, log10(pmax(1, nm)) * 0.9))
      set.seed(1L)
      jit <- function(x) x + runif(length(x), -0.01, 0.01)
      plot(jit(a), jit(cc), pch = 21, bg = col_v, col = "#00000044",
           cex = cex_v, xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02),
           xlab = sprintf("%s (annotated)", title),
           ylab = sprintf("%s (corrected)", title),
           main = sprintf("%s:  annotated vs corrected  (n=%d)", title, sum(keep)))
      abline(0, 1, col = "black", lty = 1L)
      legend("topleft", legend = names(palette), fill = palette,
             bty = "n", cex = 0.8)
      delta <- cc - a
      lost   <- sum(delta < -0.1); gained <- sum(delta >  0.1)
      unchg  <- sum(abs(delta) <= 0.1)
      mtext(sprintf("LOST=%d  GAINED=%d  UNCHANGED=%d",
                    lost, gained, unchg), side = 3, line = 0.1, cex = 0.9)
    }
  }, width = 700, height = 560)
}

# Panel 3 — TG/CA annotated vs corrected
panel_tgca <- scatter_ann_vs_corr(msa, "tgca_frac_annotated",
                                  "tgca_frac_corrected", "TG/CA frac")

# Panel 4 — TSD annotated vs corrected
panel_tsd  <- scatter_ann_vs_corr(msa, "tsd_frac_annotated",
                                  "tsd_frac_corrected", "TSD frac")

# ----------------------------------------------------------------------
# Panel 5 — per-lineage boxplots (lineage x {tgca, tsd, |shift|})
# ----------------------------------------------------------------------
panel_lineage <- render_png({
  keep <- !is.na(msa$Final_Classification)
  sub <- msa[keep, , drop = FALSE]
  tb <- table(sub$Final_Classification)
  keep_lin <- names(tb)[tb >= 3L]
  sub <- sub[sub$Final_Classification %in% keep_lin, , drop = FALSE]
  if (!nrow(sub)) {
    plot.new(); title("no lineage with >= 3 clusters")
  } else {
    short <- function(x) {
      parts <- strsplit(x, "|", fixed = TRUE)
      vapply(parts, function(p) {
        if (length(p) >= 2L) paste(tail(p, 2L), collapse = "/") else x
      }, character(1L))
    }
    sub$lin_short <- short(sub$Final_Classification)
    old <- par(mfrow = c(3L, 1L), mar = c(7, 4, 2, 1))
    on.exit(par(old), add = TRUE)
    safe_boxplot <- function(y, fill, main, ylab) {
      if (all(is.na(y))) { plot.new(); title(main); return(invisible()) }
      d <- data.frame(y = y, g = sub$lin_short)
      d <- d[!is.na(d$y), , drop = FALSE]
      if (!nrow(d) || length(unique(d$g)) < 1L) {
        plot.new(); title(main); return(invisible())
      }
      boxplot(y ~ g, data = d, las = 2L, col = fill,
              main = main, ylab = ylab, cex.axis = 0.7, xlab = "")
    }
    safe_boxplot(sub$tgca_frac_corrected, "#6baed6",
                 "TG/CA frac (corrected) per lineage", "TG/CA frac")
    safe_boxplot(sub$tsd_frac_corrected, "#fd8d3c",
                 "TSD frac (corrected) per lineage", "TSD frac")
    safe_boxplot(pmax(abs(sub$shift_5), abs(sub$shift_3)), "#cccccc",
                 "max(|shift_5|, |shift_3|) per lineage", "bp")
  }
}, width = 900, height = 900)

# ----------------------------------------------------------------------
# Panel 6 — grow-cap activity (clusters where detection fired but both
# boundaries reverted to annotated because the corrected total length
# was out of [shrink, grow] range — detection=TRUE but shift==0).
# ----------------------------------------------------------------------
panel_growcap <- render_png({
  if (!nrow(msa)) {
    plot.new(); title("no MAFFT clusters")
  } else {
    reverted <- msa$detected_5 & msa$detected_3 &
                msa$shift_5 == 0L & msa$shift_3 == 0L
    by_fam <- table(msa$family, ifelse(reverted, "grow-cap reverted",
                                       "kept as detected"))
    barplot(t(by_fam), beside = TRUE, col = c("#cccccc", "#31a354"),
            legend.text = TRUE, args.legend = list(x = "topright", bty = "n"),
            main = sprintf("grow-cap activity  (reverted=%d / detected-both=%d)",
                           sum(reverted), sum(msa$detected_5 & msa$detected_3)),
            ylab = "# clusters", las = 1L)
  }
}, width = 900, height = 420)

# ----------------------------------------------------------------------
# Panel 7 — top-N "clusters of interest" table
# Score: n_members * (|shift_5| + |shift_3|) + 200 * (|delta_tgca| + |delta_tsd|)
# ----------------------------------------------------------------------
score <- function(r) {
  abs_sh  <- abs(r$shift_5) + abs(r$shift_3)
  dtgca   <- abs(r$tgca_frac_corrected - r$tgca_frac_annotated)
  dtsd    <- abs(r$tsd_frac_corrected  - r$tsd_frac_annotated)
  dtgca[is.na(dtgca)] <- 0
  dtsd[is.na(dtsd)]   <- 0
  abs_sh[is.na(abs_sh)] <- 0
  (r$n_members * abs_sh / 10) + 200 * (dtgca + dtsd)
}
msa_with_score <- msa
msa_with_score$score <- score(msa_with_score)
msa_with_score <- msa_with_score[order(-msa_with_score$score), , drop = FALSE]
top_n_use <- min(opt$top_n, nrow(msa_with_score))
top_clusters <- if (top_n_use > 0L) msa_with_score[seq_len(top_n_use), , drop = FALSE] else msa_with_score[0, , drop = FALSE]

top_table_html <- (function() {
  if (!nrow(top_clusters)) return("<p>(no clusters)</p>")
  cols <- c("ltr_id", "Final_Classification", "n_members", "flank",
            "shift_5", "shift_3",
            "tgca_frac_annotated", "tgca_frac_corrected",
            "tsd_frac_annotated",  "tsd_frac_corrected")
  fmt_row <- function(r) {
    id    <- r[["ltr_id"]]
    cell <- function(v, is_num = FALSE) {
      if (is.null(v)) return("")
      if (is.na(v))   return('<td class="na">NA</td>')
      if (is_num && is.numeric(v)) {
        if (abs(v) < 1 && v != 0) return(sprintf('<td>%.3f</td>', v))
        return(sprintf('<td>%g</td>', v))
      }
      paste0('<td>', html_escape(v), '</td>')
    }
    link <- sprintf('<td><a href="#%s">%s</a></td>',
                    html_escape(id), html_escape(id))
    num_cols <- c("n_members", "flank", "shift_5", "shift_3",
                  "tgca_frac_annotated", "tgca_frac_corrected",
                  "tsd_frac_annotated",  "tsd_frac_corrected")
    rest <- vapply(cols[-1L], function(cn) cell(r[[cn]], cn %in% num_cols),
                   character(1L))
    paste0("<tr>", link, paste(rest, collapse = ""), "</tr>")
  }
  header <- paste0(
    "<tr>",
    paste(sprintf("<th>%s</th>", html_escape(cols)), collapse = ""),
    "</tr>"
  )
  body <- paste(vapply(seq_len(nrow(top_clusters)),
                       function(i) fmt_row(top_clusters[i, , drop = FALSE]),
                       character(1L)), collapse = "\n")
  paste0('<table class="striped">', header, body, "</table>")
})()

# ----------------------------------------------------------------------
# Panel 8 — per-cluster conservation-profile mini-panels (if alignments
# dir provided).  One plot per top-N cluster, anchored at #<ltr_id>.
# ----------------------------------------------------------------------
per_cluster_html <- ""
if (!is.null(opt$alignments) && dir.exists(opt$alignments) &&
    nrow(top_clusters) > 0L) {

  base_counts_fn <- function(mat) {
    list(A = colSums(mat == "A"), C = colSums(mat == "C"),
         G = colSums(mat == "G"), T = colSums(mat == "T"),
         gap = colSums(mat == "-"), n = nrow(mat))
  }
  conservation_profile <- function(bc) {
    non_gap_n <- bc$A + bc$C + bc$G + bc$T
    top_count <- pmax(bc$A, bc$C, bc$G, bc$T)
    tf <- ifelse(non_gap_n > 0L, top_count / non_gap_n, 0)
    tf * (1 - bc$gap / bc$n)
  }
  slide_mean <- function(x, w) {
    n <- length(x); half <- w %/% 2L
    cs <- c(0, cumsum(x))
    i <- seq_len(n)
    lo <- pmax(1L, i - half); hi <- pmin(n, i + half)
    (cs[hi + 1L] - cs[lo]) / (hi - lo + 1L)
  }

  one_panel <- function(r) {
    id <- r$ltr_id
    p  <- file.path(opt$alignments, paste0(id, ".aln.fasta"))
    if (!file.exists(p)) return(sprintf(
      '<section id="%s"><h3>%s</h3><p>(alignment missing)</p></section>',
      html_escape(id), html_escape(id)))
    aln <- suppressWarnings(readDNAStringSet(p))
    mat <- toupper(as.matrix(aln))
    n5 <- r$n_5ltr; n3 <- r$n_3ltr
    if (nrow(mat) != n5 + n3) {
      return(sprintf('<section id="%s"><h3>%s</h3><p>(row mismatch)</p></section>',
                     html_escape(id), html_escape(id)))
    }
    img <- render_png({
      old <- par(mfrow = c(2L, 1L), mar = c(4, 4, 2, 1))
      on.exit(par(old), add = TRUE)
      cons5 <- slide_mean(conservation_profile(base_counts_fn(
        mat[seq_len(n5), , drop = FALSE])), 7L)
      plot(cons5, type = "l", ylim = c(0, 1), col = "#1f77b4",
           main = sprintf("%s  5'LTR subset (n=%d)", id, n5),
           xlab = "MSA column", ylab = "conservation")
      abline(v = r$annotated_5_col, col = "red",    lty = 2L)
      abline(v = r$corrected_5_col, col = "#2ca02c", lty = 1L)
      abline(v = r$annotated_3_col, col = "red",    lty = 2L)
      abline(v = r$corrected_3_col, col = "#2ca02c", lty = 1L)
      legend("bottomright", legend = c("annotated", "corrected"),
             col = c("red", "#2ca02c"), lty = c(2L, 1L), bty = "n", cex = 0.8)
      if (n3 > 0L) {
        cons3 <- slide_mean(conservation_profile(base_counts_fn(
          mat[(n5 + 1L):(n5 + n3), , drop = FALSE])), 7L)
        plot(cons3, type = "l", ylim = c(0, 1), col = "#ff7f0e",
             main = sprintf("%s  3'LTR subset (n=%d)", id, n3),
             xlab = "MSA column", ylab = "conservation")
        abline(v = r$annotated_5_col, col = "red",    lty = 2L)
        abline(v = r$corrected_5_col, col = "#2ca02c", lty = 1L)
        abline(v = r$annotated_3_col, col = "red",    lty = 2L)
        abline(v = r$corrected_3_col, col = "#2ca02c", lty = 1L)
      } else { plot.new() }
    }, width = 900, height = 520)
    sprintf('<section id="%s"><h3>%s  (%s)</h3>
<p>n_members=%d  shift_5=%d  shift_3=%d  flank=%d<br>
TG/CA: annotated=%.3f → corrected=%.3f  |  TSD: annotated=%.3f → corrected=%.3f</p>
%s</section>',
            html_escape(id), html_escape(id),
            html_escape(r$Final_Classification),
            r$n_members, r$shift_5, r$shift_3, r$flank,
            as.numeric(r$tgca_frac_annotated), as.numeric(r$tgca_frac_corrected),
            as.numeric(r$tsd_frac_annotated),  as.numeric(r$tsd_frac_corrected),
            img)
  }
  per_cluster_html <- paste(vapply(seq_len(nrow(top_clusters)),
                                   function(i) one_panel(top_clusters[i, , drop = FALSE]),
                                   character(1L)), collapse = "\n")
}

# ----------------------------------------------------------------------
# Assemble HTML
# ----------------------------------------------------------------------
now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
css <- '
body   { font-family: system-ui, sans-serif; max-width: 960px; margin: 2em auto; padding: 0 1em; color: #222; }
h1,h2  { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
h1     { font-size: 1.5em; border-bottom: 2px solid #333; padding-bottom: 0.25em; }
h2     { font-size: 1.15em; margin-top: 2em; border-bottom: 1px solid #ccc; padding-bottom: 0.15em; }
h3     { font-family: ui-monospace, monospace; font-size: 1.05em; margin-top: 1.5em; }
img.panel { display: block; margin: 1em auto; max-width: 100%; }
table.striped { border-collapse: collapse; width: 100%; font-family: ui-monospace, monospace; font-size: 0.85em; }
table.striped th, table.striped td { padding: 4px 8px; text-align: right; border-bottom: 1px solid #eee; }
table.striped th { background: #eee; text-align: left; }
table.striped td:nth-child(1), table.striped td:nth-child(2) { text-align: left; }
table.striped tr:nth-child(even) { background: #fafafa; }
td.na { color: #bbb; }
a { color: #1f6feb; text-decoration: none; }
a:hover { text-decoration: underline; }
p.hdr { font-family: ui-monospace, monospace; color: #555; font-size: 0.9em; }
'

html <- paste0(
  '<!doctype html>\n<html lang="en"><head>\n',
  '<meta charset="utf-8">\n',
  '<title>LTR library boundary report</title>\n',
  '<style>', css, '</style>\n',
  '</head><body>\n',
  '<h1>LTR library — boundary correction report</h1>\n',
  '<p class="hdr">Input: ',  html_escape(opt$qc),
  '<br>Clusters: ', n_total, ' total, ', nrow(msa), ' via MAFFT',
  '<br>Generated: ', html_escape(now), '</p>\n',
  '<h2>1. Cluster sizes</h2>\n', panel_sizes,
  '<h2>2. Boundary shift distribution (faceted by flank)</h2>\n', panel_shifts,
  '<h2>3. TG/CA annotated vs corrected</h2>\n', panel_tgca,
  '<h2>4. TSD annotated vs corrected</h2>\n', panel_tsd,
  '<h2>5. Per-lineage breakdown</h2>\n', panel_lineage,
  '<h2>6. Grow-cap activity</h2>\n', panel_growcap,
  '<h2>7. Top ', top_n_use, ' clusters of interest</h2>\n', top_table_html,
  if (nzchar(per_cluster_html))
    paste0('<h2>8. Per-cluster conservation profiles</h2>\n', per_cluster_html)
  else '',
  '</body></html>\n'
)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
writeLines(html, opt$out)
cat(sprintf("boundary report written: %s  (%d clusters, top-N=%d)\n",
            opt$out, n_total, top_n_use))
