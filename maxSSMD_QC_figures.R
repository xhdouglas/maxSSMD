## =============================================================================
##  maxSSMD_QC_figures.R
##
##  Reproduces Figures 1-4 of
##    Zhang XD. "Quality assessment for high-content screens via maximizing the
##    SSMD of a linear combination integrating multiple imaging responses."
##
##  Input : VZVconfirmatory.QC.txt (de-identified VZV siRNA confirmatory screen,
##          tab-delimited, one row per well, 12 plates x 384 wells)
##  Output: Figure1.tiff ... Figure4.tiff (300 dpi, LZW-compressed) in the
##          output folder, and a console summary of the key numbers quoted in
##          the manuscript.
##
##  How to run
##  ----------
##  * Windows Rgui:  File > Source R code... and pick this file, or type
##        source("C:/path/to/maxSSMD_QC_figures.R")
##    If DATA_FILE below is left as "", a window opens to let you select the
##    data file. Figures are written to the folder that holds the data file
##    unless OUT_DIR is set.
##  * Command line (any OS):
##        Rscript maxSSMD_QC_figures.R [data_file] [output_dir]
##
##  In Windows paths use forward slashes ("C:/Users/me/data.txt") or doubled
##  backslashes ("C:\\Users\\me\\data.txt"); single backslashes fail in R.
##
##  Requires only base R (>= 4.0); no add-on packages.
## =============================================================================

DATA_FILE <- ""   # e.g. "C:/Users/me/Documents/VZVconfirmatory.QC.txt"
OUT_DIR   <- ""   # e.g. "C:/Users/me/Documents/figures"; "" = data-file folder

## -----------------------------------------------------------------------------
## Locate input and output
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) DATA_FILE <- args[1]
if (length(args) >= 2) OUT_DIR   <- args[2]

if (!nzchar(DATA_FILE)) {
  cand <- c("VZVconfirmatory.QC.txt", "VZVconfirmatory.txt")
  hit  <- cand[file.exists(cand)]
  if (length(hit)) {
    DATA_FILE <- hit[1]
  } else if (interactive()) {
    message("Select the data file (VZVconfirmatory.QC.txt) in the window that opens.")
    DATA_FILE <- file.choose()
  } else {
    stop("No data file given and none found in ", getwd(),
         "\nUsage: Rscript maxSSMD_QC_figures.R data_file [output_dir]")
  }
}
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", DATA_FILE, "\nWorking directory: ", getwd())
in_file <- normalizePath(DATA_FILE, winslash = "/")
out_dir <- if (nzchar(OUT_DIR)) OUT_DIR else dirname(in_file)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(out_dir, winslash = "/")
cat("Data file    :", in_file, "\nOutput folder:", out_dir, "\n")

## -----------------------------------------------------------------------------
## 0. Settings
## -----------------------------------------------------------------------------
POS_LABEL <- "Luciferase"          # uninfected: positive control (Section 2.1)
NEG_LABEL <- "Luciferase + VZV"    # infected, non-targeting: negative control
SMP_LABEL <- "Sample"

## The seven retained responses (Section 2.2), in the order used in Figure 2
RESP <- c("Nuclei No of Objects 405nm",
          "gH No Obj 488nm",
          "gH Total 2 Area 488nm",
          "gH TotTot Int",
          "gH Total 2 Peak Int 488nm",
          "gH Total 2 Peak Int Bsln 488",
          "gH Med 2 Peak Int 488nm")
SHORT <- c("Nuclei.No.Obj.405", "gH.No.Obj.488", "gH.Tot.Area.488",
           "gH.TotTot.Int", "gH.Tot.Peak.488", "gH.Tot.Peak.Bsln",
           "gH.Med.Peak.488")
NUC   <- "Nuclei No of Objects 405nm"
TOT   <- "gH TotTot Int"
AREA  <- "gH Total 2 Area 488nm"

MAD_CUT   <- 5          # dispersion-screen cutoff on MAD-scaled z (Section 2.6)
THRESH    <- c(-3, -5, -7)
CONF      <- 0.95

## Colours
COL <- list(sample = "#E9B600", neg = "#1E8C45", pos = "#7B2D9A",
            tot = "#B5382C", area = "#1E8C45", totN = "#B8961A",
            areaN = "#7F8C8D", mx = "#1F4E79", mxcv = "#5DADE2",
            lda = "#A3441B", ldacv = "#F0A04B", neg3 = "#B5382C",
            pos3 = "#1F4E79", grey = "grey45")

## -----------------------------------------------------------------------------
## 1. Data and preprocessing (Sections 2.1-2.3)
## -----------------------------------------------------------------------------
## "UTF-8-BOM" also reads plain UTF-8 and files saved from Excel with a BOM
raw <- read.delim(in_file, check.names = FALSE, stringsAsFactors = FALSE,
                  fileEncoding = "UTF-8-BOM")

## Check the required columns; tolerate names whose spaces became dots or
## underscores (e.g. "gH.TotTot.Int" for "gH TotTot Int") by matching on
## letters and digits only.
need <- c("BARCODE", "XPOS", "YPOS", "WELL_USAGE", RESP)
key  <- function(s) tolower(gsub("[^A-Za-z0-9]", "", s))
for (nm in setdiff(need, names(raw))) {
  j <- which(key(names(raw)) == key(nm))
  if (length(j) == 1) names(raw)[j] <- nm
}
miss <- setdiff(need, names(raw))
if (length(miss))
  stop("Column(s) not found in the data file: ", paste(miss, collapse = "; "),
       "\nColumns present: ", paste(names(raw), collapse = "; "))

raw$plate <- as.integer(sub("^plate", "", raw$BARCODE))
raw       <- raw[order(raw$plate, raw$XPOS, raw$YPOS), ]   # row, then column
plates    <- sort(unique(raw$plate))
P         <- length(plates)

## Per-plate floor for non-positive values, then log2 (Section 2.3)
floor_log2 <- function(x, plate) {
  out <- numeric(length(x))
  for (p in unique(plate)) {
    i  <- plate == p
    xi <- x[i]
    xi[xi <= 0] <- min(xi[xi > 0]) / 2
    out[i] <- log2(xi)
  }
  out
}
L <- sapply(RESP, function(r) floor_log2(raw[[r]], raw$plate))
colnames(L) <- RESP

grp <- ifelse(raw$WELL_USAGE == POS_LABEL, "pos",
       ifelse(raw$WELL_USAGE == NEG_LABEL, "neg",
       ifelse(raw$WELL_USAGE == SMP_LABEL, "smp", "other")))

## -----------------------------------------------------------------------------
## 2. Plate-level SSMD and its noncentral-t confidence interval (Section 2.4)
## -----------------------------------------------------------------------------
ssmd <- function(x1, x2) (mean(x1) - mean(x2)) / sqrt(var(x1) + var(x2))

## Noncentral-t CDF by adaptive quadrature over the chi-square mixing
## distribution: P(T <= t) = E_U[ Phi(t * sqrt(U / nu) - ncp) ], U ~ chi2(nu).
## Used because |t| reaches 43, beyond the accurate range of pt(ncp = ).
pnct <- function(t, nu, ncp) {
  ## S = sqrt(U / nu) has density 2 * nu * s * dchisq(nu * s^2, nu); integrate
  ## piecewise between quantiles of S (and the kink s0 = ncp / t) so that the
  ## adaptive rule never misses the region carrying the mass.
  f  <- function(s) pnorm(t * s - ncp) * 2 * nu * s * dchisq(nu * s^2, nu)
  br <- sqrt(qchisq(c(1e-16, 1e-6, 1e-3, .02, .1, .3, .5, .7, .9, .98,
                      1 - 1e-3, 1 - 1e-6, 1 - 1e-16), nu) / nu)
  s0 <- ncp / t
  if (is.finite(s0) && s0 > br[1] && s0 < br[length(br)]) br <- sort(c(br, s0))
  sum(vapply(seq_len(length(br) - 1), function(i)
    integrate(f, br[i], br[i + 1], rel.tol = 1e-12, abs.tol = 1e-16,
              subdivisions = 500L)$value, numeric(1)))
}

## Satterthwaite-approximate noncentral-t interval for the SSMD of a fixed
## score under two independent groups with unequal variances: invert the
## Welch t against the noncentral t with Welch-Satterthwaite df.
ssmd_ci <- function(x1, x2, conf = CONF) {
  n1 <- length(x1); n2 <- length(x2)
  v1 <- var(x1);    v2 <- var(x2)
  se <- sqrt(v1 / n1 + v2 / n2)
  t0 <- (mean(x1) - mean(x2)) / se
  nu <- (v1 / n1 + v2 / n2)^2 /
        ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
  k  <- sqrt(v1 + v2) / se                 # ncp = SSMD * k
  a  <- (1 - conf) / 2
  g  <- function(ncp, p) pnct(t0, nu, ncp) - p
  w  <- 10 + 2 * abs(t0)
  lamL <- uniroot(g, c(t0 - w, t0 + w), p = 1 - a, tol = 1e-10)$root
  lamU <- uniroot(g, c(t0 - w, t0 + w), p = a,     tol = 1e-10)$root
  c(est = ssmd(x1, x2), lower = lamL / k, upper = lamU / k, df = nu)
}

## SSMD of a score s on every plate (optionally with CIs)
plate_ssmd <- function(s, plate_ids = plates, ci = FALSE) {
  out <- do.call(rbind, lapply(plate_ids, function(p) {
    x1 <- s[raw$plate == p & grp == "pos"]
    x2 <- s[raw$plate == p & grp == "neg"]
    if (ci) ssmd_ci(x1, x2) else c(est = ssmd(x1, x2))
  }))
  rownames(out) <- plate_ids
  out
}

## -----------------------------------------------------------------------------
## 3. Pooled estimation of m_L and Sigma_L*, dispersion screen, directions
##    (Sections 2.5, 2.6, 2.8)
## -----------------------------------------------------------------------------
## Per-plate pieces: within-group SSCP of plate-centred residuals, df, mean diff
plate_pieces <- lapply(plates, function(p) {
  X1 <- L[raw$plate == p & grp == "pos", , drop = FALSE]
  X2 <- L[raw$plate == p & grp == "neg", , drop = FALSE]
  list(A  = crossprod(scale(X1, scale = FALSE)), d1 = nrow(X1) - 1,
       B  = crossprod(scale(X2, scale = FALSE)), d2 = nrow(X2) - 1,
       md = colMeans(X1) - colMeans(X2))
})
names(plate_pieces) <- plates

## Dispersion screen: log trace of each plate's within-plate control covariance
log_trace <- sapply(plate_pieces, function(q)
  log(sum(diag(q$A / q$d1 + q$B / q$d2))))
mad_z <- function(v) (v - median(v)) / mad(v)

## Fit a direction on a set of plates
fit_direction <- function(use, method = c("maxSSMD", "LDA")) {
  method <- match.arg(method)
  z    <- mad_z(log_trace[as.character(use)])
  keep <- use[z <= MAD_CUT]          # exclude over-dispersed plates
  pc   <- plate_pieces[as.character(keep)]
  A  <- Reduce(`+`, lapply(pc, `[[`, "A")); d1 <- sum(sapply(pc, `[[`, "d1"))
  B  <- Reduce(`+`, lapply(pc, `[[`, "B")); d2 <- sum(sapply(pc, `[[`, "d2"))
  m  <- rowMeans(sapply(pc, `[[`, "md"))
  S  <- if (method == "maxSSMD") A / d1 + B / d2 else (A + B) / (d1 + d2)
  b  <- solve(S, m)
  b  <- b / sqrt(sum(b^2))                           # unit Euclidean norm
  ## sign convention: median plate-level SSMD of the score is negative
  if (median(plate_ssmd(L %*% b, use)[, "est"]) > 0) b <- -b
  list(beta = b, Sigma = S, m = m, kept = keep)
}

## Full-data fits
fit_mx  <- fit_direction(plates, "maxSSMD")
fit_lda <- fit_direction(plates, "LDA")
score_mx  <- drop(L %*% fit_mx$beta)
score_lda <- drop(L %*% fit_lda$beta)

## Leave-one-plate-out cross-validation (Section 2.7)
lopo <- function(method) {
  res <- t(sapply(plates, function(p) {
    f  <- fit_direction(setdiff(plates, p), method)
    s  <- drop(L %*% f$beta)
    x1 <- s[raw$plate == p & grp == "pos"]; x2 <- s[raw$plate == p & grp == "neg"]
    c(ssmd_ci(x1, x2), f$beta)
  }))
  rownames(res) <- plates
  res
}
cv_mx  <- lopo("maxSSMD")
cv_lda <- lopo("LDA")

## In-sample plate-level SSMD with CIs
ci_mx  <- plate_ssmd(score_mx,  ci = TRUE)
ci_lda <- plate_ssmd(score_lda, ci = TRUE)

## Single criteria (log2) incl. the two cell-number-normalized ratios
single <- list(
  "gH.TotTot.Int"            = L[, TOT],
  "gH.Tot.Area.488"          = L[, AREA],
  "gH.TotTot.Int / Nuclei"   = L[, TOT]  - L[, NUC],
  "gH.Tot.Area.488 / Nuclei" = L[, AREA] - L[, NUC])
ssmd_single <- lapply(single, function(s) plate_ssmd(s)[, "est"])
ssmd_resp   <- sapply(RESP,   function(r) plate_ssmd(L[, r])[, "est"])

## -----------------------------------------------------------------------------
## 4. Console summary (compare with Table 1 and the text)
## -----------------------------------------------------------------------------
summ <- function(v) c(median = median(v), le5 = sum(v <= -5), le7 = sum(v <= -7))
tab1 <- rbind(t(sapply(ssmd_single, summ)),
              "maxSSMD, 7 responses"           = summ(ci_mx[, "est"]),
              "maxSSMD, 7 responses [LOPO-CV]" = summ(cv_mx[, "est"]))
cat("\n=== Table 1 (median SSMD, #plates <= -5, #plates <= -7) ===\n")
print(round(tab1, 3))
cat("\nmaxSSMD coefficients (unit norm):\n");  print(round(fit_mx$beta, 3))
cat("\nLDA coefficients (unit norm):\n");      print(round(fit_lda$beta, 3))
tr <- sapply(plate_pieces, function(q) sum(diag(q$A / q$d1 + q$B / q$d2)))
cat(sprintf("\nDispersion screen: largest MAD z = %.2f; largest trace share = %.1f%%; plates excluded: %d\n",
            max(mad_z(log_trace)), 100 * max(tr / sum(tr)), P - length(fit_mx$kept)))
cat("\nPlate-level maxSSMD, in sample and LOPO-CV:\n")
print(round(cbind(plate = plates, in_sample = ci_mx[, "est"],
                  lopo = cv_mx[, "est"]), 2))
cat(sprintf("\nMedian 95%% CI half-width: maxSSMD %.2f, LDA %.2f\n",
            median((ci_mx[, "upper"] - ci_mx[, "lower"]) / 2),
            median((ci_lda[, "upper"] - ci_lda[, "lower"]) / 2)))
cat(sprintf("Max Welch df: maxSSMD %.1f, LDA %.1f\n",
            max(ci_mx[, "df"]), max(ci_lda[, "df"])))
cat(sprintf("Cosine(maxSSMD, LDA) on the pooled-SD scale: %.2f\n\n", {
  sdp <- sqrt(diag(fit_mx$Sigma)); a <- fit_mx$beta * sdp; b <- fit_lda$beta * sdp
  sum(a * b) / sqrt(sum(a^2) * sum(b^2)) }))

## -----------------------------------------------------------------------------
## 5. Plot helpers
## -----------------------------------------------------------------------------
## TIFF device: 300 dpi, LZW compression. Cairo where available (Windows,
## macOS, Linux builds of R normally have it); otherwise the native device.
open_tiff <- function(file, w, h) {
  dev_type <- if (isTRUE(capabilities("cairo"))) "cairo" else
              if (.Platform$OS.type == "windows") "windows" else
              if (Sys.info()[["sysname"]] == "Darwin") "quartz" else "Xlib"
  tiff(file.path(out_dir, file), width = w, height = h, units = "in",
       res = 300, compression = "lzw", type = dev_type)
}
alpha <- function(col, a) adjustcolor(col, alpha.f = a)
panel_label <- function(lab, adj_x = -0.08) {
  usr <- par("usr")
  text(grconvertX(adj_x, "npc"), grconvertY(1.06, "npc"), lab,
       xpd = NA, cex = 1.1, adj = c(0, 0))
}
hlines <- function() abline(h = THRESH, lty = 3, col = COL$grey, lwd = 1.2)

## Box summary with overlaid points (Figures 1f and 4b)
box_strip <- function(vals, cols, pchs, labs, ylim) {
  k <- length(vals)
  plot(NA, xlim = c(0.5, k + 0.5), ylim = ylim, xaxt = "n", bty = "l",
       xlab = "", ylab = "SSMD", las = 1)
  hlines()
  set.seed(1)
  for (j in seq_len(k)) {
    v <- vals[[j]]; q <- quantile(v, c(.25, .5, .75), type = 7)
    rect(j - 0.25, q[1], j + 0.25, q[3], col = alpha(cols[j], 0.25),
         border = alpha(cols[j], 0.4), lwd = 1.5)
    segments(j, min(v), j, max(v), lwd = 2)
    segments(j - 0.08, c(min(v), max(v)), j + 0.08, c(min(v), max(v)), lwd = 2)
    segments(j - 0.25, q[2], j + 0.25, q[2], lwd = 2.5)
    points(j + runif(length(v), -0.15, 0.15), v, pch = pchs[j], col = cols[j],
           bg = cols[j], cex = 1.1)
  }
  text(seq_len(k), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), labs,
       srt = 35, adj = c(1, 1), xpd = NA, cex = 0.85)
}

## Plate-well series (Figures 1a-d)
series_panel <- function(y, ylab, scale_lab = NULL, vr_lab, lab) {
  idx <- grp %in% c("smp", "neg", "pos")
  ## x position: plate p occupies [p - 0.5, p + 0.5), wells ordered row/column
  x <- numeric(nrow(raw))
  for (p in plates) {
    i <- which(raw$plate == p)
    x[i] <- p - 0.5 + (seq_along(i) - 0.5) / length(i)
  }
  plot(NA, xlim = c(0.5, P + 0.5), ylim = range(y[idx]), xaxt = "n",
       bty = "l", xlab = "plate number",
       ylab = if (is.null(scale_lab)) ylab else paste0(ylab, "\n", scale_lab),
       las = 1, mgp = c(2.6, 0.7, 0))
  axis(1, at = plates, labels = plates)
  abline(v = plates[-P] + 0.5, col = "grey80", lwd = 1.2)
  s <- grp == "smp"; n <- grp == "neg"; q <- grp == "pos"
  points(x[s], y[s], pch = 16, cex = 0.25, col = COL$sample)
  points(x[n], y[n], pch = 15, cex = 0.7,  col = COL$neg)
  points(x[q], y[q], pch = 23, cex = 0.9, col = COL$pos, bg = COL$pos)
  mtext(vr_lab, side = 3, adj = 1, line = 0.3, cex = 0.7, col = "grey30")
  panel_label(lab)
}

## Median over plates of the positive/negative variance ratio
var_ratio <- function(y) median(sapply(plates, function(p)
  var(y[raw$plate == p & grp == "pos"]) / var(y[raw$plate == p & grp == "neg"])))
fmt_vr <- function(v) paste("variance ratio pos/neg: median",
                            format(signif(v, 2), scientific = FALSE))

## -----------------------------------------------------------------------------
## 6. Figure 1: plate-well series and plate-level SSMD
## -----------------------------------------------------------------------------
open_tiff("Figure1.tiff", 8.1, 10)
layout(matrix(c(1, 2, 3, 4, 5, 5, 6, 7, 8, 8), ncol = 2, byrow = TRUE),
       heights = c(1, 1, 0.18, 1.55, 0.35))
par(mar = c(3.6, 4.6, 2, 0.8), cex = 0.75)
tot_raw  <- raw[[TOT]];  area_raw <- raw[[AREA]]
series_panel(tot_raw / 1e8,  "gH.TotTot.Int", "(\u00d710\u2078)",
             fmt_vr(var_ratio(tot_raw)), "(a)")
series_panel(L[, TOT], "gH.TotTot.Int", "(log\u2082)",
             fmt_vr(var_ratio(L[, TOT])), "(b)")
series_panel(area_raw / 1e6, "gH.Total.2.Area.488nm", "(\u00d710\u2076)",
             fmt_vr(var_ratio(area_raw)), "(c)")
series_panel(L[, AREA], "gH.Total.2.Area.488nm", "(log\u2082)",
             fmt_vr(var_ratio(L[, AREA])), "(d)")
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", c("sample wells", "negative control", "positive control"),
       pch = c(16, 15, 23), col = c(COL$sample, COL$neg, COL$pos), pt.bg = COL$pos,
       pt.cex = c(0.8, 1.3, 1.8), horiz = TRUE, bty = "n", cex = 1.1,
       x.intersp = 1.2, text.width = 0.22)

## (e) plate-level SSMD vs plate number
items <- c(ssmd_single,
           list("maxSSMD, 7 responses"           = ci_mx[, "est"],
                "maxSSMD, 7 responses [LOPO-CV]" = cv_mx[, "est"]))
icol <- c(COL$tot, COL$area, COL$totN, COL$areaN, COL$mx, COL$mxcv)
ipch <- c(16, 15, 17, 25, 23, 23)
ilty <- c(1, 1, 1, 1, 1, 2)
off  <- seq(-0.15, 0.15, length.out = length(items))
yl   <- c(-12.9, 0.4)
par(mar = c(4, 4.6, 2, 0.8))
plot(NA, xlim = c(0.6, P + 0.4), ylim = yl, xaxt = "n", bty = "l",
     xlab = "plate number", ylab = "SSMD", las = 1)
axis(1, at = plates, gap.axis = 0); hlines()
for (j in seq_along(items)) {
  lines(plates + off[j], items[[j]], col = icol[j], lwd = 1.8, lty = ilty[j])
  points(plates + off[j], items[[j]], col = icol[j], bg = icol[j],
         pch = ipch[j], cex = ifelse(ipch[j] == 23, 1.3, 1.1))
}
panel_label("(e)")
## (f) box summaries
par(mar = c(6.5, 4.6, 2, 0.8))
box_strip(items, icol, ipch, names(items), yl)
panel_label("(f)")
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", names(items), col = icol, pt.bg = icol, pch = ipch,
       lty = ilty, lwd = 1.8, ncol = 3, bty = "n", cex = 1, pt.cex = 1.2)
invisible(dev.off())

## -----------------------------------------------------------------------------
## 7. Figure 2: composition of the signature
## -----------------------------------------------------------------------------
open_tiff("Figure2.tiff", 8.5, 3.6)
par(mfrow = c(1, 2), mar = c(4, 9, 2.2, 1), cex = 0.8)
ord <- rev(seq_along(RESP))                      # first response on top
bp <- barplot(fit_mx$beta[ord], horiz = TRUE, names.arg = SHORT[ord], las = 1,
              col = COL$mx, border = NA, xlim = c(-0.85, 0.85), space = 0.6,
              xlab = "coefficient (unit norm)", axes = FALSE)
axis(1, at = seq(-0.8, 0.8, 0.2)); abline(v = 0)
panel_label("(a)", 0)
par(mar = c(4, 2, 2.2, 1))
med_alone <- apply(ssmd_resp, 2, median)
barplot(med_alone[ord], horiz = TRUE, names.arg = rep("", length(RESP)),
        col = "#5D6D7E", border = NA, xlim = c(-7.3, 0.6), space = 0.6,
        xlab = "median plate-level SSMD alone", axes = FALSE)
axis(1, at = seq(-6, 0, 2)); abline(v = 0)
abline(v = THRESH, lty = 3, col = COL$grey, lwd = 1.2)
panel_label("(b)", 0)
invisible(dev.off())

## -----------------------------------------------------------------------------
## 8. Figure 3: control-well distributions on Plate 4
## -----------------------------------------------------------------------------
PLATE3 <- 4
open_tiff("Figure3.tiff", 8.2, 4.4)
layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(0.13, 1))
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", c("negative control (n = 18)", "positive control (n = 6)"),
       col = c(COL$neg3, COL$pos3), lwd = 3, horiz = TRUE, bty = "n", cex = 1.1)
dens_panel <- function(s, title, xlab) {
  i1 <- raw$plate == PLATE3 & grp == "pos"; i2 <- raw$plate == PLATE3 & grp == "neg"
  x1 <- s[i1]; x2 <- s[i2]
  rg <- range(c(x1, x2)); pad <- 0.25 * diff(rg)
  d1 <- density(x1, from = rg[1] - pad, to = rg[2] + pad)
  d2 <- density(x2, from = rg[1] - pad, to = rg[2] + pad)
  ymax <- max(d1$y, d2$y)
  par(mar = c(4, 1, 4.2, 1))
  plot(NA, xlim = range(d1$x), ylim = c(-0.3 * ymax, ymax), yaxt = "n",
       bty = "n", xlab = xlab, ylab = "")
  polygon(c(d2$x, rev(d2$x)), c(d2$y, rep(0, length(d2$y))),
          col = alpha(COL$neg3, 0.25), border = NA)
  polygon(c(d1$x, rev(d1$x)), c(d1$y, rep(0, length(d1$y))),
          col = alpha(COL$pos3, 0.25), border = NA)
  lines(d2, col = COL$neg3, lwd = 2.5); lines(d1, col = COL$pos3, lwd = 2.5)
  segments(x2, -0.05 * ymax, x2, -0.12 * ymax, col = COL$neg3, lwd = 1.5)
  segments(x1, -0.17 * ymax, x1, -0.24 * ymax, col = COL$pos3, lwd = 1.5)
  mtext(sprintf("%s\nSSMD = %.2f", title, ssmd(x1, x2)), side = 3, adj = 0,
        line = 0.5, cex = 0.85)
}
dens_panel(L[, TOT], "gH.TotTot.Int", expression(log[2] ~ "value"))
dens_panel(L[, AREA] - L[, NUC], "gH.Tot.Area.488 \u2212 Nuclei.No.Obj.405",
           expression(log[2] ~ "ratio"))
dens_panel(score_mx, "maxSSMD signature score", "score")
invisible(dev.off())

## -----------------------------------------------------------------------------
## 9. Figure 4: maxSSMD vs LDA direction, with 95% confidence intervals
## -----------------------------------------------------------------------------
it4 <- list("maxSSMD, 7 responses" = ci_mx, "maxSSMD, 7 responses [LOPO-CV]" = cv_mx,
            "LDA, 7 responses" = ci_lda,    "LDA, 7 responses [LOPO-CV]" = cv_lda)
c4 <- c(COL$mx, COL$mxcv, COL$lda, COL$ldacv)
p4 <- c(23, 23, 21, 21); l4 <- c(1, 2, 1, 2)
o4 <- c(-0.27, -0.09, 0.09, 0.27)
yl4 <- range(c(-22.5, 0, sapply(it4, function(m) range(m[, c("lower", "upper")]))))
open_tiff("Figure4.tiff", 7.8, 6.4)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.2))
par(mar = c(4, 4.6, 2.2, 0.8), cex = 0.85)
plot(NA, xlim = c(0.6, P + 0.4), ylim = yl4, xaxt = "n", bty = "l",
     xlab = "plate number", ylab = "SSMD", las = 1)
axis(1, at = plates, gap.axis = 0); hlines()
for (j in seq_along(it4)) {
  m <- it4[[j]]; x <- plates + o4[j]
  segments(x, m[, "lower"], x, m[, "upper"], col = c4[j], lwd = 1.5)
  lines(x, m[, "est"], col = c4[j], lwd = 1.8, lty = l4[j])
  points(x, m[, "est"], pch = p4[j], col = c4[j], bg = c4[j], cex = ifelse(p4[j] == 23, 1.4, 1.2))
}
panel_label("(a)")
par(mar = c(7, 4.6, 2.2, 0.8))
box_strip(lapply(it4, function(m) m[, "est"]), c4, p4, names(it4), yl4)
panel_label("(b)")
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", names(it4), col = c4, pt.bg = c4, pch = p4, lty = l4, lwd = 1.8,
       ncol = 2, bty = "n", cex = 1, pt.cex = 1.3, text.width = 0.33)
invisible(dev.off())

cat("Figures written to", normalizePath(out_dir), "\n")
