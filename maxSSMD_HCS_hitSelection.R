#!/usr/bin/env Rscript
## =============================================================================
## maxSSMD_HCS_hitSelection.R
##
## Hit selection for high-content screens via maximizing the SSMD of a linear
## combination integrating multiple imaging responses
## (VZV siRNA confirmatory screen).
##
## Reproduces, from the screen data file, the analysis, Tables 1-2 and
## Figures 1-4 of the manuscript. Base R only (written and run under R 4.3);
## no add-on packages are required. Figures are written as TIFF (LZW, 300 dpi).
##
## Usage:
##   Rscript maxSSMD_HCS_hitSelection.R VZVconfirmatory.hitS.txt results 10 5 1.2
##
## Arguments (all optional after the data file):
##   1  data file (tab-delimited export of the screen)
##   2  output directory                               [results]
##   3  VIF cutoff for variable reduction              [10]
##   4  stricter VIF cutoff shown for comparison       [5]
##   5  fold-change cutoff of the dual-flashlight rule [1.2]
##
## Section numbers in the comments refer to the manuscript.
## =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "VZVconfirmatory.hitS.txt"
outdir  <- if (length(args) >= 2) args[2] else "results"
VIF_CUT <- if (length(args) >= 3) as.numeric(args[3]) else 10
VIF_STR <- if (length(args) >= 4) as.numeric(args[4]) else 5
FC      <- if (length(args) >= 5) as.numeric(args[5]) else 1.2
LFC     <- log2(FC)
DPI     <- 300
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(outdir, f)

## ----------------------------------------------------------------------------
## 0. Constants
## ----------------------------------------------------------------------------
RESP7 <- c("Nuclei.No.of.Objects.405nm", "gH.No.Obj.488nm", "gH.Total.2.Area.488nm",
           "gH.TotTot.Int", "gH.Total.2.Peak.Int.488nm",
           "gH.Total.2.Peak.Int.Bsln.488", "gH.Med.2.Peak.Int.488nm")
SHORT <- c("Nuclei count", "gH objects", "gH area", "gH total int.", "gH peak int.",
           "gH peak int. (bsln)", "gH median peak")
names(SHORT) <- RESP7
PRIMARY <- "gH.TotTot.Int"              # primary infection endpoint
AREA    <- "gH.Total.2.Area.488nm"      # infected area
NEG_LAB <- "Luciferase + VZV"           # negative control
TGT_LAB <- "Orf51.1 + VZV"              # target control
CUTS    <- c(-3, -2)                    # SSMD cutoffs (very strong, strong)

COL <- list(blue = "#2F7ED8", orange = "#EE6A40", gold = "#E3A008", grey = "#BDBDBD",
            dgrey = "#6E6E6E", green = "#1A9850", red = "#D62728", sample = "#C9C5BD")

## ----------------------------------------------------------------------------
## 1. Helper functions
## ----------------------------------------------------------------------------
unit <- function(b) b / sqrt(sum(b^2))

## Exact interval for a paired SSMD (mean/sd of n differences) from the
## noncentral t distribution on n - 1 degrees of freedom.
ci_nct <- function(ssmd, n, level = 0.95) {
  tt <- sqrt(n) * ssmd; a <- (1 - level) / 2
  f <- function(q) uniroot(function(ncp) pt(tt, n - 1, ncp) - q,
                           c(tt - 30, tt + 30), tol = 1e-10)$root / sqrt(n)
  c(f(1 - a), f(a))
}
## Exact interval for the maximal SSMD from the one-sample Hotelling statistic:
## n D^2 (n-k) / (k (n-1)) ~ noncentral F(k, n-k; n delta^2).
ci_hotelling <- function(D2, n, k, level = 0.95) {
  Fo <- n * D2 * (n - k) / (k * (n - 1)); a <- (1 - level) / 2
  f <- function(q) {
    g <- function(ncp) pf(Fo, k, n - k, ncp) - q
    if (g(0) < 0) return(0)
    sqrt(uniroot(g, c(0, 4 * n * D2 + 50), tol = 1e-10)$root / n)
  }
  c(f(1 - a), f(a))
}
vif <- function(S) diag(solve(cov2cor(S)))
## Stepwise removal of the response with the largest VIF until all are < cut.
step_vif <- function(S, cut) {
  keep <- colnames(S); steps <- list(); removed <- character(0)
  repeat {
    v <- vif(S[keep, keep, drop = FALSE]); steps[[length(steps) + 1]] <- v
    if (max(v) < cut) break
    removed <- c(removed, names(which.max(v))); keep <- setdiff(keep, names(which.max(v)))
  }
  list(keep = keep, steps = steps, removed = removed)
}
## Two-group (unpaired) SSMD and Mann-Whitney AUROC, P(target < negative).
ssmd2 <- function(zt, zn) (mean(zt) - mean(zn)) / sqrt(var(zt) + var(zn))
mw_below <- function(zt, zn) mean(outer(zt, zn, "<")) + 0.5 * mean(outer(zt, zn, "=="))
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

## ----------------------------------------------------------------------------
## 2. Data and preprocessing (Sections 2.1-2.2)
## ----------------------------------------------------------------------------
d <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
names(d) <- make.names(iconv(names(d), "latin1", "ASCII", sub = ""))
stopifnot(all(RESP7 %in% names(d)))
d$plate  <- as.integer(sub("plate", "", d$BARCODE))
d$source <- ifelse(d$SOBARCODE == "sourcePlate1", "A", "B")
PLATES   <- sort(unique(d$plate))
SRC_OF   <- sapply(PLATES, function(p) d$source[d$plate == p][1])   # source plate of each plate

raw <- as.matrix(d[, RESP7])
## non-positive values -> half the smallest positive value of that response on
## that plate; then log2
X <- raw
for (p in PLATES) for (j in RESP7) {
  i <- d$plate == p; v <- X[i, j]; v[v <= 0] <- min(v[v > 0]) / 2; X[i, j] <- v
}
X <- log2(X)

neg <- d$WELL_USAGE == NEG_LAB
tgt <- d$WELL_USAGE == TGT_LAB
smp <- d$WELL_USAGE == "Sample" & d$SAMPLE_ID != "0"     # identified siRNAs only

## ----------------------------------------------------------------------------
## 3. Paired differences and pooled covariance (Section 2.3)
## ----------------------------------------------------------------------------
Cbar <- t(sapply(PLATES, function(p) colMeans(X[neg & d$plate == p, , drop = FALSE])))
D    <- X - Cbar[match(d$plate, PLATES), ]              # paired difference of every well

ids <- unique(d$SAMPLE_ID[smp]); m <- length(ids)
o   <- order(match(d$SAMPLE_ID[smp], ids), d$plate[smp])
Ds  <- D[smp, ][o, ]                                    # 6 rows per siRNA
gi  <- match(d$SAMPLE_ID[smp][o], ids)                  # siRNA index of each row
stopifnot(all(tabulate(gi) == 6))
NREP  <- 6
src_i <- d$source[smp][o][seq(1, length(gi), NREP)]     # source plate of each siRNA
dbar  <- rowsum(Ds, gi) / NREP
resid <- Ds - dbar[gi, ]
nu    <- nrow(Ds) - m
S7    <- crossprod(resid) / nu                          # pooled replicate covariance

## dispersion screen: MAD-scaled z-score of the log trace of the within-siRNA covariance
ztrace <- function(R) {
  lt <- log(rowsum(rowSums(resid[, R, drop = FALSE]^2), gi) / (NREP - 1))
  as.vector((lt - median(lt)) / mad(lt))
}

## ----------------------------------------------------------------------------
## 4. Stepwise variable reduction by VIF (Sections 2.4, 3.2; Table 1)
## ----------------------------------------------------------------------------
sv   <- step_vif(S7, VIF_CUT)
sv5  <- step_vif(S7, VIF_STR)
R5   <- sv$keep                 # responses retained at the working cutoff
R4   <- sv5$keep                # responses retained at the stricter cutoff

nstep <- length(sv$steps)
tab1 <- data.frame(Response = RESP7, stringsAsFactors = FALSE)
for (s in seq_len(nstep)) tab1[[paste0("Step.", s)]] <- sv$steps[[s]][RESP7]
tab1[[sprintf("Step.%d.if.using.VIF.ge.%g", nstep + 1, VIF_STR)]] <-
  sv5$steps[[length(sv5$steps)]][RESP7]
tab1$Outcome <- ifelse(RESP7 %in% R5, "retained",
                       paste("removed at step", match(RESP7, sv$removed)))
write.csv(tab1, out("Table1_vif_steps.csv"), row.names = FALSE)

## checks of the reduction on the target-control covariance
sv_T   <- step_vif(cov(D[tgt, RESP7]), VIF_CUT)
sv_TA  <- step_vif(cov(D[tgt & d$source == "A", RESP7]), VIF_CUT)
sv_TB  <- step_vif(cov(D[tgt & d$source == "B", RESP7]), VIF_CUT)

## ----------------------------------------------------------------------------
## 5. Approach 1: direction fixed by the target control (Section 2.5)
## ----------------------------------------------------------------------------
target_dir <- function(R, wells) {
  DT <- D[wells, R, drop = FALSE]; n <- nrow(DT); k <- length(R)
  mu <- colMeans(DT); ST <- cov(DT)
  D2 <- sum(mu * solve(ST, mu))
  list(beta = -solve(ST, mu),                 # oriented: target control scores negative
       ssmd = -sqrt(D2), ci = -rev(ci_hotelling(D2, n, k)),
       bc   = -sqrt(max((n - k - 2) / (n - 1) * D2 - k / n, 0)),   # bias-corrected
       n = n, k = k, single = mu / sqrt(diag(ST)), shift = mu)
}
## coefficients on the scale of the pooled SDs, normalized to unit length
coef_sd <- function(beta) unit(beta * sqrt(diag(S7))[names(beta)])

## ----------------------------------------------------------------------------
## 6. Scores, SSMD-based hit selection and sign-flip FDR (Sections 2.6-2.9)
## ----------------------------------------------------------------------------
FLIP <- as.matrix(expand.grid(rep(list(c(1, -1)), NREP)))   # 2^6 sign patterns; row 1 = identity
stopifnot(all(FLIP[1, ] == 1))

## approach 2 for one set of differences (rows = replicates): project on the
## unit-length maxSSMD direction S^-1 * mean, signed by the change in the
## primary endpoint
ap2_scores <- function(Dm, Sinv) {
  mu <- colMeans(Dm); u <- unit(as.vector(Sinv %*% mu))
  as.vector(Dm %*% u) * sign(mu[PRIMARY])
}
## Null statistics for a fixed direction: flipping the sign of a difference
## flips the sign of its projection. Column 1 is the observed (unflipped) value.
stats_fixed <- function(P) {
  mn <- (P %*% t(FLIP)) / NREP
  sdv <- sqrt((rowSums(P^2) - NREP * mn^2) / (NREP - 1))
  list(mean = mn, ssmd = mn / sdv)
}
## Null statistics for approach 2: the direction is refitted to the flipped
## differences with the same pooled covariance.
stats_ap2 <- function(R, S) {
  Sinv <- solve(S); mn <- ss <- matrix(NA_real_, m, nrow(FLIP))
  for (i in seq_len(m)) {
    Dm <- Ds[gi == i, R, drop = FALSE]
    for (f in seq_len(nrow(FLIP))) {
      z <- ap2_scores(Dm * FLIP[f, ], Sinv); mn[i, f] <- mean(z); ss[i, f] <- mean(z) / sd(z)
    }
  }
  list(mean = mn, ssmd = ss)
}
## Hit counts and sign-flip FDR estimate m F0(t) / R(t), both restricted to
## values that also pass the fold-change cutoff.
select_hits <- function(st) {
  s <- st$ssmd[, 1]; mu <- st$mean[, 1]
  one <- function(t) {
    R  <- sum(s <= t & mu <= -LFC)
    F0 <- sum(st$ssmd <= t & st$mean <= -LFC) / length(st$ssmd)
    c(hits = R, fdr = if (R > 0) m * F0 / R else NA)
  }
  pass <- mu <= -LFC
  t20  <- sort(s[pass])[min(20, sum(pass))]
  list(ssmd = s, mean = mu, cut = sapply(CUTS, one), top20 = one(t20)["fdr"],
       down = sapply(CUTS, function(t) s <= t & mu <= -LFC),
       up   = sapply(-CUTS, function(t) s >= t & mu >= LFC),
       reach = sapply(CUTS, function(t) sum(s <= t)))
}

## Full analysis on one response set R (multivariate scores only)
analyse <- function(R) {
  S  <- S7[R, R]
  tA <- target_dir(R, tgt & d$source == "A"); tB <- target_dir(R, tgt & d$source == "B")
  tAll <- target_dir(R, tgt)
  uSrc <- list(A = unit(tA$beta), B = unit(tB$beta))
  ## approach 1: projection on the source-plate target direction (unit length, log2 scale)
  P1 <- matrix(NA_real_, m, NREP)
  for (i in seq_len(m)) P1[i, ] <- Ds[gi == i, R, drop = FALSE] %*% uSrc[[src_i[i]]]
  ## approach 2: projection on each siRNA's own direction
  Sinv <- solve(S); P2 <- matrix(NA_real_, m, NREP)
  for (i in seq_len(m)) P2[i, ] <- ap2_scores(Ds[gi == i, R, drop = FALSE], Sinv)
  ## target control under approach 2 (scored exactly as an siRNA)
  zT2 <- ap2_scores(D[tgt, R, drop = FALSE], Sinv)
  list(R = R, S = S, tA = tA, tB = tB, tAll = tAll, uSrc = uSrc, P1 = P1, P2 = P2,
       h1 = select_hits(stats_fixed(P1)), h2 = select_hits(stats_ap2(R, S)),
       ctrl2 = c(mean(zT2) / sd(zT2), ci_nct(mean(zT2) / sd(zT2), length(zT2))),
       zT2 = zT2)
}

cat("Running the analysis on the retained responses ...\n")
A5 <- analyse(R5)
cat("Running the analysis on all seven responses ...\n")
A7 <- analyse(RESP7)

## single-response comparators (special case of approach 1: coordinate direction)
single <- function(r) {
  P <- matrix(Ds[, r], m, NREP, byrow = TRUE); z <- D[tgt, r]
  list(P = P, h = select_hits(stats_fixed(P)),
       ctrl = c(mean(z) / sd(z), ci_nct(mean(z) / sd(z), length(z))))
}
SG <- list(single(PRIMARY), single(AREA)); names(SG) <- c(PRIMARY, AREA)

## ----------------------------------------------------------------------------
## 7. Table 2: SSMD-based hit selection
## ----------------------------------------------------------------------------
t2row <- function(label, ctrl, h) data.frame(
  Score = label,
  Control.SSMD = ctrl[1], Control.SSMD.lower = ctrl[2], Control.SSMD.upper = ctrl[3],
  Control.AUROC = pnorm(-ctrl[1]), Control.AUROC.lower = pnorm(-ctrl[3]),
  Control.AUROC.upper = pnorm(-ctrl[2]),
  Hits.SSMD.le.m3 = h$cut["hits", 1], FDR.SSMD.le.m3 = h$cut["fdr", 1],
  Hits.SSMD.le.m2 = h$cut["hits", 2], FDR.SSMD.le.m2 = h$cut["fdr", 2],
  FDR.top20 = unname(h$top20), stringsAsFactors = FALSE)
k5 <- length(R5)
tab2 <- rbind(
  t2row(PRIMARY, SG[[PRIMARY]]$ctrl, SG[[PRIMARY]]$h),
  t2row(AREA,    SG[[AREA]]$ctrl,    SG[[AREA]]$h),
  t2row("Target direction, approach 1", c(A5$tAll$ssmd, A5$tAll$ci), A5$h1),
  t2row("maxSSMD, approach 2", A5$ctrl2, A5$h2),
  t2row("Approach 1, seven responses", c(A7$tAll$ssmd, A7$tAll$ci), A7$h1),
  t2row("Approach 2, seven responses", A7$ctrl2, A7$h2))
rownames(tab2) <- NULL
write.csv(tab2, out("Table2_hit_selection.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## 8. Supporting result files
## ----------------------------------------------------------------------------
## per-siRNA scores, means, hit calls, SSMD-derived and empirical AUROC
per_sirna <- function(A, file) {
  sc <- list(SG[[PRIMARY]]$h, SG[[AREA]]$h, A$h1, A$h2)
  Ps <- list(SG[[PRIMARY]]$P, SG[[AREA]]$P, A$P1, A$P2)
  nm <- c("TotTot", "Area", "approach1", "approach2")
  df <- data.frame(SAMPLE_ID = ids, source.plate = src_i, stringsAsFactors = FALSE)
  for (j in seq_along(sc)) {
    df[[paste0(nm[j], ".SSMD")]]  <- sc[[j]]$ssmd
    df[[paste0(nm[j], ".mean")]]  <- sc[[j]]$mean
    df[[paste0(nm[j], ".AUROC.down")]] <- pnorm(-sc[[j]]$ssmd)
    df[[paste0(nm[j], ".empirical.AUROC.down")]] <- rowMeans(Ps[[j]] < 0)
    for (c in seq_along(CUTS)) {
      df[[sprintf("%s.down.hit.SSMD.le.%g", nm[j], CUTS[c])]] <- sc[[j]]$down[, c]
      df[[sprintf("%s.up.hit.SSMD.ge.%g",   nm[j], -CUTS[c])]] <- sc[[j]]$up[, c]
    }
  }
  write.csv(df, out(file), row.names = FALSE)
}
per_sirna(A5, "per_siRNA_signed_scores_reduced.csv")
per_sirna(A7, "per_siRNA_signed_scores_full7.csv")

## coefficients of the target directions (pooled-SD scale, unit length)
coef_tab <- function(A) data.frame(Response = A$R,
  source.plate.A = coef_sd(A$tA$beta), source.plate.B = coef_sd(A$tB$beta),
  all.plates = coef_sd(A$tAll$beta), row.names = NULL)
write.csv(coef_tab(A5), out("target_direction_coefficients_reduced.csv"), row.names = FALSE)
write.csv(coef_tab(A7), out("target_direction_coefficients_full7.csv"),  row.names = FALSE)
cosAB <- function(tA, tB) sum(coef_sd(tA$beta) * coef_sd(tB$beta))
A4 <- list(R = R4, tA = target_dir(R4, tgt & d$source == "A"),
           tB = target_dir(R4, tgt & d$source == "B"), tAll = target_dir(R4, tgt))
write.csv(coef_tab(A4), out("target_direction_coefficients_strict.csv"), row.names = FALSE)

## SSMD and AUROC of the target control by source plate (Section 3.3)
src_tab <- do.call(rbind, lapply(c("A", "B"), function(sp) {
  t5 <- if (sp == "A") A5$tA else A5$tB; t7 <- if (sp == "A") A7$tA else A7$tB
  w  <- tgt & d$source == sp
  zdir <- as.vector(D[w, R5] %*% A5$uSrc[[sp]])
  rows <- data.frame(source.plate = sp, score = "target direction", SSMD = t5$ssmd,
    lower = t5$ci[1], upper = t5$ci[2], bias.corrected = t5$bc,
    AUROC = pnorm(-t5$ssmd), AUROC.lower = pnorm(-t5$ci[2]), AUROC.upper = pnorm(-t5$ci[1]),
    AUROC.bias.corrected = pnorm(-t5$bc), empirical.AUROC = mean(zdir < 0))
  rows <- rbind(rows, data.frame(source.plate = sp, score = "target direction, seven responses",
    SSMD = t7$ssmd, lower = t7$ci[1], upper = t7$ci[2], bias.corrected = t7$bc,
    AUROC = pnorm(-t7$ssmd), AUROC.lower = pnorm(-t7$ci[2]), AUROC.upper = pnorm(-t7$ci[1]),
    AUROC.bias.corrected = pnorm(-t7$bc), empirical.AUROC = NA))
  for (r in RESP7) rows <- rbind(rows, data.frame(source.plate = sp, score = r,
    SSMD = t7$single[r], lower = NA, upper = NA, bias.corrected = NA,
    AUROC = pnorm(-t7$single[r]), AUROC.lower = NA, AUROC.upper = NA,
    AUROC.bias.corrected = NA, empirical.AUROC = mean(D[w, r] < 0)))
  rows }))
rownames(src_tab) <- NULL
write.csv(src_tab, out("target_control_source_plate_SSMD_AUROC.csv"), row.names = FALSE)

## SSMD and AUROC of the target control over all plates (Table 2, Section 3.4)
zdirAll <- as.vector(D[tgt, R5] %*% unit(A5$tAll$beta))
all_tab <- data.frame(
  score = c(PRIMARY, AREA, "target direction (approach 1)", "maxSSMD (approach 2)"),
  SSMD = tab2$Control.SSMD[1:4], lower = tab2$Control.SSMD.lower[1:4],
  upper = tab2$Control.SSMD.upper[1:4], AUROC = tab2$Control.AUROC[1:4],
  AUROC.lower = tab2$Control.AUROC.lower[1:4], AUROC.upper = tab2$Control.AUROC.upper[1:4],
  empirical.AUROC = c(mean(D[tgt, PRIMARY] < 0), mean(D[tgt, AREA] < 0),
                      mean(zdirAll < 0), mean(A5$zT2 < 0)))
write.csv(all_tab, out("target_control_SSMD_AUROC.csv"), row.names = FALSE)

## plate-level contrasts: target-control vs negative-control wells of a plate
## as two groups (Figure 3b)
plate_tab <- do.call(rbind, lapply(PLATES, function(p) {
  tt <- tgt & d$plate == p; nn <- neg & d$plate == p; u <- A5$uSrc[[SRC_OF[match(p, PLATES)]]]
  sc <- list(as.vector(X[tt, R5] %*% u), as.vector(X[nn, R5] %*% u),
             X[tt, PRIMARY], X[nn, PRIMARY], X[tt, AREA], X[nn, AREA])
  s <- c(ssmd2(sc[[1]], sc[[2]]), ssmd2(sc[[3]], sc[[4]]), ssmd2(sc[[5]], sc[[6]]))
  mw <- c(mw_below(sc[[1]], sc[[2]]), mw_below(sc[[3]], sc[[4]]), mw_below(sc[[5]], sc[[6]]))
  data.frame(plate = p, source.plate = SRC_OF[match(p, PLATES)],
    SSMD.direction = s[1], SSMD.TotTot = s[2], SSMD.Area = s[3],
    AUROC.direction = pnorm(-s[1]), AUROC.TotTot = pnorm(-s[2]), AUROC.Area = pnorm(-s[3]),
    MannWhitney.direction = mw[1], MannWhitney.TotTot = mw[2], MannWhitney.Area = mw[3]) }))
write.csv(plate_tab, out("target_control_plate_SSMD_reduced.csv"), row.names = FALSE)

## plate-to-plate variation of the negative control (Sections 2.8, 3.1)
anova_tab <- do.call(rbind, lapply(RESP7, function(r) do.call(rbind, lapply(c("original", "log2"), function(sc) {
  y <- if (sc == "original") raw[neg, r] else X[neg, r]; g <- factor(d$plate[neg])
  a <- anova(lm(y ~ g)); pm <- tapply(y, g, mean)
  data.frame(Response = r, scale = sc, F = a$`F value`[1], df1 = a$Df[1], df2 = a$Df[2],
    p = a$`Pr(>F)`[1], min.plate.mean = min(pm), max.plate.mean = max(pm),
    min.plate = PLATES[which.min(pm)], max.plate = PLATES[which.max(pm)],
    CV.plate.means = sd(pm) / mean(pm), SD.plate.means = sd(pm),
    SD.expected.within = sqrt(a$`Mean Sq`[2] / mean(table(g)))) }))))
write.csv(anova_tab, out("plate_variation_anova.csv"), row.names = FALSE)

## AUROC values and distribution-free lower bounds for the SSMD cutoffs (Section 2.10)
cut_tab <- data.frame(abs.SSMD = c(1, 2, 3))
cut_tab$AUROC.normal <- pnorm(cut_tab$abs.SSMD)
cut_tab$AUROC.lower.bound.symmetric.unimodal <- 1 - 2 / (9 * cut_tab$abs.SSMD^2)
cut_tab$AUROC.lower.bound.unimodal <- 1 - 4 / (9 * cut_tab$abs.SSMD^2)
write.csv(cut_tab, out("SSMD_cutoffs_AUROC.csv"), row.names = FALSE)

## empirical AUROC = 1 (all six differences below zero) versus SSMD-based selection (Section 3.4)
auroc_sum <- do.call(rbind, lapply(1:4, function(j) {
  h <- list(SG[[PRIMARY]]$h, SG[[AREA]]$h, A5$h1, A5$h2)[[j]]
  P <- list(SG[[PRIMARY]]$P, SG[[AREA]]$P, A5$P1, A5$P2)[[j]]
  alln <- rowSums(P < 0) == NREP & h$mean <= -LFC; nh <- alln & !h$down[, 1]
  data.frame(score = c(PRIMARY, AREA, "approach 1", "approach 2")[j],
    pass.FC.all.six.below.zero = sum(alln), hits.SSMD.le.m3 = sum(h$down[, 1]),
    hits.SSMD.le.m2 = sum(h$down[, 2]), all.six.below.zero.not.hit.m3 = sum(nh),
    median.AUROC.not.hit.m3 = median(pnorm(-h$ssmd[nh]))) }))
write.csv(auroc_sum, out("hit_selection_AUROC_summary.csv"), row.names = FALSE)

## approach 2 directions against the target direction (Section 3.3)
Sinv5 <- solve(A5$S)
cos2 <- sapply(seq_len(m), function(i) {
  b <- -as.vector(Sinv5 %*% dbar[i, R5]); names(b) <- R5        # oriented as the target direction
  sum(coef_sd(b) * coef_sd(list(A = A5$tA$beta, B = A5$tB$beta)[[src_i[i]]])) })

## key numbers quoted in the text
kn <- list()
add <- function(name, value) kn[[length(kn) + 1]] <<- data.frame(quantity = name, value = value)
add("siRNAs", m); add("sample wells", nrow(Ds)); add("pooled covariance df", nu)
add("max dispersion z, seven responses", max(ztrace(RESP7)))
add("max dispersion z, retained responses", max(ztrace(R5)))
add("condition number, seven responses", kappa(S7, exact = TRUE))
add("condition number, retained responses", kappa(S7[R5, R5], exact = TRUE))
add("retained responses", paste(R5, collapse = "; "))
add("removed responses (in order)", paste(sv$removed, collapse = "; "))
add("retained at stricter cutoff", paste(R4, collapse = "; "))
add("target-control covariance: removed (in order)", paste(sv_T$removed, collapse = "; "))
add("target-control covariance, source plate A: retained", paste(sv_TA$keep, collapse = "; "))
add("target-control covariance, source plate B: retained", paste(sv_TB$keep, collapse = "; "))
add("cosine A vs B, retained responses", cosAB(A5$tA, A5$tB))
add("cosine A vs B, seven responses", cosAB(A7$tA, A7$tB))
add("cosine A vs B, stricter cutoff", cosAB(A4$tA, A4$tB))
add("cosine A vs all plates", sum(coef_sd(A5$tA$beta) * coef_sd(A5$tAll$beta)))
add("cosine B vs all plates", sum(coef_sd(A5$tB$beta) * coef_sd(A5$tAll$beta)))
add("approach 2 vs target direction: min cosine", min(cos2))
add("approach 2 vs target direction: max cosine", max(cos2))
add("approach 2 vs target direction: median cosine", median(cos2))
add("approach 2 vs target direction: n cosine > 0.9", sum(cos2 > 0.9))
add("approach 2 vs target direction: n cosine < 0", sum(cos2 < 0))
add("hits shared by approaches 1 and 2, SSMD <= -3", sum(A5$h1$down[, 1] & A5$h2$down[, 1]))
add("hits shared by approaches 1 and 2, SSMD <= -2", sum(A5$h1$down[, 2] & A5$h2$down[, 2]))
add("up-regulated hits, approach 1", sum(A5$h1$up[, 1])); add("up-regulated hits, approach 2", sum(A5$h2$up[, 1]))
for (c in 1:2) {
  add(sprintf("approach 1, retained: reach SSMD <= %g", CUTS[c]), A5$h1$reach[c])
  add(sprintf("approach 1, retained: removed by fold-change cutoff at %g", CUTS[c]), A5$h1$reach[c] - A5$h1$cut["hits", c])
  add(sprintf("approach 1, seven: reach SSMD <= %g", CUTS[c]), A7$h1$reach[c])
  add(sprintf("approach 1, seven: removed by fold-change cutoff at %g", CUTS[c]), A7$h1$reach[c] - A7$h1$cut["hits", c])
  add(sprintf("approach 1: seven-response hits retained at %g", CUTS[c]), sum(A5$h1$down[, c] & A7$h1$down[, c]))
}
add("Spearman correlation of approach 1 SSMD, retained vs seven", cor(A5$h1$ssmd, A7$h1$ssmd, method = "spearman"))
proj_T <- function(A) mean(sapply(which(tgt), function(w) sum(D[w, A$R] * A$uSrc[[d$source[w]]])))
add("mean projected difference of target control, retained", proj_T(A5))
add("mean projected difference of target control, seven", proj_T(A7))
medsh <- apply(sweep(dbar[, R5], 2, sqrt(diag(A5$S)), "/"), 2, median)
for (r in R5) add(paste("median siRNA shift in pooled SDs:", r), medsh[r])
for (sp in c("A", "B")) for (r in R5)
  add(sprintf("target-control shift in pooled SDs, source plate %s: %s", sp, r),
      (if (sp == "A") A5$tA else A5$tB)$shift[r] / sqrt(S7[r, r]))
kn <- do.call(rbind, kn)
write.csv(kn, out("key_numbers.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## 9. Figures (TIFF)
## ----------------------------------------------------------------------------
open_tiff <- function(file, w, h)
  tiff(out(file), width = w, height = h, units = "in", res = DPI, compression = "lzw")
plab <- function(lab, fx, fy) {               # panel letter in device (NDC) coordinates
  op <- par(xpd = NA, new = FALSE)
  text(grconvertX(fx, "ndc", "user"), grconvertY(fy, "ndc", "user"), lab, font = 2, cex = 1.7)
  par(op)
}

## ---- Figure 1: plate-to-plate variation of the infected area ----------------
fig1 <- function() {
  open_tiff("Figure1_plate_variation.tiff", 11, 8.4)
  on.exit(dev.off())
  par(oma = c(0.3, 0, 1.2, 0))
  layout(matrix(1:3, 3, 1), heights = c(1, 0.13, 1))
  shown <- d$WELL_USAGE %in% c("Sample", NEG_LAB, TGT_LAB)
  xpos <- numeric(nrow(d))
  for (p in PLATES) {
    i <- which(shown & d$plate == p); oo <- order(d$XPOS[i], d$YPOS[i])
    xpos[i[oo]] <- (match(p, PLATES) - 1) + (seq_along(i) - 0.5) / length(i)
  }
  draw <- function(y, ylim, ylab, lab) {
    par(mar = c(1.2, 5.6, 3.2, 0.8), mgp = c(3.4, 0.8, 0), las = 1)
    plot(NA, xlim = c(0, length(PLATES)), ylim = ylim, xaxs = "i", yaxs = "i", xaxt = "n",
         bty = "l", xlab = "", ylab = ylab, cex.axis = 1.5, cex.lab = 1.6, col.axis = "grey25")
    abline(v = seq_len(length(PLATES) - 1), col = "grey85")
    w <- shown & d$WELL_USAGE == "Sample"; points(xpos[w], y[w], pch = 16, cex = 0.45, col = COL$sample)
    w <- neg; points(xpos[w], y[w], pch = 16, cex = 0.9, col = COL$blue)
    w <- tgt; points(xpos[w], y[w], pch = 16, cex = 0.9, col = COL$orange)
    for (p in PLATES) {
      k <- match(p, PLATES); segments(k - 1, mean(y[neg & d$plate == p]), k, lwd = 2.4, col = COL$blue)
      mtext(paste0(p, SRC_OF[k]), side = 3, at = k - 0.5, line = 0.5, cex = 1.15, col = "grey30")
    }
    plab(lab, 0.035, grconvertY(1, "npc", "ndc") + 0.04)
  }
  draw(raw[, AREA] / 1e6, c(0, 1.45), bquote(.(AREA) ~ "(" * "" %*% 10^6 * ")"), "a")
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", cex = 1.5, text.width = c(0.085, 0.16, 0.14, 0.26),
         legend = c("Sample", "Negative control", "Target control", "plate mean, negative control"),
         pch = c(16, 16, 16, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2.4), pt.cex = 1.6,
         col = c(COL$sample, COL$blue, COL$orange, COL$blue))
  draw(X[, AREA], c(15.5, 20.65), bquote(log[2] ~ .(AREA)), "b")
}

## ---- Figure 2: collinearity and stepwise VIF reduction ----------------------
fig2 <- function() {
  open_tiff("Figure2_collinearity_VIF.tiff", 10.5, 4.6)
  on.exit(dev.off())
  layout(matrix(1:3, 1, 3), widths = c(1, 0.17, 1.05))
  Rm <- cov2cor(S7); k <- length(RESP7)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
  ## (a) correlation matrix
  par(mar = c(9.5, 10.5, 3, 0.5), las = 1)
  image(1:k, 1:k, t(Rm[k:1, ]), zlim = c(-1, 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(2, at = k:1, labels = SHORT, cex.axis = 1.25, lwd = 0, lwd.ticks = 1)
  axis(1, at = 1:k, labels = FALSE, lwd = 0, lwd.ticks = 1)
  text(1:k + 0.15, par("usr")[3] - 0.35, SHORT, srt = 40, adj = 1, xpd = NA, cex = 1.25)
  for (i in 1:k) for (j in 1:k)
    text(j, k + 1 - i, sprintf("%.2f", Rm[i, j]), cex = 1.05, col = if (abs(Rm[i, j]) > 0.55) "white" else "black")
  plab("a", 0.03, 0.95)
  ## colour bar
  par(mar = c(9.5, 0.6, 3, 5.2))
  image(1, seq(-1, 1, length.out = 201), matrix(seq(-1, 1, length.out = 201), 1), col = pal,
        axes = FALSE, xlab = "", ylab = ""); box()
  axis(4, at = seq(-1, 1, 0.25), labels = sprintf("%.2f", seq(-1, 1, 0.25)), cex.axis = 1.2)
  mtext("Pooled replicate correlation", side = 4, line = 3.9, las = 0, cex = 0.95)
  ## (b) VIF by step
  par(mar = c(9.5, 6, 3, 1), las = 1, mgp = c(3.6, 0.8, 0))
  allv <- unlist(sv$steps); yl <- c(1, 10^ceiling(log10(max(allv))))
  plot(NA, xlim = c(0.5, k + 0.5), ylim = yl, log = "y", xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "VIF from the pooled replicate covariance", cex.lab = 1.3)
  ticks <- 10^(0:log10(yl[2])); abline(h = ticks, col = "grey85")
  axis(2, at = ticks, labels = parse(text = paste0("10^", log10(ticks))), cex.axis = 1.3)
  axis(1, at = 1:k, labels = FALSE)
  text(1:k + 0.15, 10^(par("usr")[3] - 0.07 * diff(par("usr")[3:4])), SHORT, srt = 40, adj = 1, xpd = NA, cex = 1.25)
  abline(h = VIF_CUT, lty = 2, lwd = 1.5); abline(h = VIF_STR, lty = 3, lwd = 2)
  text(4.9, VIF_CUT * 1.35, sprintf("VIF = %g", VIF_CUT), cex = 1.1)
  text(4.9, VIF_STR / 1.35, sprintf("VIF = %g", VIF_STR), cex = 1.1)
  pch <- c(16, 15, 17, 18); cols <- c(COL$dgrey, COL$orange, COL$gold, COL$blue)
  off <- seq(-0.2, 0.2, length.out = max(nstep, 2))[seq_len(nstep)]
  for (s in seq_len(nstep)) {
    v <- sv$steps[[s]]; points(match(names(v), RESP7) + off[s], v, pch = pch[s], col = cols[s], cex = 2)
  }
  for (s in seq_along(sv$removed)) {
    r <- sv$removed[s]
    text(match(r, RESP7) + off[s] + 0.2, sv$steps[[s]][r] * 1.25, sprintf("removed at step %d", s), adj = 0, cex = 1.05, xpd = NA)
  }
  legend("topleft", bty = "n", pch = pch[seq_len(nstep)], col = cols[seq_len(nstep)], pt.cex = 1.7, cex = 1.15,
         legend = sprintf("step %d (%d responses)", seq_len(nstep), sapply(sv$steps, length)), inset = c(0.02, 0.0))
  plab("b", 0.575, 0.95)
}

## ---- Figure 3: target direction and its coefficients ------------------------
fig3 <- function() {
  open_tiff("Figure3_target_direction.tiff", 11, 8.8)
  on.exit(dev.off())
  layout(rbind(c(1, 1, 1, 2, 2, 2), c(3, 3, 4, 4, 5, 5)), heights = c(1.12, 1))
  cA <- COL$blue; cB <- COL$orange
  ## (a) single-response SSMD of the target control and the maximal SSMD
  k <- length(R5); ypos <- c(rev(seq_len(k)) + 1, 1)          # rows top to bottom; last = direction
  par(mar = c(14, 12.5, 3.5, 2), las = 1, mgp = c(3, 0.8, 0))
  xr <- range(c(A5$tA$single, A5$tB$single, A5$tA$ci, A5$tB$ci, A7$tA$ssmd, A7$tB$ssmd))
  xl <- c(floor(xr[1]), ceiling(xr[2] * 2) / 2 + 0.3)
  plot(NA, xlim = xl, ylim = c(0.45, k + 1.6), yaxt = "n", bty = "n", xlab = "", ylab = "", cex.axis = 1.3, xaxt = "n")
  abline(v = pretty(xl), col = "grey88"); abline(v = 0, lwd = 1.3)
  axis(1, at = pretty(xl), cex.axis = 1.3)
  axis(2, at = ypos, labels = c(SHORT[R5], "target direction\n(95% CI)"), cex.axis = 1.3, lwd = 1)
  for (j in seq_len(k)) {
    rect(0, ypos[j], A5$tA$single[R5[j]], ypos[j] + 0.36, col = cA, border = NA)
    rect(0, ypos[j] - 0.4, A5$tB$single[R5[j]], ypos[j] - 0.04, col = cB, border = NA)
  }
  dia <- function(t5, t7, y, col) {
    segments(t5$ci[1], y, t5$ci[2], y, col = col, lwd = 1.8)
    points(t7$ssmd, y, pch = 23, bg = "white", col = col, cex = 2.1, lwd = 1.4)
    points(t5$ssmd, y, pch = 23, bg = col, col = col, cex = 2.4)
  }
  dia(A5$tA, A7$tA, 1.2, cA); dia(A5$tB, A7$tB, 0.8, cB)
  mtext("SSMD, target control vs\nnegative control (paired)", side = 1, line = 4.6, cex = 0.95)
  legend(grconvertX(0.27, "npc"), grconvertY(-0.42, "npc"), xpd = NA, bty = "n", cex = 1.3,
         legend = c("source plate A", "source plate B", "seven-response direction"),
         pch = c(15, 15, 23), col = c(cA, cB, "grey30"), pt.bg = c(NA, NA, "white"), pt.cex = c(1.8, 1.8, 1.7))
  plab("a", 0.04, 0.975)
  ## (b) plate-level SSMD
  par(mar = c(14, 8, 3.5, 1.5), las = 1, mgp = c(3, 0.8, 0))
  yy <- as.matrix(plate_tab[, c("SSMD.TotTot", "SSMD.Area", "SSMD.direction")])
  yl <- c(floor(min(yy) * 2) / 2, max(0.3, ceiling(max(yy) * 2) / 2))
  plot(NA, xlim = c(0.5, length(PLATES) + 0.5), ylim = yl, xaxt = "n", bty = "n", xlab = "", ylab = "", cex.axis = 1.3)
  abline(h = pretty(yl), col = "grey88"); abline(h = 0, lwd = 1.3)
  axis(1, at = seq_along(PLATES), labels = PLATES, cex.axis = 1.2)
  points(seq_along(PLATES) - 0.2, yy[, 1], pch = 16, col = cA, cex = 2)
  points(seq_along(PLATES),       yy[, 2], pch = 15, col = cB, cex = 1.8)
  points(seq_along(PLATES) + 0.2, yy[, 3], pch = 18, col = COL$gold, cex = 2.6)
  mtext("Plate", side = 1, line = 3, cex = 0.95)
  mtext("Plate-level SSMD, target control\nvs negative control", side = 2, line = 4, las = 0, cex = 0.95)
  legend(grconvertX(0.3, "npc"), grconvertY(-0.42, "npc"), xpd = NA, bty = "n", cex = 1.3,
         legend = c(PRIMARY, AREA, "target direction"), pch = c(16, 15, 18),
         col = c(cA, cB, COL$gold), pt.cex = c(1.8, 1.7, 2.3))
  plab("b", 0.545, 0.975)
  ## (c-e) coefficients on three response sets
  coefpanel <- function(A, title, lab, fx) {
    R <- A$R; k <- length(R); a <- coef_sd(A$tA$beta); b <- coef_sd(A$tB$beta); yp <- rev(seq_len(k))
    par(mar = c(9.5, 11, 5, 1.5), las = 1, mgp = c(3, 0.8, 0))
    xl <- c(-0.7, 1.05)
    plot(NA, xlim = xl, ylim = c(0.4, k + 0.6), yaxt = "n", bty = "n", xlab = "", ylab = "", cex.axis = 1.3, xaxt = "n")
    abline(v = seq(-0.5, 1, 0.5), col = "grey88"); abline(v = 0, lwd = 1.3)
    axis(1, at = seq(-0.5, 1, 0.5), labels = sprintf("%.1f", seq(-0.5, 1, 0.5)), cex.axis = 1.3)
    axis(2, at = yp, labels = SHORT[R], cex.axis = 1.3)
    for (j in seq_len(k)) {
      rect(0, yp[j], a[j], yp[j] + 0.36, col = cA, border = NA)
      rect(0, yp[j] - 0.4, b[j], yp[j] - 0.04, col = cB, border = NA)
    }
    mtext(title, side = 3, line = 0.8, cex = 0.95)
    mtext("Coefficient (unit length,\npooled-SD scale)", side = 1, line = 4.6, cex = 0.95)
    mtext(sprintf("cosine, A vs B = %.3f", cosAB(A$tA, A$tB)), side = 1, line = 7.3, cex = 0.95, font = 3)
    plab(lab, fx, 0.43)
  }
  coefpanel(A5, sprintf("%d responses retained\n(VIF cutoff %g)", length(R5), VIF_CUT), "c", 0.012)
  coefpanel(A7, sprintf("All %d responses\n(no reduction)", length(RESP7)), "d", 0.347)
  coefpanel(A4, sprintf("%d responses retained\n(VIF cutoff %g)", length(R4), VIF_STR), "e", 0.68)
}

## ---- Figure 4: dual-flashlight plots ---------------------------------------
fig4 <- function() {
  open_tiff("Figure4_dual_flashlight.tiff", 10.5, 9)
  on.exit(dev.off())
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  fcl <- gsub("\\.?0+$", "", format(FC))
  one <- function(h, title, xlab, ylab, lab, fx, fy, xline = 3.2) {
    par(mar = c(7, 6.2, 4.2, 1.2), las = 1, mgp = c(3, 0.8, 0))
    dn <- h$down[, 1]; up <- h$up[, 1]
    yl <- range(h$ssmd) + c(-0.6, 0.6) * 1; yl[2] <- yl[2] + 0.34 * diff(yl)
    plot(NA, xlim = range(h$mean) + c(-0.06, 0.06) * diff(range(h$mean)), ylim = yl, bty = "l",
         xlab = "", ylab = "", cex.axis = 1.25)
    grid(col = "grey90", lty = 1); abline(h = 0, v = 0, col = "grey50")
    abline(h = c(-3, 3), lty = 2, lwd = 1.3); abline(v = c(-LFC, LFC), lty = 3, lwd = 2)
    points(h$mean[!dn & !up], h$ssmd[!dn & !up], pch = 16, col = COL$grey, cex = 1)
    points(h$mean[up], h$ssmd[up], pch = 16, col = COL$red, cex = 1.5)
    points(h$mean[dn], h$ssmd[dn], pch = 16, col = COL$green, cex = 1.5)
    mtext(title, side = 3, adj = 0, line = 0.6, cex = 1.15)
    mtext(xlab, side = 1, line = xline, cex = 0.95)
    mtext(ylab, side = 2, line = 3.4, las = 0, cex = 0.95)
    lg <- c(expression("siRNA"),
      bquote("up-regulated, SSMD" >= 3 ~ "and mean" >= log[2] * .(fcl) ~ "(" * .(sum(up)) * ")"),
      bquote("down-regulated, SSMD" <= -3 ~ "and mean" <= -log[2] * .(fcl) ~ "(" * .(sum(dn)) * ")"))
    legend("topleft", legend = as.expression(lg), pch = 16, col = c(COL$grey, COL$red, COL$green),
           pt.cex = c(1, 1.4, 1.4), bty = "n", cex = 1.05, inset = c(-0.01, -0.01))
    plab(lab, fx, fy)
  }
  xl <- expression("Mean paired " * log[2] * " difference")
  one(SG[[PRIMARY]]$h, PRIMARY, xl, "SSMD", "a", 0.022, 0.975)
  one(SG[[AREA]]$h,    AREA,    xl, "SSMD", "b", 0.522, 0.975)
  one(A5$h1, "Approach 1: target direction",
      expression(atop("Mean paired difference along the", "target unit direction (" * log[2] * " units)")),
      "SSMD along target direction", "c", 0.022, 0.475, xline = 5)
  one(A5$h2, "Approach 2: siRNA-specific maxSSMD",
      expression(atop("Mean paired difference along the siRNA-specific", "unit direction (" * log[2] * " units, signed)")),
      "SSMD along siRNA-specific direction\n(signed)", "d", 0.522, 0.475, xline = 5)
}

cat("Drawing figures ...\n")
fig1(); fig2(); fig3(); fig4()

## ----------------------------------------------------------------------------
## 10. Console summary
## ----------------------------------------------------------------------------
cat("\n=== Table 1. Stepwise variable reduction by VIF ===\n")
print(cbind(tab1[1], signif(tab1[sapply(tab1, is.numeric)], 3), tab1["Outcome"]), row.names = FALSE)
cat("\n=== Target control by source plate (retained responses) ===\n")
for (sp in c("A", "B")) { t <- if (sp == "A") A5$tA else A5$tB
  cat(sprintf("source plate %s: SSMD %.2f (%.2f to %.2f; bias-corrected %.2f); AUROC %.3f (%.3f to %.3f; bias-corrected %.3f)\n",
      sp, t$ssmd, t$ci[1], t$ci[2], t$bc, pnorm(-t$ssmd), pnorm(-t$ci[2]), pnorm(-t$ci[1]), pnorm(-t$bc))) }
cat(sprintf("cosine A vs B: %.3f (retained), %.3f (seven), %.3f (stricter cutoff)\n",
            cosAB(A5$tA, A5$tB), cosAB(A7$tA, A7$tB), cosAB(A4$tA, A4$tB)))
cat("\n=== Table 2. SSMD-based hit selection ===\n")
print(data.frame(Score = tab2$Score,
  Control.SSMD = sprintf("%.3f (%.3f, %.3f)", tab2$Control.SSMD, tab2$Control.SSMD.lower, tab2$Control.SSMD.upper),
  Control.AUROC = sprintf("%.3f (%.3f, %.3f)", tab2$Control.AUROC, tab2$Control.AUROC.lower, tab2$Control.AUROC.upper),
  Hits.m3 = tab2$Hits.SSMD.le.m3, FDR.m3 = fmt_pct(tab2$FDR.SSMD.le.m3),
  Hits.m2 = tab2$Hits.SSMD.le.m2, FDR.m2 = fmt_pct(tab2$FDR.SSMD.le.m2),
  FDR.top20 = fmt_pct(tab2$FDR.top20)), row.names = FALSE)
cat("\nOutput written to:", normalizePath(outdir), "\n")
