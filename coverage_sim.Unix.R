# Coverage and bias of estimators / intervals for the maximal SSMD
#   theta = sqrt( m' Sigma_*^{-1} m ),  Sigma_* = Sigma1 + Sigma2 - Sigma12 - Sigma12'
# Situations: (a) independent, equal covariance; (b) independent, unequal covariance;
#             (c) paired / correlated.
# Intervals: exact noncentral-F (a, c), Nel-van der Merwe F (b), bootstrap percentile and BCa (all).
suppressMessages(library(MASS))
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 1000
B     <- if (length(args) >= 2) as.integer(args[2]) else 500
alpha <- 0.05
theta0 <- 2                      # true maximal SSMD in every setting

ar1 <- function(k, rho) rho^abs(outer(1:k, 1:k, "-"))
msqrt <- function(A) { e <- eigen(A, symmetric = TRUE); e$vectors %*% diag(sqrt(e$values), nrow(A)) %*% t(e$vectors) }

make_setting <- function(sit, k) {
  S1 <- ar1(k, 0.5)
  if (sit == "a") { S2 <- S1; S12 <- matrix(0, k, k) }
  if (sit == "b") { S2 <- 3 * ar1(k, 0.2); S12 <- matrix(0, k, k) }
  if (sit == "c") { S2 <- 2 * ar1(k, 0.3); S12 <- 0.6 * msqrt(S1) %*% msqrt(S2) }
  Sst <- S1 + S2 - S12 - t(S12)
  m0 <- rep(1, k)
  m <- m0 * theta0 / sqrt(drop(t(m0) %*% solve(Sst, m0)))   # scale so that theta = theta0
  list(sit = sit, k = k, S1 = S1, S2 = S2, S12 = S12, m = m,
       J = rbind(cbind(S1, S12), cbind(t(S12), S2)))
}

# ---- plug-in maximal SSMD for a data set --------------------------------------
plugin <- function(sit, X1, X2) {
  if (sit == "c") { D <- X1 - X2; d <- colMeans(D); S <- cov(D) }
  else {
    d <- colMeans(X1) - colMeans(X2)
    if (sit == "a") { n1 <- nrow(X1); n2 <- nrow(X2)
      S <- 2 * ((n1 - 1) * cov(X1) + (n2 - 1) * cov(X2)) / (n1 + n2 - 2) }
    else S <- cov(X1) + cov(X2)
  }
  q <- tryCatch(drop(crossprod(d, solve(S, d))), error = function(e) NA_real_)
  sqrt(max(q, 0))
}

# ---- invert noncentral F for lambda -------------------------------------------
ncp_ci <- function(Fobs, df1, df2) {
  f <- function(l, p) pf(Fobs, df1, df2, ncp = l) - p
  hi_lim <- 10; while (pf(Fobs, df1, df2, ncp = hi_lim) > alpha / 2) hi_lim <- hi_lim * 2
  lo <- if (pf(Fobs, df1, df2, ncp = 0) <= 1 - alpha / 2) 0 else uniroot(f, c(0, hi_lim), p = 1 - alpha / 2)$root
  hi <- if (pf(Fobs, df1, df2, ncp = 0) <= alpha / 2) 0 else uniroot(f, c(0, hi_lim), p = alpha / 2)$root
  c(lo, hi)
}

analytic <- function(sit, X1, X2) {
  k <- ncol(X1); n1 <- nrow(X1); n2 <- nrow(X2)
  if (sit == "a") {
    nu <- n1 + n2 - 2
    Sp <- ((n1 - 1) * cov(X1) + (n2 - 1) * cov(X2)) / nu
    d <- colMeans(X1) - colMeans(X2); D2 <- drop(crossprod(d, solve(Sp, d)))
    D2u <- (nu - k - 1) / nu * D2 - k * (n1 + n2) / (n1 * n2)
    bc <- sqrt(max(D2u, 0) / 2)
    T2 <- n1 * n2 / (n1 + n2) * D2; Fobs <- (nu - k + 1) / (k * nu) * T2
    lam <- ncp_ci(Fobs, k, nu - k + 1)
    ci <- sqrt((n1 + n2) / (2 * n1 * n2) * lam)
  } else if (sit == "c") {
    D <- X1 - X2; nu <- n1 - 1
    d <- colMeans(D); D2 <- drop(crossprod(d, solve(cov(D), d)))
    D2u <- (nu - k - 1) / nu * D2 - k / n1
    bc <- sqrt(max(D2u, 0))
    T2 <- n1 * D2; Fobs <- (n1 - k) / (k * (n1 - 1)) * T2
    lam <- ncp_ci(Fobs, k, n1 - k)
    ci <- sqrt(lam / n1)
  } else {                                   # (b) Nel-van der Merwe
    A1 <- cov(X1) / n1; A2 <- cov(X2) / n2; Se <- A1 + A2
    d <- colMeans(X1) - colMeans(X2)
    T2 <- drop(crossprod(d, solve(Se, d)))
    tr <- function(M) sum(diag(M))
    nuNM <- (tr(Se %*% Se) + tr(Se)^2) /
      ((tr(A1 %*% A1) + tr(A1)^2) / (n1 - 1) + (tr(A2 %*% A2) + tr(A2)^2) / (n2 - 1))
    Fobs <- (nuNM - k + 1) / (k * nuNM) * T2
    lam <- ncp_ci(Fobs, k, nuNM - k + 1)
    ci <- sqrt(lam / n1)                        # exact map when n1 = n2
    bc <- NA
  }
  c(bc = bc, lo = ci[1], hi = ci[2])
}

boot_ci <- function(sit, X1, X2, est) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  bs <- numeric(B)
  for (b in 1:B) {
    if (sit == "c") { i <- sample.int(n1, n1, TRUE); bs[b] <- plugin(sit, X1[i, , drop = FALSE], X2[i, , drop = FALSE]) }
    else bs[b] <- plugin(sit, X1[sample.int(n1, n1, TRUE), , drop = FALSE], X2[sample.int(n2, n2, TRUE), , drop = FALSE])
  }
  fail <- mean(is.na(bs)); bs <- bs[!is.na(bs)]
  pct <- quantile(bs, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  # BCa: bias correction z0 and jackknife acceleration
  z0 <- qnorm(mean(bs < est) + 0.5 * mean(bs == est))
  if (sit == "c") jk <- sapply(1:n1, function(i) plugin(sit, X1[-i, , drop = FALSE], X2[-i, , drop = FALSE]))
  else jk <- c(sapply(1:n1, function(i) plugin(sit, X1[-i, , drop = FALSE], X2)),
               sapply(1:n2, function(i) plugin(sit, X1, X2[-i, , drop = FALSE])))
  jk <- jk[!is.na(jk)]; u <- mean(jk) - jk
  a <- sum(u^3) / (6 * sum(u^2)^1.5)
  zq <- qnorm(c(alpha / 2, 1 - alpha / 2))
  p <- pnorm(z0 + (z0 + zq) / (1 - a * (z0 + zq)))
  bca <- quantile(bs, p, names = FALSE)
  c(pct, bca, fail)
}

one_rep <- function(st, n) {
  k <- st$k
  if (st$sit == "c") {
    Z <- mvrnorm(n, c(st$m, rep(0, k)), st$J)
    X1 <- Z[, 1:k, drop = FALSE]; X2 <- Z[, k + 1:k, drop = FALSE]
  } else {
    X1 <- mvrnorm(n, st$m, st$S1); X2 <- mvrnorm(n, rep(0, k), st$S2)
  }
  est <- plugin(st$sit, X1, X2)
  an <- analytic(st$sit, X1, X2)
  bt <- boot_ci(st$sit, X1, X2, est)
  c(plug = est, an, pct_lo = bt[1], pct_hi = bt[2], bca_lo = bt[3], bca_hi = bt[4], bfail = bt[5])
}

grid <- expand.grid(n = c(20, 50), k = c(2, 5, 10), sit = c("a", "b", "c"), stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG"); set.seed(20261025)
res <- list()
for (g in seq_len(nrow(grid))) {
  st <- make_setting(grid$sit[g], grid$k[g]); n <- grid$n[g]
  out <- do.call(rbind, mclapply(1:R_REP, function(r) one_rep(st, n), mc.cores = detectCores(), mc.set.seed = TRUE))
  cov_ <- function(lo, hi) mean(!is.na(lo) & !is.na(hi) & lo <= theta0 & theta0 <= hi)
  res[[g]] <- data.frame(sit = grid$sit[g], k = grid$k[g], n = n, theta = theta0,
    plugin_mean = mean(out[, "plug"]), bc_mean = mean(out[, "bc"]),
    cov_exact = cov_(out[, "lo"], out[, "hi"]), width_exact = mean(out[, "hi"] - out[, "lo"]),
    cov_pct = cov_(out[, "pct_lo"], out[, "pct_hi"]),
    cov_bca = cov_(out[, "bca_lo"], out[, "bca_hi"]), width_bca = median(out[, "bca_hi"] - out[, "bca_lo"], na.rm = TRUE),
    boot_fail = mean(out[, "bfail"]))
  print(res[[g]], digits = 3); flush.console()
}
tab <- do.call(rbind, res)
write.csv(tab, "coverage_results.csv", row.names = FALSE)
cat("R_REP =", R_REP, " B =", B, "\n")
print(tab, digits = 3)
