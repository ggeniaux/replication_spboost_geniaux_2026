## Cross-fitted CFE mini-experiment (angle "correct the score, not the smoother")
## S2xC2, rho0 in {0, 0.6}, reps 1:10.
## Arms:
##  F0-PL / F0-CF   : oracle, plugin vs cross-fitted (sanity check)
##  F2-PL / F2-CF   : full + free s(u,v), plugin vs cross-fitted  <- central test
##  K100-PL / K100-CF : spx surface kspec=100 mu=1e4, plugin vs cross-fitted
## Cross-fitting: out-of-fold residuals of the 3 auxiliary regressions (Y, WY, W2Y),
## folds = 50 scattered kmeans blocks assigned round-robin to K=5 folds; same quadratic
## a rho^2 - b rho + c = 0; final f_hat = full refit on YS = Y - rho WY.
suppressMessages({library(spboost); library(Matrix); library(mgcv)})

K <- 5L; NCELL <- 50L
scatter_folds <- function(coords, K, ncell, seed) {
  set.seed(seed)
  cl <- kmeans(coords, ncell, nstart = 3, iter.max = 50)$cluster
  ((cl - 1L) %% K) + 1L
}

cfe_root <- function(e0, e1, e2) {
  a <- sum(e1 * e2); b <- sum(e0 * e2) + sum(e1 * e1); c <- sum(e0 * e1)
  D <- b^2 - 4 * a * c
  if (is.finite(D) && D >= 0 && abs(a) > 1e-12) {
    list(rho = (b - sqrt(D)) / (2 * a), method = "cfe_exact")
  } else {
    list(rho = c / sum(e1 * e1), method = "cfe_ratio")
  }
}

## plugin residuals (same 3 regressions, fitted on the full sample)
resid_plugin <- function(fml, d) {
  f0 <- update.formula(fml, Y ~ .); f1 <- update.formula(fml, WY ~ .)
  f2 <- update.formula(fml, W2Y ~ .)
  list(e0 = residuals(gam(f0, data = d)),
       e1 = residuals(gam(f1, data = d)),
       e2 = residuals(gam(f2, data = d)))
}

## out-of-fold residuals
resid_cf <- function(fml, d, folds) {
  f0 <- update.formula(fml, Y ~ .); f1 <- update.formula(fml, WY ~ .)
  f2 <- update.formula(fml, W2Y ~ .)
  n <- nrow(d); e0 <- e1 <- e2 <- rep(NA_real_, n)
  for (k in sort(unique(folds))) {
    tr <- folds != k; te <- !tr
    dtr <- d[tr, , drop = FALSE]
    ## for the spx arms: W restricted to the training fold
    if (!is.null(d_env$Wfull)) {
      assign("W", d_env$Wfull[tr, tr], envir = d_env)
    }
    m0 <- gam(f0, data = dtr); m1 <- gam(f1, data = dtr); m2 <- gam(f2, data = dtr)
    e0[te] <- d$Y[te]   - as.numeric(predict(m0, d[te, , drop = FALSE]))
    e1[te] <- d$WY[te]  - as.numeric(predict(m1, d[te, , drop = FALSE]))
    e2[te] <- d$W2Y[te] - as.numeric(predict(m2, d[te, , drop = FALSE]))
  }
  if (!is.null(d_env$Wfull)) assign("W", d_env$Wfull, envir = d_env)
  list(e0 = e0, e1 = e1, e2 = e2)
}

d_env <- new.env()

run_arm <- function(arm, fml, d, Wfull, folds, eigs25, truth2) {
  t0 <- proc.time()[3]
  out <- tryCatch({
    rs <- if (grepl("-CF$", arm)) resid_cf(fml, d, folds) else resid_plugin(fml, d)
    rt <- cfe_root(rs$e0, rs$e1, rs$e2)
    rho <- max(min(rt$rho, 0.99), -0.99)
    ## final interpretable refit on YS
    d$YS <- d$Y - rho * d$WY
    fs <- update.formula(fml, YS ~ .)
    assign("W", Wfull, envir = d_env)
    fit <- gam(fs, data = d)
    pk <- pi_k_diagnostic(fit, Wfull, k = 25, eigs = eigs25)$pi_k
    tp <- predict(fit, type = "terms")
    c2 <- grep("s\\(X2", colnames(tp), value = TRUE)[1]
    m2v <- if (!is.na(c2)) {
      fh <- as.numeric(tp[, c2]); fh <- fh - mean(fh)
      mean((fh - truth2)^2)
    } else NA_real_
    data.frame(arm = arm, rho_hat = rho, method = rt$method,
               pibar = mean(pk), mise2 = m2v, err = NA_character_)
  }, error = function(e) data.frame(arm = arm, rho_hat = NA, method = NA,
                                    pibar = NA, mise2 = NA,
                                    err = conditionMessage(e)))
  out$time_s <- round(proc.time()[3] - t0, 1)
  out
}

res <- list()
for (rho0 in c(0, 0.6)) for (r in 1:10) {
  sim <- dgp(n = 1000, rho = rho0, sigma2 = 1, model = "SAR", nonlin = TRUE,
             myseed = r, SNR = 0.7, f_sp = TRUE, f_corsp = TRUE)
  d <- sim$data; Wfull <- sim$W
  d$WY <- as.numeric(Wfull %*% d$Y); d$W2Y <- as.numeric(Wfull %*% d$WY)
  coords <- cbind(d$x, d$y)
  folds <- scatter_folds(coords, K, NCELL, seed = r)
  Ws <- (Wfull + Matrix::t(Wfull)) / 2
  eg <- RSpectra::eigs_sym(Ws, k = 100, which = "LA")
  ord <- order(eg$values, decreasing = TRUE)
  e100 <- list(vectors = eg$vectors[, ord], values = eg$values[ord])
  eigs25 <- list(vectors = e100$vectors[, 1:25], values = e100$values[1:25])
  truth2 <- d$f2 - mean(d$f2)

  assign("Wfull", Wfull, envir = d_env)
  assign("W", Wfull, envir = d_env)

  base7 <- paste(sprintf("s(X%d)", 1:7), collapse = " + ")
  oracle <- paste(sprintf("s(X%d)", 1:3), collapse = " + ")
  ## formulas evaluated in d_env so that W refers to the current W (fold or full)
  mkf <- function(rhs) as.formula(paste("Y ~", rhs), env = d_env)
  spx100 <- "s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = 100), sp = c(-1, 1e4))"

  arms <- list(
    `F0-PL`   = mkf(oracle),
    `F0-CF`   = mkf(oracle),
    `F2-PL`   = mkf(paste(base7, "+ s(x, y, k = 60)")),
    `F2-CF`   = mkf(paste(base7, "+ s(x, y, k = 60)")),
    `K100-PL` = mkf(paste(base7, "+", spx100)),
    `K100-CF` = mkf(paste(base7, "+", spx100))
  )
  for (a in names(arms)) {
    row <- run_arm(a, arms[[a]], d, Wfull, folds, eigs25, truth2)
    row$rho0 <- rho0; row$rep <- r
    res[[paste(rho0, r, a)]] <- row
    }
}
out <- do.call(rbind, res)
saveRDS(out, "mini_crossfit.rds")
for (rh in c(0, 0.6)) {
  cat("\n== rho0 =", rh, "==\n")
  o <- out[out$rho0 == rh, ]
  agg <- aggregate(cbind(rho_hat, pibar, mise2, time_s) ~ arm, o, mean,
                   na.action = na.pass)
  agg$bias <- agg$rho_hat - rh
  print(agg[order(match(agg$arm, c("F0-PL","F0-CF","F2-PL","F2-CF",
                                   "K100-PL","K100-CF"))), ],
        digits = 3, row.names = FALSE)
}
errs <- out$err[!is.na(out$err)]
if (length(errs)) print(table(errs))
