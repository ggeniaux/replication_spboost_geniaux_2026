## E3' design mini-experiment: where does the spectral leakage live?
## S2xC2, rho0=0.6, reps 1:10. Arms:
##  F2            : free surface (trap reference)
##  K25/K50/K100  : spx surface only, fixed mu=1e4, varying kspec
##  T25/T100      : spx surface + spx on s(X2) and s(X4..X6) (contaminated terms)
suppressMessages({library(spboost); library(Matrix)})
KMAX <- 100L
res <- list()
for (r in 1:10) {
  sim <- dgp(n = 1000, rho = 0.6, sigma2 = 1, model = "SAR", nonlin = TRUE,
             myseed = r, SNR = 0.7, f_sp = TRUE, f_corsp = TRUE)
  d <- sim$data; W <- sim$W
  Ws <- (W + Matrix::t(W)) / 2
  eg <- RSpectra::eigs_sym(Ws, k = KMAX, which = "LA")
  ord <- order(eg$values, decreasing = TRUE)
  EG <- list(vectors = eg$vectors[, ord], values = eg$values[ord])
  sub <- function(k) list(vectors = EG$vectors[, 1:k, drop = FALSE],
                          values = EG$values[1:k])
  e25 <- sub(25L); e50 <- sub(50L); e100 <- sub(100L)

  sxy <- function(eig, k) sprintf(
    "s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = %d, eigs = %s), sp = c(-1, 1e4))", k, eig)
  sx1 <- function(v, eig, k) sprintf(
    "s(%s, bs = 'spx', xt = list(W = W, kspec = %d, eigs = %s), sp = c(-1, 1e4))", v, k, eig)
  base7 <- paste(sprintf("s(X%d)", 1:7), collapse = " + ")
  base_no_cont <- paste(c("s(X1)", "s(X3)", "s(X7)"), collapse = " + ")

  arms <- list(
    F2   = paste(base7, "+ s(x, y, k = 60)"),
    K25  = paste(base7, "+", sxy("e25", 25)),
    K50  = paste(base7, "+", sxy("e50", 50)),
    K100 = paste(base7, "+", sxy("e100", 100)),
    T25  = paste(base_no_cont, "+",
                 paste(vapply(c("X2","X4","X5","X6"), sx1, "", eig = "e25", k = 25),
                       collapse = " + "), "+", sxy("e25", 25)),
    T100 = paste(base_no_cont, "+",
                 paste(vapply(c("X2","X4","X5","X6"), sx1, "", eig = "e100", k = 100),
                       collapse = " + "), "+", sxy("e100", 100))
  )
  for (a in names(arms)) {
    fml <- as.formula(paste("Y ~", arms[[a]]), env = environment())
    row <- tryCatch({
      fit <- GAM_SAR_CFE(fml, data = d, W = W)
      pk <- pi_k_diagnostic(fit, W, k = 25, eigs = e25)$pi_k
      tp <- mgcv::predict.gam(fit, type = "terms")
      c2 <- grep("s\\(X2", colnames(tp), value = TRUE)[1]
      fh <- as.numeric(tp[, c2]); fh <- fh - mean(fh)
      tr <- d$f2 - mean(d$f2)
      data.frame(rep = r, arm = a, rho_hat = as.numeric(fit$rho)[1],
                 pibar = mean(pk), mise2 = mean((fh - tr)^2),
                 err = NA_character_)
    }, error = function(e) data.frame(rep = r, arm = a, rho_hat = NA,
                                      pibar = NA, mise2 = NA,
                                      err = conditionMessage(e)))
    res[[paste(r, a)]] <- row
  }
}
out <- do.call(rbind, res)
saveRDS(out, "mini_kspec.rds")
agg <- aggregate(cbind(rho_hat, pibar, mise2) ~ arm, out, mean, na.action = na.pass)
agg$bias <- agg$rho_hat - 0.6
print(agg[order(match(agg$arm, c("F2","K25","K50","K100","T25","T100"))), ],
      digits = 3, row.names = FALSE)
errs <- out$err[!is.na(out$err)]
if (length(errs)) print(unique(errs))
