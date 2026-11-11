## Quick fit: knn4 only, spx surface with mu_s by REML (full n).
## Output: out_lucas/lucas_knn4_reml.rds
suppressMessages({ library(spboost); library(Matrix) })
prep <- readRDS("lucas_prep.rds"); d <- prep$d
wname <- Sys.getenv("WNAME", "knn4")
W <- prep$Wlist[[wname]]
out_rds <- sprintf("out_lucas/lucas_%s_reml.rds", wname)
KSPEC <- 100L
dir.create("out_lucas", FALSE, TRUE)

S <- (W + Matrix::t(W)) / 2
eg <- RSpectra::eigs_sym(S, k = KSPEC, which = "LA")
EIGS <- list(vectors = eg$vectors, values = eg$values)
cat(sprintf("eigs ok (%d modes)\n", KSPEC))

cont   <- c("yrbuilt", "TLA", "lotsize", "frontage", "depth",
            "garagesqft", "dcbd")
counts <- c("rooms", "beds", "baths", "halfbaths")
fact   <- c("wall", "garage", "stories", "syear")
rhs_core <- paste(c(sprintf("s(%s)", cont), counts, fact), collapse = " + ")
f_reml <- as.formula(paste(
  "lprice ~", rhs_core,
  sprintf("+ s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = %d, eigs = EIGS))",
          KSPEC)))

t0 <- proc.time()[3]
## engine="gam": bam builds the smooths on a subsample
## (discretisation or chunk.size), incompatible with the spx constructor,
## which requires the n rows aligned with W
fit <- GAM_SAR_CFE(f_reml, data = d, W = W, engine = "gam")
saveRDS(list(partial = TRUE), out_rds)
el <- proc.time()[3] - t0
res <- list(W = wname, n = nrow(d), arm = "spx_REML",
            rho = as.numeric(fit$rho)[1],
            var_rho  = as.numeric(fit$var_rho)[1],
            var_hc0  = as.numeric(fit$var_rho_hc0)[1],
            var_shac = as.numeric(fit$var_rho_shac)[1],
            rho_method = fit$rho_method, rmse = fit$rmse,
            engine = fit$engine, elapsed = el)
saveRDS(res, out_rds)   # result without pi_k first
cat(sprintf("REML %s : rho=%.4f  se: chf=%.4f hc0=%.4f shac=%.4f  (%s, %.0fs)\n",
            wname, res$rho, sqrt(res$var_rho), sqrt(res$var_hc0),
            sqrt(res$var_shac), res$engine, el))
res$pik <- tryCatch(
  pi_k_diagnostic(fit, W, k = 25L,
                  eigs = list(vectors = EIGS$vectors[, 1:25],
                              values = EIGS$values[1:25])),
  error = function(e) conditionMessage(e))
saveRDS(res, out_rds)
if (is.data.frame(res$pik))
  cat("pi_k (first 5 modes):",
      paste(sprintf("%.3f", res$pik$pi_k[1:5]), collapse = " "), "\n")
