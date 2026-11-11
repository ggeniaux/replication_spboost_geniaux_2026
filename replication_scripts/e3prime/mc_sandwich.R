## MC validation of the robust variances (rem:sandwich, appendix D).
## Design: oracle GAM_SAR_CFE, n=1000, rho0 in {0.2, 0.6},
## two error regimes with the same mean variance sigma2 (the one calibrated
## at SNR=0.7 by dgp):
##   homo   : eps_i ~ N(0, sigma2)
##   hetero : sd_i proportional to exp(1.5*(X1-0.5)) (ratio ~4.5:1),
##            normalized so that E[sd_i^2] = sigma2
## Outputs: 95 % CI coverage for var_rho (closed form), HC0, SHAC.
## Env : NR=1000 NCORES=4 OUT=mc_sandwich.rds
suppressMessages({library(spboost); library(Matrix); library(parallel)})
env <- Sys.getenv
NR <- as.integer(env("NR", "1000"))
ncores <- max(1L, as.integer(env("NCORES", "4")))
out_f <- env("OUT", "mc_sandwich.rds")
n <- 1000L

one <- function(r) {
  rows <- list()
  for (rho0 in c(0.2, 0.6)) {
    sim <- dgp(n = n, rho = rho0, sigma2 = 1, model = "SAR", nonlin = TRUE,
               myseed = r, SNR = 0.7)
    d <- sim$data; W <- sim$W
    s2 <- as.numeric(sim$sigma2)[1]
    f_tot <- d$f0 + d$f1 + d$f2 + d$f3
    IA <- Matrix::Diagonal(n) - rho0 * W
    set.seed(10000 + r)
    z <- rnorm(n)
    sdh <- exp(1.5 * (d$X1 - 0.5))
    sdh <- sdh * sqrt(s2 / mean(sdh^2))
    for (reg in c("homo", "hetero")) {
      eps <- if (reg == "homo") sqrt(s2) * z else sdh * z
      d$Y <- as.numeric(Matrix::solve(IA, f_tot + eps))
      fit <- tryCatch(
        GAM_SAR_CFE(Y ~ s(X1) + s(X2) + s(X3), data = d, W = W),
        error = function(e) NULL)
      if (is.null(fit) || inherits(fit, "error")) next
      rows[[paste(rho0, reg)]] <- data.frame(
        rep = r, rho0 = rho0, regime = reg,
        rho_hat = as.numeric(fit$rho)[1],
        v_cf = as.numeric(fit$var_rho)[1],
        v_hc0 = fit$var_rho_hc0, v_shac = fit$var_rho_shac)
    }
  }
  do.call(rbind, rows)
}

res <- do.call(rbind, mclapply(seq_len(NR), one, mc.cores = ncores,
                               mc.preschedule = FALSE))
saveRDS(res, out_f)

z975 <- qnorm(0.975)
for (v in c("v_cf", "v_hc0", "v_shac")) {
  res[[paste0("cov_", v)]] <-
    as.numeric(abs(res$rho_hat - res$rho0) / sqrt(pmax(res[[v]], 0)) <= z975)
}
agg <- aggregate(cbind(rho_hat, cov_v_cf, cov_v_hc0, cov_v_shac) ~
                   rho0 + regime, res, mean)
agg$bias <- agg$rho_hat - agg$rho0
sdv <- aggregate(rho_hat ~ rho0 + regime, res, sd)
names(sdv)[3] <- "sd_emp"
vm <- aggregate(cbind(v_cf, v_hc0, v_shac) ~ rho0 + regime, res,
                function(x) sqrt(mean(x)))
names(vm)[3:5] <- c("sd_cf", "sd_hc0", "sd_shac")
agg <- Reduce(merge, list(agg, sdv, vm))
print(agg[, c("rho0", "regime", "bias", "sd_emp", "sd_cf", "sd_hc0",
              "sd_shac", "cov_v_cf", "cov_v_hc0", "cov_v_shac")],
      digits = 3, row.names = FALSE)
