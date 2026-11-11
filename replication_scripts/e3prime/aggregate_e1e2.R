## Aggregation of the E1/E2(i) rerun (see DESIGN.md).
## Produces e1e2_all.rds and e1e2_cells.csv: per (exp, model, n, snr, rho0,
## estimator) -- bias/RMSE of the spatial parameter and of f0, RMSE of f1-f3,
## 95 % CI coverage (when var_par is available).
## Usage: OUT=out_e1e2 Rscript aggregate_e1e2.R

out_dir <- Sys.getenv("OUT", "out_e1e2")
fs <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
cat("files:", length(fs), "\n")
all <- do.call(rbind, lapply(fs, readRDS))
saveRDS(all, "e1e2_all.rds")

all$cover <- as.numeric(abs(all$par_hat - all$rho0) /
                          sqrt(pmax(all$var_par, 0)) <= qnorm(0.975))
key <- c("exp", "model", "n", "snr", "rho0", "est")
sp <- split(all, all[key], drop = TRUE)
agg <- do.call(rbind, lapply(sp, function(g) {
  m <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  s <- function(x) if (all(is.na(x))) NA_real_ else sd(x, na.rm = TRUE)
  data.frame(g[1, key],
             nrep = nrow(g),
             bias_par = m(g$par_hat) - g$rho0[1],
             rmse_par = sqrt((m(g$par_hat) - g$rho0[1])^2 + s(g$par_hat)^2),
             cover = m(g$cover),
             bias_f0 = m(g$f0_hat),
             rmse_f0 = sqrt(m(g$f0_hat)^2 + s(g$f0_hat)^2),
             rmse_f1 = sqrt(m(g$mise1)), rmse_f2 = sqrt(m(g$mise2)),
             rmse_f3 = sqrt(m(g$mise3)),
             n_err = sum(!is.na(g$err)),
             stringsAsFactors = FALSE, row.names = NULL)
}))
agg <- agg[order(agg$exp, agg$model, agg$n, agg$snr, agg$rho0, agg$est), ]
write.csv(agg, "e1e2_cells.csv", row.names = FALSE)

cat("\n-- E1 detail SAR (n=1000, SNR=0.9):\n")
print(agg[agg$exp == "E1detail" & agg$model == "SAR",
          c("rho0", "est", "bias_par", "rmse_par", "cover",
            "rmse_f1", "rmse_f2", "rmse_f3", "n_err")],
      digits = 3, row.names = FALSE)
cat("\n-- E2(i):\n")
print(agg[agg$exp == "E2i",
          c("rho0", "est", "bias_par", "rmse_par", "cover", "rmse_f1",
            "n_err")],
      digits = 3, row.names = FALSE)
