## Aggregation of the mstop-CV sub-experiment (see DESIGN.md).
## Reads OUT (one RDS per cell x rho0 x rep task) and produces:
## - e3prime_mstopcv_all.rds: raw table
## - e3prime_mstopcv_cells.csv: per (cell, rho0, arm) -- bias, sd, RMSE,
##   size of the test of rho=0 (rej0) or power, selected mstop (median, IQR)
## Usage: OUT=out_mstopcv Rscript aggregate_mstopcv.R

out_dir <- Sys.getenv("OUT", "out_mstopcv")
fs <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
cat("files:", length(fs), "\n")
all <- do.call(rbind, lapply(fs, readRDS))
saveRDS(all, "e3prime_mstopcv_all.rds")

all$z0 <- with(all, rho_hat / sqrt(pmax(var_rho, 0)))
all$rej0 <- as.numeric(abs(all$z0) > qnorm(0.975))
all$cover <- as.numeric(abs(all$rho_hat - all$rho0) /
                          sqrt(pmax(all$var_rho, 0)) <= qnorm(0.975))

mstat <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
agg <- do.call(rbind, lapply(
  split(all, list(all$cell, all$rho0, all$arm), drop = TRUE),
  function(g) data.frame(
    cell = g$cell[1], rho0 = g$rho0[1], arm = g$arm[1], n = nrow(g),
    bias = mstat(g$rho_hat, mean) - g$rho0[1],
    sd_rho = mstat(g$rho_hat, sd),
    rej0 = mstat(g$rej0, mean), cover = mstat(g$cover, mean),
    mstop_med = mstat(g$best_mstop, median),
    mstop_q1 = if (all(is.na(g$best_mstop))) NA_real_ else
      quantile(g$best_mstop, 0.25, na.rm = TRUE),
    mstop_q3 = if (all(is.na(g$best_mstop))) NA_real_ else
      quantile(g$best_mstop, 0.75, na.rm = TRUE),
    n_err = sum(!is.na(g$err)),
    stringsAsFactors = FALSE)))
agg$rmse <- sqrt(agg$bias^2 + agg$sd_rho^2)
ord_arm <- c("m100", "m400", "m1000", "cvrand", "cvblock")
agg <- agg[order(agg$cell, agg$rho0, match(agg$arm, ord_arm)), ]
write.csv(agg, "e3prime_mstopcv_cells.csv", row.names = FALSE)
print(agg[, c("cell", "rho0", "arm", "n", "bias", "rmse", "rej0",
              "mstop_med", "mstop_q1", "mstop_q3", "n_err")],
      digits = 3, row.names = FALSE)
