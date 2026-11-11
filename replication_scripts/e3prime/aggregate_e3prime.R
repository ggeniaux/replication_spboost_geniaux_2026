## Aggregation of the E3' pilot: reads OUT (one RDS per task), produces
## - e3prime_pilot_all.rds: stacked raw table
## - e3prime_pilot_cells.csv: summary per cell (S,C,rho0,est,arm)
## - e3prime_pilot_traj.csv: path rho(mu_s) (S2xC2, reps 1-10)
## Usage: OUT=out_pilot Rscript aggregate_e3prime.R

out_dir <- Sys.getenv("OUT", "out_pilot")
fs <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
cat("files:", length(fs), "\n")
all <- do.call(rbind, lapply(fs, readRDS))
all$rho_hat <- as.numeric(all$rho_hat)
saveRDS(all, "e3prime_pilot_all.rds")

## the fixed-mu_s arms run on all cells in the campaign:
## they appear in the cell table AND feed the path table
is_traj <- grepl("^F3mu", all$arm)
main <- all
agg <- aggregate(
  cbind(rho_hat, pibar, edf_sxy, mu_s, mise1, mise2, mise3, time_s) ~
    S + C + rho0 + est + arm,
  data = main, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
sdr <- aggregate(rho_hat ~ S + C + rho0 + est + arm, data = main,
                 FUN = function(x) sd(x, na.rm = TRUE), na.action = na.pass)
names(sdr)[6] <- "sd_rho"
nerr <- aggregate(err ~ S + C + rho0 + est + arm, data = main,
                  FUN = function(x) sum(!is.na(x)), na.action = na.pass)
names(nerr)[6] <- "n_err"
nr <- aggregate(rep ~ S + C + rho0 + est + arm, data = main,
                FUN = length, na.action = na.pass)
names(nr)[6] <- "n"
## test of H0: rho = 0 (z = rho_hat / sqrt(var_rho), threshold 1.96) and
## coverage of the 95 % CI around rho0 -- var_rho is the closed-form
## CFE variance, known to be conservative
main$z0 <- with(main, rho_hat / sqrt(pmax(var_rho, 0)))
main$rej0 <- as.numeric(abs(main$z0) > qnorm(0.975))
main$cover <- as.numeric(abs(main$rho_hat - main$rho0) /
                           sqrt(pmax(main$var_rho, 0)) <= qnorm(0.975))
tst <- aggregate(cbind(rej0, cover) ~ S + C + rho0 + est + arm, data = main,
                 FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
agg <- Reduce(function(a, b) merge(a, b, all.x = TRUE),
              list(agg, sdr, nerr, nr, tst))
agg$bias <- agg$rho_hat - agg$rho0
agg$rmse <- sqrt(agg$bias^2 + agg$sd_rho^2)
agg <- agg[order(agg$rho0, agg$S, agg$C, agg$est, agg$arm), ]
write.csv(agg, "e3prime_pilot_cells.csv", row.names = FALSE)

traj <- all[is_traj, ]
if (nrow(traj)) {
  traj$mu <- as.numeric(sub("^F3mu", "", traj$arm))
  tagg <- aggregate(cbind(rho_hat, pibar, edf_sxy) ~ rho0 + mu, data = traj,
                    FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
  write.csv(tagg, "e3prime_pilot_traj.csv", row.names = FALSE)
}

## console preview: the trap corner and its baseline
sel <- agg[agg$est == "GAM" &
             agg$arm %in% c("F0", "F1", "F2", "F3r", "F3mu10000",
                            "F0-CF", "F2-CF", "F3-CF", "F4"),
           c("S", "C", "rho0", "arm", "n", "bias", "rmse", "pibar",
             "rej0", "cover", "n_err")]
print(sel, digits = 3, row.names = FALSE)
cat("\nBSPA (rho):\n")
selb <- agg[agg$est == "BSPA",
            c("S", "C", "rho0", "arm", "n", "bias", "rmse", "rej0",
              "cover", "n_err")]
print(selb, digits = 3, row.names = FALSE)
