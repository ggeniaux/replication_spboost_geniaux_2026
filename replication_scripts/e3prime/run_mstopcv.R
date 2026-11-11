## mstop-CV sub-experiment (outside the factorial, see DESIGN.md):
## mstop is the complexity dial of boosting (the analogue of mu_s).
## Question: does a PREDICTIVE tuning of the dial (CV) work in the trap,
## as REML does for mu_s? Cells S0xC0 (control) and S2xC2 (trap),
## rho0 in {0, 0.6}. BSPA-F2 arms:
##   m100 / m400 / m1000 : fixed mstop
##   cvrand  : mstop by CV (grid 1..2000), random folds (common practice)
##   cvblock : mstop by CV (grid 1..2000), spatial folds (blockCV)
## Env : REPS="1:250" NCORES=1 OUT=out_mstopcv N=1000 SNR=0.7
suppressMessages({
  library(spboost); library(Matrix); library(mboost); library(parallel)
})
env   <- Sys.getenv
reps  <- eval(parse(text = paste0("c(", env("REPS", "1:2"), ")")))
ncores <- max(1L, as.integer(env("NCORES", "1")))
out_dir <- env("OUT", "out_mstopcv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
n    <- as.integer(env("N", "1000"))
SNRv <- as.numeric(env("SNR", "0.7"))

cells <- list(S0_C0 = list(), S2_C2 = list(f_sp = TRUE, f_corsp = TRUE))
couples <- expand.grid(rep = reps, rho0 = c(0, 0.6), cell = names(cells),
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
cat(sprintf("[mstopcv] %d tasks, %d cores, out=%s\n",
            nrow(couples), ncores, out_dir))

run_couple <- function(i) {
  cp <- couples[i, ]
  id <- sprintf("%s_rho%02.0f_rep%04d", cp$cell, 100 * cp$rho0, cp$rep)
  f_out <- file.path(out_dir, paste0(id, ".rds"))
  if (file.exists(f_out)) return(invisible(f_out))

  rows <- list()
  withCallingHandlers({
    args <- c(list(n = n, rho = cp$rho0, sigma2 = 1, model = "SAR",
                   nonlin = TRUE, myseed = cp$rep, SNR = SNRv),
              cells[[cp$cell]])
    sim <- do.call(dgp, args)
    d <- sim$data; W <- sim$W
    cand <- if (cp$cell == "S0_C0") paste0("X", 1:3) else paste0("X", 1:7)
    fml <- as.formula(paste(
      "Y ~", paste(sprintf("bbs(%s)", cand), collapse = " + "),
      "+ bspatial(x, y)"), env = environment())

    one <- function(arm, fitter) {
      t0 <- proc.time()[3]
      row <- data.frame(cell = cp$cell, rho0 = cp$rho0, rep = cp$rep,
                        arm = arm, rho_hat = NA_real_, var_rho = NA_real_,
                        best_mstop = NA_real_, cv_mode_used = NA_character_,
                        sel = NA_character_, time_s = NA_real_,
                        err = NA_character_, stringsAsFactors = FALSE)
      out <- tryCatch({
        fit <- fitter()
        if (inherits(fit, "error")) stop(as.character(fit))
        row$rho_hat <- as.numeric(fit$rho)[1]
        vr <- fit$var_rho; if (!is.null(vr)) row$var_rho <- as.numeric(vr)[1]
        ms <- fit$mstop_selection
        if (!is.null(ms)) {
          row$best_mstop <- as.numeric(ms$best_mstop)[1]
          row$cv_mode_used <- as.character(ms$cv_mode_used)[1]
        }
        row$sel <- tryCatch({
          bl <- names(fit$baselearner)
          tb <- table(bl[mboost::selected(fit)])
          paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = ";")
        }, error = function(e) NA_character_)
        row
      }, error = function(e) { row$err <- conditionMessage(e); row })
      out$time_s <- round(proc.time()[3] - t0, 2)
      out
    }

    for (m in c(100L, 400L, 1000L)) {
      rows[[paste0("m", m)]] <- one(paste0("m", m), function()
        BSPA_SAR_CFE(fml, data = d, W = W,
                     control = mboost::boost_control(mstop = m)))
    }
    rows$cvrand <- one("cvrand", function()
      spbgam(fml, data = d, W = W, DGP = "SAR", method = "BSPA_SAR_CFE",
             control = list(mstop_criterion = "CV", cv_mode = "random",
                            nfold = 5L,
                            control_gamboost = mboost::boost_control(mstop = 2000L),
                            mstop_init = 2000L)))
    rows$cvblock <- one("cvblock", function()
      spbgam(fml, data = d, W = W, DGP = "SAR", method = "BSPA_SAR_CFE",
             control = list(mstop_criterion = "CV", cv_mode = "spatial_block",
                            nfold = 5L,
                            control_gamboost = mboost::boost_control(mstop = 2000L),
                            mstop_init = 2000L)))
  }, warning = function(w) invokeRestart("muffleWarning"))

  saveRDS(do.call(rbind, rows), f_out)
  invisible(f_out)
}

t0 <- proc.time()[3]
if (ncores > 1L) {
  invisible(mclapply(seq_len(nrow(couples)), run_couple,
                     mc.cores = ncores, mc.preschedule = FALSE))
} else {
  for (i in seq_len(nrow(couples))) run_couple(i)
}
cat(sprintf("[mstopcv] done: %d tasks in %.0fs\n",
            nrow(couples), proc.time()[3] - t0))
