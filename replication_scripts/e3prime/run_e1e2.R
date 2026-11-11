## Rerun of E1 / E2(i) -- scope restricted to what the article tables
## show (see DESIGN.md, section "Rerun of E1 / E2(i)"):
##  E1detail : n=1000, SNR=0.9, rho/lambda in {0,0.2,0.6,0.9}
##     SAR : gam, gamboost, lagsarlm, BSPA_SAR_ML, BSPA_SAR_CFE,
##           GAM_SAR_ML, GAM_SAR_CFE          (tab_E1_1)
##     SEM : gam, gamboost, errorsarlm, BSPA_SEM_ML, BSPA_SEM_CFE,
##           GAM_SEM_CFE                      (tab_E1_2)
##  E1synth  : n in {1000,2000,5000} x SNR in {0.7,0.9} x 4 rho values,
##     CFE only (BSPA_*_CFE, GAM_*_CFE), SAR and SEM
##     (tab_E1_synth, tab_E1_synth_SEM, tab_coverage_*)
##  E2i      : n=2000, SNR=0.7, X_het=TRUE, SAR, 4 rho values,
##     lagsarlm, BSPA_SAR_CFE, GAM_SAR_CFE    (tab_E2_1)
## mstop of the boosting arms: selected by spatial-block CV (grid
## 1..2000, 5 folds) -- the protocol recommended by the E3' mstop-CV
## sub-experiment ("practice what we preach"); the median selected mstop
## is reported in the tables.
## One task = one replication (all cells, all estimators) -> one RDS.
## Env : REPS="1:1000" NCORES=1 OUT=out_e1e2 EXPS="E1detail,E1synth,E2i"

suppressMessages({
  library(spboost); library(Matrix); library(mboost); library(parallel)
})
env   <- Sys.getenv
reps  <- eval(parse(text = paste0("c(", env("REPS", "1"), ")")))
ncores <- max(1L, as.integer(env("NCORES", "1")))
out_dir <- env("OUT", "out_e1e2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
exps <- strsplit(env("EXPS", "E1detail,E1synth,E2i"), ",")[[1]]
RHOS <- c(0, 0.2, 0.6, 0.9)
MSTOP_MAX <- 2000L
GAM_SEM_CFE_ <- utils::getFromNamespace("GAM_SEM_CFE", "spboost")
scatter_folds_ <- utils::getFromNamespace("spb_scatter_folds", "spboost")
CV_CTRL <- list(mstop_criterion = "CV", cv_mode = "spatial_block",
                nfold = 5L,
                control_gamboost = mboost::boost_control(mstop = MSTOP_MAX),
                mstop_init = MSTOP_MAX)

cells <- do.call(rbind, list(
  if ("E1detail" %in% exps)
    expand.grid(exp = "E1detail", model = c("SAR", "SEM"), n = 1000L,
                snr = 0.9, rho0 = RHOS, stringsAsFactors = FALSE),
  if ("E1synth" %in% exps)
    expand.grid(exp = "E1synth", model = c("SAR", "SEM"),
                n = c(1000L, 2000L, 5000L), snr = c(0.7, 0.9), rho0 = RHOS,
                stringsAsFactors = FALSE),
  if ("E2i" %in% exps)
    expand.grid(exp = "E2i", model = "SAR", n = 2000L, snr = 0.7,
                rho0 = RHOS, stringsAsFactors = FALSE)))
cat(sprintf("[e1e2] %d cells x %d reps, %d cores, out=%s\n",
            nrow(cells), length(reps), ncores, out_dir))

## ---- metric extraction ----------------------------------------------
## truths: columns f1..f3 of dgp (centred); true f0 = 0
mise_from_mat <- function(fh_list, d) {
  out <- rep(NA_real_, 3)
  for (j in 1:3) {
    fh <- fh_list[[j]]
    if (!is.null(fh)) {
      fh <- as.numeric(fh); fh <- fh - mean(fh)
      tr <- d[[paste0("f", j)]]; tr <- tr - mean(tr)
      out[j] <- mean((fh - tr)^2)
    }
  }
  out
}

metrics_spboost <- function(fit, d) {
  dec <- fitted_decomp_spboost(fit)
  m <- mise_from_mat(list(dec$X1, dec$X2, dec$X3), d)
  ms <- tryCatch(as.numeric(fit$mstop_selection$best_mstop)[1],
                 error = function(e) NA_real_)
  if (is.null(ms) || !length(ms)) ms <- NA_real_
  list(par = as.numeric(fit$rho)[1],
       var_par = tryCatch(as.numeric(fit$var_rho)[1], error = function(e) NA_real_),
       f0 = mean(dec$Intercept), mise = m, mstop_sel = ms)
}
metrics_gam <- function(fit, d) {
  tp <- mgcv::predict.gam(fit, type = "terms")
  g <- function(j) if (paste0("s(X", j, ")") %in% colnames(tp))
    tp[, paste0("s(X", j, ")")] else NULL
  list(par = NA_real_, var_par = NA_real_,
       f0 = unname(coef(fit)[1]),
       mise = mise_from_mat(list(g(1), g(2), g(3)), d))
}
metrics_gamboost <- function(fit, d) {
  g <- function(j) tryCatch(
    as.numeric(predict(fit, which = paste0("bbs(X", j, ")"))),
    error = function(e) NULL)
  list(par = NA_real_, var_par = NA_real_,
       f0 = tryCatch(as.numeric(fit$offset)[1], error = function(e) NA_real_),
       mise = mise_from_mat(list(g(1), g(2), g(3)), d))
}
metrics_linear <- function(fit, d, par, separ) {
  cf <- coef(fit)
  g <- function(j) {
    nm <- paste0("X", j)
    if (nm %in% names(cf)) cf[[nm]] * d[[nm]] else NULL
  }
  list(par = as.numeric(par)[1],
       var_par = if (is.null(separ) || !is.finite(separ)) NA_real_ else separ^2,
       f0 = unname(cf[["(Intercept)"]]),
       mise = mise_from_mat(list(g(1), g(2), g(3)), d))
}

run_rep <- function(r) {
  f_out <- file.path(out_dir, sprintf("rep%04d.rds", r))
  if (file.exists(f_out)) return(invisible(f_out))
  rows <- list()

  for (ic in seq_len(nrow(cells))) tryCatch(withCallingHandlers({
    ce <- cells[ic, ]
    args <- list(n = ce$n, rho = ce$rho0, sigma2 = 1, model = ce$model,
                 nonlin = TRUE, myseed = r, SNR = ce$snr)
    if (ce$exp == "E2i") { args$X_het <- TRUE; args$zeta <- 0.5 }
    ## article E2(i): error shift 0.5*sqrt(|u-0.35|+|v-0.65|) -> zeta=0.5
    sim <- do.call(dgp, args)
    d <- sim$data; W <- sim$W
    f_gam <- Y ~ s(X1) + s(X2) + s(X3)
    f_bbs <- Y ~ bbs(X1) + bbs(X2) + bbs(X3)
    f_lin <- Y ~ X1 + X2 + X3
    lw <- NULL
    mk_lw <- function() {
      if (is.null(lw)) lw <<- spdep::mat2listw(W, style = "W")
      lw
    }
    ## non-spatial gamboost: scattered-block CV via cvrisk
    cv_gamboost <- function() {
      fit <- mboost::gamboost(f_bbs, data = d,
                              control = mboost::boost_control(mstop = MSTOP_MAX))
      fl <- scatter_folds_(cbind(d$x, d$y), nfold = 5L, seed = r)
      fw <- sapply(sort(unique(fl)), function(k) as.numeric(fl != k))
      cvr <- mboost::cvrisk(fit, folds = fw,
                            grid = seq(10L, MSTOP_MAX, by = 10L),
                            papply = lapply)
      ms <- mboost::mstop(cvr)
      m <- metrics_gamboost(fit[ms], d)
      m$mstop_sel <- ms
      m
    }

    arms <- list()
    cfe_arms_sar <- list(
      BSPA_SAR_CFE = function() metrics_spboost(
        spbgam(f_bbs, data = d, W = W, DGP = "SAR", method = "BSPA_SAR_CFE",
               control = CV_CTRL), d),
      GAM_SAR_CFE = function() metrics_spboost(
        GAM_SAR_CFE(f_gam, data = d, W = W), d))
    cfe_arms_sem <- list(
      BSPA_SEM_CFE = function() metrics_spboost(
        spbgam(f_bbs, data = d, W = W, DGP = "SEM", method = "BSPA_SEM_CFE",
               control = CV_CTRL), d),
      GAM_SEM_CFE = function() metrics_spboost(
        GAM_SEM_CFE_(f_gam, data = d, W = W), d))

    if (ce$exp %in% c("E1detail", "E2i") && ce$model == "SAR") {
      if (ce$exp == "E1detail") {
        arms$gam <- function() metrics_gam(mgcv::gam(f_gam, data = d), d)
        arms$gamboost <- cv_gamboost
        arms$BSPA_SAR_ML <- function() metrics_spboost(
          spbgam(f_bbs, data = d, W = W, DGP = "SAR",
                 method = "BSPA_SAR_ML", control = CV_CTRL), d)
        arms$GAM_SAR_ML <- function() metrics_spboost(
          GAM_SAR_ML(f_gam, data = d, W = W), d)
      }
      arms$lagsarlm <- function() {
        fit <- spatialreg::lagsarlm(f_lin, data = d, listw = mk_lw(),
                                    method = "LU")
        metrics_linear(fit, d, fit$rho, fit$rho.se)
      }
      arms <- c(arms, cfe_arms_sar)
    } else if (ce$exp == "E1detail" && ce$model == "SEM") {
      arms$gam <- function() metrics_gam(mgcv::gam(f_gam, data = d), d)
      arms$gamboost <- cv_gamboost
      arms$errorsarlm <- function() {
        fit <- spatialreg::errorsarlm(f_lin, data = d, listw = mk_lw(),
                                      method = "LU")
        metrics_linear(fit, d, fit$lambda, fit$lambda.se)
      }
      arms$BSPA_SEM_ML <- function() metrics_spboost(
        spbgam(f_bbs, data = d, W = W, DGP = "SEM",
               method = "BSPA_SEM_ML", control = CV_CTRL), d)
      arms <- c(arms, cfe_arms_sem)
    } else if (ce$exp == "E1synth") {
      arms <- if (ce$model == "SAR") cfe_arms_sar else cfe_arms_sem
    }

    for (a in names(arms)) {
      t0 <- proc.time()[3]
      row <- data.frame(exp = ce$exp, model = ce$model, n = ce$n,
                        snr = ce$snr, rho0 = ce$rho0, rep = r, est = a,
                        par_hat = NA_real_, var_par = NA_real_,
                        f0_hat = NA_real_, mise1 = NA_real_,
                        mise2 = NA_real_, mise3 = NA_real_,
                        mstop_sel = NA_real_,
                        time_s = NA_real_, err = NA_character_,
                        stringsAsFactors = FALSE)
      res <- tryCatch(arms[[a]](), error = function(e)
        list(err = conditionMessage(e)))
      if (is.null(res$err)) {
        row$par_hat <- res$par; row$var_par <- res$var_par
        row$f0_hat <- res$f0
        if (!is.null(res$mstop_sel)) row$mstop_sel <- res$mstop_sel
        row$mise1 <- res$mise[1]; row$mise2 <- res$mise[2]
        row$mise3 <- res$mise[3]
      } else row$err <- res$err
      row$time_s <- round(proc.time()[3] - t0, 2)
      rows[[paste(ic, a)]] <- row
    }
  }, warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) {
      ce <- cells[ic, ]
      rows[[paste(ic, "CELL")]] <<- data.frame(
        exp = ce$exp, model = ce$model, n = ce$n, snr = ce$snr,
        rho0 = ce$rho0, rep = r, est = "CELL",
        par_hat = NA_real_, var_par = NA_real_, f0_hat = NA_real_,
        mise1 = NA_real_, mise2 = NA_real_, mise3 = NA_real_,
        mstop_sel = NA_real_,
        time_s = NA_real_, err = conditionMessage(e),
        stringsAsFactors = FALSE)
    })

  saveRDS(do.call(rbind, rows), f_out)
  invisible(f_out)
}

t0 <- proc.time()[3]
if (ncores > 1L) {
  invisible(mclapply(reps, run_rep, mc.cores = ncores,
                     mc.preschedule = FALSE))
} else {
  for (r in reps) run_rep(r)
}
cat(sprintf("[e1e2] done: %d reps in %.0fs\n",
            length(reps), proc.time()[3] - t0))
