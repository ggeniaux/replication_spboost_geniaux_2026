## E3' campaign -- "selection and spatial confounding" factorial (frozen design,
## see DESIGN.md). Axes:
##   S : S0 iid | S1 X1=(I-0.7W)^-1 X1 | S2 X2=f(u,v)
##   C : C0 none | C1 X6,X7 correlated (0.2/0.7)+noise | C2 X4-X6 distances+noise
##   rho0 in {0, 0.2, 0.6} ; n=1000, SNR=0.7, W 4-NN.
## GAM arms (GAM_SAR_CFE):
##   F0 oracle | F1 full (C!=C0) | F2 full+free s(u,v)
##   F3r  full+spx(kspec=100), mu_s by REML        (documents the REML failure)
##   F3mu{0.01,0.1,1,10,100,1e4} : fixed mu_s grid (path/plateau/delta)
##   F0-CF, F2-CF, F3-CF (mu=1e4) : cross-fitted residuals (scattered-block folds)
##   F4   (S2 cells, plugin only) : spx PER TERM (RSR-type foil)
## BSPA arms (BSPA_SAR_CFE): F0 | F1 (C!=C0) | F2 (+bspatial), mstop=100;
##   on S2xC2: F2 with mstop 400 and 1000 (is the immunity due to underfitting?)
## Env : REPS="1:1000" NCORES=1 OUT=out_full N=1000 SNR=0.7 RHOS="0,0.2,0.6"
## One task = one tuple (rho0,S,C,rep) -> one RDS (skipped if present).

suppressMessages({
  library(spboost); library(Matrix); library(mboost); library(parallel)
})

env   <- Sys.getenv
reps  <- eval(parse(text = paste0("c(", env("REPS", "1:2"), ")")))
ncores <- max(1L, as.integer(env("NCORES", "1")))
out_dir <- env("OUT", "out_full")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
n     <- as.integer(env("N", "1000"))
SNRv  <- as.numeric(env("SNR", "0.7"))
rhos  <- as.numeric(strsplit(env("RHOS", "0,0.2,0.6"), ",")[[1]])
KSPEC <- 100L
KPI   <- 25L
MU_GRID <- c(1e-2, 1e-1, 1, 10, 100, 1e4)

S_flags <- list(S0 = list(), S1 = list(X_sp = TRUE), S2 = list(f_sp = TRUE))
C_flags <- list(C0 = list(), C1 = list(X_cor = TRUE), C2 = list(f_corsp = TRUE))

couples <- expand.grid(rep = reps, rho0 = rhos,
                       S = names(S_flags), C = names(C_flags),
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
## optional cell filter, e.g. CELLS="S0_C0,S2_C2"
cells <- env("CELLS", "")
if (nzchar(cells)) {
  keep <- strsplit(cells, ",")[[1]]
  couples <- couples[paste(couples$S, couples$C, sep = "_") %in% keep, ,
                     drop = FALSE]
}
cat(sprintf("[e3prime] %d tasks, %d cores, out=%s\n",
            nrow(couples), ncores, out_dir))

empty_row <- function(cp, est, arm) {
  data.frame(S = cp$S, C = cp$C, rho0 = cp$rho0, rep = cp$rep,
             est = est, arm = arm,
             rho_hat = NA_real_, rho_method = NA_character_,
             var_rho = NA_real_,
             pibar = NA_real_, pimax = NA_real_,
             mise1 = NA_real_, mise2 = NA_real_, mise3 = NA_real_,
             edf_sxy = NA_real_, mu_s = NA_real_,
             sel = NA_character_, snr = NA_real_,
             time_s = NA_real_, err = NA_character_,
             stringsAsFactors = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

run_couple <- function(i) {
  cp <- couples[i, ]
  id <- sprintf("%s_%s_rho%02.0f_rep%04d", cp$S, cp$C, 100 * cp$rho0, cp$rep)
  f_out <- file.path(out_dir, paste0(id, ".rds"))
  if (file.exists(f_out)) return(invisible(f_out))

  rows <- list()
  withCallingHandlers({
    args <- c(list(n = n, rho = cp$rho0, sigma2 = 1, model = "SAR",
                   nonlin = TRUE, myseed = cp$rep, SNR = SNRv),
              S_flags[[cp$S]], C_flags[[cp$C]])
    sim <- do.call(dgp, args)
    d <- sim$data; W <- sim$W
    snr_real <- tryCatch(as.numeric(sim$snr$target %||% sim$snr)[1],
                         error = function(e) NA_real_)

    Ws <- (W + Matrix::t(W)) / 2
    eg <- RSpectra::eigs_sym(Ws, k = KSPEC, which = "LA")
    ord <- order(eg$values, decreasing = TRUE)
    eigs100 <- list(vectors = eg$vectors[, ord], values = eg$values[ord])
    eigs25 <- list(vectors = eigs100$vectors[, 1:KPI],
                   values = eigs100$values[1:KPI])

    cand <- if (cp$C == "C0") paste0("X", 1:3) else paste0("X", 1:7)
    oracle_rhs <- paste(sprintf("s(%s)", paste0("X", 1:3)), collapse = " + ")
    full_rhs   <- paste(sprintf("s(%s)", cand), collapse = " + ")
    sxy_free <- "s(x, y, k = 60)"
    spx_term <- function(sp_str) paste0(
      "s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = ", KSPEC,
      ", eigs = eigs100)", sp_str, ")")
    spx1d <- function(v) paste0(
      "s(", v, ", bs = 'spx', xt = list(W = W, kspec = ", KSPEC,
      ", eigs = eigs100), sp = c(-1, 1e4))")
    mkf <- function(rhs) as.formula(paste("Y ~", rhs), env = environment())

    gam_metrics <- function(fit, arm) {
      r <- list(rho_hat = as.numeric(fit$rho)[1],
                rho_method = paste(fit$rho_method, collapse = "/"),
                var_rho = as.numeric(fit$var_rho %||% NA_real_)[1],
                pibar = NA_real_, pimax = NA_real_,
                mise1 = NA_real_, mise2 = NA_real_, mise3 = NA_real_,
                edf_sxy = NA_real_, mu_s = NA_real_, sel = NA_character_)
      pk <- tryCatch(pi_k_diagnostic(fit, W, k = KPI, eigs = eigs25)$pi_k,
                     error = function(e) NULL)
      if (!is.null(pk)) { r$pibar <- mean(pk); r$pimax <- max(pk) }
      tp <- tryCatch(mgcv::predict.gam(fit, type = "terms"),
                     error = function(e) NULL)
      if (!is.null(tp)) for (j in 1:3) {
        cj <- paste0("s(X", j, ")")
        if (cj %in% colnames(tp)) {
          fh <- as.numeric(tp[, cj]); fh <- fh - mean(fh)
          tr <- d[[paste0("f", j)]];  tr <- tr - mean(tr)
          r[[paste0("mise", j)]] <- mean((fh - tr)^2)
        }
      }
      sms <- tryCatch(fit$smooth, error = function(e) NULL)
      if (!is.null(sms)) for (sm in sms) {
        if (identical(sm$term, c("x", "y"))) {
          r$edf_sxy <- sum(fit$edf[sm$first.para:sm$last.para])
        }
      }
      if (identical(arm, "F3r")) r$mu_s <- unname(tail(fit$sp, 1))
      r
    }

    bspa_metrics <- function(fit, arm) {
      r <- list(rho_hat = as.numeric(fit$rho)[1],
                rho_method = paste(fit$rho_method, collapse = "/"),
                var_rho = as.numeric(fit$var_rho %||% NA_real_)[1],
                pibar = NA_real_, pimax = NA_real_,
                mise1 = NA_real_, mise2 = NA_real_, mise3 = NA_real_,
                edf_sxy = NA_real_, mu_s = NA_real_, sel = NA_character_)
      r$sel <- tryCatch({
        bl <- names(fit$baselearner)
        tb <- table(bl[mboost::selected(fit)])
        paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = ";")
      }, error = function(e) NA_character_)
      r
    }

    do_gam <- function(arm, fml, cf = FALSE, mu_fixed = NA_real_) {
      t0 <- proc.time()[3]
      row <- empty_row(cp, "GAM", arm)
      res <- tryCatch({
        fit <- if (cf) GAM_SAR_CFE(fml, data = d, W = W, crossfit = TRUE,
                                   cf_seed = cp$rep)
               else GAM_SAR_CFE(fml, data = d, W = W)
        if (inherits(fit, "error")) stop(as.character(fit))
        m <- gam_metrics(fit, arm)
        for (nm in names(m)) row[[nm]] <- m[[nm]]
        row
      }, error = function(e) { row$err <- conditionMessage(e); row })
      if (is.finite(mu_fixed)) res$mu_s <- mu_fixed
      res$time_s <- round(proc.time()[3] - t0, 2)
      res$snr <- snr_real
      res
    }
    do_bspa <- function(arm, fml, mstop = NULL) {
      t0 <- proc.time()[3]
      row <- empty_row(cp, "BSPA", arm)
      res <- tryCatch({
        ctrl <- if (is.null(mstop)) mboost::boost_control()
                else mboost::boost_control(mstop = mstop)
        fit <- BSPA_SAR_CFE(fml, data = d, W = W, control = ctrl)
        if (inherits(fit, "error")) stop(as.character(fit))
        m <- bspa_metrics(fit, arm)
        for (nm in names(m)) row[[nm]] <- m[[nm]]
        row
      }, error = function(e) { row$err <- conditionMessage(e); row })
      res$time_s <- round(proc.time()[3] - t0, 2)
      res$snr <- snr_real
      res
    }

    ## --- GAM, plugin ---
    rows$GAM_F0 <- do_gam("F0", mkf(oracle_rhs))
    if (cp$C != "C0") rows$GAM_F1 <- do_gam("F1", mkf(full_rhs))
    rows$GAM_F2 <- do_gam("F2", mkf(paste(full_rhs, "+", sxy_free)))
    rows$GAM_F3r <- do_gam("F3r", mkf(paste(full_rhs, "+", spx_term(""))))
    for (mu in MU_GRID) {
      arm <- sprintf("F3mu%g", mu)
      sp_str <- paste0(", sp = c(-1, ", format(mu, scientific = FALSE), ")")
      rows[[paste0("GAM_", arm)]] <-
        do_gam(arm, mkf(paste(full_rhs, "+", spx_term(sp_str))), mu_fixed = mu)
    }
    ## --- GAM, cross-fit ---
    rows$GAM_F0CF <- do_gam("F0-CF", mkf(oracle_rhs), cf = TRUE)
    rows$GAM_F2CF <- do_gam("F2-CF", mkf(paste(full_rhs, "+", sxy_free)),
                            cf = TRUE)
    rows$GAM_F3CF <- do_gam(
      "F3-CF", mkf(paste(full_rhs, "+", spx_term(", sp = c(-1, 10000)"))),
      cf = TRUE, mu_fixed = 1e4)
    ## --- F4 foil (S2 only): spx per term ---
    if (cp$S == "S2") {
      vterms <- if (cp$C == "C0") "X2" else c("X2", "X4", "X5", "X6")
      others <- setdiff(cand, vterms)
      rhs4 <- paste(c(sprintf("s(%s)", others),
                      vapply(vterms, spx1d, ""),
                      spx_term(", sp = c(-1, 10000)")), collapse = " + ")
      rows$GAM_F4 <- do_gam("F4", mkf(rhs4), mu_fixed = 1e4)
    }
    ## --- BSPA ---
    bb_rhs <- function(v) paste(sprintf("bbs(%s)", v), collapse = " + ")
    rows$BSPA_F0 <- do_bspa("F0", mkf(bb_rhs(paste0("X", 1:3))))
    if (cp$C != "C0") rows$BSPA_F1 <- do_bspa("F1", mkf(bb_rhs(cand)))
    f2b <- mkf(paste(bb_rhs(cand), "+ bspatial(x, y)"))
    rows$BSPA_F2 <- do_bspa("F2", f2b)
    if (cp$S == "S2" && cp$C == "C2") {
      rows$BSPA_F2m400  <- do_bspa("F2m400", f2b, mstop = 400L)
      rows$BSPA_F2m1000 <- do_bspa("F2m1000", f2b, mstop = 1000L)
    }
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
cat(sprintf("[e3prime] done: %d tasks in %.0fs\n",
            nrow(couples), proc.time()[3] - t0))
