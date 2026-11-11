## Lucas County estimation: full paper-1 workflow over the W menu.
##   A. GAM_SAR_CFE, free surface s(u,v)        (trap exposure), +/- crossfit
##   B. GAM_SAR_CFE, mu_s path (spx, kspec)    -> rho_hat plateau; + REML
##   C. GAM_SAR_CFE, retained plateau + crossfit   (final specification)
##   D. BSPA_SAR_CFE, selection, mstop by spatial-block CV (cap MSTOP_MAX)
##   E. pi_k (25 modes) on the free and the penalized surface
## Env: WSEL=knn4,knn10,dband  MSTOP_MAX=2000  KSPEC=100
##       MU_GRID=0.01,1,100,1000,10000  SMOKE=0  OUT=out_lucas  SKIP_BSPA=0
## Usage: Rscript 02_fit_lucas.R   (after 01_prep_lucas.R, same directory)

suppressMessages({ library(spboost); library(mboost); library(Matrix) })
env <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

prep <- readRDS("lucas_prep.rds")
d <- prep$d
wsel      <- strsplit(env("WSEL", "knn4,knn10,dband"), ",")[[1]]
MSTOP_MAX <- as.integer(env("MSTOP_MAX", "2000"))
KSPEC     <- as.integer(env("KSPEC", "100"))
mu_grid   <- as.numeric(strsplit(env("MU_GRID", "0.01,1,100,1000,10000"), ",")[[1]])
out_dir   <- env("OUT", "out_lucas"); dir.create(out_dir, FALSE, TRUE)
smoke     <- as.integer(env("SMOKE", "0"))
skip_bspa <- as.integer(env("SKIP_BSPA", "0")) == 1L

if (smoke > 0L) {               # pipeline check on a subsample
  ## NB: the W menu is REBUILT on the subsample (a subset
  ## of the full W = perforated graph -> rho_hat ~ 0 mechanically).
  set.seed(1); idx <- sort(sample.int(nrow(d), smoke))
  d <- d[idx, ]
  cc <- as.matrix(d[, c("x", "y")]) * 1000
  lw2sparse <- function(lw) as(spatialreg::as_dgRMatrix_listw(lw), "CsparseMatrix")
  knn_W <- function(k) lw2sparse(spdep::nb2listw(
    spdep::knn2nb(spdep::knearneigh(cc, k = k)), style = "W"))
  d1 <- unlist(spdep::nbdists(spdep::knn2nb(spdep::knearneigh(cc, k = 1L)), cc))
  prep$Wlist <- list(
    knn4 = knn_W(4L), knn10 = knn_W(10L),
    dband = lw2sparse(spdep::nb2listw(
      spdep::dnearneigh(cc, -1, max(d1) * 1.001), style = "W")))
}
n <- nrow(d)

cont   <- c("yrbuilt", "TLA", "lotsize", "frontage", "depth",
            "garagesqft", "dcbd")
counts <- c("rooms", "beds", "baths", "halfbaths")
fact   <- c("wall", "garage", "stories", "syear")

f_bbs <- as.formula(paste("lprice ~",
  paste(c(sprintf("bols(%s)", c(cont, counts)),
          sprintf("bbs(%s, df = 4)", cont),
          sprintf("bols(%s)", fact),
          "bspatial(x, y)"), collapse = " + ")))
rhs_core <- paste(c(sprintf("s(%s)", cont), counts, fact), collapse = " + ")
f_free <- as.formula(paste("lprice ~", rhs_core, "+ s(x, y, k = 60)"))
mk_spx <- function(sp_str) as.formula(paste(
  "lprice ~", rhs_core,
  sprintf("+ s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = %d, eigs = EIGS)%s)",
          KSPEC, sp_str)))

CV_CTRL <- list(mstop_criterion = "CV", cv_mode = "spatial_block", nfold = 5L,
                control_gamboost = mboost::boost_control(mstop = MSTOP_MAX),
                mstop_init = MSTOP_MAX)

grab <- function(fit) {
  if (!inherits(fit, "spboost")) return(list(error = as.character(fit)))
  list(rho = as.numeric(fit$rho)[1],
       var_rho   = tryCatch(as.numeric(fit$var_rho)[1],      error = function(e) NA_real_),
       var_hc0   = tryCatch(as.numeric(fit$var_rho_hc0)[1],  error = function(e) NA_real_),
       var_shac  = tryCatch(as.numeric(fit$var_rho_shac)[1], error = function(e) NA_real_),
       rho_method = tryCatch(fit$rho_method, error = function(e) NA_character_),
       mstop = tryCatch(as.numeric(fit$mstop_selection$best_mstop)[1],
                        error = function(e) NA_real_),
       rmse = tryCatch(fit$rmse, error = function(e) NA_real_),
       elapsed = tryCatch(fit$elapsed, error = function(e) NA_real_))
}
safely <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))
nz <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]

for (wn in wsel) {
  W <- prep$Wlist[[wn]]
  cat(sprintf("== W = %s (n=%d) ==\n", wn, n))
  ## leading eigenvectors of (W+W')/2, shared by spx and pi_k
  S <- (W + Matrix::t(W)) / 2
  eg <- RSpectra::eigs_sym(S, k = KSPEC, which = "LA")
  EIGS <- list(vectors = eg$vectors, values = eg$values)
  res <- list(W = wn, n = n)

  ## A. free surface, without then with crossfit
  t0 <- proc.time()[3]
  fit_free <- safely(GAM_SAR_CFE(f_free, data = d, W = W))
  res$free <- if (is.character(fit_free)) list(error = fit_free) else grab(fit_free)
  res$free_cf <- grab(safely(GAM_SAR_CFE(f_free, data = d, W = W,
                                         crossfit = TRUE,
                                         cf_coord_vars = c("x", "y"))))
  cat(sprintf("  free: rho=%.4f (cf %.4f)  [%.0fs]\n",
              nz(res$free$rho), nz(res$free_cf$rho), proc.time()[3] - t0))

  ## B. mu_s path + REML
  ## engine="gam" for the spx fits: bam builds the smooths on a
  ## subsample (discretisation / chunk.size), incompatible with the
  ## spx constructor, which requires the n rows aligned with W
  res$mu_path <- lapply(mu_grid, function(mu) {
    r <- grab(safely(GAM_SAR_CFE(mk_spx(sprintf(", sp = c(-1, %g)", mu)),
                                 data = d, W = W, engine = "gam")))
    r$mu <- mu; r
  })
  res$mu_reml <- grab(safely(GAM_SAR_CFE(mk_spx(""), data = d, W = W,
                                         engine = "gam")))
  rhos <- vapply(res$mu_path, function(r)
    if (is.null(r$rho)) NA_real_ else r$rho, 0)
  cat(sprintf("  mu-path rho: %s | REML rho=%.4f\n",
              paste(sprintf("%.4f", rhos), collapse = " "), nz(res$mu_reml$rho)))

  ## plateau: smallest mu such that |rho(mu) - rho(mu_max)| <= 0.1 * sd(path)
  ok <- which(is.finite(rhos))
  mu_star <- if (length(ok) >= 2) {
    tol <- max(0.002, 0.1 * stats::sd(rhos[ok]))
    mu_grid[ok][which(abs(rhos[ok] - rhos[ok][length(ok)]) <= tol)[1]]
  } else NA_real_
  res$mu_star <- mu_star

  ## C. final specification: plateau + crossfit
  if (is.finite(mu_star)) {
    fit_star <- safely(GAM_SAR_CFE(mk_spx(sprintf(", sp = c(-1, %g)", mu_star)),
                                   data = d, W = W, engine = "gam",
                                   crossfit = TRUE,
                                   cf_coord_vars = c("x", "y")))
    res$star_cf <- if (is.character(fit_star)) list(error = fit_star) else grab(fit_star)
    res$star_cf$mu <- mu_star
    if (!is.null(res$star_cf$error)) {
      cat(sprintf("  plateau mu*=%g : ERROR %s\n", mu_star, res$star_cf$error))
    } else {
      cat(sprintf("  plateau mu*=%g : rho=%.4f (se chf=%.4f hc0=%.4f shac=%.4f)\n",
                  mu_star, res$star_cf$rho, sqrt(res$star_cf$var_rho),
                  sqrt(res$star_cf$var_hc0), sqrt(res$star_cf$var_shac)))
    }
    saveRDS(res, file.path(out_dir, sprintf("lucas_%s%s.rds", wn,
                                            if (smoke > 0) "_smoke" else "")))
  }

  ## D. BSPA selection (mstop by spatial-block CV)
  if (!skip_bspa) {
    t0 <- proc.time()[3]
    fit_b <- safely(spbgam(f_bbs, data = d, W = W, DGP = "SAR",
                           method = "BSPA_SAR_CFE", control = CV_CTRL))
    res$bspa <- if (is.character(fit_b)) list(error = fit_b) else grab(fit_b)
    if (!is.character(fit_b)) {
      sel <- tryCatch(names(sort(table(fit_b$xselect()), decreasing = TRUE)),
                      error = function(e) NULL)
      res$bspa$selected <- sel
    }
    cat(sprintf("  BSPA: rho=%.4f mstop=%s [%.0fs]\n",
                nz(res$bspa$rho), nz(res$bspa$mstop), proc.time()[3] - t0))
  }

  ## E. pi_k on the free and the penalized surface
  eigs25 <- list(vectors = EIGS$vectors[, 1:25], values = EIGS$values[1:25])
  if (!is.character(fit_free))
    res$pik_free <- safely(pi_k_diagnostic(fit_free, W, k = 25L, eigs = eigs25))
  if (is.finite(mu_star) && exists("fit_star") && !is.character(fit_star))
    res$pik_star <- safely(pi_k_diagnostic(fit_star, W, k = 25L, eigs = eigs25))

  saveRDS(res, file.path(out_dir, sprintf("lucas_%s%s.rds", wn,
                                          if (smoke > 0) "_smoke" else "")))
}
cat("done.\n")
