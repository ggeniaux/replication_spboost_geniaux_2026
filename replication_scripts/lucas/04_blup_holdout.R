## Out-of-sample RMSE, spatially dispersed test set (Lucas, knn4).
## Sampling: kmeans of coordinates into NTEST cells (10 % of n), then ONE
## point drawn at random per cell at each replication (test points
## spread over space, never clustered).
## Arms (BSPA_SAR_CFE, BLUP prediction of type "BPN"):
##   nos_rand : without s(u,v), standard-CV mstop (random)
##   nos_spat : without s(u,v), spatial-block-CV mstop   (controlled)
##   s_rand   : with bspatial(x,y), standard-CV mstop (exposed to the trap)
##   s_spat   : with bspatial(x,y), spatial-block-CV mstop (controlled)
## Env: NREP=10  MSTOP_MAX=3000  FRAC=0.10  OUT=out_blup  SEED0=100
##       REPS=<i[,j,...]> : runs only these replications (slurm array mode)
##       ARMS=<a[,b,...]> : (re)runs only these arms; arms that already succeeded
##                          in an existing RDS are kept and skipped
## Usage: Rscript 04_blup_holdout.R   (after 01_prep_lucas.R)

suppressMessages({ library(spboost); library(mboost); library(Matrix)
  library(spdep) })
env <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
NREP  <- as.integer(env("NREP", "10"))
MMAX  <- as.integer(env("MSTOP_MAX", "3000"))
FRAC  <- as.numeric(env("FRAC", "0.10"))
SEED0 <- as.integer(env("SEED0", "100"))
out_dir <- env("OUT", "out_blup"); dir.create(out_dir, FALSE, TRUE)

prep <- readRDS("lucas_prep.rds"); d <- prep$d
n <- nrow(d); ntest <- round(FRAC * n)
cc <- as.matrix(d[, c("x", "y")])

## fixed kmeans cells (reproducible), 1 test point per cell
set.seed(1)
km <- kmeans(cc, centers = ntest, iter.max = 60, nstart = 1,
             algorithm = "Lloyd")
cl <- km$cluster
cat(sprintf("kmeans: %d cells, sizes %d-%d\n",
            ntest, min(table(cl)), max(table(cl))))

knn_W <- function(xy, k = 4L) {
  nb <- spdep::knn2nb(spdep::knearneigh(xy, k = k))
  as(spatialreg::as_dgRMatrix_listw(spdep::nb2listw(nb, style = "W")),
     "CsparseMatrix")
}

cont   <- c("yrbuilt", "TLA", "lotsize", "frontage", "depth",
            "garagesqft", "dcbd")
counts <- c("rooms", "beds", "baths", "halfbaths")
fact   <- c("wall", "garage", "stories", "syear")
rhs0 <- paste(c(sprintf("bols(%s)", c(cont, counts)),
                sprintf("bbs(%s, df = 4)", cont),
                sprintf("bols(%s)", fact)), collapse = " + ")
f_nos <- as.formula(paste("lprice ~", rhs0))
f_s   <- as.formula(paste("lprice ~", rhs0, "+ bspatial(x, y)"))
ctrl <- function(mode) list(mstop_criterion = "CV", cv_mode = mode,
                            nfold = 5L,
                            control_gamboost = mboost::boost_control(mstop = MMAX),
                            mstop_init = MMAX)
arms <- list(nos_rand = list(f = f_nos, cv = "random"),
             nos_spat = list(f = f_nos, cv = "spatial_block"),
             s_rand   = list(f = f_s,   cv = "random"),
             s_spat   = list(f = f_s,   cv = "spatial_block"))
ARMS <- strsplit(env("ARMS", paste(names(arms), collapse = ",")), ",")[[1]]

reps <- if (nzchar(Sys.getenv("REPS"))) {
  as.integer(strsplit(Sys.getenv("REPS"), ",")[[1]])
} else seq_len(NREP)
for (r in reps) {
  fout <- file.path(out_dir, sprintf("blup_rep%02d.rds", r))
  res0 <- if (file.exists(fout)) readRDS(fout) else NULL
  ok0 <- vapply(ARMS, function(a) !is.null(res0[[a]]) &&
                  is.null(res0[[a]]$error) && is.null(res0[[a]]$error_pred),
                TRUE)
  if (!is.null(res0) && all(ok0)) { cat("rep", r, "already done\n"); next }
  set.seed(SEED0 + r)
  te <- vapply(split(seq_len(n), cl), function(ix)
    if (length(ix) == 1L) ix else sample(ix, 1L), 0L)
  te <- sort(unname(te)); tr <- setdiff(seq_len(n), te)
  d_tr <- d[tr, ]; d_te <- d[te, ]
  W_tr   <- knn_W(cc[tr, ])                    # 4-NN graph of the training set only
  W_full <- knn_W(cc[c(tr, te), ])             # full 4-NN graph, train-then-test order
  res <- if (is.null(res0)) list(rep = r, ntest = length(te)) else res0
  for (a in ARMS) {
    if (isTRUE(ok0[a])) next
    t0 <- proc.time()[3]
    fit <- tryCatch(
      spbgam(arms[[a]]$f, data = d_tr, W = W_tr, DGP = "SAR",
             method = "BSPA_SAR_CFE", control = ctrl(arms[[a]]$cv)),
      error = function(e) conditionMessage(e))
    if (is.character(fit)) { res[[a]] <- list(error = fit); next }
    pr <- tryCatch(
      as.numeric(predict(fit, newdata = d_te, data = d_tr, W = W_full,
                         type = "BPN")),
      error = function(e) conditionMessage(e))
    res[[a]] <- if (is.character(pr)) {
      list(error_pred = pr, rho = as.numeric(fit$rho)[1])
    } else {
      list(rho = as.numeric(fit$rho)[1],
           mstop = tryCatch(as.numeric(fit$mstop_selection$best_mstop)[1],
                            error = function(e) NA_real_),
           rmse = sqrt(mean((d_te$lprice - pr)^2)),
           mae = mean(abs(d_te$lprice - pr)),
           elapsed = proc.time()[3] - t0)
    }
    cat(sprintf("rep %d %-8s rho=%.4f mstop=%s rmse=%.4f [%.0fs]\n", r, a,
                res[[a]]$rho,
                if (is.null(res[[a]]$mstop)) "-" else res[[a]]$mstop,
                if (is.null(res[[a]]$rmse)) NA else res[[a]]$rmse,
                proc.time()[3] - t0))
  }
  saveRDS(res, fout)
}
cat("done.\n")
