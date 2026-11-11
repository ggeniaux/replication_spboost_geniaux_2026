## Supporting computations for the extended Lucas section (knn4):
##  a) linear lagsarlm (rho, se, linear age coefficient, residual Moran)
##  b) grids of the free vs penalized s(u,v) surface (map figure)
##  c) impact multipliers: total 1/(1-rho), direct tr((I-rho W)^-1)/n
##     (Hutchinson, 20 probes), indirect share; decomposed AME of TLA
##  d) residual Moran of the penalized specification
## Output: out_lucas_ext/lucas_extras.rds
## Usage: Rscript 05_extras_lucas.R   (from replication_scripts/lucas)

suppressMessages({ library(spboost); library(Matrix); library(spdep);
  library(spatialreg) })
prep <- readRDS("lucas_prep.rds"); d <- prep$d
cc <- as.matrix(d[, c("x", "y")]) * 1000
W <- prep$Wlist$knn4
n <- nrow(d)
KSPEC <- 100L
res <- list()

## ---- a) linear SAR (lagsarlm) ----------------------------------------
nb <- spdep::knn2nb(spdep::knearneigh(cc, k = 4L))
lw <- spdep::nb2listw(nb, style = "W")
f_lin <- lprice ~ yrbuilt + TLA + lotsize + frontage + depth + garagesqft +
  dcbd + rooms + beds + baths + halfbaths + wall + garage + stories + syear
t0 <- proc.time()[3]
ml <- spatialreg::lagsarlm(f_lin, data = d, listw = lw, method = "LU")
res$lagsarlm <- list(
  rho = as.numeric(ml$rho), se_rho = as.numeric(ml$rho.se),
  coef_yrbuilt = as.numeric(coef(ml)["yrbuilt"]),
  rmse = sqrt(mean(residuals(ml)^2)),
  moran_z = as.numeric(spdep::moran.test(residuals(ml), lw)$statistic),
  elapsed = proc.time()[3] - t0)
cat(sprintf("lagsarlm: rho=%.4f (%.4f) coef_yrbuilt=%.5f moran_z=%.1f [%.0fs]\n",
            res$lagsarlm$rho, res$lagsarlm$se_rho, res$lagsarlm$coef_yrbuilt,
            res$lagsarlm$moran_z, res$lagsarlm$elapsed))

## ---- free and penalized GAM fits (reused by b, c, d) -------------------
S <- (W + Matrix::t(W)) / 2
eg <- RSpectra::eigs_sym(S, k = KSPEC, which = "LA")
EIGS <- list(vectors = eg$vectors, values = eg$values)
cont   <- c("yrbuilt", "TLA", "lotsize", "frontage", "depth",
            "garagesqft", "dcbd")
counts <- c("rooms", "beds", "baths", "halfbaths")
fact   <- c("wall", "garage", "stories", "syear")
rhs <- paste(c(sprintf("s(%s)", cont), counts, fact), collapse = " + ")
f_free <- as.formula(paste("lprice ~", rhs, "+ s(x, y, k = 60)"))
f_star <- as.formula(paste("lprice ~", rhs,
  sprintf("+ s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = %d, eigs = EIGS), sp = c(-1, 10000))", KSPEC)))
fit_free <- GAM_SAR_CFE(f_free, data = d, W = W)
fit_star <- GAM_SAR_CFE(f_star, data = d, W = W, engine = "gam")
res$rho_free <- as.numeric(fit_free$rho); res$rho_star <- as.numeric(fit_star$rho)
cat(sprintf("GAM: rho free=%.4f, penalized=%.4f\n", res$rho_free, res$rho_star))

## ---- b) surface grids ------------------------------------------------
gx <- seq(min(d$x), max(d$x), length.out = 160)
gy <- seq(min(d$y), max(d$y), length.out = 160)
grid <- expand.grid(x = gx, y = gy)
## mask: keep cells within 1.5 km of a sale
nn1 <- spdep::knearneigh(as.matrix(grid) * 1000, k = 1)  # placeholder (unused)
dist1 <- FNN::get.knnx(cc / 1000, as.matrix(grid), k = 1)$nn.dist[, 1]
keep <- dist1 <= 1.5
nd <- d[rep(1L, nrow(grid)), ]
nd$x <- grid$x; nd$y <- grid$y
surf_of <- function(fit) {
  tt <- mgcv::predict.gam(fit, newdata = nd, type = "terms")
  js <- grep("s\\(x,y\\)", colnames(tt))
  v <- as.numeric(tt[, js]); v[!keep] <- NA; v
}
res$surf <- list(gx = gx, gy = gy, keep = keep,
                 free = surf_of(fit_free), star = surf_of(fit_star))
cat("surfaces ok\n")

## ---- c) impact multipliers ------------------------------------------
diag_mean_inv <- function(rho, nprobe = 20L) {
  set.seed(42); acc <- 0
  A <- Matrix::Diagonal(n) - rho * W
  for (p in seq_len(nprobe)) {
    z <- sample(c(-1, 1), n, replace = TRUE)
    acc <- acc + sum(z * as.numeric(Matrix::solve(A, z))) / n
  }
  acc / nprobe
}
mult <- function(rho) {
  dir <- diag_mean_inv(rho)
  tot <- 1 / (1 - rho)
  c(rho = rho, direct = dir, total = tot, indirect = tot - dir,
    ind_share = (tot - dir) / tot)
}
res$mult_free <- mult(res$rho_free)
res$mult_star <- mult(res$rho_star)
res$mult_lin  <- mult(res$lagsarlm$rho)
## multipliers at the retained rho (F2pc: penalized + cross-fit, run 02)
f_cf <- "out_lucas_ext/lucas_knn4.rds"
if (file.exists(f_cf))
  res$mult_star_cf <- mult(as.numeric(readRDS(f_cf)$star_cf$rho))
print(round(rbind(free = res$mult_free, star = res$mult_star,
                  star_cf = res$mult_star_cf, lagsarlm = res$mult_lin), 4))

## AME of TLA (mean derivative of the smooth, per 100 sq ft) on the penalized fit
h <- 50
nd2 <- d; nd2$TLA <- d$TLA + h
t1 <- mgcv::predict.gam(fit_star, newdata = d,  type = "terms")
t2 <- mgcv::predict.gam(fit_star, newdata = nd2, type = "terms")
jT <- grep("s\\(TLA\\)", colnames(t1))
res$ame_TLA_per100 <- mean((t2[, jT] - t1[, jT]) / h) * 100
cat(sprintf("AME TLA (penalized): %.5f log-point / 100 sqft\n", res$ame_TLA_per100))

## ---- d) residual Moran of the penalized fit ---------------------------------------
res$moran_star_z <- as.numeric(spdep::moran.test(residuals(fit_star), lw)$statistic)
res$moran_free_z <- as.numeric(spdep::moran.test(residuals(fit_free), lw)$statistic)
cat(sprintf("Moran z: lagsarlm=%.1f, GAM free=%.1f, GAM penalized=%.1f\n",
            res$lagsarlm$moran_z, res$moran_free_z, res$moran_star_z))

saveRDS(res, "out_lucas_ext/lucas_extras.rds")
cat("written: out_lucas_ext/lucas_extras.rds\n")
