## ---------------------------------------------------------------------
## Figure: spectral projection diagnostic pi_k, design E1 vs design E2
## Reproduces fig_pi_k_E1_E2.pdf (Figure B.10, Appendix B of the paper)
## Usage: Rscript fig_pi_k_E1_E2.R   (from replication_scripts/); output in
## $OUT_DIR (default "output/")
##
## pi_k = phi_k' S_mu phi_k, computed with spboost::pi_k_diagnostic() on
## the auxiliary additive regression gam(Y ~ s(X1)+s(X2)+s(X3)) fitted on
## each simulated dataset. E1: baseline nonlinear design (iid covariates).
## E2: spatially autocorrelated covariate X1 = (I - 0.7 W)^{-1} nu
## (dgp(..., X_sp = TRUE)), the design where the spectral separation
## condition (HSS) is most at risk.
## ---------------------------------------------------------------------
library(spboost)

n     <- 1000L
K     <- 25L     # number of dominant eigenvectors
nrep  <- 10L     # replications averaged per design
rho0  <- 0.6
snr   <- 0.7

one_run <- function(seed, X_sp, geo = FALSE) {
  sim <- dgp(n = n, rho = rho0, sigma2 = 1, model = "SAR", nonlin = TRUE,
             X_sp = X_sp, myseed = seed, SNR = snr)
  fml <- if (geo) Y ~ s(X1) + s(X2) + s(X3) + s(x, y) else
                  Y ~ s(X1) + s(X2) + s(X3)
  fit <- mgcv::gam(fml, data = sim$data)
  pi_k_diagnostic(fit, sim$W, k = K)$pi_k
}

set.seed(1)
pi_E1  <- rowMeans(vapply(seq_len(nrep), one_run, numeric(K), X_sp = FALSE))
pi_E2  <- rowMeans(vapply(seq_len(nrep), one_run, numeric(K), X_sp = TRUE))
pi_geo <- rowMeans(vapply(seq_len(nrep), one_run, numeric(K),
                          X_sp = FALSE, geo = TRUE))

cat("mean pi_k over the 10 dominant modes: E1 =",
    round(mean(pi_E1[1:10]), 3), "| E2 =", round(mean(pi_E2[1:10]), 3),
    "| E1+s(x,y) =", round(mean(pi_geo[1:10]), 3), "\n")
cat("max pi_k: E1 =", round(max(pi_E1), 3),
    "| E2 =", round(max(pi_E2), 3),
    "| E1+s(x,y) =", round(max(pi_geo), 3), "\n")

## ---------------------------------------------------------------------
## Data-level diagnostic: spectral_overlap on the same designs
## (overlap of the spectral energy profile of each variable with that of
## the spatial lag component rho*Wy, on the K dominant modes).
## Variables screened: the three covariates and the coordinates (u,v),
## stored as columns x,y in dgp() output.
## ---------------------------------------------------------------------
ov_run <- function(seed, X_sp) {
  sim <- dgp(n = n, rho = rho0, sigma2 = 1, model = "SAR", nonlin = TRUE,
             X_sp = X_sp, myseed = seed, SNR = snr)
  res <- spboost:::spb_compute_spectral_overlap(
    W = sim$W, data = sim$data, formula = Y ~ X1 + X2 + X3,
    rho = rho0, k = K, vars_pattern = "^(X[123]|x|y)$")
  bv <- res$by_var
  bv <- bv[match(c("X1", "X2", "X3", "x", "y"), bv$variable), ]
  bv$overlap
}

ov_E1 <- rowMeans(vapply(seq_len(nrep), ov_run, numeric(5), X_sp = FALSE))
ov_E2 <- rowMeans(vapply(seq_len(nrep), ov_run, numeric(5), X_sp = TRUE))
ov_tab <- rbind(E1 = ov_E1, E2 = ov_E2)
colnames(ov_tab) <- c("X1", "X2", "X3", "u", "v")
cat("\nspectral_overlap (mean over", nrep, "replications):\n")
print(round(ov_tab, 3))

## ---------------------------------------------------------------------
## Dictionary approximation error of the true components
## (Remark rem:approx of the paper): ||(I - P_B) f_j|| / ||f_j||,
## with P_B the projection onto the mgcv spline design used by the
## auxiliary regressions (default basis sizes) and f_j the true
## centred components stored by dgp() in columns f1, f2, f3.
## ---------------------------------------------------------------------
sim0 <- dgp(n = n, rho = rho0, sigma2 = 1, model = "SAR", nonlin = TRUE,
            myseed = 1, SNR = snr)
d0   <- sim0$data
fit0 <- mgcv::gam(Y ~ s(X1) + s(X2) + s(X3), data = d0)
Bmat <- stats::model.matrix(fit0)
Qb   <- qr.Q(qr(Bmat))
eps_j <- vapply(c("f1", "f2", "f3"), function(v) {
  fv <- as.numeric(d0[[v]])
  fv <- fv - mean(fv)
  rv <- fv - Qb %*% crossprod(Qb, fv)
  sqrt(sum(rv^2) / sum(fv^2))
}, numeric(1))
cat("\ndictionary approximation error ||(I-P_B)f_j||/||f_j||:\n")
print(signif(eps_j, 3))

od  <- Sys.getenv("OUT_DIR", "output"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
out <- file.path(od, "fig_pi_k_E1_E2.pdf")
pdf(out, width = 7, height = 4.2, pointsize = 11)
op <- par(mar = c(4.2, 4.4, 1.2, 0.6))
plot(seq_len(K), pi_geo, type = "b", pch = 17, col = "#B2182B",
     ylim = c(0, 1),
     xlab = expression(paste("rank ", k, " of the eigenvector of the symmetrised ", W)),
     ylab = expression(pi[k] == phi[k]^T * S[hat(mu)] * phi[k]),
     las = 1)
lines(seq_len(K), pi_E2, type = "b", pch = 19, col = "#E8820C")
lines(seq_len(K), pi_E1, type = "b", pch = 1, col = "grey35")
legend("topright", bty = "n",
       legend = c("E1 + geographic smoother s(x,y)",
                  expression(paste("E2: ", X[1], " spatially autocorrelated")),
                  "E1: baseline (iid covariates)"),
       col = c("#B2182B", "#E8820C", "grey35"),
       pch = c(17, 19, 1), lty = 1)
par(op)
dev.off()
cat("written:", out, "\n")
