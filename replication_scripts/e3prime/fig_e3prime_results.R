## Figures of the E3' + mstop-CV results (base R, no dependency).
## Inputs (produced by fetch_and_aggregate.sh):
##   <PREFIX>_cells.csv, <PREFIX>_all.rds   (PREFIX=e3prime_full by default)
##   e3prime_mstopcv_all.rds                (if present)
## Outputs (in $OUT_DIR, default "output/"): fig_e3prime_heatmaps.pdf, fig_e3prime_traj.pdf,
##           fig_e3prime_size.pdf, fig_mstopcv.pdf
## Usage: PREFIX=e3prime_full Rscript fig_e3prime_results.R
od <- Sys.getenv("OUT_DIR", "output"); dir.create(od, showWarnings = FALSE, recursive = TRUE)

prefix <- Sys.getenv("PREFIX", "e3prime_full")
cells <- read.csv(paste0(prefix, "_cells.csv"))
Svals <- c("S0", "S1", "S2"); Cvals <- c("C0", "C1", "C2")

## Arm codes of the paper (axis F): F2p = penalized s(u,v),
## suffix c = cross-fitted, F3 = over-penalized (RSR-type)
arm_lab <- c(
  "F0"        = "F0",
  "F1"        = "F1",
  "F2"        = "F2",
  "F3r"       = "F2p (REML)",
  "F3mu10000" = "F2p",
  "F0-CF"     = "F0c",
  "F2-CF"     = "F2c",
  "F3-CF"     = "F2pc",
  "F4"        = "F3")
## Column-header gloss for the heatmap figure
arm_gloss <- c(
  "F2"        = "F1 + s(u,v)",
  "F3mu10000" = "F1 + penalized s(u,v)",
  "F2-CF"     = "F1 + s(u,v) + cross-fit",
  "F3-CF"     = "F1 + penalized s(u,v) + cross-fit",
  "F4"        = "F1 + over-penalized s(u,v)")
lab_of <- function(a) ifelse(a %in% names(arm_lab), arm_lab[a], a)

## ---- 1. trap checkerboard: bias per S x C cell ---------------------
heat_panel <- function(df, main, zlim) {
  z <- matrix(NA_real_, 3, 3, dimnames = list(Svals, Cvals))
  for (i in seq_len(nrow(df))) z[df$S[i], df$C[i]] <- df$bias[i]
  pal <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(101)
  ## grey background for cells without the arm (e.g. F4 outside S2)
  image(1:3, 1:3, matrix(0, 3, 3), col = "grey88", axes = FALSE,
        xlab = "", ylab = "", main = main, cex.main = 0.95)
  image(1:3, 1:3, t(z[3:1, ]), zlim = zlim, col = pal, axes = FALSE,
        add = TRUE)
  axis(1, 1:3, Cvals, tick = FALSE, line = -0.7)
  axis(2, 1:3, rev(Svals), tick = FALSE, las = 1, line = -0.7)
  for (i in 1:3) for (j in 1:3) {
    v <- z[4 - i, j]
    if (is.finite(v)) text(j, i, sprintf("%.3f", v), cex = 0.85)
  }
  box()
}
arms_hm <- c("F2", "F3mu10000", "F2-CF", "F3-CF", "F4")
rhos_hm <- sort(unique(cells$rho0))
nc <- length(arms_hm)
pdf(file.path(od, "fig_e3prime_heatmaps.pdf"), width = 2.6 * nc,
    height = 2.7 * length(rhos_hm) + 0.55)
## header row: one title box per column (fitted model),
## then the panels (subtitle = value of rho_0 only)
layout(matrix(seq_len((1 + length(rhos_hm)) * nc), ncol = nc,
              byrow = TRUE),
       heights = c(0.22, rep(1, length(rhos_hm))))
op <- par(mar = c(2.2, 2.6, 2.0, 0.6))
for (a in arms_hm) {
  par(mar = c(0, 2.6, 0.2, 0.6)); plot.new()
  text(0.5, 0.45, paste0(lab_of(a), " : ", arm_gloss[a]),
       font = 2, cex = 1.0)
}
par(mar = c(2.2, 2.6, 2.0, 0.6))
g0 <- cells[cells$est == "GAM" & cells$arm %in% arms_hm, ]
zl <- max(abs(g0$bias), na.rm = TRUE) * c(-1, 1)
for (r in rhos_hm) for (a in arms_hm) {
  heat_panel(g0[g0$rho0 == r & g0$arm == a, ],
             bquote(rho[0] == .(r)), zl)
}
par(op); dev.off()

## ---- 2. paths rho_hat(mu_s) and pi_bar(mu_s) ----------------------
all <- readRDS(paste0(prefix, "_all.rds"))
tr <- all[grepl("^F3mu", all$arm), ]
if (nrow(tr)) {
  tr$mu <- as.numeric(sub("^F3mu", "", tr$arm))
  pdf(file.path(od, "fig_e3prime_traj.pdf"), width = 9, height = 4.2)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  for (what in c("rho_hat", "pibar")) {
    plot(NA, xlim = range(tr$mu), ylim = range(tr[[what]], na.rm = TRUE),
         log = "x", xlab = expression(mu[s]),
         ylab = if (what == "rho_hat") expression(hat(rho)) else
           expression(bar(pi)),
         main = if (what == "rho_hat") "Path of the spatial parameter"
         else "Path of the spectral diagnostic")
    cells_tr <- unique(tr[, c("S", "C", "rho0")])
    ## one curve per cell x rho0 (mean over reps), S2C2 in bold
    for (k in seq_len(nrow(cells_tr))) {
      ck <- cells_tr[k, ]
      g <- tr[tr$S == ck$S & tr$C == ck$C & tr$rho0 == ck$rho0, ]
      m <- aggregate(g[[what]], list(mu = g$mu), mean, na.rm = TRUE)
      trap <- ck$S == "S2" & ck$C == "C2"
      lines(m$mu, m$x, col = if (trap) "#B2182B" else "grey70",
            lwd = if (trap) 2.5 else 0.8)
      if (what == "rho_hat" && trap)
        abline(h = ck$rho0, lty = 3, col = "#B2182B")
    }
  }
  par(op); dev.off()
}

## ---- 3. size of the test of rho = 0 --------------------------------
sz <- cells[cells$rho0 == 0 & cells$est == "GAM" &
              cells$arm %in% c("F0", "F1", "F2", "F3mu10000",
                               "F0-CF", "F2-CF", "F3-CF", "F4"), ]
if (nrow(sz)) {
  pdf(file.path(od, "fig_e3prime_size.pdf"), width = 8, height = 4.5)
  op <- par(mar = c(6, 4, 2.5, 1))
  sz$cell <- paste(sz$S, sz$C, sep = "x")
  arms <- c("F0", "F0-CF", "F1", "F2", "F3mu10000", "F2-CF", "F3-CF",
            "F4")
  cls <- c("grey40", "#92C5DE", "grey60", "#B2182B", "#E08214",
           "#2166AC", "#1B7837", "#542788")
  m <- tapply(sz$rej0, list(sz$arm, sz$cell), mean)
  m <- m[intersect(arms, rownames(m)), , drop = FALSE]
  barplot(m, beside = TRUE, col = cls[match(rownames(m), arms)],
          las = 2,
          ylab = expression(paste("rejection rate of ", H[0], ": ",
                                  rho == 0, " (nominal 5%)")),
          main = expression(paste("Empirical size of the test, ",
                                  rho[0], " = 0")))
  abline(h = 0.05, lty = 2)
  legend("topleft", lab_of(rownames(m)),
         fill = cls[match(rownames(m), arms)],
         bty = "n", cex = 0.75, ncol = 2)
  par(op); dev.off()
}

## ---- 4. mstop dose-response + CV points -----------------------------
if (file.exists("e3prime_mstopcv_all.rds")) {
  mc <- readRDS("e3prime_mstopcv_all.rds")
  pdf(file.path(od, "fig_mstopcv.pdf"), width = 9, height = 4.5)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  for (ce in unique(mc$cell)) {
    g <- mc[mc$cell == ce, ]
    fixed <- g[g$arm %in% c("m100", "m400", "m1000"), ]
    fixed$mstop <- as.numeric(sub("^m", "", fixed$arm))
    plot(NA, xlim = c(80, 2200), log = "x",
         ylim = range(g$rho_hat - g$rho0, na.rm = TRUE),
         xlab = "mstop", ylab = expression(paste("bias of ", hat(rho))),
         main = paste0("BSPA_SAR_CFE, free s(u,v), ",
                       sub("_", "\u00d7", ce)))
    abline(h = 0, lty = 3)
    for (r0 in unique(g$rho0)) {
      f <- fixed[fixed$rho0 == r0, ]
      m <- aggregate(f$rho_hat - r0, list(mstop = f$mstop), mean,
                     na.rm = TRUE)
      lines(m$mstop, m$x, type = "b", pch = 19,
            col = if (r0 == 0) "#2166AC" else "#B2182B")
      for (cv in c("cvrand", "cvblock")) {
        gc <- g[g$arm == cv & g$rho0 == r0, ]
        if (nrow(gc)) points(median(gc$best_mstop, na.rm = TRUE),
                             mean(gc$rho_hat - r0, na.rm = TRUE),
                             pch = if (cv == "cvrand") 17 else 15,
                             cex = 1.4,
                             col = if (r0 == 0) "#2166AC" else "#B2182B")
      }
    }
    legend("bottomleft",
           expression(rho[0] == 0, rho[0] == 0.6,
                      "random CV (median mstop)", "spatial-block CV"),
           col = c("#2166AC", "#B2182B", "grey30", "grey30"),
           pch = c(19, 19, 17, 15), lty = c(1, 1, NA, NA), bty = "n",
           cex = 0.75)
  }
  par(op); dev.off()
}
cat("figures written\n")
