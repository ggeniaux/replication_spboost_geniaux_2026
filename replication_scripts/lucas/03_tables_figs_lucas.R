## Paper tables and figures for the Lucas County empirical section.
## Inputs: out_lucas_ext/lucas_knn4.rds, lucas_knn10.rds (extended run),
##         out_lucas/lucas_knn4_reml.rds, lucas_knn10_reml.rds (REML pi_k),
##         lucas_prep.rds (vintage-curve refit).
## Outputs (in $OUT_DIR, default "output/"): tab_lucas_rho.tex,
##   fig_lucas_mu_path.pdf, fig_lucas_pik.pdf, fig_lucas_vintage.pdf
## Usage: Rscript 03_tables_figs_lucas.R   (from replication_scripts/lucas)

suppressMessages({ library(spboost); library(Matrix) })
pdir <- Sys.getenv("OUT_DIR", "output"); dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
WS <- c(knn4 = "4-NN", knn10 = "10-NN", dband = "dist. band")
WS <- WS[file.exists(sprintf("out_lucas_ext/lucas_%s.rds", names(WS)))]
res <- lapply(names(WS), function(w)
  readRDS(sprintf("out_lucas_ext/lucas_%s.rds", w)))
names(res) <- names(WS)
extras <- if (file.exists("out_lucas_ext/lucas_extras.rds"))
  readRDS("out_lucas_ext/lucas_extras.rds") else NULL
reml <- lapply(names(WS), function(w) {
  f <- sprintf("out_lucas/lucas_%s_reml.rds", w)
  if (file.exists(f)) readRDS(f) else NULL })
names(reml) <- names(WS)

fm  <- function(x, d = 3) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "--")
se2 <- function(a) sprintf("(%s, %s)", fm(sqrt(a$var_rho), 4), fm(sqrt(a$var_shac), 4))
pib <- function(p) if (is.data.frame(p)) fm(mean(p$pi_k), 3) else "--"

## ---- table rho x specification x W ------------------------------------
star_mu <- function(r) r$star_cf$mu
path_at <- function(r, mu) {
  i <- which(vapply(r$mu_path, function(z) isTRUE(z$mu == mu), TRUE))
  r$mu_path[[i]]
}
row4 <- function(est, lab, get_a, get_p = function(r) NULL,
                 bold_w = NULL) {
  cells <- unlist(lapply(names(res), function(w) {
    r <- res[[w]]
    a <- get_a(r)
    v <- if (is.null(a$rho)) c("--", "--", "--")
         else c(fm(a$rho), se2(a), pib(get_p(r)))
    if (w %in% bold_w) v <- sprintf("\\textbf{%s}", v)
    v
  }))
  paste(c(est, lab, cells), collapse = " & ")
}
lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Lucas County: estimates of the spatial parameter across",
  "specifications and spatial weight matrices ($n=25{,}357$). Standard",
  "errors in parentheses: closed-form, then \\texttt{SHAC}. $\\bar\\pi$:",
  "mean spectral projection of the fitted surface on the 25 dominant",
  "eigenvectors of $W$. Configuration codes as in axis F of the Monte",
  "Carlo design. The retained specification (F2pc) combines the spectral",
  "penalty at the plateau value $\\mu_s^\\star$ ($10^4$ for the kNN",
  "matrices, $10^5$ for the distance band) with spatial cross-fitting. \\texttt{BSPA\\_SAR\\_CFE}: free $s(u,v)$ surface,",
  sprintf("stopping iteration by spatial-block cross-validation (%s).",
          paste(vapply(res, function(r) as.integer(r$bspa$mstop), 0L),
                collapse = ", ")),
  "The distance-band matrix uses the smallest radius leaving no",
  "observation isolated (median of 488 neighbours per observation).",
  "}",
  "\\label{tab:lucas_rho}",
  "\\resizebox{\\linewidth}{!}{%",
  sprintf("\\begin{tabular}{ll%s}", strrep("ccc", length(res))),
  "\\toprule",
  paste0(" & ", paste(ifelse(names(WS) == "knn4",
                    sprintf("& \\multicolumn{3}{c}{\\textbf{$W$: %s}}", WS),
                    sprintf("& \\multicolumn{3}{c}{$W$: %s}", WS)),
                    collapse = " "), "\\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}", 3 * seq_along(res),
                3 * seq_along(res) + 2), collapse = ""),
  paste0("Estimator & Configuration ",
         strrep("& $\\hat\\rho$ & (se, se$_{\\mathrm{SHAC}}$) & $\\bar\\pi$ ",
                length(res)), "\\\\"),
  "\\midrule",
  paste0(row4("\\multirow[t]{5}{*}{\\texttt{GAM\\_SAR\\_CFE}}",
              "F2 (free $s(u,v)$)", function(r) r$free,
              function(r) r$pik_free), "\\\\"),
  paste0(row4("", "F2c (F2 + cross-fit)", function(r) r$free_cf), "\\\\"),
  paste0(row4("", "F2p, REML $\\mu_s$", function(r) r$mu_reml,
              function(r) if (is.null(reml[[r$W]])) NULL else reml[[r$W]]$pik), "\\\\"),
  paste0(row4("", "F2p, $\\mu_s=\\mu_s^\\star$",
              function(r) path_at(r, star_mu(r))), "\\\\"),
  paste0(row4("", "\\textbf{F2pc (retained)}",
              function(r) r$star_cf,
              function(r) r$pik_star, bold_w = "knn4"), "\\\\"),
  "\\addlinespace[0.6ex]",
  paste0(row4("\\texttt{BSPA\\_SAR\\_CFE}",
              "F2, spatial-CV \\texttt{mstop}",
              function(r) r$bspa), "\\\\"),
  if (!is.null(extras)) paste0(
    "\\addlinespace[0.6ex]\n",
    paste(c("SAR (ML)", "linear covariates",
            fm(extras$lagsarlm$rho),
            sprintf("(%s, --)", fm(extras$lagsarlm$se_rho, 4)), "--",
            rep("--", 3 * (length(res) - 1))), collapse = " & "),
    "\\\\"),
  "\\bottomrule",
  "\\end{tabular}}",
  "\\end{table}")
writeLines(lines, file.path(pdir, "tab_lucas_rho.tex"))
cat("written: tab_lucas_rho.tex\n")

## ---- mu_s path figure: one panel per W, direct annotations ------------
pdf(file.path(pdir, "fig_lucas_mu_path.pdf"), width = 10.5, height = 3.8)
par(mfrow = c(1, length(WS)), mar = c(4.2, 4.2, 2.4, 0.8))
mus <- vapply(res$knn4$mu_path, `[[`, 0, "mu")
lmu <- log10(mus)
ylim <- range(unlist(lapply(res, function(r) c(
  vapply(r$mu_path, `[[`, 0, "rho"), r$free$rho, r$star_cf$rho))))
ylim <- ylim + c(-0.012, 0.012)
first <- TRUE
for (w in names(WS)) {
  r <- res[[w]]
  rho_p <- vapply(r$mu_path, `[[`, 0, "rho")
  lstar <- log10(r$star_cf$mu)
  plot(lmu, rho_p, type = "n", ylim = ylim,
       xlab = expression(log[10](mu[s])), ylab = expression(hat(rho)),
       main = sprintf("W: %s", WS[w]), cex.main = 0.95)
  ## shaded plateau zone (from the retained mu* to the right)
  usr <- par("usr")
  rect(lstar, usr[3], usr[2], usr[4], col = "#0000000A", border = NA)
  grid(col = "grey88")
  abline(h = r$free$rho, col = "grey40", lty = 3)
  lines(lmu, rho_p, col = "#2166AC", lwd = 2)
  points(lmu, rho_p, col = "#2166AC", pch = 16, cex = 0.8)
  points(lmu[1], r$mu_reml$rho, pch = 1, cex = 1.5, lwd = 1.6)
  points(lstar, r$star_cf$rho, pch = 17, cex = 1.5, col = "#B2182B")
  if (first) {   # explanatory annotations on the first panel
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    text(lmu[1] + 0.07 * dx, r$free$rho - 0.05 * dy,
         "free surface", col = "grey30", cex = 0.8, adj = 0)
    arrows(lmu[1] + 0.11 * dx, r$mu_reml$rho + 0.14 * dy,
           lmu[1] + 0.015 * dx, r$mu_reml$rho + 0.025 * dy,
           length = 0.07, col = "grey30")
    text(lmu[1] + 0.115 * dx, r$mu_reml$rho + 0.17 * dy,
         "REML lands here", col = "grey30", cex = 0.8, adj = 0)
    arrows(lstar - 0.16 * dx, r$star_cf$rho + 0.09 * dy,
           lstar - 0.025 * dx, r$star_cf$rho + 0.015 * dy,
           length = 0.07, col = "#B2182B")
    text(lstar - 0.16 * dx, r$star_cf$rho + 0.13 * dy,
         "retained:\nplateau + cross-fit", col = "#B2182B", cex = 0.8,
         adj = c(1, 0.5))
    text(mean(c(lstar, usr[2])), usr[3] + 0.05 * dy, "plateau",
         col = "grey45", cex = 0.8)
    first <- FALSE
  }
}
dev.off()
cat("written: fig_lucas_mu_path.pdf\n")

## ---- pi_k figure -------------------------------------------------------
pdf(file.path(pdir, "fig_lucas_pik.pdf"), width = 3.5 * length(WS), height = 3.6)
par(mfrow = c(1, length(WS)), mar = c(4.2, 4.2, 2, 0.6))
for (w in names(WS)) {
  r <- res[[w]]
  plot(r$pik_free$mode, r$pik_free$pi_k, type = "h", lwd = 3,
       col = "grey55", ylim = c(0, max(r$pik_free$pi_k) * 1.05),
       xlab = "eigenvector rank k", ylab = expression(pi[k]),
       main = sprintf("W: %s", WS[w]), cex.main = 0.95)
  points(r$pik_star$mode + 0.35, r$pik_star$pi_k, type = "h", lwd = 3,
         col = "#2166AC")
  legend("topright", bty = "n", cex = 0.8, lwd = 3,
         col = c("grey55", "#2166AC"),
         legend = c("free surface", "spectral penalty"))
}
dev.off()
cat("written: fig_lucas_pik.pdf\n")

## ---- vintage curve (refit of the retained specification, 4-NN) --------
prep <- readRDS("lucas_prep.rds"); d <- prep$d
W <- prep$Wlist$knn4
KSPEC <- 100L
S <- (W + Matrix::t(W)) / 2
eg <- RSpectra::eigs_sym(S, k = KSPEC, which = "LA")
EIGS <- list(vectors = eg$vectors, values = eg$values)
cont   <- c("yrbuilt", "TLA", "lotsize", "frontage", "depth",
            "garagesqft", "dcbd")
counts <- c("rooms", "beds", "baths", "halfbaths")
fact   <- c("wall", "garage", "stories", "syear")
rhs <- paste(c(sprintf("s(%s)", cont), counts, fact), collapse = " + ")
f_star <- as.formula(paste("lprice ~", rhs,
  sprintf("+ s(x, y, k = 60, bs = 'spx', xt = list(W = W, kspec = %d, eigs = EIGS), sp = c(-1, %g))",
          KSPEC, res$knn4$star_cf$mu)))
fit <- GAM_SAR_CFE(f_star, data = d, W = W, engine = "gam")
f_free2 <- as.formula(paste("lprice ~", rhs, "+ s(x, y, k = 60)"))
fit_free2 <- GAM_SAR_CFE(f_free2, data = d, W = W)
cat(sprintf("refit vintage: rho=%.4f (expected %.4f without crossfit)\n",
            fit$rho, path_at(res$knn4, star_mu(res$knn4))$rho))

term_curve <- function(var, grid, f = fit) {
  nd <- d[rep(1L, length(grid)), ]
  nd[[var]] <- grid
  tt <- mgcv::predict.gam(f, newdata = nd, type = "terms", se.fit = TRUE)
  jc <- grep(var, colnames(tt$fit))
  list(est = tt$fit[, jc], se = tt$se.fit[, jc])
}
pdf(file.path(pdir, "fig_lucas_effects.pdf"), width = 10, height = 4.2)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2, 0.6))
## (a) vintage curve + linear effect from the linear SAR
grid_yr <- seq(min(d$yrbuilt), max(d$yrbuilt), by = 1)
cv <- term_curve("yrbuilt", grid_yr)
cvf <- term_curve("yrbuilt", grid_yr, f = fit_free2)
lin <- if (!is.null(extras))
  extras$lagsarlm$coef_yrbuilt * (grid_yr - mean(d$yrbuilt)) else NULL
cat(sprintf("vintage: range controlled %.3f | free %.3f\n",
            diff(range(cv$est)), diff(range(cvf$est))))
plot(grid_yr, cv$est, type = "n",
     ylim = range(cv$est - 2 * cv$se, cv$est + 2 * cv$se, cvf$est, lin),
     xlab = "year built", ylab = "partial effect on log price",
     main = "(a) vintage", cex.main = 0.95)
grid(col = "grey85")
polygon(c(grid_yr, rev(grid_yr)),
        c(cv$est - 2 * cv$se, rev(cv$est + 2 * cv$se)),
        col = "#2166AC22", border = NA)
lines(grid_yr, cvf$est, col = "grey45", lwd = 2, lty = 4)
lines(grid_yr, cv$est, col = "#2166AC", lwd = 2)
if (!is.null(extras)) {
  lines(grid_yr, lin, col = "#B2182B", lwd = 2, lty = 2)
  legend("bottomright", bty = "n", cex = 0.85, lwd = 2, lty = c(1, 4, 2),
         col = c("#2166AC", "grey45", "#B2182B"),
         legend = c("retained specification", "free-surface specification",
                    "linear SAR (ML)"))
}
rug(d$yrbuilt[sample.int(nrow(d), 3000)], col = "#00000022")
## (b) distance-to-CBD gradient
grid_db <- seq(quantile(d$dcbd, 0.002), quantile(d$dcbd, 0.998),
               length.out = 200)
cb <- term_curve("dcbd", grid_db)
cbf0 <- term_curve("dcbd", grid_db, f = fit_free2)
lin_db0 <- if (!is.null(extras) && !is.null(extras$lagsarlm$coefs))
  extras$lagsarlm$coefs["dcbd"] * (grid_db - mean(d$dcbd)) else NULL
plot(grid_db, cb$est, type = "n",
     ylim = range(cb$est - 2 * cb$se, cb$est + 2 * cb$se, cbf0$est, lin_db0),
     xlab = "distance to CBD (km)", ylab = "partial effect on log price",
     main = "(b) distance to the CBD", cex.main = 0.95)
grid(col = "grey85")
polygon(c(grid_db, rev(grid_db)),
        c(cb$est - 2 * cb$se, rev(cb$est + 2 * cb$se)),
        col = "#2166AC22", border = NA)
lines(grid_db, cb$est, col = "#2166AC", lwd = 2)
cbf <- term_curve("dcbd", grid_db, f = fit_free2)
lines(grid_db, cbf$est, col = "grey45", lwd = 2, lty = 4)
if (!is.null(extras) && !is.null(extras$lagsarlm$coefs)) {
  lin_db <- extras$lagsarlm$coefs["dcbd"] * (grid_db - mean(d$dcbd))
  lines(grid_db, lin_db, col = "#B2182B", lwd = 2, lty = 2)
  legend("bottomleft", bty = "n", cex = 0.85, lwd = 2, lty = c(1, 4, 2),
         col = c("#2166AC", "grey45", "#B2182B"),
         legend = c("retained specification", "free-surface specification",
                    "linear SAR (ML)"))
}
rug(d$dcbd[sample.int(nrow(d), 3000)], col = "#00000022")
dev.off()
cat("written: fig_lucas_effects.pdf\n")


## ---- impact multiplier table (knn4) ------------------------------------
if (!is.null(extras)) {
  mrow <- function(lab, m) sprintf(
    "%s & %.3f & %.3f & %.3f & %.3f & %.0f\\%%\\\\",
    lab, m["rho"], m["direct"], m["indirect"], m["total"],
    100 * m["ind_share"])
  writeLines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Lucas County: average impact multipliers implied by",
    "$\\hat\\rho$ under the 4-NN matrix. Direct: $n^{-1}\\mathrm{tr}\\,(I-\\hat\\rho W)^{-1}$",
    "(Hutchinson estimator, 20 probes); total: $1/(1-\\hat\\rho)$;",
    "indirect: total minus direct. Configuration codes as in axis F.",
    "Any covariate effect scales by these",
    "factors: the average marginal effect of living area in the retained",
    sprintf("specification (F2pc), %.4f log point per 100 sq.\\,ft., splits into", extras$ame_TLA_per100),
    sprintf("%.4f direct and %.4f indirect.}",
            extras$ame_TLA_per100 * extras$mult_star_cf["direct"],
            extras$ame_TLA_per100 * extras$mult_star_cf["indirect"]),
    "\\label{tab:lucas_impacts}",
    "\\begin{tabular}{lccccc}",
    "\\toprule",
    "Configuration & $\\hat\\rho$ & Direct & Indirect & Total & Indirect share\\\\",
    "\\midrule",
    mrow("F2 (free $s(u,v)$)", extras$mult_free),
    mrow("F2pc (retained)", extras$mult_star_cf),
    mrow("linear SAR (ML)", extras$mult_lin),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"), file.path(pdir, "tab_lucas_impacts.tex"))
  cat("written: tab_lucas_impacts.tex\n")

  ## ---- surface maps ---------------------------------------------
  sf <- extras$surf
  bm <- if (file.exists("lucas_basemap.rds")) readRDS("lucas_basemap.rds") else NULL
  if (!is.null(bm)) suppressMessages(library(sf))
  zmax <- max(abs(sf$free), na.rm = TRUE)
  pal <- colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(101)
  panel <- function(v, zlim, main, xlab = "x (km)", ylab = "y (km)") {
    z <- matrix(v, length(sf$gx), length(sf$gy))
    image(sf$gx, sf$gy, z, zlim = zlim, col = pal, asp = 1,
          xlab = xlab, ylab = ylab, main = main, cex.main = 0.95,
          useRaster = TRUE)
    if (!is.null(bm)) {
      plot(bm$water, add = TRUE, col = "#C7DEED", border = "#9FC4DB",
           lwd = 0.4)
      plot(bm$boundary, add = TRUE, col = NA, border = "grey25", lwd = 1.1)
    }
  }
  zs <- max(abs(sf$star), na.rm = TRUE)
  kzoom <- round(zmax / zs)
  pdf(file.path(pdir, "fig_lucas_surfaces.pdf"), width = 10.5, height = 3.9)
  layout(matrix(1:3, 1), widths = c(1, 1, 0.62))
  par(mar = c(4, 4, 2.4, 0.8))
  panel(sf$free, c(-zmax, zmax),
        sprintf("(a) free surface (range %.2f)", 2 * zmax))
  panel(sf$star, c(-zmax, zmax), "(b) penalized surface, same scale")
  ## reduced panel (c): same surface as (b), colour scale x kzoom
  par(mar = c(7.2, 5.2, 5.6, 0.4), cex.main = 0.8, cex.axis = 0.8,
      cex.lab = 0.8)
  panel(sf$star, c(-zs, zs),
        sprintf("(c) = (b), colour scale \u00d7%d", kzoom), ylab = "")
  u <- par("usr"); par(xpd = NA)
  dx <- diff(u[1:2]); ya <- u[3] + 0.8 * diff(u[3:4])
  arrows(u[1] - 0.33 * dx, ya, u[1] - 0.14 * dx, ya,
         length = 0.09, lwd = 1.8, col = "grey30")
  text(u[1] - 0.28 * dx, ya + 0.09 * diff(u[3:4]),
       sprintf("\u00d7%d", kzoom), col = "grey30", cex = 0.95, font = 2)
  par(xpd = FALSE)
  dev.off()
  cat("written: fig_lucas_surfaces.pdf\n")
}
