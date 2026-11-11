## Generates the compact E3' table of the article (tab_E3p_main.tex) from the
## frozen aggregates: diagnostic cells S0xC0 (control) and S2xC2 (trap),
## rho0 in {0, 0.6}. Columns: bias, RMSE, size of the test rho=0 (at rho0=0)
## or 95% CI coverage (at rho0=0.6), RMSE of f2, pi_bar.
## Usage: Rscript gen_tab_e3p.R  (from replication_scripts/e3prime)

a <- read.csv("e3prime_full_cells.csv")
out_dir <- Sys.getenv("OUT_DIR", "output"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(out_dir, "tab_E3p_main.tex")  # output: $OUT_DIR/tab_E3p_main.tex

arm_order <- c("F0", "F1", "F2", "F3r", "F3mu10000", "F2-CF", "F3-CF", "F4",
               "BSPA-F2")
arm_lab <- c(
  "F0"        = "F0 (oracle)",
  "F1"        = "F1 (full specification)",
  "F2"        = "F2 (F1 + $s(u,v)$)",
  "F3r"       = "F2p, REML $\\mu_s$",
  "F3mu10000" = "F2p, $\\mu_s=10^4$",
  "F2-CF"     = "F2c (F2 + cross-fit)",
  "F3-CF"     = "F2pc (F2p + cross-fit)",
  "F4"        = "F3 (over-penalized)",
  "BSPA-F2"   = "F2, boosted (\\texttt{BSPA\\_SAR\\_CFE})")

g <- a[(a$est == "GAM" & a$arm %in% arm_order) |
         (a$est == "BSPA" & a$arm == "F2"), ]
g$arm[g$est == "BSPA"] <- "BSPA-F2"

fmt <- function(x, d = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f"), x))

lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Experiment E3: the geoadditive trap and its remedies in the",
  "control cell (S0$\\times$C0) and the adversarial cell (S2$\\times$C2).",
  "Monte Carlo results over 1000 replications, $n=1000$, SNR$=0.7$.",
  "Size: empirical rejection rate of the 5\\% $z$-test of $H_0:\\rho=0$",
  "($\\rho_0=0$ block); coverage: empirical coverage of the 95\\%",
  "confidence interval ($\\rho_0=0.6$ block), both based on the closed-form",
  "variance. $\\bar\\pi$: mean spectral projection over the 25 dominant",
  "eigenvectors of $W$ (Appendix~\\ref{app:orthogonality} and",
  "Figure~\\ref{fig:pi_k}). Fitted-model codes as in axis F of the design;",
  "the over-penalized arm F3 runs in S2 cells only.}",
  "\\label{tab:E3p_main}",
  "\\resizebox{0.9\\linewidth}{!}{%",
  "\\begin{tabular}{llccccc}",
  "\\toprule",
  "Cell & Fitted model & Bias & RMSE & Size / Cov. & RMSE $f_2$ & $\\bar\\pi$\\\\")

for (r0 in c(0, 0.6)) {
  lines <- c(lines, "\\midrule",
             sprintf("\\multicolumn{7}{l}{\\textbf{$\\rho_0 = %g$ --- %s}}\\\\",
                     r0, if (r0 == 0) "size of the test of $\\rho=0$"
                     else "coverage of the 95\\% CI"))
  for (cell in c("S0 C0", "S2 C2")) {
    sc <- strsplit(cell, " ")[[1]]
    first <- TRUE
    for (arm in arm_order) {
      row <- g[g$S == sc[1] & g$C == sc[2] & g$rho0 == r0 & g$arm == arm, ]
      if (!nrow(row)) next
      tc <- if (r0 == 0) row$rej0 else row$cover
      lines <- c(lines, sprintf(
        "%s & %s & %s & %s & %s & %s & %s\\\\",
        if (first) sprintf("S%s$\\times$C%s",
                           substr(sc[1], 2, 2), substr(sc[2], 2, 2)) else "",
        arm_lab[arm], fmt(row$bias), fmt(row$rmse), fmt(tc),
        fmt(sqrt(row$mise2)), fmt(row$pibar, 2)))
      first <- FALSE
    }
    lines <- c(lines, "\\addlinespace[0.6ex]")
  }
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
writeLines(lines, out)
cat("written:", out, "\n")
