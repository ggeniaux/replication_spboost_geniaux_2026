## Regenerates the E1/E2 tables of the article from the e1e2cv rerun
## (mstop of the boosting arms by spatial-block CV).
## Input: e1e2_cells.csv + e1e2_all.rds (aggregate_e1e2.R on OUT=out_e1e2cv)
## Outputs (in $OUT_DIR, default "output/"): tab_E1_1.tex, tab_E1_2.tex,
## tab_E1_synth.tex, tab_E1_synth_SEM.tex, tab_coverage_SAR.tex,
## tab_coverage_SEM.tex, tab_E2_1.tex
## Usage: Rscript gen_tab_e1e2.R   (from replication_scripts/e3prime)

dest <- Sys.getenv("OUT_DIR", "output"); dir.create(dest, showWarnings = FALSE, recursive = TRUE)
cells <- read.csv("e1e2_cells.csv")
all <- readRDS("e1e2_all.rds")

## median of the selected mstop + variance ratio per group
key <- c("exp", "model", "n", "snr", "rho0", "est")
ms <- aggregate(mstop_sel ~ exp + model + n + snr + rho0 + est, all,
                function(x) median(x, na.rm = TRUE), na.action = na.pass)
names(ms)[7] <- "mstop_med"
vr <- do.call(rbind, lapply(split(all, all[key], drop = TRUE), function(g) {
  data.frame(g[1, key],
             ratio = stats::var(g$par_hat) / mean(g$var_par, na.rm = TRUE),
             row.names = NULL)
}))
cells <- Reduce(function(a, b) merge(a, b, all.x = TRUE), list(cells, ms, vr))

fm <- function(x, d = 4) ifelse(is.na(x), "-", sprintf(paste0("%.", d, "f"), x))
fm0 <- function(x) ifelse(is.na(x), "-", sprintf("%.0f", x))

## ---- detail tables (tab_E1_1 / tab_E1_2) ---------------------------
detail_tab <- function(model, ests, par_sym, file, label) {
  g <- cells[cells$exp == "E1detail" & cells$model == model, ]
  L <- c("\\begin{table}[!h]", "\\centering",
    sprintf("\\caption{Results for experiment E1, nonlinear %s, SNR$=0.9$, $n=1000$. Monte Carlo results over 1000 replications; \\texttt{mstop} of the boosting arms selected by spatial-block cross-validation (median reported).}", model),
    sprintf("\\label{%s}", label), "\\centering",
    "\\resizebox{0.85\\linewidth}{!}{%", "\\begin{tabular}[t]{llcccccccc}",
    "\\toprule",
    "\\multicolumn{2}{c}{ } & \\multicolumn{2}{c}{Spatial parameter} & \\multicolumn{2}{c}{Intercept} & \\multicolumn{3}{c}{Nonlinear components} & \\multicolumn{1}{c}{Boosting} \\\\",
    "\\cmidrule(l{3pt}r{3pt}){3-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){7-9} \\cmidrule(l{3pt}r{3pt}){10-10}",
    sprintf("\\multicolumn{2}{c}{ } & \\multicolumn{2}{c}{$%s$} & \\multicolumn{2}{c}{$f_0$} & \\multicolumn{1}{c}{$f_1$} & \\multicolumn{1}{c}{$f_2$} & \\multicolumn{1}{c}{$f_3$} & \\multicolumn{1}{c}{$m_{stop}$} \\\\", par_sym),
    "Estimator &   & Bias & RMSE & Bias & RMSE & RMSE & RMSE & RMSE & Median\\\\",
    "\\midrule")
  for (e in ests) {
    ge <- g[g$est == e, ]; ge <- ge[order(ge$rho0), ]
    for (i in seq_len(nrow(ge))) {
      r <- ge[i, ]
      lab <- if (i == nrow(ge))
        sprintf("\\multirow{-%d}{*}{\\raggedright\\arraybackslash %s}",
                nrow(ge), gsub("_", "\\\\_", e)) else ""
      L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\",
        lab, r$rho0, fm(r$bias_par), fm(r$rmse_par), fm(r$bias_f0),
        fm(r$rmse_f0), fm(r$rmse_f1), fm(r$rmse_f2), fm(r$rmse_f3),
        fm0(r$mstop_med)))
    }
    if (e != ests[length(ests)]) L <- c(L, "\\cmidrule{1-10}")
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(L, file.path(dest, file))
}
detail_tab("SAR", c("gam", "gamboost", "lagsarlm", "BSPA_SAR_ML",
                    "BSPA_SAR_CFE", "GAM_SAR_ML", "GAM_SAR_CFE"),
           "\\rho", "tab_E1_1.tex", "tab:SAR_E1_1000_09")
detail_tab("SEM", c("gam", "gamboost", "errorsarlm", "BSPA_SEM_ML",
                    "BSPA_SEM_CFE", "GAM_SEM_CFE"),
           "\\lambda", "tab_E1_2.tex", "tab:SEM_E1_1000_09")

## ---- synthesis tables (CFE, n x SNR) --------------------------------
synth_tab <- function(model, ests, par_sym, file, label) {
  g <- cells[cells$exp == "E1synth" & cells$model == model, ]
  L <- c("\\begin{table}[htbp]", "\\centering", "\\scriptsize",
    sprintf("\\caption{Results for experiment E1, nonlinear %s, varying size and SNR (CFE estimators, 1000 replications; \\texttt{mstop} by spatial-block cross-validation).}", model),
    sprintf("\\label{%s}", label),
    "\\resizebox{0.65\\linewidth}{!}{%", "\\begin{tabular}{llccccccc}",
    "\\toprule",
    "\\multicolumn{2}{c}{ } & \\multicolumn{2}{c}{Spatial parameter} & \\multicolumn{2}{c}{Intercept} & \\multicolumn{3}{c}{Nonlinear components}\\\\",
    sprintf("\\multicolumn{2}{c}{ } & \\multicolumn{2}{c}{$%s$} & \\multicolumn{2}{c}{$f_0$} & \\multicolumn{1}{c}{$f_1$} & \\multicolumn{1}{c}{$f_2$} & \\multicolumn{1}{c}{$f_3$}\\\\", par_sym),
    "Estimator &   & Bias & RMSE & Bias & RMSE & RMSE & RMSE & RMSE\\\\",
    "\\midrule")
  for (nn in sort(unique(g$n))) for (sn in sort(unique(g$snr))) {
    L <- c(L, sprintf("\\multicolumn{9}{l}{\\textbf{$n=%d,\\;SNR=%.1f$}}\\\\\\addlinespace[1ex]", nn, sn))
    for (e in ests) {
      ge <- g[g$est == e & g$n == nn & g$snr == sn, ]
      ge <- ge[order(ge$rho0), ]
      for (i in seq_len(nrow(ge))) {
        r <- ge[i, ]
        lab <- if (i == 1)
          sprintf("\\multirow{4}{*}{ %s}", gsub("_", "\\\\_", e)) else ""
        L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\",
          lab, r$rho0, fm(r$bias_par), fm(r$rmse_par), fm(r$bias_f0),
          fm(r$rmse_f0), fm(r$rmse_f1), fm(r$rmse_f2), fm(r$rmse_f3)))
      }
      L <- c(L, "\\addlinespace[1ex]")
    }
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(L, file.path(dest, file))
}
synth_tab("SAR", c("BSPA_SAR_CFE", "GAM_SAR_CFE"), "\\rho",
          "tab_E1_synth.tex", "tab:SAR_E1_synt")
synth_tab("SEM", c("BSPA_SEM_CFE", "GAM_SEM_CFE"), "\\lambda",
          "tab_E1_synth_SEM.tex", "tab:SEM_E1_synth")

## ---- coverage (CFE, SNR=0.7) ------------------------------------------
cov_tab2 <- function(model, ests, file, label, par_sym) {
  g <- cells[cells$exp == "E1synth" & cells$model == model &
               cells$snr == 0.7, ]
  L <- c("\\begin{table}[htbp]", "\\centering", "\\small",
    sprintf("\\caption{Empirical coverage rates and variance ratios for $\\hat%s_{CFE}\\pm 1.96\\,\\hat\\sigma(%s)$ in the nonlinear %s design (E1, SNR$=0.7$, 1000 replications). Ratio: empirical variance of $\\hat%s$ over the mean closed-form variance; values below 1 indicate a conservative formula.}",
            par_sym, par_sym, model, par_sym),
    sprintf("\\label{%s}", label),
    "\\resizebox{\\linewidth}{!}{%", "\\begin{tabular}{lcccccccccccc}",
    "\\toprule",
    sprintf("$%s$ & \\multicolumn{6}{c}{\\texttt{%s}} & \\multicolumn{6}{c}{\\texttt{%s}} \\\\",
            par_sym, gsub("_", "\\\\_", ests[1]), gsub("_", "\\\\_", ests[2])),
    "\\cmidrule(lr){2-7} \\cmidrule(lr){8-13}",
    " & \\multicolumn{2}{c}{$n=1000$} & \\multicolumn{2}{c}{$n=2000$} & \\multicolumn{2}{c}{$n=5000$} & \\multicolumn{2}{c}{$n=1000$} & \\multicolumn{2}{c}{$n=2000$} & \\multicolumn{2}{c}{$n=5000$} \\\\",
    " & Cov. & Ratio & Cov. & Ratio & Cov. & Ratio & Cov. & Ratio & Cov. & Ratio & Cov. & Ratio \\\\",
    "\\midrule")
  for (r0 in sort(unique(g$rho0))) {
    vals <- character(0)
    for (e in ests) for (nn in c(1000, 2000, 5000)) {
      r <- g[g$est == e & g$n == nn & g$rho0 == r0, ]
      vals <- c(vals,
                if (nrow(r)) sprintf("%.3f", r$cover) else "-",
                if (nrow(r)) sprintf("%.2f", r$ratio) else "-")
    }
    L <- c(L, paste0(paste(c(r0, vals), collapse = " & "), " \\\\"))
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(L, file.path(dest, file))
}
cov_tab2("SAR", c("BSPA_SAR_CFE", "GAM_SAR_CFE"), "tab_coverage_SAR.tex",
         "tab:coverage_SAR", "\\rho")
cov_tab2("SEM", c("BSPA_SEM_CFE", "GAM_SEM_CFE"), "tab_coverage_SEM.tex",
         "tab:coverage_SEM", "\\lambda")

## ---- E2(i) -------------------------------------------------------------
{
  g <- cells[cells$exp == "E2i", ]
  ests <- c("lagsarlm", "BSPA_SAR_CFE", "GAM_SAR_CFE")
  L <- c("\\begin{table}[htbp]", "\\centering",
    "\\caption{Results for experiment E2 (unobserved spatial heterogeneity in the disturbances), nonlinear SAR, $n=2000$, SNR$=0.7$, 1000 replications; \\texttt{mstop} by spatial-block cross-validation.}",
    "\\label{tab:SAR_E2}",
    "\\resizebox{0.8\\linewidth}{!}{%", "\\begin{tabular}{llcccccccc}",
    "\\toprule",
    "\\multicolumn{2}{c}{ } & \\multicolumn{2}{c}{$\\rho$} & \\multicolumn{2}{c}{$f_0$} & \\multicolumn{1}{c}{$f_1$} & \\multicolumn{1}{c}{$f_2$} & \\multicolumn{1}{c}{$f_3$} & Coverage\\\\",
    "Estimator &   & Bias & RMSE & Bias & RMSE & RMSE & RMSE & RMSE & 95\\% CI\\\\",
    "\\midrule")
  for (e in ests) {
    ge <- g[g$est == e, ]; ge <- ge[order(ge$rho0), ]
    for (i in seq_len(nrow(ge))) {
      r <- ge[i, ]
      lab <- if (i == nrow(ge))
        sprintf("\\multirow{-%d}{*}{\\raggedright\\arraybackslash %s}",
                nrow(ge), gsub("_", "\\\\_", e)) else ""
      L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\",
        lab, r$rho0, fm(r$bias_par), fm(r$rmse_par), fm(r$bias_f0),
        fm(r$rmse_f0), fm(r$rmse_f1), fm(r$rmse_f2), fm(r$rmse_f3),
        fm(r$cover, 3)))
    }
    if (e != ests[length(ests)]) L <- c(L, "\\cmidrule{1-10}")
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(L, file.path(dest, "tab_E2_1.tex"))
}
cat("tables written to", dest, "\n")
