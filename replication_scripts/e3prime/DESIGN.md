# E3′ — frozen design (2026-11-05)

"Selection and spatial confounding" factorial of paper 1 (Geographical
Analysis). Frozen after a pilot with nr=100 (the HPC cluster, 1800 pairs, 0 errors) and two
local mini-experiments (kspec / per-term spx; cross-fitted CFE), both
on S2×C2, 10 repl., ρ₀ ∈ {0 ; 0.6}.

## Axes

| Axis | Levels |
|---|---|
| S — spatial structure of the covariates | S0 iid · S1 x̃₁=(I−0.7W)⁻¹x₁ (`X_sp`) · S2 x̃₂=f(u,v) distance (`f_sp`) |
| C — contamination of the candidates | C0 none · C1 X6,X7 correlated 0.2/0.7 + noise (`X_cor`) · C2 X4–X6 distances + X7 noise (`f_corsp`) |
| ρ₀ | 0 ; 0.2 ; 0.6 (0 tests the **size**: spurious ρ̂) |

n = 1000, SNR = 0.7, W 4-NN (dgp), nr = 1000, seeds = replication number
(CRN across cells and arms). Pruning: under C0, F1 ≡ F0.

## GAM arms (`GAM_SAR_CFE`)

| Arm | Specification | Role |
|---|---|---|
| F0 | oracle s(X1..X3) | benchmark |
| F1 | full s(X1..X7) (C≠C0) | implicit selection |
| F2 | full + **free** s(u,v) | the trap |
| F3r | full + spx(kspec=100), μ_s **REML** | documents the failure of automatic tuning |
| F3muX | μ_s fixed ∈ {10⁻²,…,10⁴} (6 pts) | trajectory, plateau rule, Δπ̄≤δ post hoc |
| F0-CF, F2-CF, F3-CF(μ=10⁴) | **cross-fitted** residuals | estimation-error channel; F3-CF = complete remedy |
| F4 (S2 only) | **per-term** spx (X2, X4–X6) + surface | counter-example: RSR-type pathology (spurious ρ̂ +0.42 at ρ₀=0) |

Cross-fitting: "scattered blocks" folds (50 kmeans cells → 5 folds
round-robin), `crossfit` option of `GAM_SAR_CFE`; the final interpretable
refit remains full-sample.

## BSPA arms (`BSPA_SAR_CFE`)

F0, F1 (C≠C0), F2 (+`bspatial(x,y)`), mstop=100 (default); on S2×C2:
F2 at mstop 400 and 1000 (is the immunity observed in the pilot
underfitting?).

## Outputs per fit

ρ̂, root-finding method, **var_rho** (closed form, for size/power of the
test of ρ=0), π̄ and π_max (k=25), MISE of f̂₁–f̂₃ (`dgp` truths),
edf of s(u,v), μ_s, BSPA selection, realized SNR, time, error.

## Decisions from the pilot and the mini-experiments

- kspec = 100 (25 leaves mid-frequency leakage: bias −0.089 →
  −0.037); justification 25/50/100 in the appendix.
- Absolute δ rule abandoned (π̄ floor ≈ 0.23 carried by the smoothers of
  legitimate covariates); replaced by the μ_s grid + plateau rule
  ρ̂(μ_s), evaluated post hoc per replication.
- Over-constraining the covariates is a counter-example, not a remedy (arm F4):
  π̄ is not an objective function.
- Cross-fitted CFE: the partialling-out (Robinson) score is Neyman-orthogonal;
  out-of-fold residuals remove the "own-observation" channel. Mini-exp:
  F2 bias −0.122 → −0.001 at ρ₀=0; combined spx+CF ≈ oracle at both ρ₀,
  MISE f̂₂ 0.06. Theory under spatial dependence → paper 2.
- Headline result to highlight: **the test of ρ=0 is invalid under
  free s(u,v)** and cross-fitted CFE restores the size.

## Agreed extensions (2026-11-05, after freeze)

- **ρ₀ = 0.9**: grid completed for continuity with E1–E4 (job `e3rho9`,
  same 33 cells, nr=1000).
- **Sub-axis n ∈ {1000, 5000}**: diagnostic cells S0×C0 and S2×C2,
  ρ₀ ∈ {0 ; 0.6}, nr=250 (job `e3n5000`, `CELLS=`/`N=`/`OUT=` of the sbatch) —
  checks the decay of the trap in K/n predicted by the fixed-dictionary
  regime.
- **mstop-CV: dedicated sub-experiment, outside the factorial** (special
  treatment — mstop is the complexity dial of boosting, the analogue
  of μ_s: tuning by a predictive criterion is structurally
  indifferent to the ρWy ↔ smoothers allocation, same pathology expected as
  REML for μ_s). Cells S0×C0 and S2×C2, ρ₀ ∈ {0 ; 0.6}, arm
  BSPA-F2: mstop ∈ {100, 400, 1000, CV random folds, CV spatial
  folds}, nr=250. If CV erodes the immunity, the principled remedy is
  cross-fitted BSPA (same "correct the score" angle), to be implemented
  then.

## Rerun of E1 / E2(i) — restricted scope (decided 2026-11-05)

Objective: a single version stamp (spboost 0.7.0, the HPC cluster pipeline, archived
seeds) and harmonized metrics (rej0, coverage) on EVERYTHING the
paper shows — and nothing more. The former E2(ii)-(iii)/E3/E4 are
superseded by E3′ and are not rerun.

- **E1 detail** (tab_E1_1/_2): n=1000, SNR=0.9, ρ/λ ∈ {0; 0.2; 0.6; 0.9},
  7 SAR estimators (gam, gamboost, lagsarlm, BSPA_SAR_ML/CFE,
  GAM_SAR_ML/CFE) + SEM analogues. **The ML estimators run only at n=1000**:
  their absence at large n is the paper's message (cost), documented by
  the comptime table.
- **E1 summary + coverage** (tab_E1_synth/_SEM, tab_coverage_*):
  n ∈ {1000, 2000, 5000} × SNR ∈ {0.7; 0.9}, **CFE only**
  (BSPA_*_CFE, GAM_*_CFE), var_rho in output for coverage.
- **E2(i)** (tab_E2_1): n=2000, SNR=0.7, X_het=TRUE, 3 estimators
  (BSPA_SAR_CFE, GAM_SAR_CFE, lagsarlm).
- nr=1000 throughout; total cost ≈ 60–80 core-h (~1 h on the HPC cluster).

## Execution

- the HPC cluster (slurm): `sbatch run_e3prime.sbatch` (array 1-1000%100, single-core; ≤120 cores).
- the HPC cluster (SGE): `qsub run_e3prime.sge` (backup).
- Aggregation: `OUT=out_full Rscript aggregate_e3prime.R`.
- Restart: already computed pairs (RDS present) are skipped.
