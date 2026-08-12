# Regenerated supplementary tables (2026-08-12)

Produced by `regenerate_S21_S28_S29.R` (run from the repo root) with the
**consistently built G**: diagonal = univariate GEMMA LMM PVE, off-diagonals =
pairwise mvLMM genetic correlations r_g = Vg21/√(Vg11·Vg22) from the 66 logs in
`multivarGWAS/all_impmGWAS_logoutput/`, projected to the nearest valid
correlation matrix (`Matrix::nearPD`, max entry change 0.11), then
G = D_PVE^½ · R_g · D_PVE^½. Decision logged in
`../../CHANGELOG.md` (2026-08-12). The deposited
`data/Gmatrix_imp_multivar_corr_matrix.rda` mixed raw covariances with the PVE
diagonal, which produced the negative eigenvalues in the old Table S29.

## Status

| File | Status | Notes |
|---|---|---|
| `G_regenerated.rda/.csv` | **final** | PSD (two zero eigenvalues from the nearPD projection; none negative) |
| `S29_regenerated.csv` | **final** | eigendecomposition of the new G; replaces old S29 (no negative eigenvalues, no sampling-error footnote needed) |
| `S28_regenerated_mli.csv` | **final** | evolvability/respondability/angle/cosine from the **published S13 betas** (env `mli`) + new G. Survival 54.3°, Fitness 50.9°, Fecundity 65.9° |
| `S28_regenerated_allenvs.csv` | **final** | same, all four environments — pick which env(s) the supplement should report |
| `S21_regenerated_PROVISIONAL.csv` | **provisional** | Z = Gβ with 500-rep bootstrap, but computed from the repo's phenotype CSV, which is a *different export* than the Nov-2022 file behind published S13/S14 (betas correlate r ≈ 0.89, not 1.0) |
| `S13_validation.csv` | diagnostic | regenerated betas vs deposited S13/S14, per trait × fitness × env |

## To finalize S21

Fetch from the cluster: `~/moilab/projects/natvar/data/atlas1001_phenotype_matrix_imputed_withID.csv`
(**space-separated** — the repo's comma-separated file of the same name is a
different, later export). Re-run the script pointing at it; if the fresh betas
then match S13 exactly, the S21 output is final (drop the PROVISIONAL suffix).

## Notes for the supplement/workbook

- Old S28's "Genetic Variance Direction" column is dropped (it duplicated
  Evolvability); "Percent Selection Alignment" is renamed "Cosine Alignment"
  (it never was a percent). cos θ = Evolvability/Respondability identically.
- Old S28 published values (survival 0.291/0.412/45.0°) were never reproducible
  from deposited files (script was in Megan's folder) and rested on the
  inconsistent G; these regenerated values supersede them.
- S29 Top_Traits: consider adding signed loadings (audit note S29).
- Standardisation caveat for the supplement text: these are variance-standardised
  (SD-scale) quantities, not Hansen & Houle's mean-standardised evolvabilities —
  δ¹³C's negative mean makes mean-scaling meaningless.
