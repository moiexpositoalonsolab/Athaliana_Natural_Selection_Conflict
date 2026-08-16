# Data folder — provenance and completeness

Updated 2026-08-12. Files needed by `AtlasofPhenotypes.Rmd` (paths in the code are
relative to the repo root, e.g. `./data/...`). Note: files with an `.rda`
extension in this folder are **RDS-format** (read with `readRDS()`, not `load()`),
matching how the code reads them. CSV siblings of the small analysis objects are
provided for readability; the `.rda` files are the ones the code uses.

## Analysis objects (exported 2026-08-12 from `resubmissiona-nalyses/` rstudio exports)

| File | Contents |
|---|---|
| `Gmatrix_imp_multivar_corr_matrix.rda/.csv` | 12×12 "G matrix" used for Tables S21/S28/S29 and Fig. 2D. **Known construction issue** — off-diagonals are pairwise GEMMA mvLMM genetic covariances (log line 26, field 1 = Vg[2,1]) while the diagonal was overwritten with univariate PVE; the matrix is not positive semi-definite (two negative eigenvalues, cf. Table S29). See `TableS28_evolvability_explainer.Rmd` in the parent folder for full forensics. Built in `AtlasofPhenotypes.Rmd` (~lines 2930–3010, 5815–5870). |
| `multivariate_SelectionResults_targettraits_1122122.rda/.csv` | Direct (β) and total (s) selection estimates, 12 traits × 3 fitness measures × 4 environments (= Tables S13/S14). |
| `TargetPhenoMatrix.rda/.csv` | 515 accessions × 65 phenotype matrix for the field-experiment accessions. |
| `imputed_mutant_dC13.rda/.csv` | δ¹³C for FLC mutant/wt lines (306 × 5), CRISPR experiment. |
| `flc_droughtexp_data.rda/.csv` | FLC drought experiment data (47 × 5). |

## Uli Lutz line lists (added 2026-08-12, from Downloads / email; keep for posterity)

| File | Contents |
|---|---|
| `dLFC_T3_Cas9_transgene_homoz_UlrichLutz_01112019.xlsx` | Uli's Nov 2019 T3 list: 214 homozygous-Cas9 T3 *flc* lines across **80 accessions** — the full transformation effort, of which the paper's 62-accession panel is the sequenced subset. Sheet `ID_verification` holds the 2019 SNPmatch check of parent stocks (~20 label mismatches, e.g. 6680→6990, 7307→7316, 9536→9736, 9705→7186, 9928→9908). |
| `dFLC_list_lines_ulutz_042021.xlsx` | Uli's Apr 2021 curated list: 113 lines / 63 accessions with corrected accession names = the paper's 112-line panel (`lines_info_annot.csv`) plus ID9897-03-R04 (class 4, dropped). |

Why they matter: the growth-chamber drought experiment (Fig. 3C/D; sown Oct 2020,
`field_prep_year1/crispr/growth chamber/wrangled_data_crispr.csv`) uses **pre-correction**
accession labels. Its 47 *flc* entries = 36 panel accessions (5 of them under old names:
6680, 7307, 9536, 9705, 9928) + 9897 + 10 unsequenced T3 accessions (6040, 6177, 7008,
7130, 7353, 7416, 9525, 9706, 9826, 9935). Verify genotype-keyed joins (imputed δ¹³C,
Fig. 3D GWA) used corrected ids.

## Phenotype matrices (pre-existing)

`atlas1001_phenotype_matrix_imputed_withID.csv`, `atlas1001_phenotypes_matrix_MR.csv`,
`atlas_phenotype_matrix_imputedwithgeneticpcs.csv`, `atlas_phenotype_matrix_withid.csv`
(copied 2026-08-12 from the project root), `2029gaccession.csv`, `2029gclimate.csv`,
`imputed_wt&mut_dc13.csv`.

## GEMMA mvLMM provenance (in `../multivarGWAS/`)

Log files + run scripts for the two locally available pairwise runs
(`mGWAS_FT16_Delta_13C`, `mGWAS_FT16_Growth_rate`), with and without PC
covariates. The REMLE Vg blocks in these logs are the source format for the
G-matrix off-diagonals. The full `.assoc.txt` outputs (48–75 MB each) and the
genotype files (`1001gbi.bed` 3.3 GB, `.bim`, `.sXX.txt`) are **not** in the repo
— they remain in `natvar/resubmissiona-nalyses/` on the lab Shared Drive.

## Known missing (cluster only, `~/safedata/natvar/` at Berkeley) — fetch list

Needed to fully re-run or knit `AtlasofPhenotypes.Rmd`:

- `analyses/phenoselection_multi_FUNCTIONS-copy.R` — the `PHENOSELECTION()`
  function (breeder's-equation predictions behind Table S21 / Fig. 2D)
- `data/d4.rda` — field experiment data (515 accessions), read throughout
- `data/pheno_fromgoogle.tsv`, `data/atlas1001_phenotypes_matrix.csv`,
  `data/atlas1001_rawPheno_Quantile.tsv`, `data/atlas1001_imputedPheno_Quantile.tsv`
- `data/all_gwa_stats_table.rda` — univariate GWA summary (PVE diagonal source)
- `multivarGWAS/all_impmGWAS_logoutput/` — **all 66 pairwise imputed mvLMM logs;
  required to rebuild a valid G matrix** (extract full Vg per pair → r_g →
  assemble → scale by one consistent h² set)
- `figs/` and `figs/tmpobjects/` caches — required for the Rmd to knit as-is
  (~150 saved ggplot objects)
