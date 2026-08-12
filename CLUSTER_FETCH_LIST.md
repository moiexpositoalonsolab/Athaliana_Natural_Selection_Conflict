# Tier-1 fetch list for cluster Claude — 2026-08-12

All source paths are relative to `~/moilab/projects/natvar/` on the cluster.
Goal: tar these up (keeping relative paths) so they can be dropped into the
`Athaliana_Natural_Selection_Conflict` repo with the same structure. Everything
here should be small-to-medium (KB–tens of MB). Do NOT include `*.assoc.txt`,
`*.bed/.bim/.fam`, or `*.sXX.txt` — too large, not wanted.

Suggested command once verified:
`tar czvf natvar_tier1.tgz -C ~/moilab/projects/natvar -T filelist.txt` (plus the two
`~/moilab/projects/`-level items separately).

## 1. Highest priority — G-matrix rebuild (needed for Tables S21/S28/S29 fix)

- `multivarGWAS/all_impmGWAS_logoutput/` — whole folder (66 pairwise GEMMA mvLMM
  `.log.txt` files, imputed phenotypes). THE critical item.
- `multivarGWAS/all_mGWAS_logOutput/` — whole folder (non-imputed counterparts,
  for comparison), if present.

## 2. Analysis functions (needed to reproduce S21 / Fig 2D)

- `analyses/phenoselection_multi_FUNCTIONS-copy.R`  ← contains PHENOSELECTION()
- `analyses/QQPlotbyMatthewFlickinger.R`
- `analysesphenotypeselection/R/phenoselection_multi_FUNCTIONS.R` (if it exists)

## 3. data/ objects read by AtlasofPhenotypes.Rmd

- `data/d4.rda`                      ← field experiment data, used everywhere
- `data/pheno_fromgoogle.tsv`        ← trait metadata, used everywhere
- `data/phenotypenames.rda`
- `data/all_gwa_stats_table.rda`     ← univariate PVE table (G diagonal source)
- `data/atlas1001_phenotypes_matrix.csv`
- `data/atlas1001_rawPheno_Quantile.tsv`
- `data/atlas1001_imputedPheno_Quantile.tsv`   ← mvLMM input phenotypes
- `data/atlas_phenotype_matrix_imputedwithpcs.csv`
- `data/atlas_dir.tsv`
- `data/atlas_pheno_NewCategors.tsv`
- `data/atlas_phenotype_names_numaccessions_willannotatemanual_strategies_MR.tsv`
- `data/importanttraits.csv`
- `data/idex515_2029.rda`
- `data/mypalette.rda`
- `data/pcaTarget.rda`
- `data/preds.rda`
- `data/FTcor.rda`
- `data/FTvariance_df.rda`
- `data/phenofit.rda`
- `data/allphenotypes.rda`
- `data/imp_multivar_corr_matrix.Rda`          ← pre-diagonal-overwrite matrix!
- `data/multivar_corr_matrix.Rda`              ← non-imputed counterpart
- `data/target_Zcors.rda`
- `data/Gmatrix_Zcors.rda`
- `data/multivariate_SelectionResults_targettraits.rda`
- `data/multivariate_SelectionResults_NoGrowthRate.rda`
- `data/SelectionRepResults_AllPhenos_112122_100bsreps.rda`
- `data/SelectionRepResults_OtherPhenos_112122_100bsreps.rda`
- `data/PhenotypeSelectionSummary.rda`
- `data/stats_table_raw.rda`
- `data/stats_table_norm.rda`
- `data/stats_table_bslmm.rda`
- `data/bslmm.filt.dat.rda`
- `data/Zcor.rda`
- `data/Zcor_all.rda`
- `data/GenPhenCor_101822.rda`
- `data/GenPhenCor_all.rda`
- `data/List_PhenoCorrelations.rda`
- `data/Phenotype_Imp&Raw_TargetCorrelations.rda`   (note the & in the name)
- `data/fdr_growthchamber_exp.rda`
- `data/imputed_mutant_dC13_withSeeds.rda`
- `data/deltaC_extended_rep3.csv`
- `data/LD_peak1.rda`
- `data/gene_pos_raw.tsv`
- `data/mGWA_FT&dC13_andLFitnessGWA_top05hits.rda`  (note the &)
- `data/mGWA_FT&dC13_andLFitnessGWA_all.rda`        (note the &)
- `data/Arabidopsis_thaliana_world_accessions_list.tsv`

NOT needed (already in the repo): atlas1001_phenotype_matrix_imputed_withID.csv,
atlas1001_phenotypes_matrix_MR.csv, atlas_phenotype_matrix_withid.csv,
Gmatrix_imp_multivar_corr_matrix.rda, TargetPhenoMatrix.rda,
multivariate_SelectionResults_targettraits_1122122.rda, imputed_mutant_dC13.rda,
flc_droughtexp_data.rda, 2029gaccession.csv, 2029gclimate.csv.

## 4. tables/ — whole folder

- `tables/` — all TSV/CSV summary outputs (~40 small files referenced by the Rmd).

## 5. figs/ caches — required for the Rmd to knit

- `figs/**/*.rda` — ALL .rda under figs/ including `figs/tmpobjects/`
  (~150 saved ggplot objects). Only `.rda` files — skip pdf/png if easier.

## 6. Misc singles

- `hqsnps515nature.txt`
- `data-raw/allarabidopsisaccessions.csv`
- `data-raw/atlas_phenotype_names.csv`
- `climate/worldclim2/2029g_climate_accessions.csv`
- `MeganAnalysis/GwaCor_NormvImput.rda`
- `MeganAnalysis/GwaCor_RawvImput.rda`
- `MeganAnalysis/GwaCor_RawvNorm.rda`
- `LauraAnalysis/field_prep_year1/crispr/growth chamber/seed_weights.csv`
- `LauraAnalysis/field_prep_year1/crispr/growth chamber/wrangled_data_crispr.csv`
- `multivarGWAS/mGWAS_FT16_Delta_13C/peak1snps.txt` (have a copy; grab if easy)

One level UP from natvar (`~/moilab/projects/`):
- `~/moilab/projects/tables/mega_hypparams_table.rda`
- `~/moilab/projects/tables/mega_hypparams_table_2.rda`
- `~/moilab/projects/data/atlas1001_phenotype_matrix_imputed_onlypheno.rda`
- `~/moilab/projects/phenotypes/phenotypes_list.txt`

## 7. Bonus hunts while you're there (5 min, high value)

1. **Find the Table S28 script** — the evolvability/respondability/angle
   calculation. It is NOT in the repo or the deposited exports. Try:
   `grep -rln "respondability\|acos\|evolvab" ~/moilab/projects/natvar --include="*.R" --include="*.Rmd"`
   Also check `~/moilab/projects/natvar/analyses/` and shell history if accessible.
2. **Find what wrote Table S21** — any script calling PHENOSELECTION and writing
   `SelectionResponse`/`deltaZ`-style outputs, to pin down the exact call.
3. If found, include those scripts in the tarball.

## Notes for the cluster Claude

- Filenames contain `&` and spaces ("growth chamber") — quote paths.
- `.rda` files in this project are often RDS-format (readRDS, not load).
- If a listed file is missing, note it and move on — don't hunt substitutes.
- Report back: what was found, what was missing, tarball size.
