# Athaliana_Natural_Selection_Conflict

Code and data repository for **Ruffley, Lutz et al.** — natural selection on
flowering time and water-use efficiency (δ¹³C) in *Arabidopsis thaliana*
(Nature Ecology & Evolution resubmission).

The master analysis file is **`AtlasofPhenotypes.Rmd`**, which contains R code
for every analysis and figure referenced in the main text and the supplement.
All paths in the code are relative to this folder (e.g. `./data/...`).

## Top-level files

| File | Contents |
|---|---|
| `AtlasofPhenotypes.Rmd` | Master R Markdown: all analyses and figures (main text + supplement). |
| `AtlasofPhenotypes_selection-differentials.html` | Rendered HTML output of the selection-differential sections of the Rmd. |
| `Review.Rmd` | Statistics recomputed for the NEE revision (flowering time × δ¹³C correlations, slopes, covariances) — written to answer Reviewer 2 comment 1 and Reviewer 1's points on Fig. 3B/3C. |
| `SuppTables_RuffleyLutz.xlsx` | Current supplementary tables workbook. |
| `Ruffley_etal_2023.pdf` | PDF of the earlier (2023) version of the manuscript. |
| `Rplots.pdf` | Stray plot output from script runs; safe to ignore/regenerate. |

## Folders

| Folder | Contents |
|---|---|
| `data/` | All input files required by `AtlasofPhenotypes.Rmd`: phenotype matrices, selection-analysis results, G matrix, FLC mutant/drought experiment data, etc. **Has its own `README.md`** with per-file provenance. Note: `.rda` files here are RDS-format — read with `readRDS()`, not `load()`. |
| `data-raw/` | Raw lookup tables: the full *A. thaliana* accession list and the phenotype-name mapping (study ID → trait name) for the phenotype atlas. |
| `analysesphenotypeselection/R/` | `phenoselection_multi_FUNCTIONS.R` — the canonical function library for the multivariate phenotype-selection analyses, sourced by the Rmd. |
| `analyses/` | Standalone helper scripts: a copy of the phenotype-selection functions and QQ-plot code (by Matthew Flickinger) used for GWAS diagnostics. |
| `climate/` | Climate data for the accessions: `2029gclimate.csv` and the WorldClim2 per-accession extraction (`worldclim2/`). |
| `multivarGWAS/` | GEMMA multivariate GWAS (mvLMM) runs. `all_impmGWAS_logoutput/` = log files for all 66 pairwise runs on imputed phenotypes (source of the genetic correlations in the G matrix); `all_mGWAS_logOutput/` = 35 non-imputed runs; `mGWAS_FT16_Delta_13C/` and `mGWAS_FT16_Growth_rate/` = full run directories with launch scripts. Large `.assoc.txt` outputs and genotype files are gitignored (see `.gitignore` and `data/README.md`). |
| `figs/` | Saved ggplot objects (`.rda`) used to assemble the composite figures; `tmpobjects/` holds intermediate plot objects. |
| `tables/` | Generated tables underlying the supplement: phenotype info, missingness/imputation accuracy, selection estimates per trait, mGWAS top hits and GO annotations, BSLMM parameter estimates, etc. |
| `regenerated-tables/` | 2026 regeneration of supplementary Tables S21/S28/S29 with a consistently built G matrix (fixes the non-PSD issue in the deposited G). **Has its own `README.md`** with the method and per-file status. |
| `Data_code_Ruffley_Lutz_et_al/` | Self-contained data + Jupyter notebook (`Figures_5_6_...ipynb`) for Figures 5–6: the FLC CRISPR/mutant experiment (δ¹³C measurements, line annotations, FLC coverage analysis). |
| `Ruffleyetal2023_old/` | Supplementary tables/figures from the previous (2023) version of the manuscript, which did not include the *flc* mutant-line analyses. |
