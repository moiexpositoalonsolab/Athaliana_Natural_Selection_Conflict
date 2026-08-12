#!/bin/bash
# Fetch the files AtlasofPhenotypes.Rmd needs from the cluster (~/moilab/projects/natvar)
# to make this repo self-contained. Generated 2026-08-12 from the Rmd's file
# references (see data/README.md).
#
# Usage:  bash fetch_from_cluster.sh <user@cluster-host>
# Runs rsync into the repo folder this script lives in. Re-runnable (rsync -a).
# TIER 2 (multi-GB GWAS outputs) is OFF by default — set FETCH_LARGE=1 to include.

set -euo pipefail
HOST="${1:?usage: bash fetch_from_cluster.sh user@cluster-host}"
BASE="~/moilab/projects/natvar"
DEST="$(cd "$(dirname "$0")" && pwd)"
RS="rsync -av --progress"

############ TIER 1 — essential, small-to-medium ############

# The PHENOSELECTION() function (breeder's equation, S21/Fig 2D) + helpers
$RS "$HOST:$BASE/analyses/phenoselection_multi_FUNCTIONS-copy.R" \
    "$HOST:$BASE/analyses/QQPlotbyMatthewFlickinger.R" \
    "$DEST/analyses/"
$RS "$HOST:$BASE/analysesphenotypeselection/R/phenoselection_multi_FUNCTIONS.R" \
    "$DEST/analysesphenotypeselection/R/" || true

# ALL 66 pairwise mvLMM logs — REQUIRED to rebuild a valid G matrix
$RS "$HOST:$BASE/multivarGWAS/all_impmGWAS_logoutput/" "$DEST/multivarGWAS/all_impmGWAS_logoutput/"
$RS "$HOST:$BASE/multivarGWAS/all_mGWAS_logOutput/"    "$DEST/multivarGWAS/all_mGWAS_logOutput/" || true

# data/ objects referenced by the Rmd (skip the ones already in the repo)
$RS --exclude="atlas1001_phenotype_matrix_imputed_withID.csv" \
    --exclude="atlas1001_phenotypes_matrix_MR.csv" \
    "$HOST:$BASE/data/d4.rda" \
    "$HOST:$BASE/data/pheno_fromgoogle.tsv" \
    "$HOST:$BASE/data/phenotypenames.rda" \
    "$HOST:$BASE/data/all_gwa_stats_table.rda" \
    "$HOST:$BASE/data/atlas1001_phenotypes_matrix.csv" \
    "$HOST:$BASE/data/atlas1001_rawPheno_Quantile.tsv" \
    "$HOST:$BASE/data/atlas1001_imputedPheno_Quantile.tsv" \
    "$HOST:$BASE/data/atlas_phenotype_matrix_imputedwithpcs.csv" \
    "$HOST:$BASE/data/atlas_dir.tsv" \
    "$HOST:$BASE/data/atlas_pheno_NewCategors.tsv" \
    "$HOST:$BASE/data/atlas_phenotype_names_numaccessions_willannotatemanual_strategies_MR.tsv" \
    "$HOST:$BASE/data/importanttraits.csv" \
    "$HOST:$BASE/data/idex515_2029.rda" \
    "$HOST:$BASE/data/mypalette.rda" \
    "$HOST:$BASE/data/pcaTarget.rda" \
    "$HOST:$BASE/data/preds.rda" \
    "$HOST:$BASE/data/FTcor.rda" \
    "$HOST:$BASE/data/FTvariance_df.rda" \
    "$HOST:$BASE/data/phenofit.rda" \
    "$HOST:$BASE/data/allphenotypes.rda" \
    "$HOST:$BASE/data/imp_multivar_corr_matrix.Rda" \
    "$HOST:$BASE/data/multivar_corr_matrix.Rda" \
    "$HOST:$BASE/data/target_Zcors.rda" \
    "$HOST:$BASE/data/Gmatrix_Zcors.rda" \
    "$HOST:$BASE/data/multivariate_SelectionResults_targettraits.rda" \
    "$HOST:$BASE/data/multivariate_SelectionResults_NoGrowthRate.rda" \
    "$HOST:$BASE/data/SelectionRepResults_AllPhenos_112122_100bsreps.rda" \
    "$HOST:$BASE/data/SelectionRepResults_OtherPhenos_112122_100bsreps.rda" \
    "$HOST:$BASE/data/PhenotypeSelectionSummary.rda" \
    "$HOST:$BASE/data/stats_table_raw.rda" \
    "$HOST:$BASE/data/stats_table_norm.rda" \
    "$HOST:$BASE/data/stats_table_bslmm.rda" \
    "$HOST:$BASE/data/bslmm.filt.dat.rda" \
    "$HOST:$BASE/data/Zcor.rda" \
    "$HOST:$BASE/data/Zcor_all.rda" \
    "$HOST:$BASE/data/GenPhenCor_101822.rda" \
    "$HOST:$BASE/data/GenPhenCor_all.rda" \
    "$HOST:$BASE/data/List_PhenoCorrelations.rda" \
    "$HOST:$BASE/data/Phenotype_Imp&Raw_TargetCorrelations.rda" \
    "$HOST:$BASE/data/fdr_growthchamber_exp.rda" \
    "$HOST:$BASE/data/imputed_mutant_dC13_withSeeds.rda" \
    "$HOST:$BASE/data/deltaC_extended_rep3.csv" \
    "$HOST:$BASE/data/LD_peak1.rda" \
    "$HOST:$BASE/data/gene_pos_raw.tsv" \
    "$HOST:$BASE/data/mGWA_FT&dC13_andLFitnessGWA_top05hits.rda" \
    "$HOST:$BASE/data/mGWA_FT&dC13_andLFitnessGWA_all.rda" \
    "$HOST:$BASE/data/Arabidopsis_thaliana_world_accessions_list.tsv" \
    "$DEST/data/"

# tables/ — small TSV/CSV outputs, referenced ~40x; take the whole folder
$RS "$HOST:$BASE/tables/" "$DEST/tables/"

# figs/ caches — REQUIRED for the Rmd to knit (saved ggplot .rda objects)
$RS --include="*/" --include="*.rda" --exclude="*" "$HOST:$BASE/figs/" "$DEST/figs/"

# misc single files
$RS "$HOST:$BASE/hqsnps515nature.txt" "$DEST/" || true
$RS "$HOST:$BASE/data-raw/allarabidopsisaccessions.csv" \
    "$HOST:$BASE/data-raw/atlas_phenotype_names.csv" "$DEST/data-raw/" || true
$RS "$HOST:$BASE/climate/worldclim2/2029g_climate_accessions.csv" "$DEST/climate/worldclim2/" || true
$RS "$HOST:$BASE/MeganAnalysis/GwaCor_NormvImput.rda" \
    "$HOST:$BASE/MeganAnalysis/GwaCor_RawvImput.rda" \
    "$HOST:$BASE/MeganAnalysis/GwaCor_RawvNorm.rda" "$DEST/MeganAnalysis/" || true
$RS "$HOST:$BASE/LauraAnalysis/field_prep_year1/crispr/growth chamber/seed_weights.csv" \
    "$HOST:$BASE/LauraAnalysis/field_prep_year1/crispr/growth chamber/wrangled_data_crispr.csv" \
    "$DEST/LauraAnalysis/field_prep_year1/crispr/growth chamber/" || true
# parent-relative references (../ from the cluster natvar dir = ~/safedata)
$RS "$HOST:~/moilab/projects/tables/mega_hypparams_table.rda" \
    "$HOST:~/moilab/projects/tables/mega_hypparams_table_2.rda" "$DEST/../tables-safedata/" || true
$RS "$HOST:~/moilab/projects/data/atlas1001_phenotype_matrix_imputed_onlypheno.rda" "$DEST/data/" || true
$RS "$HOST:~/moilab/projects/phenotypes/phenotypes_list.txt" "$DEST/phenotypes/" || true

############ TIER 2 — LARGE GWAS outputs (assoc.txt, ~50MB–2GB each) ############
# Only needed to re-run Manhattan plots / GWA-dependent chunks from scratch.
# Too large for GitHub — if fetched, keep out of git (see .gitignore).
if [ "${FETCH_LARGE:-0}" = "1" ]; then
  $RS "$HOST:$BASE/multivarGWAS/mGWAS_FT16_Delta_13C/output/"      "$DEST/multivarGWAS/mGWAS_FT16_Delta_13C/output/"
  $RS "$HOST:$BASE/multivarGWAS/mGWAS_FT16_Growth_rate/output/"    "$DEST/multivarGWAS/mGWAS_FT16_Growth_rate/output/"
  $RS "$HOST:$BASE/multivarGWAS/imp_mGWAS_FT16_Delta_13C/output/"  "$DEST/multivarGWAS/imp_mGWAS_FT16_Delta_13C/output/"
  $RS "$HOST:$BASE/multivarGWAS/imp_mGWAS_FT16_Growth_rate/output/" "$DEST/multivarGWAS/imp_mGWAS_FT16_Growth_rate/output/"
  $RS "$HOST:$BASE/multivarGWAS/covary_output/"                    "$DEST/multivarGWAS/covary_output/"
  $RS "$HOST:$BASE/gwaresults/"                                    "$DEST/gwaresults/"
  $RS "$HOST:$BASE/FitnessResiduals/output/"                       "$DEST/FitnessResiduals/output/"
  $RS "$HOST:$BASE/FLC_koData/output/"                             "$DEST/FLC_koData/output/"
  $RS "$HOST:$BASE/phenotypes/1001_Consortium_Cell_2016_PID_27293186/" "$DEST/phenotypes/1001_Consortium_Cell_2016_PID_27293186/"
  $RS "$HOST:$BASE/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/" "$DEST/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/"
fi

echo "Done. Check data/README.md for what each file is."
