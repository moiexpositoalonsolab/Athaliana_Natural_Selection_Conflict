## Regenerate Supplementary Tables S21, S28, S29 with a consistently built G.
##
## Decision (Moi, 2026-08-12): keep the original sourcing — diagonal = univariate
## GEMMA LMM PVE (larger per-trait n), off-diagonals = pairwise mvLMM — but the
## mvLMM enters as the genetic CORRELATION r_g = Vg21/sqrt(Vg11*Vg22), rescaled
## by the univariate PVE:   G = D_PVE^1/2 %*% nearPD(R_g) %*% D_PVE^1/2.
## (The deposited matrix pasted the raw mvLMM covariance next to the PVE
## diagonal, which is what produced the negative eigenvalues in old Table S29.)
##
## Run from the repo root:  Rscript regenerated-tables/regenerate_S21_S28_S29.R
## Inputs (all in this repo): multivarGWAS/all_impmGWAS_logoutput/ (66 logs),
##   data/Gmatrix_imp_multivar_corr_matrix.rda (for the PVE diagonal + trait order),
##   data/atlas1001_phenotype_matrix_imputed_withID.csv, data/d4.rda,
##   analyses/phenoselection_multi_FUNCTIONS-copy.R
## Outputs: regenerated-tables/{G_regenerated.rda,G_regenerated.csv,
##   S29_regenerated.csv, S28_regenerated_mli.csv, S28_regenerated_allenvs.csv,
##   S21_regenerated.csv, S13_validation.csv}

suppressMessages({ library(Matrix); library(boot) })
outdir <- "regenerated-tables"

## ---- 1. Build G from the 66 pairwise logs + univariate PVE diagonal ----
G_dep  <- readRDS("data/Gmatrix_imp_multivar_corr_matrix.rda")
traits <- rownames(G_dep)
pve    <- diag(G_dep)                       # univariate LMM PVE (decided diagonal)

logs <- list.files("multivarGWAS/all_impmGWAS_logoutput", pattern = "log.txt$",
                   full.names = TRUE)
stopifnot(length(logs) == 66)
Rg <- matrix(NA_real_, 12, 12, dimnames = list(traits, traits)); diag(Rg) <- 1
for (f in logs) {
  base <- sub("^imp_mGWAS_", "", sub("\\.log\\.txt$", "", basename(f)))
  t1 <- traits[startsWith(base, paste0(traits, "_"))][1]
  t2 <- substring(base, nchar(t1) + 2)
  stopifnot(t2 %in% traits)
  ll  <- readLines(f, warn = FALSE)
  i   <- grep("REMLE estimate for Vg", ll)[1]
  v11 <- as.numeric(strsplit(trimws(ll[i + 1]), "\\s+")[[1]][1])
  v2  <- as.numeric(strsplit(trimws(ll[i + 2]), "\\s+")[[1]])
  stopifnot(abs(G_dep[t1, t2] - v2[1]) < 1e-6)   # provenance guard: log == deposited
  Rg[t1, t2] <- Rg[t2, t1] <- v2[1] / sqrt(v11 * v2[2])
}
Rg_pd <- as.matrix(nearPD(Rg, corr = TRUE)$mat)   # 66 pairwise r_g are jointly non-PSD
D <- diag(sqrt(pve))
Gmatrix <- D %*% Rg_pd %*% D                      # global: responselinear() finds it here
dimnames(Gmatrix) <- list(traits, traits)
saveRDS(Gmatrix, file.path(outdir, "G_regenerated.rda"))
write.csv(as.data.frame(Gmatrix), file.path(outdir, "G_regenerated.csv"), row.names = TRUE)
cat("G eigenvalues:", round(eigen(Gmatrix, symmetric = TRUE)$values, 4), "\n")

## ---- 2. Table S29: eigendecomposition of G ----
ev <- eigen(Gmatrix, symmetric = TRUE)
top_traits <- apply(ev$vectors, 2, function(v)
  paste(traits[order(abs(v), decreasing = TRUE)][1:4], collapse = ", "))
S29 <- data.frame(Axis = paste0("g", 1:12),
                  Eigenvalue = round(ev$values, 8),
                  Prop_Variance = round(ev$values / sum(ev$values), 3),
                  Top_Traits = top_traits)
write.csv(S29, file.path(outdir, "S29_regenerated.csv"), row.names = FALSE)

## ---- 3. PHENOSELECTION with the new G (S21 + fresh betas) ----
source("analyses/phenoselection_multi_FUNCTIONS-copy.R")
fn <- function(x) as.numeric(as.character(x))          # moiR::fn equivalent

pheno <- read.csv("data/atlas1001_phenotype_matrix_imputed_withID.csv")
load("data/d4.rda")                                    # field experiment -> d4
idex515 <- which(pheno$id %in% unique(d4$id))
cat("field accessions:", length(idex515), "\n")

df <- pheno[idex515, traits]                           # columns in G's trait order
df <- apply(df, 2, fn); df <- apply(df, 2, scale)

fits <- c("rFitness", "rSurvival_fruit", "rSeeds")
envs <- c("mli", "mlp", "thi", "thp")
set.seed(0)
lsel <- NULL
for (i in fits) for (j in envs) {
  w <- pheno[[paste0(i, "_", j)]][idex515]
  res <- PHENOSELECTION(Variables = df, Fitness = w, Gmatrix = Gmatrix, replicates = 500)
  for (p in c("gradient_linear", "coefficient_linear", "response_linear")) {
    tmp <- do.call(cbind, parseformatted(res[[p]]))
    tmp <- data.frame(mean = fn(tmp[, 1]), se = fn(tmp[, 2]), signi = tmp[, 3],
                      param = c(gradient_linear = "beta", coefficient_linear = "s",
                                response_linear = "z")[p],
                      trait = traits, fitness = i, env = j)
    lsel <- rbind(lsel, tmp)
  }
  cat("done:", i, j, "\n")
}
lsel$signi[lsel$signi == ""] <- "ns"
lsel$lower <- lsel$mean - 1.96 * lsel$se
lsel$upper <- lsel$mean + 1.96 * lsel$se

## PROVISIONAL: the repo's phenotype CSV is a different export than the
## Nov-2022 space-separated file the published S13/S14 betas came from
## (point-estimate betas correlate r ~ 0.89 with S13, not 1.0). The exact
## original lives only on the cluster: data/atlas1001_phenotype_matrix_imputed_withID.csv
## (sep = " "). Fetch it, re-run, and if betas then match S13 exactly, drop
## the PROVISIONAL suffix.
S21 <- lsel[lsel$param == "z",
            c("trait", "fitness", "env", "mean", "se", "signi", "lower", "upper")]
names(S21) <- c("trait", "fitness", "env", "Z", "SE", "significance",
                "lower est", "upper est")
write.csv(S21, file.path(outdir, "S21_regenerated_PROVISIONAL.csv"), row.names = FALSE)

## ---- 4. Validation: fresh betas vs deposited Table S13/S14 ----
sel_dep <- readRDS("data/multivariate_SelectionResults_targettraits_1122122.rda")
map <- c(ABA = "ABA_96h_low_water_potential", Delta13C = "Delta_13C",
         Dormancy = "DSDS10", FloweringTime = "FT16", GerminationPerc = "d8_10C_perc",
         GrowthRate = "Growth_rate", RGR = "RGR",
         RootRGR = "Relative_root_growth_rate_day002.day003",
         RootHorizIndex = "Root_horizontal_index_day001",
         StomataDensity = "stomata_density", StomataSize = "Stomatal_index_in_first_leaf",
         Vernalization = "X72_Vern_Growth")
sel_dep$g_trait <- map[sel_dep$trait]
val <- merge(lsel[lsel$param == "beta", c("trait", "fitness", "env", "mean")],
             sel_dep[, c("g_trait", "fitness", "env", "mean.x")],
             by.x = c("trait", "fitness", "env"), by.y = c("g_trait", "fitness", "env"))
names(val)[4:5] <- c("beta_regenerated", "beta_S13_deposited")
write.csv(val, file.path(outdir, "S13_validation.csv"), row.names = FALSE)
cat(sprintf("beta validation vs S13/S14: r = %.4f, max |diff| = %.3f (n = %d)\n",
            cor(val$beta_regenerated, val$beta_S13_deposited),
            max(abs(val$beta_regenerated - val$beta_S13_deposited)), nrow(val)))

## ---- 5. Table S28: evolvability / respondability / alignment ----
## Uses the DEPOSITED S13/S14 betas (the published selection gradients), so S28
## is final regardless of the S21 provisional status above.
geom <- function(b) {
  b  <- b / sqrt(sum(b^2))
  dz <- as.vector(Gmatrix %*% b)
  e  <- sum(b * dz); r <- sqrt(sum(dz^2))
  c(Evolvability = e, Respondability = r,
    `Angle Degrees` = acos(e / r) * 180 / pi, `Cosine Alignment` = e / r)
}
S28 <- NULL
for (i in fits) for (j in envs) {
  bsub <- sel_dep[sel_dep$fitness == i & sel_dep$env == j, ]
  b <- setNames(bsub$mean.x, bsub$g_trait)[traits]
  S28 <- rbind(S28, data.frame(`Fitness Trait` = i, env = j, t(round(geom(b), 5)),
                               check.names = FALSE))
}
lab <- c(rFitness = "Fitness", rSurvival_fruit = "Survival", rSeeds = "Fecundity (seeds set)")
S28$`Fitness Trait` <- lab[S28$`Fitness Trait`]
write.csv(S28, file.path(outdir, "S28_regenerated_allenvs.csv"), row.names = FALSE)
S28_mli <- S28[S28$env == "mli", setdiff(names(S28), "env")]
S28_mli <- S28_mli[match(c("Survival", "Fitness", "Fecundity (seeds set)"),
                         S28_mli$`Fitness Trait`), ]
write.csv(S28_mli, file.path(outdir, "S28_regenerated_mli.csv"), row.names = FALSE)
print(S28_mli, row.names = FALSE)
cat("done.\n")
