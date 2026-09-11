# =============================================================================
# Genome-Wide Association Study (GWAS) on the HapMap_3_r3_1 tutorial dataset
# R + PLINK pipeline.  Companion Python twin: gwas_hapmap.ipynb
# =============================================================================
# Question: which SNPs associate with case/control status, after standard quality
# control and correcting for population structure (PCA) and multiple testing?
# Data: classic HapMap tutorial genotypes (HapMap_3_r3_1, Marees et al. 2018).
#
# Requires PLINK 1.9 (https://www.cog-genomics.org/plink/). Put `plink` on your PATH,
# OR set an env var with the full path, e.g. Sys.setenv(PLINK="C:/tools/plink/plink.exe").
# Run in RStudio: Session -> Set Working Directory -> To Source File Location -> Source.
# =============================================================================

## ---- Setup ----
setwd("D:/GIT HUB PROJECTS/gwas-from-genotypes")
getwd()
#if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
#  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")
library(qqman)
dir.create("data_R", showWarnings = FALSE); dir.create("results_R", showWarnings = FALSE)
PLINK <- Sys.getenv("PLINK", unset = "plink")     # PLINK 1.9 on PATH, or set full path
plink <- function(...) system2(PLINK, c(...))
stopifnot(nchar(Sys.which(PLINK)) > 0 || file.exists(PLINK))

## ---- 1. Get the dataset (PLINK binary set: .bed/.bim/.fam) ----
# .fam column 6 = phenotype (1 = control, 2 = case, -9 = missing).
if (!file.exists("data_R/HapMap_3_r3_1.bed")) {
  options(timeout = 600)
  zip <- "data_R/1_QC_GWAS.zip"
  download.file("https://github.com/MareesAT/GWA_tutorial/raw/master/1_QC_GWAS.zip", zip, mode = "wb")
  unzip(zip, exdir = "data_R")
  bed <- list.files("data_R", recursive = TRUE, pattern = "HapMap_3_r3_1\\.bed$", full.names = TRUE)[1]
  src <- dirname(bed)
  if (src != "data_R") file.copy(list.files(src, full.names = TRUE), "data_R", overwrite = TRUE)
}
fam <- read.table("data_R/HapMap_3_r3_1.fam")
cat("samples:", nrow(fam), "| phenotype table:\n"); print(table(fam$V6))

## ---- 2. Quality control ----
# Remove low-call-rate samples/SNPs, rare variants, and HWE-violating SNPs (genotyping errors).
plink("--bfile","data_R/HapMap_3_r3_1","--mind","0.02","--geno","0.02",
      "--maf","0.01","--hwe","1e-6","--make-bed","--out","data_R/qc")

## ---- 3. LD pruning -> PCA (ancestry covariates) ----
# Prune correlated SNPs so PCs reflect genome-wide ancestry; keep the top 10 PCs as covariates.
plink("--bfile","data_R/qc","--indep-pairwise","50","5","0.2","--out","data_R/prune")
plink("--bfile","data_R/qc","--extract","data_R/prune.prune.in","--pca","10","--out","data_R/pca")

## ---- 4. Association: logistic regression, WITH PC covariates ----
plink("--bfile","data_R/qc","--logistic","--covar","data_R/pca.eigenvec",
      "--hide-covar","--out","results_R/assoc")

## ---- 5. Association WITHOUT covariates (to expose the confounding) ----
plink("--bfile","data_R/qc","--logistic","--out","results_R/assoc_nocov")

## ---- 6. Genomic inflation lambda, Manhattan, QQ ----
# lambda = median(chi-square)/0.456; ~1 = well calibrated, >1 = residual structure.
lambda <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(1 - p, 1)) / qchisq(0.5, 1) }
res  <- read.table("results_R/assoc.assoc.logistic", header = TRUE); res  <- res[res$TEST == "ADD", ]
res0 <- read.table("results_R/assoc_nocov.assoc.logistic", header = TRUE); res0 <- res0[res0$TEST == "ADD", ]
cat(sprintf("lambda WITH PCs: %.3f | lambda NO PCs: %.3f\n", lambda(res$P), lambda(res0$P)))

res <- res[!is.na(res$P), ]
top <- head(res[order(res$P), ], 50)                       # commit only the top hits, not 1M rows
write.csv(top, "results_R/top_hits.csv", row.names = FALSE)

png("results_R/manhattan.png", width = 1600, height = 700, res = 150)
manhattan(res, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
          main = "GWAS Manhattan (HapMap, PC-corrected)",
          suggestiveline = FALSE, genomewideline = -log10(5e-8))
dev.off()

graphics.off()
pv <- res$P[is.finite(res$P) & res$P > 0]                  # qq() breaks on p==0
png("results_R/qq.png", width = 700, height = 700, res = 150)
qq(pv, main = sprintf("QQ (lambda = %.3f)", lambda(pv)))
dev.off()

## ---- Interpretation ----
# After QC + PC covariates the test is well-calibrated (lambda near 1); comparing lambda WITH vs
# WITHOUT the PCs shows population structure inflating the uncorrected analysis. On this teaching
# dataset the Manhattan is largely null (limited real signal) - the value is the correct method:
# clean QC, structure controlled, significance judged against 5e-8. A hit would tag an LD region,
# not a proven causal gene.
