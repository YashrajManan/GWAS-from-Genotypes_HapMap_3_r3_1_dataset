# GWAS from Genotypes — HapMap_3_r3_1

A full genome-wide association study on the classic HapMap tutorial dataset: quality control, population-structure
correction with PCA, and per-SNP association — in **R + PLINK**, with a **Python twin** that re-implements the
association by hand.

## Research question

Which SNPs are associated with case/control status, once the genotypes are cleaned and the two things that would
otherwise give false answers — **population structure** and **multiple testing** — are properly handled?

## Biological background

A **SNP** is a genome position where individuals carry different alleles; GWAS tests common variants (MAF above
~1%). Each person is coded **0/1/2** by their count of one allele (the *additive* model), and the trait is
regressed on that dosage — logistic regression for a case/control phenotype. Two ideas dominate: **population
structure** (ancestry differs in both allele frequency and disease rate, so it confounds association — corrected
by adding genotype **PCA** components as covariates), and **multiple testing** (testing ~10⁵–10⁶ SNPs forces the
genome-wide significance threshold to **5 × 10⁻⁸**). Because of **linkage disequilibrium**, a hit tags a *region*,
not a proven causal gene.

## Data

Classic HapMap tutorial genotypes — `HapMap_3_r3_1` from Marees et al. (2018),
[*A tutorial on conducting genome-wide association studies*](https://github.com/MareesAT/GWA_tutorial) — a PLINK
binary set (`.bed/.bim/.fam`) with a case/control phenotype. Fetched automatically by the scripts.

## Methods

Load PLINK binary set → **QC** (sample/SNP missingness, MAF, HWE) → **LD-prune → PCA** (ancestry covariates) →
**logistic association** with the PCs → repeat **without** the PCs to expose the confounding → **Manhattan + QQ**
plots and the **genomic inflation factor λ**.

| Stage | R + PLINK | Python twin |
|---|---|---|
| QC / prune / PCA | PLINK 1.9 | (reuses R output) |
| Association | PLINK `--logistic --covar` | `statsmodels` logistic, per SNP |
| Calibration | λ from `qchisq` | λ from `scipy.stats.chi2` |
| Plots | `qqman` (Manhattan/QQ) | `matplotlib` (chr22) |

## Key result

After QC and PC covariates the test is **well-calibrated (λ ≈ 1)**, and comparing λ **with vs without** the PCs
shows population structure inflating the uncorrected analysis — the central demonstration of why GWAS must
control for ancestry.

![Manhattan](results_R/manhattan.png)
![QQ](results_R/qq.png)

## Interpretation

On a teaching dataset the genome-wide signal is limited, so the value is the **method**: clean QC, ancestry
controlled, significance judged against 5 × 10⁻⁸. Any hit would tag an LD region rather than a single causal
variant. The Python twin reproduces the chr22 result independently, cross-validating the pipeline.

## Limitations

Teaching-scale dataset (limited real signal); common variants only; a hit is a locus, not a gene (LD); PCs
correct structure but not cryptic relatedness or batch effects. Rigorous next steps: mixed-model association
(e.g. BOLT-LMM/GCTA) and fine-mapping of any hits.

## Repository

```
gwas_hapmap.R        # R + PLINK pipeline (primary)
gwas_hapmap.ipynb    # Python twin (per-SNP logistic on chr22)
results_R/           # Manhattan + QQ + top_hits.csv
results_py/          # chr22 Manhattan + QQ + results
```

## Run

**R (RStudio):** install PLINK 1.9 (https://www.cog-genomics.org/plink/) and put `plink` on your PATH — or set
`Sys.setenv(PLINK = "C:/path/to/plink.exe")`. Open `gwas_hapmap.R`, set the working directory to the file
location, and source it. It downloads the data, runs QC → PCA → association, and writes the figures.

**Python:** run `gwas_hapmap.R` first (it produces `data_R/qc.*` and `data_R/pca.eigenvec`), then run
`gwas_hapmap.ipynb` top to bottom.
