# ============================================================
# CRP–VaD Genetic Correlation Using LDSC
# ============================================================
#
# Purpose:
#   Estimate the genetic correlation between genetically
#   predicted C-reactive protein (CRP) and a vascular dementia
#   (VaD) phenotype using LD score regression.
#
# Current analysis:
#   CRP <-> subcortical vascular dementia (VD_SC)
#
# Notes:
#   - CRP is treated as a continuous trait.
#   - VaD is treated as a binary trait and converted to the
#     liability scale using sample and population prevalence.
#   - HapMap3 SNPs and European LD scores are used.
#   - Cross-trait LDSC intercept is retained for QC.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(GenomicSEM)
library(here)


# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

# CRP sample size
CRP_N <- 575531

# Subcortical VaD sample information
VD_CASES <- 1031
VD_CONTROLS <- 474702

VD_SAMPLE_PREVALENCE <- VD_CASES /
  (VD_CASES + VD_CONTROLS)

# Prespecified population prevalence for subcortical VaD.
#
# IMPORTANT:
# This value must be justified by the study protocol,
# epidemiological evidence, or a prespecified assumption.
VD_POPULATION_PREVALENCE <- 0.002


# ------------------------------------------------------------
# 3. Project paths
# ------------------------------------------------------------

# Munged LDSC summary statistics
crp_sumstats <- here::here(
  "results",
  "LDSC",
  "munged",
  "CRP.sumstats.gz"
)

vd_sumstats <- here::here(
  "results",
  "LDSC",
  "munged",
  "VD_SC.sumstats.gz"
)

# European LD scores
ld_ref_dir <- here::here(
  "data",
  "reference",
  "eur_w_ld_chr"
)

# Output directory
output_dir <- here::here(
  "results",
  "LDSC"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Input validation
# ------------------------------------------------------------

required_files <- c(
  crp_sumstats,
  vd_sumstats
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following LDSC input files were not found:\n",
    paste(missing_files, collapse = "\n"),
    "\nPlease run the corresponding munging workflow first."
  )
}

if (!dir.exists(ld_ref_dir)) {
  stop(
    "LD score reference directory was not found: ",
    ld_ref_dir
  )
}


# ------------------------------------------------------------
# 5. Basic parameter checks
# ------------------------------------------------------------

if (VD_CASES <= 0 || VD_CONTROLS <= 0) {
  stop(
    "VD_CASES and VD_CONTROLS must both be positive."
  )
}

if (
  VD_POPULATION_PREVALENCE <= 0 ||
  VD_POPULATION_PREVALENCE >= 1
) {
  stop(
    "VD_POPULATION_PREVALENCE must be between 0 and 1."
  )
}

message(
  "Subcortical VaD sample prevalence: ",
  round(VD_SAMPLE_PREVALENCE, 6)
)

message(
  "Subcortical VaD population prevalence: ",
  VD_POPULATION_PREVALENCE
)


# ------------------------------------------------------------
# 6. LDSC
# ------------------------------------------------------------

res <- GenomicSEM::ldsc(
  
  traits = c(
    crp_sumstats,
    vd_sumstats
  ),
  
  sample.prev = c(
    NA,
    VD_SAMPLE_PREVALENCE
  ),
  
  population.prev = c(
    NA,
    VD_POPULATION_PREVALENCE
  ),
  
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  
  trait.names = c(
    "CRP",
    "VD_SC"
  ),
  
  stand = TRUE
)


# ------------------------------------------------------------
# 7. Extract genetic correlation
# ------------------------------------------------------------

rg_est <- res$S_Stand[1, 2]

rg_se <- sqrt(
  res$V_Stand[2, 2]
)

rg_z <- rg_est / rg_se

rg_p <- 2 * pnorm(
  -abs(rg_z)
)


# ------------------------------------------------------------
# 8. Extract cross-trait LDSC intercept
# ------------------------------------------------------------

cross_trait_intercept <- res$I[1, 2]

cross_trait_intercept_se <- sqrt(
  res$V[2, 2]
)


# ------------------------------------------------------------
# 9. Extract number of SNPs
# ------------------------------------------------------------

n_snp <- res$m


# ------------------------------------------------------------
# 10. Create results table
# ------------------------------------------------------------

result <- data.frame(
  
  Trait1 = "CRP",
  
  Trait2 = "VD_SC",
  
  rg = rg_est,
  
  SE = rg_se,
  
  Z = rg_z,
  
  P = rg_p,
  
  cross_trait_intercept =
    cross_trait_intercept,
  
  cross_trait_intercept_SE =
    cross_trait_intercept_se,
  
  N_cases = VD_CASES,
  
  N_controls = VD_CONTROLS,
  
  sample_prevalence =
    VD_SAMPLE_PREVALENCE,
  
  population_prevalence =
    VD_POPULATION_PREVALENCE,
  
  n_SNP = n_snp,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11. Save results
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "CRP_VD_SC_LDSC.csv"
)

write.csv(
  result,
  output_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Print results
# ------------------------------------------------------------

print(result)

message(
  "CRP–subcortical VaD LDSC analysis completed."
)

message(
  "Results saved to: ",
  output_file
)