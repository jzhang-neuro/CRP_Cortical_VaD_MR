############################################################
# MVMR analyses
#
# Each model includes:
#   Exposure 1: CRP
#   Exposure 2: one prespecified covariate
#   Outcome: vascular dementia
#
# Covariates:
#   1. Type 2 diabetes
#   2. Hyperlipidemia
#   3. Obesity
#   4. SBP
#   5. DBP
#   6. Education
#   7. Physical activity
#   8. Smoking
############################################################

library(here)
library(vroom)
library(dplyr)
library(TwoSampleMR)
library(MVMR)

############################################################
# 1. Project configuration
############################################################

crp_raw_file <- here(
  "data", "raw", "CRP",
  "CRP_with_EAF.tsv"
)

crp_instrument_file <- here(
  "data", "processed", "CRP",
  "CRP_Ffiltered.csv"
)

outcome_file <- here(
  "data", "raw", "Outcome_VD",
  "finngen_R12_VD_SC.gz"
)

############################################################
# 2. Prespecified covariates
############################################################

covariate_info <- tibble::tribble(
  
  ~name,                ~raw_file, ~instrument_file,
  
  "T2D",
  "...",
  "...",
  
  "Hyperlipidemia",
  "...",
  "...",
  
  "Obesity",
  "...",
  "...",
  
  "SBP",
  "...",
  "...",
  
  "DBP",
  "...",
  "...",
  
  "Education",
  "...",
  "...",
  
  "Physical_activity",
  "...",
  "...",
  
  "Smoking",
  "...",
  "..."
)

############################################################
# 3. Output directories
############################################################

output_dir <- here(
  "results",
  "Multivariable_MR"
)

harmonised_dir <- file.path(
  output_dir,
  "harmonised"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  harmonised_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 4. Read CRP data
############################################################

crp_raw <- vroom(
  crp_raw_file,
  show_col_types = FALSE
)

crp_instruments <- read.csv(
  crp_instrument_file
)

############################################################
# 5. Read VaD outcome
############################################################

outcome_raw <- vroom(
  outcome_file,
  show_col_types = FALSE
)

############################################################
# 6. Format CRP outcome/exposure information
############################################################

crp_instruments <- crp_instruments %>%
  mutate(
    exposure = "CRP",
    id.exposure = "CRP"
  )

############################################################
# 7. Loop over the eight covariates
############################################################

mvmr_results_all <- list()
strength_results_all <- list()

for (i in seq_len(nrow(covariate_info))) {
  
  covariate_name <- covariate_info$name[i]
  
  message(
    "\n========================================\n",
    "Running MVMR: CRP + ",
    covariate_name,
    " -> VD_SC\n",
    "========================================"
  )
  
  ##########################################################
  # 7.1 Read covariate GWAS
  ##########################################################
  
  covariate_raw <- vroom(
    covariate_info$raw_file[i],
    show_col_types = FALSE
  )
  
  covariate_instruments <- read.csv(
    covariate_info$instrument_file[i]
  )
  
  covariate_instruments <- covariate_instruments %>%
    mutate(
      exposure = covariate_name,
      id.exposure = covariate_name
    )
  
  ##########################################################
  # 7.2 Construct union of CRP and covariate instruments
  ##########################################################
  
  exposure_combined <- bind_rows(
    crp_instruments,
    covariate_instruments
  )
  
  union_snps <- exposure_combined %>%
    select(
      SNP,
      pval.exposure
    ) %>%
    group_by(SNP) %>%
    summarise(
      pval.exposure = min(
        pval.exposure,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      id.exposure = "MVMR_union"
    )
  
  ##########################################################
  # 7.3 LD clumping
  ##########################################################
  
  union_clumped <- clump_data(
    union_snps,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 5e-8,
    pop = "EUR"
  )
  
  clumped_snps <- union_clumped$SNP
  
  message(
    "Union SNPs after LD clumping: ",
    length(clumped_snps)
  )
  
  ##########################################################
  # 7.4 Extract CRP SNPs
  ##########################################################
  
  exp_crp <- crp_raw %>%
    filter(
      SNP %in% clumped_snps
    ) %>%
    transmute(
      SNP,
      CHR,
      effect_allele.exposure =
        toupper(effect_allele.exposure),
      other_allele.exposure =
        toupper(other_allele.exposure),
      eaf.exposure,
      beta.exposure,
      se.exposure,
      pval.exposure,
      samplesize.exposure,
      exposure = "CRP",
      id.exposure = "CRP"
    ) %>%
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  ##########################################################
  # 7.5 Extract covariate SNPs
  ##########################################################
  
  exp_covariate <- covariate_raw %>%
    
    # Adapt these variable names to the specific GWAS
    mutate(
      SNP = rsid,
      beta_calc = log(as.numeric(odds_ratio)),
      se_calc =
        (
          log(as.numeric(ci_upper)) -
            log(as.numeric(ci_lower))
        ) / (2 * 1.96)
    ) %>%
    
    filter(
      SNP %in% clumped_snps,
      !is.na(beta_calc),
      !is.na(se_calc),
      is.finite(beta_calc),
      is.finite(se_calc)
    ) %>%
    
    transmute(
      SNP,
      CHR = chromosome,
      effect_allele.exposure =
        toupper(effect_allele),
      other_allele.exposure =
        toupper(other_allele),
      eaf.exposure =
        effect_allele_frequency,
      beta.exposure =
        beta_calc,
      se.exposure =
        se_calc,
      pval.exposure =
        p_value,
      samplesize.exposure =
        n,
      exposure =
        covariate_name,
      id.exposure =
        covariate_name
    ) %>%
    
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  ##########################################################
  # 7.6 Extract VaD outcome
  ##########################################################
  
  outcome_dat <- outcome_raw %>%
    mutate(
      SNP = rsids
    ) %>%
    filter(
      SNP %in% clumped_snps
    ) %>%
    transmute(
      SNP,
      beta.outcome = beta,
      se.outcome = sebeta,
      effect_allele.outcome =
        toupper(alt),
      other_allele.outcome =
        toupper(ref),
      eaf.outcome =
        af_alt,
      pval.outcome =
        pval,
      outcome = "VD_SC",
      id.outcome =
        "FinnGen_R12_VD_SC"
    ) %>%
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  ##########################################################
  # 7.7 Keep SNPs available for all three datasets
  ##########################################################
  
  common_snps <- Reduce(
    intersect,
    list(
      exp_crp$SNP,
      exp_covariate$SNP,
      outcome_dat$SNP
    )
  )
  
  exp_crp <- exp_crp %>%
    filter(SNP %in% common_snps)
  
  exp_covariate <- exp_covariate %>%
    filter(SNP %in% common_snps)
  
  outcome_dat <- outcome_dat %>%
    filter(SNP %in% common_snps)
  
  message(
    "Common SNPs: ",
    length(common_snps)
  )
  
  ##########################################################
  # 7.8 Format exposure data
  ##########################################################
  
  exposure_dat <- bind_rows(
    exp_crp,
    exp_covariate
  )
  
  exposure_dat <- format_data(
    exposure_dat,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col =
      "effect_allele.exposure",
    other_allele_col =
      "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col =
      "samplesize.exposure",
    phenotype_col = "exposure",
    id_col = "id.exposure"
  )
  
  ##########################################################
  # 7.9 Format outcome data
  ##########################################################
  
  outcome_fmt <- format_data(
    outcome_dat,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col =
      "effect_allele.outcome",
    other_allele_col =
      "other_allele.outcome",
    eaf_col = "eaf.outcome",
    pval_col = "pval.outcome",
    phenotype_col = "outcome",
    id_col = "id.outcome"
  )
  
  ##########################################################
  # 7.10 Harmonisation
  ##########################################################
  
  mvdat <- mv_harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_fmt
  )
  
  ##########################################################
  # 7.11 Save harmonised data
  ##########################################################
  
  saveRDS(
    mvdat,
    file = file.path(
      harmonised_dir,
      paste0(
        "CRP_plus_",
        covariate_name,
        "_VD_SC_MVMR_harmonised.rds"
      )
    )
  )
  
  ##########################################################
  # 7.12 MVMR input
  ##########################################################
  
  MRMVInputObject <- mr_mvinput(
    bx = mvdat$exposure_beta,
    bxse = mvdat$exposure_se,
    by = mvdat$outcome_beta,
    byse = mvdat$outcome_se
  )
  
  ##########################################################
  # 7.13 IVW
  ##########################################################
  
  MRMVObject.IVW <- mr_mvivw(
    MRMVInputObject,
    model = "default",
    correl = FALSE,
    distribution = "normal",
    alpha = 0.05
  )
  
  ##########################################################
  # 7.14 MR-Egger
  ##########################################################
  
  MRMVObject.Egger <- mr_mvegger(
    MRMVInputObject,
    orientate = 1,
    correl = FALSE,
    distribution = "normal",
    alpha = 0.05
  )
  
  ##########################################################
  # 7.15 Weighted median
  ##########################################################
  
  MRMVObject.Median <- mr_mvmedian(
    MRMVInputObject,
    distribution = "normal",
    alpha = 0.05,
    iterations = 10000,
    seed = 314159265
  )
  
  ##########################################################
  # 7.16 Extract results
  ##########################################################
  
  extract_mvmr_result <- function(
    obj,
    method,
    exposure_name,
    covariate_name) {
    
    data.frame(
      Model = method,
      Exposure1 = exposure_name,
      Exposure2 = covariate_name,
      Outcome = "VD_SC",
      Estimate = as.numeric(obj@Estimate),
      StdError = as.numeric(obj@StdError),
      CILower = as.numeric(obj@CILower),
      CIUpper = as.numeric(obj@CIUpper),
      Pvalue = as.numeric(obj@Pvalue),
      SNPs = as.numeric(obj@SNPs),
      CondFstat = as.numeric(obj@CondFstat),
      OR = exp(as.numeric(obj@Estimate)),
      OR_Lower = exp(as.numeric(obj@CILower)),
      OR_Upper = exp(as.numeric(obj@CIUpper))
    )
  }
  
  ivw_result <- extract_mvmr_result(
    MRMVObject.IVW,
    "IVW",
    "CRP",
    covariate_name
  )
  
  egger_result <- extract_mvmr_result(
    MRMVObject.Egger,
    "MR-Egger",
    "CRP",
    covariate_name
  )
  
  median_result <- extract_mvmr_result(
    MRMVObject.Median,
    "Weighted median",
    "CRP",
    covariate_name
  )
  
  combined_results <- bind_rows(
    ivw_result,
    egger_result,
    median_result
  )
  
  ##########################################################
  # 7.17 Save results for this covariate
  ##########################################################
  
  write.csv(
    combined_results,
    file = file.path(
      output_dir,
      paste0(
        "CRP_plus_",
        covariate_name,
        "_VD_SC_MVMR_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  mvmr_results_all[[covariate_name]] <-
    combined_results
  
  ##########################################################
  # 7.18 Conditional F-statistics
  ##########################################################
  
  F_data <- format_mvmr(
    BXGs = mvdat$exposure_beta,
    BYG = mvdat$outcome_beta,
    seBXGs = mvdat$exposure_se,
    seBYG = mvdat$outcome_se,
    RSID = mvdat$exp_dat$SNP
  )
  
  strength_result <- strength_mvmr(
    r_input = F_data,
    gencov = 0
  )
  
  strength_result$Exposure1 <- "CRP"
  strength_result$Exposure2 <- covariate_name
  
  write.csv(
    strength_result,
    file = file.path(
      output_dir,
      paste0(
        "CRP_plus_",
        covariate_name,
        "_Conditional_F_statistics.csv"
      )
    ),
    row.names = FALSE
  )
  
  strength_results_all[[covariate_name]] <-
    strength_result
}

############################################################
# 8. Combine all MVMR results
############################################################

all_mvmr_results <- bind_rows(
  mvmr_results_all
)

write.csv(
  all_mvmr_results,
  file = file.path(
    output_dir,
    "CRP_plus_covariates_VD_SC_MVMR_all_results.csv"
  ),
  row.names = FALSE
)

############################################################
# 9. Combine all conditional F-statistics
############################################################

all_strength_results <- bind_rows(
  strength_results_all
)

write.csv(
  all_strength_results,
  file = file.path(
    output_dir,
    "CRP_plus_covariates_VD_SC_Conditional_F_statistics_all.csv"
  ),
  row.names = FALSE
)

message(
  "\nAll MVMR analyses completed."
)
