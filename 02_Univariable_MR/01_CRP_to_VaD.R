# ============================================================
# CRP -> Vascular Dementia
# Two-sample Mendelian randomization analysis
#
# Primary exposure:
#   C-reactive protein (CRP)
#
# Outcomes:
#   1. Overall vascular dementia
#   2. Multi-infarct vascular dementia
#   3. Mixed vascular dementia
#   4. Subcortical vascular dementia
#   5. Undefined vascular dementia
#
# Primary MR method:
#   Random-effects inverse-variance weighted MR
#
# Complementary MR methods:
#   Weighted median
#   Weighted mode
#   MR-Egger
#   MR-RAPS
#
# Additional analysis:
#   Steiger directionality test
#   Overall R2 and F statistic
#
# Sensitivity analyses including MR-PRESSO, heterogeneity,
# horizontal pleiotropy, leave-one-out analysis, and related
# plots are handled separately in the sensitivity-analysis
# module.
#
# Input:
#   data/processed/CRP/CRP_instruments_Ffiltered.csv
#   data/raw/Outcome_VD/
#
# Output:
#   results/CRP_to_VaD/
# ============================================================

library(here)
library(vroom)
library(dplyr)
library(TwoSampleMR)
library(mr.raps)

# ------------------------------------------------------------
# 1. Define directories
# ------------------------------------------------------------

exposure_file <- here::here(
  "data",
  "processed",
  "CRP",
  "CRP_instruments_Ffiltered.csv"
)

outcome_dir <- here::here(
  "data",
  "raw",
  "Outcome_VD"
)

result_dir <- here::here(
  "results",
  "CRP_to_VaD"
)

harmonised_dir <- file.path(
  result_dir,
  "harmonised_data"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  harmonised_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Analysis parameters
# ------------------------------------------------------------

p_threshold <- 5e-8
F_threshold <- 10

# Five prespecified CRP -> VaD outcome tests
bonferroni_threshold <- 0.05 / 5

# ------------------------------------------------------------
# 3. Read final CRP instruments
# ------------------------------------------------------------

if (!file.exists(exposure_file)) {
  stop(
    "CRP instrument file not found: ",
    exposure_file
  )
}

exposure_dat <- read.csv(
  exposure_file,
  stringsAsFactors = FALSE
)

if (nrow(exposure_dat) < 2) {
  stop(
    "Fewer than 2 CRP instruments are available."
  )
}

# ------------------------------------------------------------
# 4. Extract exposure sample size
# ------------------------------------------------------------

N_values <- unique(
  exposure_dat$samplesize.exposure[
    !is.na(exposure_dat$samplesize.exposure)
  ]
)

if (length(N_values) == 0) {
  
  exposure_N <- NA_real_
  
} else {
  
  exposure_N <- as.numeric(
    N_values[1]
  )
}

message(
  "CRP exposure sample size: ",
  exposure_N
)

# ------------------------------------------------------------
# 5. Define vascular dementia outcomes
# ------------------------------------------------------------

outcome_traits <- data.frame(
  
  Trait = c(
    "F5_VASCDEM",
    "VD_MI",
    "VD_MX",
    "VD_SC",
    "VD_U"
  ),
  
  File = c(
    "finngen_R12_F5_VASCDEM.gz",
    "finngen_R12_VD_MI.gz",
    "finngen_R12_VD_MX.gz",
    "finngen_R12_VD_SC.gz",
    "finngen_R12_VD_U.gz"
  ),
  
  Case = c(
    3624,
    620,
    531,
    1031,
    184
  ),
  
  Control = c(
    475484,
    474702,
    474702,
    474702,
    474702
  ),
  
  stringsAsFactors = FALSE
) %>%
  mutate(
    SampleSize = Case + Control
  )

# ------------------------------------------------------------
# 6. Load outcome GWAS data
# ------------------------------------------------------------

outcome_list <- vector(
  "list",
  nrow(outcome_traits)
)

names(outcome_list) <- outcome_traits$Trait

for (j in seq_len(nrow(outcome_traits))) {
  
  trait <- outcome_traits$Trait[j]
  
  outcome_file <- file.path(
    outcome_dir,
    outcome_traits$File[j]
  )
  
  if (!file.exists(outcome_file)) {
    stop(
      "Outcome file not found: ",
      outcome_file
    )
  }
  
  outcome_list[[trait]] <- vroom(
    outcome_file,
    col_names = TRUE,
    show_col_types = FALSE
  ) %>%
    mutate(
      samplesize =
        outcome_traits$SampleSize[j]
    )
  
  message(
    "Loaded outcome: ",
    trait,
    " | N = ",
    outcome_traits$SampleSize[j]
  )
}

# ------------------------------------------------------------
# 7. Initialize result containers
# ------------------------------------------------------------

mr_results_all <- list()
steiger_results_all <- list()

# ------------------------------------------------------------
# 8. Main CRP -> VaD MR loop
# ------------------------------------------------------------

for (trait in names(outcome_list)) {
  
  message(
    "\n===================================================="
  )
  
  message(
    "CRP -> ",
    trait
  )
  
  message(
    "===================================================="
  )
  
  # ----------------------------------------------------------
  # 8.1 Format outcome data
  # ----------------------------------------------------------
  
  out_dat <- tryCatch(
    
    format_data(
      outcome_list[[trait]],
      type = "outcome",
      snps = exposure_dat$SNP,
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "sebeta",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      eaf_col = "af_alt",
      pval_col = "pval",
      samplesize_col = "samplesize"
    ),
    
    error = function(e) {
      
      warning(
        "Outcome formatting failed for ",
        trait,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(out_dat) || nrow(out_dat) < 2) {
    
    message(
      "Not enough outcome SNPs for ",
      trait
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 8.2 Harmonise exposure and outcome
  # ----------------------------------------------------------
  
  dat <- tryCatch(
    
    harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = out_dat,
      action = 2
    ),
    
    error = function(e) {
      
      warning(
        "Harmonisation failed for ",
        trait,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(dat)) {
    next
  }
  
  # ----------------------------------------------------------
  # 8.3 Save harmonised dataset
  # ----------------------------------------------------------
  
  harmonised_file <- file.path(
    harmonised_dir,
    paste0(
      "harmonised_CRP_to_",
      trait,
      ".csv"
    )
  )
  
  write.csv(
    dat,
    harmonised_file,
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 8.4 Retain valid instruments
  # ----------------------------------------------------------
  
  valid_dat <- dat %>%
    filter(
      mr_keep == TRUE
    )
  
  if (nrow(valid_dat) < 2) {
    
    message(
      "Fewer than 2 valid instruments after harmonisation: ",
      trait
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 8.5 Steiger directionality test
  #
  # The result is reported but is not used as a hard
  # exclusion criterion.
  # ----------------------------------------------------------
  
  steiger_result <- tryCatch(
    
    directionality_test(valid_dat),
    
    error = function(e) {
      
      warning(
        "Steiger test failed for ",
        trait,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (!is.null(steiger_result)) {
    
    steiger_df <- data.frame(
      
      exposure = "CRP",
      
      outcome = trait,
      
      nsnp = nrow(valid_dat),
      
      snp_r2_exposure =
        steiger_result$snp_r2.exposure,
      
      snp_r2_outcome =
        steiger_result$snp_r2.outcome,
      
      correct_causal_direction =
        steiger_result$correct_causal_direction,
      
      steiger_pval =
        steiger_result$steiger_pval,
      
      stringsAsFactors = FALSE
    )
    
    steiger_results_all[[
      length(steiger_results_all) + 1
    ]] <- steiger_df
  }
  
  # ----------------------------------------------------------
  # 8.6 Primary and complementary MR
  # ----------------------------------------------------------
  
  method_list <- "mr_ivw"
  
  if (nrow(valid_dat) >= 3) {
    
    method_list <- c(
      method_list,
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_egger_regression"
    )
  }
  
  res <- tryCatch(
    
    mr(
      valid_dat,
      method_list = method_list
    ),
    
    error = function(e) {
      
      warning(
        "MR analysis failed for ",
        trait,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(res)) {
    next
  }
  
  # ----------------------------------------------------------
  # 8.7 MR-RAPS
  # ----------------------------------------------------------
  
  if (nrow(valid_dat) >= 3) {
    
    raps_result <- tryCatch(
      
      mr.raps(
        b_exp = valid_dat$beta.exposure,
        b_out = valid_dat$beta.outcome,
        se_exp = valid_dat$se.exposure,
        se_out = valid_dat$se.outcome
      ),
      
      error = function(e) {
        
        warning(
          "MR-RAPS failed for ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (!is.null(raps_result)) {
      
      raps_df <- data.frame(
        
        id.exposure =
          unique(valid_dat$id.exposure)[1],
        
        id.outcome =
          unique(valid_dat$id.outcome)[1],
        
        exposure = "CRP",
        
        outcome = trait,
        
        method = "MR-RAPS",
        
        nsnp = nrow(valid_dat),
        
        b = raps_result$beta.hat,
        
        se = raps_result$beta.se,
        
        pval = raps_result$beta.p.value,
        
        stringsAsFactors = FALSE
      )
      
      res <- bind_rows(
        res,
        raps_df
      )
    }
  }
  
  # ----------------------------------------------------------
  # 8.8 Generate odds ratios
  # ----------------------------------------------------------
  
  res <- tryCatch(
    
    generate_odds_ratios(res),
    
    error = function(e) {
      
      warning(
        "Odds-ratio generation failed for ",
        trait,
        ": ",
        conditionMessage(e)
      )
      
      res
    }
  )
  
  # ----------------------------------------------------------
  # 8.9 Extract IVW result
  # ----------------------------------------------------------
  
  ivw_result <- res %>%
    filter(
      method == "Inverse variance weighted"
    )
  
  if (nrow(ivw_result) == 0) {
    
    warning(
      "IVW result unavailable for ",
      trait
    )
    
    next
  }
  
  ivw_pval <- ivw_result$pval[1]
  
  # ----------------------------------------------------------
  # 8.10 Calculate overall R2 and F statistic
  # ----------------------------------------------------------
  
  overall_R2 <- NA_real_
  overall_F <- NA_real_
  
  if (
    "R2" %in% names(valid_dat) &&
    !is.na(exposure_N)
  ) {
    
    overall_R2 <- sum(
      valid_dat$R2,
      na.rm = TRUE
    )
    
    k <- nrow(valid_dat)
    
    if (
      is.finite(overall_R2) &&
      overall_R2 > 0 &&
      overall_R2 < 1
    ) {
      
      overall_F <- (
        overall_R2 *
          (exposure_N - k - 1)
      ) /
        (
          k *
            (1 - overall_R2)
        )
    }
  }
  
  # ----------------------------------------------------------
  # 8.11 Add analysis metadata
  # ----------------------------------------------------------
  
  res <- res %>%
    mutate(
      
      exposure = "CRP",
      
      outcome = trait,
      
      overall_R2 =
        overall_R2,
      
      overall_F_statistic =
        overall_F,
      
      nominally_significant =
        pval < 0.05,
      
      bonferroni_significant =
        pval < bonferroni_threshold
    )
  
  # ----------------------------------------------------------
  # 8.12 Save individual MR result
  # ----------------------------------------------------------
  
  result_file <- file.path(
    result_dir,
    paste0(
      "CRP_to_",
      trait,
      "_MR.csv"
    )
  )
  
  write.csv(
    res,
    result_file,
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 8.13 Store result
  # ----------------------------------------------------------
  
  result_key <- paste(
    "CRP",
    trait,
    sep = "__"
  )
  
  mr_results_all[[result_key]] <- res
  
  message(
    "IVW P = ",
    format.pval(
      ivw_pval,
      digits = 3
    ),
    " | Bonferroni threshold = ",
    format(
      bonferroni_threshold,
      scientific = TRUE,
      digits = 3
    )
  )
  
  # ----------------------------------------------------------
  # 8.14 Release temporary objects
  # ----------------------------------------------------------
  
  rm(
    out_dat,
    dat,
    valid_dat,
    res
  )
  
  gc()
}

# ------------------------------------------------------------
# 9. Combine MR results
# ------------------------------------------------------------

if (length(mr_results_all) > 0) {
  
  mr_results_combined <- bind_rows(
    mr_results_all
  )
  
  write.csv(
    mr_results_combined,
    file.path(
      result_dir,
      "CRP_to_VaD_all_MR_results.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 9.1 IVW results
  # ----------------------------------------------------------
  
  ivw_results <- mr_results_combined %>%
    filter(
      method == "Inverse variance weighted"
    ) %>%
    arrange(
      pval
    )
  
  write.csv(
    ivw_results,
    file.path(
      result_dir,
      "CRP_to_VaD_IVW_results.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 9.2 Nominally significant results
  # ----------------------------------------------------------
  
  nominal_results <- ivw_results %>%
    filter(
      pval < 0.05
    )
  
  write.csv(
    nominal_results,
    file.path(
      result_dir,
      "CRP_to_VaD_nominally_significant.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 9.3 Bonferroni-significant results
  # ----------------------------------------------------------
  
  bonferroni_results <- ivw_results %>%
    filter(
      pval < bonferroni_threshold
    )
  
  write.csv(
    bonferroni_results,
    file.path(
      result_dir,
      "CRP_to_VaD_Bonferroni_significant.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 10. Save Steiger results
# ------------------------------------------------------------

if (length(steiger_results_all) > 0) {
  
  steiger_results <- bind_rows(
    steiger_results_all
  )
  
  write.csv(
    steiger_results,
    file.path(
      result_dir,
      "CRP_to_VaD_Steiger_results.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 11. Save analysis parameters
# ------------------------------------------------------------

analysis_parameters <- data.frame(
  
  parameter = c(
    "genome_wide_significance",
    "LD_r2",
    "LD_window_kb",
    "F_threshold",
    "number_of_outcomes",
    "Bonferroni_threshold"
  ),
  
  value = c(
    p_threshold,
    0.001,
    10000,
    F_threshold,
    nrow(outcome_traits),
    bonferroni_threshold
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file.path(
    result_dir,
    "CRP_to_VaD_analysis_parameters.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Save session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    result_dir,
    "sessionInfo.txt"
  )
)

message(
  "\n===================================================="
)

message(
  "CRP -> vascular dementia MR analysis completed."
)

message(
  "===================================================="
)