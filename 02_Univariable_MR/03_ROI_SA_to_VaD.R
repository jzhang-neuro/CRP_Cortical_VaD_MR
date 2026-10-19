# ============================================================
# MR analysis of regional cortical surface area and vascular
# dementia
#
# Purpose:
#   Test the associations between 180 regional cortical
#   surface-area phenotypes and vascular dementia outcomes.
#
# Exposure:
#   Regional cortical surface area
#
# Outcomes:
#   1. Overall vascular dementia
#   2. Multi-infarct vascular dementia
#   3. Mixed vascular dementia
#   4. Subcortical vascular dementia
#   5. Undefined vascular dementia
#
# Main MR method:
#   Random-effects inverse-variance weighted MR
#
# Additional MR methods:
#   Weighted median
#   Weighted mode
#   MR-Egger
#   MR-RAPS
#
# Additional analyses:
#   Steiger directionality test
#   Overall R2 and F statistic
#
# Input:
#   data/processed/exp_dat_eaf_Ffilter_ROI_SA_CAMIDE/
#   data/raw/Outcome_VD/
#
# Output:
#   results/ROI_SA_to_VaD/
# ============================================================

library(here)
library(vroom)
library(dplyr)
library(TwoSampleMR)
library(mr.raps)

# ------------------------------------------------------------
# 1. Define directories
# ------------------------------------------------------------

exposure_dir <- here::here(
  "data",
  "processed",
  "exp_dat_eaf_Ffilter_ROI_SA_CAMIDE"
)

outcome_dir <- here::here(
  "data",
  "raw",
  "Outcome_VD"
)

result_dir <- here::here(
  "results",
  "ROI_SA_to_VaD"
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
# 2. Define outcome traits
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
)

outcome_traits <- outcome_traits %>%
  mutate(
    SampleSize = Case + Control
  )

# ------------------------------------------------------------
# 3. Identify exposure files
# ------------------------------------------------------------

exposure_files <- list.files(
  exposure_dir,
  pattern = "^exp_dat_.*_Ffiltered\\.csv$",
  full.names = TRUE
)

if (length(exposure_files) == 0) {
  stop(
    "No F-filtered exposure files were found."
  )
}

exposure_files <- sort(exposure_files)

message(
  "Number of ROI exposure files: ",
  length(exposure_files)
)

# ------------------------------------------------------------
# 4. Define Bonferroni threshold
# ------------------------------------------------------------

bonferroni_threshold <- 0.05 /(n) # Defined by your framework

message(
  "Bonferroni threshold: ",
  format(
    bonferroni_threshold,
    scientific = TRUE,
    digits = 4
  )
)

# ------------------------------------------------------------
# 5. Read outcome data once
# ------------------------------------------------------------

outcome_raw_list <- vector(
  "list",
  nrow(outcome_traits)
)

names(outcome_raw_list) <- outcome_traits$Trait

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
  
  message(
    "Loading outcome: ",
    trait
  )
  
  outcome_raw_list[[trait]] <- vroom(
    outcome_file,
    col_names = TRUE,
    show_col_types = FALSE
  ) %>%
    mutate(
      samplesize = outcome_traits$SampleSize[j]
    )
}

# ------------------------------------------------------------
# 6. Initialize result containers
# ------------------------------------------------------------

mr_results_all <- list()
steiger_results_all <- list()

# ------------------------------------------------------------
# 7. Main MR loop
# ------------------------------------------------------------

for (i in seq_along(exposure_files)) {
  
  exposure_file <- exposure_files[i]
  
  idp_name <- sub(
    "^exp_dat_",
    "",
    basename(exposure_file)
  )
  
  idp_name <- sub(
    "_Ffiltered\\.csv$",
    "",
    idp_name
  )
  
  message(
    "\n===================================================="
  )
  
  message(
    "Exposure ",
    i,
    "/",
    length(exposure_files),
    ": ",
    idp_name
  )
  
  message(
    "===================================================="
  )
  
  # ----------------------------------------------------------
  # 7.1 Read exposure data
  # ----------------------------------------------------------
  
  exposure_dat <- read.csv(
    exposure_file,
    stringsAsFactors = FALSE
  )
  
  if (nrow(exposure_dat) < 2) {
    
    warning(
      "Fewer than 2 instruments for ",
      idp_name
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 7.2 Obtain exposure sample size
  # ----------------------------------------------------------
  
  N_values <- unique(
    exposure_dat$samplesize.exposure[
      !is.na(exposure_dat$samplesize.exposure)
    ]
  )
  
  if (length(N_values) == 0) {
    
    warning(
      "Exposure sample size unavailable for ",
      idp_name
    )
    
    N_exposure <- NA_real_
    
  } else {
    
    N_exposure <- as.numeric(
      N_values[1]
    )
  }
  
  message(
    "Exposure sample size: ",
    N_exposure
  )
  
  # ----------------------------------------------------------
  # 7.3 Loop over vascular dementia outcomes
  # ----------------------------------------------------------
  
  for (j in seq_len(nrow(outcome_traits))) {
    
    trait <- outcome_traits$Trait[j]
    
    message(
      "\n  Outcome: ",
      trait
    )
    
    VD_dat <- outcome_raw_list[[trait]]
    
    # --------------------------------------------------------
    # 7.4 Format outcome data
    # --------------------------------------------------------
    
    out_dat <- tryCatch(
      
      format_data(
        VD_dat,
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
          idp_name,
          " -> ",
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
        idp_name,
        " -> ",
        trait
      )
      
      next
    }
    
    # --------------------------------------------------------
    # 7.5 Harmonise exposure and outcome
    # --------------------------------------------------------
    
    dat <- tryCatch(
      
      harmonise_data(
        exposure_dat = exposure_dat,
        outcome_dat = out_dat,
        action = 2
      ),
      
      error = function(e) {
        
        warning(
          "Harmonisation failed for ",
          idp_name,
          " -> ",
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
    
    # --------------------------------------------------------
    # 7.6 Save harmonised data
    # --------------------------------------------------------
    
    harmonised_file <- file.path(
      harmonised_dir,
      paste0(
        idp_name,
        "_on_",
        trait,
        "_harmonised.csv"
      )
    )
    
    write.csv(
      dat,
      harmonised_file,
      row.names = FALSE
    )
    
    # --------------------------------------------------------
    # 7.7 Keep valid instruments
    # --------------------------------------------------------
    
    valid_dat <- dat %>%
      filter(
        mr_keep == TRUE
      )
    
    if (nrow(valid_dat) < 2) {
      
      message(
        "Fewer than 2 valid instruments after harmonisation."
      )
      
      next
    }
    
    # --------------------------------------------------------
    # 7.8 Steiger directionality test
    #
    # Important:
    #   The Steiger result is reported but is NOT used as a
    #   hard exclusion criterion for the MR analysis.
    # --------------------------------------------------------
    
    steiger_result <- tryCatch(
      
      directionality_test(valid_dat),
      
      error = function(e) {
        
        warning(
          "Steiger test failed for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (!is.null(steiger_result)) {
      
      steiger_df <- data.frame(
        exposure = idp_name,
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
    
    # --------------------------------------------------------
    # 7.9 Main MR analyses
    #
    # IVW is the primary method.
    # MR-Egger, weighted median and weighted mode require
    # multiple instruments and are used as complementary
    # sensitivity estimators.
    # --------------------------------------------------------
    
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
          idp_name,
          " -> ",
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
    
    # --------------------------------------------------------
    # 7.10 MR-RAPS
    # --------------------------------------------------------
    
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
            idp_name,
            " -> ",
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
          exposure = idp_name,
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
    
    # --------------------------------------------------------
    # 7.11 Add odds ratios
    #
    #   VaD outcomes are binary traits.
    # --------------------------------------------------------
    
    res <- tryCatch(
      
      generate_odds_ratios(res),
      
      error = function(e) {
        
        warning(
          "Odds-ratio generation failed for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        res
      }
    )
    
    # --------------------------------------------------------
    # 7.12 Extract IVW p-value
    # --------------------------------------------------------
    
    ivw_rows <- res %>%
      filter(
        method == "Inverse variance weighted"
      )
    
    if (nrow(ivw_rows) == 0) {
      
      warning(
        "IVW result unavailable for ",
        idp_name,
        " -> ",
        trait
      )
      
      next
    }
    
    ivw_pval <- ivw_rows$pval[1]
    
    # --------------------------------------------------------
    # 7.13 Calculate overall R2 and F statistic
    # --------------------------------------------------------
    
    overall_R2 <- NA_real_
    overall_F <- NA_real_
    
    if ("R2" %in% names(valid_dat) &&
        !is.na(N_exposure)) {
      
      overall_R2 <- sum(
        valid_dat$R2,
        na.rm = TRUE
      )
      
      k <- nrow(valid_dat)
      
      if (is.finite(overall_R2) &&
          overall_R2 > 0 &&
          overall_R2 < 1) {
        
        overall_F <- (
          overall_R2 *
            (N_exposure - k - 1)
        ) /
          (
            k *
              (1 - overall_R2)
          )
      }
    }
    
    # --------------------------------------------------------
    # 7.14 Add metadata
    # --------------------------------------------------------
    
    res <- res %>%
      mutate(
        exposure = idp_name,
        outcome = trait,
        overall_R2 = overall_R2,
        overall_F_statistic = overall_F,
        nominally_significant =
          pval < 0.05,
        bonferroni_significant =
          pval < bonferroni_threshold
      )
    
    # --------------------------------------------------------
    # 7.15 Store result
    # --------------------------------------------------------
    
    result_key <- paste(
      idp_name,
      trait,
      sep = "__"
    )
    
    mr_results_all[[result_key]] <- res
    
    # --------------------------------------------------------
    # 7.16 Save individual result
    # --------------------------------------------------------
    
    result_file <- file.path(
      result_dir,
      paste0(
        idp_name,
        "_on_",
        trait,
        "_MR.csv"
      )
    )
    
    write.csv(
      res,
      result_file,
      row.names = FALSE
    )
    
    message(
      "Finished: ",
      idp_name,
      " -> ",
      trait,
      " | IVW P = ",
      format.pval(ivw_pval, digits = 3),
      " | IVW Bonferroni = ",
      ifelse(
        ivw_pval < bonferroni_threshold,
        "YES",
        "NO"
      )
    )
    
    # --------------------------------------------------------
    # 7.17 Release temporary objects
    # --------------------------------------------------------
    
    rm(
      VD_dat,
      out_dat,
      dat,
      valid_dat,
      res
    )
    
    gc()
  }
}

# ------------------------------------------------------------
# 8. Combine all MR results
# ------------------------------------------------------------

if (length(mr_results_all) > 0) {
  
  mr_results_combined <- bind_rows(
    mr_results_all
  )
  
  write.csv(
    mr_results_combined,
    file.path(
      result_dir,
      "ROI_SA_to_VaD_all_MR_results.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 8.1 IVW summary
  # ----------------------------------------------------------
  
  ivw_summary <- mr_results_combined %>%
    filter(
      method == "Inverse variance weighted"
    ) %>%
    arrange(
      pval
    )
  
  write.csv(
    ivw_summary,
    file.path(
      result_dir,
      "ROI_SA_to_VaD_IVW_results.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 8.2 Bonferroni-significant results
  # ----------------------------------------------------------
  
  bonferroni_results <- ivw_summary %>%
    filter(
      pval < bonferroni_threshold
    )
  
  write.csv(
    bonferroni_results,
    file.path(
      result_dir,
      "ROI_SA_to_VaD_Bonferroni_significant.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 8.3 Nominally significant results
  # ----------------------------------------------------------
  
  nominal_results <- ivw_summary %>%
    filter(
      pval < 0.05
    )
  
  write.csv(
    nominal_results,
    file.path(
      result_dir,
      "ROI_SA_to_VaD_nominally_significant.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 9. Save Steiger results
# ------------------------------------------------------------

if (length(steiger_results_all) > 0) {
  
  steiger_results <- bind_rows(
    steiger_results_all
  )
  
  write.csv(
    steiger_results,
    file.path(
      result_dir,
      "ROI_SA_to_VaD_Steiger_results.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 10. Save analysis parameters
# ------------------------------------------------------------

analysis_parameters <- data.frame(
  parameter = c(
    "number_of_exposures",
    "number_of_outcomes",
    "total_MR_tests",
    "genome_wide_significance",
    "LD_r2_threshold",
    "LD_window_kb",
    "F_threshold",
    "Bonferroni_alpha"
  ),
  
  value = c(
    n_exposures,
    n_outcomes,
    n_exposures * n_outcomes,
    5e-8,
    0.001,
    10000,
    10,
    bonferroni_threshold
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file.path(
    result_dir,
    "ROI_SA_to_VaD_analysis_parameters.csv"
  ),
  row.names = FALSE
)

message(
  "\n===================================================="
)

message(
  "ROI surface-area -> vascular dementia MR analysis completed."
)

message(
  "Bonferroni threshold = ",
  format(
    bonferroni_threshold,
    scientific = TRUE,
    digits = 4
  )
)

message(
  "===================================================="
)
