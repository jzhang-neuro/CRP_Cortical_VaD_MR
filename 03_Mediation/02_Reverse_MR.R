# ============================================================
# Reverse Mendelian randomization:
# Candidate cortical surface area IDPs -> CRP
#
# Exposures:
#   IDP_surface_area
#
# Outcome:
#   C-reactive protein (CRP)
#
# Purpose:
#   Assess evidence for reverse causality between candidate
#   cortical surface area IDPs and CRP.
#
# ============================================================

library(here)
library(vroom)
library(dplyr)
library(TwoSampleMR)
library(mr.raps)


# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_root <- here::here()

idp_instrument_dir <- here(
  "data",
  "processed",
  "IDP_ROI_SA_CAMIDE"
)

crp_outcome_file <- here(
  "data",
  "raw",
  "CRP",
  "CRP_with_EAF.tsv"
)

output_dir <- here(
  "results",
  "Reverse_MR_IDP_to_CRP"
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


# ------------------------------------------------------------
# 2. Lock candidate IDPs
# ------------------------------------------------------------
#
# These are the same candidate IDPs identified in the
# preceding regional screening analysis.
#
# Do not expand the candidate set automatically here.

candidate_idps <- c(
  "surface_area-SA_plinkmeta101",
  "surface_area-SA_plinkmeta105"
)

exposure_files <- file.path(
  idp_instrument_dir,
  paste0(
    candidate_idps,
    "_Ffiltered.csv"
  )
)

missing_files <- exposure_files[
  !file.exists(exposure_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following IDP instrument files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 3. Read CRP outcome GWAS
# ------------------------------------------------------------

crp_raw <- vroom(
  crp_outcome_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 4. Format CRP outcome
# ------------------------------------------------------------
#
# CRP_with_EAF.tsv should contain standardized GWAS columns:
#
#   SNP
#   BETA
#   SE
#   A1
#   A2
#   EAF
#   P
#   SAMPLESIZE
#
# If the current CRP file still uses .exposure column names,
# standardize it in the CRP data-preparation script rather
# than introducing a special format here.

crp_outcome_dat <- format_data(
  crp_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  pval_col = "P",
  samplesize_col = "SAMPLESIZE"
)


# ------------------------------------------------------------
# 5. Containers for results
# ------------------------------------------------------------

mr_results_all <- list()
directionality_results <- list()


# ------------------------------------------------------------
# 6. Loop over candidate IDPs
# ------------------------------------------------------------

for (i in seq_along(candidate_idps)) {
  
  exposure_name <- candidate_idps[i]
  exposure_file <- exposure_files[i]
  
  message("----------------------------------------")
  message(
    "Reverse MR: ",
    exposure_name,
    " -> CRP"
  )
  message("----------------------------------------")
  
  
  # ----------------------------------------------------------
  # 6.1 Read IDP exposure instruments
  # ----------------------------------------------------------
  
  exposure_dat <- read.csv(
    exposure_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  if (nrow(exposure_dat) < 2) {
    message(
      "Fewer than 2 IDP instruments available. Skipping."
    )
    next
  }
  
  
  # ----------------------------------------------------------
  # 6.2 Exposure sample size
  # ----------------------------------------------------------
  
  N_exposure <- unique(
    exposure_dat$samplesize.exposure
  )
  
  if (
    length(N_exposure) != 1 ||
    is.na(N_exposure)
  ) {
    stop(
      "A unique valid sample size is required for: ",
      exposure_name
    )
  }
  
  message(
    "Exposure sample size: ",
    N_exposure
  )
  
  message(
    "Number of exposure instruments: ",
    nrow(exposure_dat)
  )
  
  
  # ----------------------------------------------------------
  # 6.3 Harmonise IDP exposure and CRP outcome
  # ----------------------------------------------------------
  
  dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = crp_outcome_dat,
    action = 2
  )
  
  
  if (is.null(dat) || nrow(dat) < 2) {
    
    message(
      "Fewer than 2 SNPs after harmonisation. Skipping."
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # 6.4 Retain valid instruments
  # ----------------------------------------------------------
  
  valid_dat <- dat %>%
    filter(mr_keep == TRUE)
  
  
  message(
    "Valid SNPs after harmonisation: ",
    nrow(valid_dat)
  )
  
  
  if (nrow(valid_dat) < 2) {
    
    message(
      "Fewer than 2 valid SNPs for MR. Skipping."
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # 6.5 Save harmonised data
  # ----------------------------------------------------------
  
  harmonised_file <- file.path(
    harmonised_dir,
    paste0(
      exposure_name,
      "_to_CRP_harmonised.csv"
    )
  )
  
  write.csv(
    valid_dat,
    harmonised_file,
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.6 Steiger directionality test
  #
  # Steiger is used as a directionality diagnostic.
  # It is NOT used as a hard exclusion criterion.
  # ----------------------------------------------------------
  
  steiger_result <- tryCatch(
    directionality_test(valid_dat),
    error = function(e) {
      
      message(
        "Steiger test failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (!is.null(steiger_result)) {
    
    steiger_df <- data.frame(
      exposure = exposure_name,
      outcome = "CRP",
      nsnp = nrow(valid_dat),
      snp_r2.exposure =
        steiger_result$snp_r2.exposure,
      snp_r2.outcome =
        steiger_result$snp_r2.outcome,
      correct_causal_direction =
        steiger_result$correct_causal_direction,
      steiger_pval =
        steiger_result$steiger_pval
    )
    
    directionality_results[[exposure_name]] <-
      steiger_df
  }
  
  
  # ----------------------------------------------------------
  # 6.7 Main MR
  #
  # Random-effects IVW is the primary method.
  # Additional methods are conditional on SNP number.
  # ----------------------------------------------------------
  
  method_list <- c(
    "mr_ivw"
  )
  
  if (nrow(valid_dat) >= 3) {
    
    method_list <- c(
      method_list,
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_egger_regression"
    )
  }
  
  
  res <- mr(
    valid_dat,
    method_list = method_list
  )
  
  
  # ----------------------------------------------------------
  # 6.8 MR-RAPS
  # ----------------------------------------------------------
  
  if (nrow(valid_dat) >= 3) {
    
    raps_result <- tryCatch(
      
      mr.raps::mr.raps(
        b_exp = valid_dat$beta.exposure,
        b_out = valid_dat$beta.outcome,
        se_exp = valid_dat$se.exposure,
        se_out = valid_dat$se.outcome
      ),
      
      error = function(e) {
        
        message(
          "MR-RAPS failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    
    if (!is.null(raps_result)) {
      
      raps_df <- data.frame(
        
        id.exposure =
          unique(valid_dat$id.exposure),
        
        id.outcome =
          unique(valid_dat$id.outcome),
        
        exposure = exposure_name,
        
        outcome = "CRP",
        
        method = "MR-RAPS",
        
        nsnp = nrow(valid_dat),
        
        b = raps_result$beta.hat,
        
        se = raps_result$beta.se,
        
        pval = raps_result$beta.p.value
      )
      
      res <- bind_rows(
        res,
        raps_df
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # 6.9 Add metadata
  # ----------------------------------------------------------
  
  res$exposure <- exposure_name
  res$outcome <- "CRP"
  
  
  # ----------------------------------------------------------
  # 6.10 Overall R2 and F statistic
  #
  # R2 is obtained from the formal calculate_R2_F()
  # instrument-preparation pipeline.
  # ----------------------------------------------------------
  
  if ("R2" %in% names(valid_dat)) {
    
    overall_R2 <- sum(
      valid_dat$R2,
      na.rm = TRUE
    )
    
    k <- nrow(valid_dat)
    
    overall_F <- (
      overall_R2 *
        (N_exposure - k - 1)
    ) /
      (
        k *
          (1 - overall_R2)
      )
    
    res$overall_R2 <- overall_R2
    res$overall_F_statistic <- overall_F
    
  } else {
    
    warning(
      "R2 column not found for ",
      exposure_name
    )
    
    res$overall_R2 <- NA_real_
    res$overall_F_statistic <- NA_real_
  }
  
  
  # ----------------------------------------------------------
  # 6.11 Save MR results
  # ----------------------------------------------------------
  
  output_file <- file.path(
    output_dir,
    paste0(
      exposure_name,
      "_to_CRP_MR_results.csv"
    )
  )
  
  write.csv(
    res,
    output_file,
    row.names = FALSE
  )
  
  
  mr_results_all[[exposure_name]] <-
    res
  
  
  message(
    "Completed: ",
    exposure_name,
    " -> CRP"
  )
}


# ------------------------------------------------------------
# 7. Combine MR results
# ------------------------------------------------------------

if (length(mr_results_all) > 0) {
  
  all_mr_results <- bind_rows(
    mr_results_all
  )
  
  write.csv(
    all_mr_results,
    file.path(
      output_dir,
      "reverse_MR_IDP_to_CRP_all_results.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  warning(
    "No reverse MR analyses were successfully completed."
  )
}


# ------------------------------------------------------------
# 8. Combine Steiger results
# ------------------------------------------------------------

if (length(directionality_results) > 0) {
  
  all_directionality_results <- bind_rows(
    directionality_results
  )
  
  write.csv(
    all_directionality_results,
    file.path(
      output_dir,
      "reverse_MR_IDP_to_CRP_Steiger_results.csv"
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 9. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


message(
  "Reverse MR analysis completed."
)