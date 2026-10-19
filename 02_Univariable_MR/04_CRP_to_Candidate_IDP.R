# ============================================================
# CRP -> Candidate Cortical IDPs
# Two-sample Mendelian randomization
#
# Exposure:
#   Genetically predicted C-reactive protein (CRP)
#
# Outcomes:
#   Candidate cortical surface area IDPs identified in
#   the preceding MR screening analysis
#
# ============================================================

library(here)
library(vroom)
library(dplyr)
library(TwoSampleMR)

# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_root <- here::here()

crp_instrument_file <- here(
  "data", "processed", "CRP",
  "CRP_instruments_F_filtered.csv"
)

idp_dir <- here(
  "data", "processed",
  "IDP_ROI_SA_CAMIDE"
)

output_dir <- here(
  "results",
  "CRP_to_Candidate_IDP"
)

harmonised_dir <- file.path(output_dir, "harmonised")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(harmonised_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Candidate IDPs
# ------------------------------------------------------------
#
# These IDPs should be the significant ROI phenotypes
# identified in the preceding screening analysis.
#
# Keep the original GWAS/IDP identifiers here rather than
# replacing them with generic labels such as IDP1 or IDP2.

candidate_idps <- c(
  "surface_area-SA_plinkmeta101",
  "surface_area-SA_plinkmeta105"
)

candidate_files <- file.path(
  idp_dir,
  paste0(candidate_idps, "_rsID_QC.tsv.gz")
)

# Check whether all candidate IDP files exist
missing_files <- candidate_files[!file.exists(candidate_files)]

if (length(missing_files) > 0) {
  stop(
    "The following candidate IDP files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 3. Read final CRP instruments
# ------------------------------------------------------------

exposure_dat <- read.csv(
  crp_instrument_file,
  header = TRUE
)

if (nrow(exposure_dat) < 2) {
  stop("Fewer than 2 CRP instruments are available.")
}

message(
  "Number of CRP instruments: ",
  nrow(exposure_dat)
)

N_exposure <- unique(exposure_dat$samplesize.exposure)

if (length(N_exposure) != 1 || is.na(N_exposure)) {
  stop("A unique valid exposure sample size is required.")
}

message(
  "CRP exposure sample size: ",
  N_exposure
)


# ------------------------------------------------------------
# 4. Containers for results
# ------------------------------------------------------------

mr_results_all <- list()
directionality_results <- list()


# ------------------------------------------------------------
# 5. CRP -> candidate IDP MR
# ------------------------------------------------------------

for (j in seq_along(candidate_idps)) {
  
  idp_name <- candidate_idps[j]
  idp_file <- candidate_files[j]
  
  message("----------------------------------------")
  message("Processing outcome: ", idp_name)
  message("----------------------------------------")
  
  
  # ----------------------------------------------------------
  # 5.1 Read candidate IDP GWAS
  # ----------------------------------------------------------
  
  outcome_raw <- vroom(
    idp_file,
    show_col_types = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 5.2 Format outcome data
  # ----------------------------------------------------------
  
  outcome_dat <- format_data(
    outcome_raw,
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
  
  
  # ----------------------------------------------------------
  # 5.3 Harmonise CRP and candidate IDP
  # ----------------------------------------------------------
  
  dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat,
    action = 2
  )
  
  if (is.null(dat) || nrow(dat) < 2) {
    message(
      "Fewer than 2 SNPs available after harmonisation: ",
      idp_name
    )
    next
  }
  
  
  # ----------------------------------------------------------
  # 5.4 Retain valid harmonised SNPs
  # ----------------------------------------------------------
  
  valid_dat <- dat %>%
    filter(mr_keep == TRUE)
  
  message(
    "Valid SNPs after harmonisation: ",
    nrow(valid_dat)
  )
  
  if (nrow(valid_dat) < 2) {
    message(
      "Fewer than 2 valid SNPs for MR: ",
      idp_name
    )
    next
  }
  
  
  # ----------------------------------------------------------
  # 5.5 Save harmonised data
  # ----------------------------------------------------------
  
  harmonised_file <- file.path(
    harmonised_dir,
    paste0(
      "CRP_to_",
      idp_name,
      "_harmonised.csv"
    )
  )
  
  write.csv(
    valid_dat,
    harmonised_file,
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 5.6 Steiger directionality test
  #
  # Steiger is treated as a diagnostic rather than a
  # hard exclusion criterion.
  # ----------------------------------------------------------
  
  steiger_result <- tryCatch(
    directionality_test(valid_dat),
    error = function(e) {
      message(
        "Steiger test failed for ",
        idp_name,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (!is.null(steiger_result)) {
    
    steiger_df <- data.frame(
      exposure = "CRP",
      outcome = idp_name,
      nsnp = nrow(valid_dat),
      snp_r2.exposure = steiger_result$snp_r2.exposure,
      snp_r2.outcome = steiger_result$snp_r2.outcome,
      correct_causal_direction =
        steiger_result$correct_causal_direction,
      steiger_pval =
        steiger_result$steiger_pval
    )
    
    directionality_results[[idp_name]] <- steiger_df
  }
  
  
  # ----------------------------------------------------------
  # 5.7 MR analysis
  #
  # Random-effects IVW is the primary analysis.
  # Additional methods are included according to SNP count.
  # ----------------------------------------------------------
  
  method_list <- c("mr_ivw")
  
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
  # 5.8 MR-RAPS
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
          "MR-RAPS failed for ",
          idp_name,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(raps_result)) {
      
      raps_df <- data.frame(
        id.exposure = unique(valid_dat$id.exposure),
        id.outcome = unique(valid_dat$id.outcome),
        outcome = idp_name,
        exposure = "CRP",
        method = "MR-RAPS",
        nsnp = nrow(valid_dat),
        b = raps_result$beta.hat,
        se = raps_result$beta.se,
        pval = raps_result$beta.p.value
      )
      
      res <- bind_rows(res, raps_df)
    }
  }
  
  
  # ----------------------------------------------------------
  # 5.9 Odds ratios
  #
  # Not meaningful for continuous IDP outcomes, therefore
  # no OR transformation is performed here.
  # ----------------------------------------------------------
  
  # ----------------------------------------------------------
  # 5.10 Add metadata
  # ----------------------------------------------------------
  
  res$exposure <- "CRP"
  res$outcome <- idp_name
  res$nsnp <- ifelse(
    is.na(res$nsnp),
    nrow(valid_dat),
    res$nsnp
  )
  
  res$bonferroni_threshold <- 0.05 / length(candidate_idps)
  
  res$bonferroni_significant <-
    res$pval < res$bonferroni_threshold
  
  
  # ----------------------------------------------------------
  # 5.11 Overall R2 and F statistic
  #
  # The formal analysis uses calculate_R2_F.R only.
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
      "R2 column not found in CRP instrument file: ",
      crp_instrument_file
    )
    
    res$overall_R2 <- NA_real_
    res$overall_F_statistic <- NA_real_
  }
  
  
  # ----------------------------------------------------------
  # 5.12 Save MR results
  # ----------------------------------------------------------
  
  output_file <- file.path(
    output_dir,
    paste0(
      "CRP_to_",
      idp_name,
      "_MR_results.csv"
    )
  )
  
  write.csv(
    res,
    output_file,
    row.names = FALSE
  )
  
  mr_results_all[[idp_name]] <- res
  
  message(
    "Finished: ",
    idp_name
  )
}


# ------------------------------------------------------------
# 6. Combine all MR results
# ------------------------------------------------------------

if (length(mr_results_all) > 0) {
  
  all_mr_results <- bind_rows(
    mr_results_all
  )
  
  write.csv(
    all_mr_results,
    file.path(
      output_dir,
      "CRP_to_Candidate_IDP_MR_all_results.csv"
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 7. Combine Steiger results
# ------------------------------------------------------------

if (length(directionality_results) > 0) {
  
  all_directionality_results <- bind_rows(
    directionality_results
  )
  
  write.csv(
    all_directionality_results,
    file.path(
      output_dir,
      "CRP_to_Candidate_IDP_Steiger_results.csv"
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 8. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

message("CRP -> Candidate IDP MR analysis completed.")