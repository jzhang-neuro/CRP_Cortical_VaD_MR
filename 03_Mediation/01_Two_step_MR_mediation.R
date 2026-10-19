# ============================================================
# Two-step MR mediation:
# CRP -> Candidate cortical surface area IDP -> VD
#
#
# Estimation:
#   Indirect effect = beta_A * beta_B
#   Delta-method standard error
#   Mediation proportion = indirect effect / total effect
#
# ============================================================

library(here)
library(dplyr)
library(openxlsx)


# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_root <- here::here()

crp_to_idp_dir <- here(
  "results",
  "CRP_to_Candidate_IDP"
)

idp_to_vd_dir <- here(
  "results",
  "Regional_IDP_to_VaD"
)

crp_to_vd_dir <- here(
  "results",
  "CRP_to_VaD"
)

output_dir <- here(
  "results",
  "Two_step_MR_mediation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Lock candidate mediators
# ------------------------------------------------------------
#
# These are the candidate IDPs identified in the preceding
# screening analysis. Do not automatically search for or
# expand the candidate set here.

candidate_idps <- c(
  "surface_area-SA_plinkmeta101",
  "surface_area-SA_plinkmeta105"
)


# ------------------------------------------------------------
# 3. Lock VD outcomes
# ------------------------------------------------------------

vd_types <- c(
  "VD_SC"
)

# ------------------------------------------------------------
# 4. Helper function
# ------------------------------------------------------------

read_ivw_result <- function(
    file,
    exposure_name = NULL,
    outcome_name = NULL
) {
  
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  dat <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Identify the IVW result
  ivw <- dat[
    dat$method == "Inverse variance weighted",
    ,
    drop = FALSE
  ]
  
  if (nrow(ivw) == 0) {
    stop(
      "No 'Inverse variance weighted' result found in: ",
      basename(file)
    )
  }
  
  if (nrow(ivw) > 1) {
    warning(
      "More than one IVW row found in: ",
      basename(file),
      ". The first row will be used."
    )
  }
  
  ivw <- ivw[1, ]
  
  result <- list(
    beta = as.numeric(ivw$b),
    se = as.numeric(ivw$se),
    pval = as.numeric(ivw$pval),
    nsnp = if ("nsnp" %in% names(ivw)) {
      as.numeric(ivw$nsnp)
    } else {
      NA_real_
    }
  )
  
  return(result)
}


# ------------------------------------------------------------
# 5. Two-step MR mediation function
# ------------------------------------------------------------

mediation_analysis <- function(
    file_A,
    file_B,
    file_total,
    mediator,
    vd
) {
  
  # ----------------------------------------------------------
  # Path A:
  # CRP -> Candidate IDP
  # ----------------------------------------------------------
  
  path_A <- read_ivw_result(
    file_A,
    exposure_name = "CRP",
    outcome_name = mediator
  )
  
  beta_A <- path_A$beta
  se_A <- path_A$se
  p_A <- path_A$pval
  
  
  # ----------------------------------------------------------
  # Path B:
  # Candidate IDP -> VD
  # ----------------------------------------------------------
  
  path_B <- read_ivw_result(
    file_B,
    exposure_name = mediator,
    outcome_name = vd
  )
  
  beta_B <- path_B$beta
  se_B <- path_B$se
  p_B <- path_B$pval
  
  
  # ----------------------------------------------------------
  # Total effect:
  # CRP -> VD
  # ----------------------------------------------------------
  
  total_effect <- read_ivw_result(
    file_total,
    exposure_name = "CRP",
    outcome_name = vd
  )
  
  beta_total <- total_effect$beta
  se_total <- total_effect$se
  p_total <- total_effect$pval
  
  
  # ----------------------------------------------------------
  # Indirect effect
  # ----------------------------------------------------------
  
  beta_indirect <- beta_A * beta_B
  
  se_indirect <- sqrt(
    beta_A^2 * se_B^2 +
      beta_B^2 * se_A^2
  )
  
  if (
    is.na(se_indirect) ||
    !is.finite(se_indirect) ||
    se_indirect <= 0
  ) {
    
    z_indirect <- NA_real_
    p_indirect <- NA_real_
    indirect_lower <- NA_real_
    indirect_upper <- NA_real_
    
  } else {
    
    z_indirect <- beta_indirect / se_indirect
    
    p_indirect <- 2 * pnorm(
      abs(z_indirect),
      lower.tail = FALSE
    )
    
    indirect_lower <-
      beta_indirect - 1.96 * se_indirect
    
    indirect_upper <-
      beta_indirect + 1.96 * se_indirect
  }
  
  
  # ----------------------------------------------------------
  # Total-effect 95% CI
  # ----------------------------------------------------------
  
  total_lower <-
    beta_total - 1.96 * se_total
  
  total_upper <-
    beta_total + 1.96 * se_total
  
  
  # ----------------------------------------------------------
  # Mediation proportion
  #
  # Avoid unstable estimates when total effect is zero or
  # extremely close to zero.
  # ----------------------------------------------------------
  
  total_effect_epsilon <- 1e-12
  
  if (
    is.na(beta_total) ||
    !is.finite(beta_total) ||
    abs(beta_total) < total_effect_epsilon
  ) {
    
    mediation_proportion <- NA_real_
    mediation_proportion_se <- NA_real_
    mediation_proportion_lower <- NA_real_
    mediation_proportion_upper <- NA_real_
    
  } else {
    
    mediation_proportion <-
      beta_indirect / beta_total
    
    # Delta-method SE for:
    # indirect / total
    
    mediation_proportion_se <- sqrt(
      (se_indirect / beta_total)^2 +
        (
          beta_indirect *
            se_total /
            beta_total^2
        )^2
    )
    
    mediation_proportion_lower <-
      mediation_proportion -
      1.96 * mediation_proportion_se
    
    mediation_proportion_upper <-
      mediation_proportion +
      1.96 * mediation_proportion_se
  }
  
  
  # ----------------------------------------------------------
  # Assemble results
  # ----------------------------------------------------------
  
  results <- data.frame(
    
    mediator = mediator,
    
    outcome = vd,
    
    path = c(
      "Total effect",
      "Path A: CRP -> IDP",
      "Path B: IDP -> VD",
      "Indirect effect",
      "Mediation proportion"
    ),
    
    beta = c(
      beta_total,
      beta_A,
      beta_B,
      beta_indirect,
      mediation_proportion
    ),
    
    se = c(
      se_total,
      se_A,
      se_B,
      se_indirect,
      mediation_proportion_se
    ),
    
    lower_95_CI = c(
      total_lower,
      beta_A - 1.96 * se_A,
      beta_B - 1.96 * se_B,
      indirect_lower,
      mediation_proportion_lower
    ),
    
    upper_95_CI = c(
      total_upper,
      beta_A + 1.96 * se_A,
      beta_B + 1.96 * se_B,
      indirect_upper,
      mediation_proportion_upper
    ),
    
    p_value = c(
      p_total,
      p_A,
      p_B,
      p_indirect,
      NA_real_
    ),
    
    z_statistic = c(
      NA_real_,
      NA_real_,
      NA_real_,
      z_indirect,
      NA_real_
    ),
    
    nsnp = c(
      total_effect$nsnp,
      path_A$nsnp,
      path_B$nsnp,
      NA_real_,
      NA_real_
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Add overall path information
  # ----------------------------------------------------------
  
  results$beta_A <- beta_A
  results$se_A <- se_A
  results$p_A <- p_A
  
  results$beta_B <- beta_B
  results$se_B <- se_B
  results$p_B <- p_B
  
  results$beta_total <- beta_total
  results$se_total <- se_total
  results$p_total <- p_total
  
  
  # ----------------------------------------------------------
  # Save individual mediation result
  # ----------------------------------------------------------
  
  output_file <- file.path(
    output_dir,
    paste0(
      "mediation_CRP_",
      mediator,
      "_on_",
      vd,
      ".csv"
    )
  )
  
  write.csv(
    results,
    output_file,
    row.names = FALSE
  )
  
  
  return(results)
}


# ------------------------------------------------------------
# 6. Locate the three required MR results
# ------------------------------------------------------------
#
# Important:
# The candidate IDP and VD names are explicitly matched.
# No broad pattern matching over all surface-area files is used.

all_results <- list()


for (mediator in candidate_idps) {
  
  for (vd in vd_types) {
    
    message("----------------------------------------")
    message(
      "Mediation path: CRP -> ",
      mediator,
      " -> ",
      vd
    )
    message("----------------------------------------")
    
    
    # --------------------------------------------------------
    # Path A:
    # CRP -> Candidate IDP
    # --------------------------------------------------------
    
    file_A <- file.path(
      crp_to_idp_dir,
      paste0(
        "CRP_to_",
        mediator,
        "_MR_results.csv"
      )
    )
    
    
    # --------------------------------------------------------
    # Path B:
    # Candidate IDP -> VD
    # --------------------------------------------------------
    
    file_B <- file.path(
      idp_to_vd_dir,
      paste0(
        mediator,
        "_to_",
        vd,
        "_MR_results.csv"
      )
    )
    
    
    # --------------------------------------------------------
    # Total effect:
    # CRP -> VD
    # --------------------------------------------------------
    
    file_total <- file.path(
      crp_to_vd_dir,
      paste0(
        "CRP_to_",
        vd,
        "_MR_results.csv"
      )
    )
    
    
    # --------------------------------------------------------
    # Check the three required files
    # --------------------------------------------------------
    
    required_files <- c(
      file_A,
      file_B,
      file_total
    )
    
    missing_files <- required_files[
      !file.exists(required_files)
    ]
    
    if (length(missing_files) > 0) {
      
      warning(
        "Skipping mediation path because required files ",
        "are missing:\n",
        paste(missing_files, collapse = "\n")
      )
      
      next
    }
    
    
    # --------------------------------------------------------
    # Run two-step mediation
    # --------------------------------------------------------
    
    result <- mediation_analysis(
      file_A = file_A,
      file_B = file_B,
      file_total = file_total,
      mediator = mediator,
      vd = vd
    )
    
    all_results[[paste(mediator, vd, sep = "__")]] <-
      result
  }
}


# ------------------------------------------------------------
# 7. Combine all mediation results
# ------------------------------------------------------------

if (length(all_results) == 0) {
  
  stop(
    "No valid mediation paths were completed."
  )
  
}

mediation_summary <- bind_rows(
  all_results
)


# ------------------------------------------------------------
# 8. Add multiple-testing information
# ------------------------------------------------------------
#
# Four mediation pathways are tested:
#   2 candidate IDPs x 1 VD outcome
#
# The correction is applied to the indirect-effect P values.

n_mediation_tests <-
  length(candidate_idps) * length(vd_types)

bonferroni_threshold <-
  0.05 / n_mediation_tests


mediation_summary$bonferroni_threshold <-
  bonferroni_threshold


mediation_summary$bonferroni_significant <-
  mediation_summary$p_value <
  bonferroni_threshold


# ------------------------------------------------------------
# 9. Save complete mediation summary
# ------------------------------------------------------------

write.csv(
  mediation_summary,
  file.path(
    output_dir,
    "two_step_MR_mediation_all_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Create a compact summary for reporting
# ------------------------------------------------------------

compact_summary <- mediation_summary %>%
  filter(
    path %in% c(
      "Total effect",
      "Indirect effect",
      "Mediation proportion"
    )
  ) %>%
  select(
    mediator,
    outcome,
    path,
    beta,
    se,
    lower_95_CI,
    upper_95_CI,
    p_value,
    bonferroni_threshold,
    bonferroni_significant
  )


write.csv(
  compact_summary,
  file.path(
    output_dir,
    "two_step_MR_mediation_compact_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Export Excel summary
# ------------------------------------------------------------

write.xlsx(
  list(
    All_results = mediation_summary,
    Compact_summary = compact_summary
  ),
  file = file.path(
    output_dir,
    "two_step_MR_mediation_summary.xlsx"
  ),
  rowNames = FALSE
)


# ------------------------------------------------------------
# 12. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


message(
  "Two-step MR mediation analysis completed."
)