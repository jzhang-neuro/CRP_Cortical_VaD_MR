# ============================================================
# Global cortical IDP -> vascular dementia MR
#
# Purpose:
#   To evaluate the associations between genetically predicted
#   global cortical imaging-derived phenotypes (IDPs) and
#   vascular dementia (VaD) and its clinically defined subtypes.
#
# Exposure:
#   Global cortical IDPs from the CAMIDE neuroimaging GWAS
#
# Outcomes:
#   1. Overall vascular dementia
#   2. Multi-infarct vascular dementia
#   3. Mixed vascular dementia
#   4. Subcortical vascular dementia
#   5. Undefined vascular dementia
#
# Primary method:
#   Random-effects inverse-variance weighted (IVW) MR
#
# Additional analyses:
#   MR-RAPS
#   Weighted median
#   Weighted mode
#   MR-Egger
#   MR-PRESSO
#   Cochran's Q test
#   MR-Egger intercept
#   Leave-one-out analysis
#   Steiger directionality test
#
# Input:
#   results/instruments/global_idp/F_filtered/
#   data/external/Outcome_VD/
#
# Output:
#   results/mr/global_idp_to_vad/
#
# Notes:
#   This script uses only summary-level GWAS data.
#   No individual-level participant data are processed.
# ============================================================


# ------------------------------
# 1. Load packages
# ------------------------------

library(vroom)
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(mr.raps)
library(MRPRESSO)
library(ggplot2)


# ------------------------------
# 2. Define analysis parameters
# ------------------------------

f_threshold <- 10

p_threshold <- 0.05


# ------------------------------
# 3. Define project directories
# ------------------------------

exposure_dir <- file.path(
  "results",
  "instruments",
  "global_idp",
  "F_filtered"
)

outcome_dir <- file.path(
  "data",
  "external",
  "Outcome_VD"
)

result_dir <- file.path(
  "results",
  "mr",
  "global_idp_to_vad"
)

data_dir <- file.path(
  result_dir,
  "harmonized_data"
)

plot_dir <- file.path(
  result_dir,
  "plots"
)

presso_dir <- file.path(
  result_dir,
  "MR_PRESSO"
)


# Create output directories

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  presso_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------
# 4. Check input directories
# ------------------------------

if (!dir.exists(exposure_dir)) {
  stop(
    "Exposure instrument directory not found: ",
    exposure_dir
  )
}

if (!dir.exists(outcome_dir)) {
  stop(
    "VaD outcome directory not found: ",
    outcome_dir
  )
}


# ------------------------------
# 5. Define VaD outcomes
# ------------------------------

outcome_traits <- data.frame(
  Trait = c(
    "Overall_VaD",
    "Multi_infarct_VaD",
    "Mixed_VaD",
    "Subcortical_VaD",
    "Undefined_VaD"
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

outcome_traits$Sample_size <-
  outcome_traits$Case +
  outcome_traits$Control


# ------------------------------
# 6. Identify exposure files
# ------------------------------

exposure_files <- list.files(
  exposure_dir,
  pattern = "^exp_dat_.*_Ffiltered\\.csv$",
  full.names = TRUE
)

if (length(exposure_files) == 0) {
  stop(
    "No F-filtered global cortical IDP instruments were found in: ",
    exposure_dir
  )
}

message(
  "Number of global cortical IDPs identified: ",
  length(exposure_files)
)


# ------------------------------
# 7. Initialize result containers
# ------------------------------

all_mr_results <- list()

significant_results <- data.frame()

heterogeneity_results <- data.frame()

pleiotropy_results <- data.frame()

directionality_results <- data.frame()


# ------------------------------
# 8. Main IDP -> VaD analysis
# ------------------------------

for (i in seq_along(exposure_files)) {
  
  exposure_file <- exposure_files[i]
  
  exposure_filename <- basename(
    exposure_file
  )
  
  idp_name <- sub(
    "^exp_dat_(.*)_Ffiltered\\.csv$",
    "\\1",
    exposure_filename
  )
  
  message("\n========================================")
  message("Processing IDP: ", idp_name)
  message("========================================")
  
  
  # ------------------------------
  # 8.1 Read exposure data
  # ------------------------------
  
  exposure_dat <- fread(
    exposure_file,
    data.table = FALSE
  )
  
  
  # ------------------------------
  # 8.2 Check exposure sample size
  # ------------------------------
  
  if (!"samplesize.exposure" %in% names(exposure_dat)) {
    
    stop(
      "Column 'samplesize.exposure' is missing for ",
      idp_name
    )
  }
  
  N <- unique(
    exposure_dat$samplesize.exposure[
      !is.na(exposure_dat$samplesize.exposure)
    ]
  )
  
  if (length(N) != 1) {
    
    stop(
      "A unique exposure sample size could not be identified for ",
      idp_name
    )
  }
  
  message(
    "Exposure sample size: ",
    N
  )
  
  
  # ------------------------------
  # 8.3 Process each VaD outcome
  # ------------------------------
  
  for (j in seq_len(nrow(outcome_traits))) {
    
    trait <- outcome_traits$Trait[j]
    
    outcome_file <- file.path(
      outcome_dir,
      outcome_traits$File[j]
    )
    
    message(
      "\nAnalyzing: ",
      idp_name,
      " -> ",
      trait
    )
    
    
    # ------------------------------
    # 8.3.1 Check outcome file
    # ------------------------------
    
    if (!file.exists(outcome_file)) {
      
      warning(
        "Outcome file not found: ",
        outcome_file
      )
      
      next
    }
    
    
    # ------------------------------
    # 8.3.2 Read outcome GWAS
    # ------------------------------
    
    VD_dat <- vroom(
      outcome_file,
      col_names = TRUE,
      show_col_types = FALSE
    )
    
    
    # ------------------------------
    # 8.3.3 Add outcome sample size
    # ------------------------------
    
    VD_dat <- VD_dat %>%
      mutate(
        samplesize =
          outcome_traits$Sample_size[j]
      )
    
    
    # ------------------------------
    # 8.3.4 Format outcome data
    # ------------------------------
    
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
        pval_col = "pval"
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
      
      warning(
        "Insufficient outcome SNPs for ",
        idp_name,
        " -> ",
        trait
      )
      
      next
    }
    
    
    # ------------------------------
    # 8.3.5 Harmonize exposure and outcome
    # ------------------------------
    
    dat <- tryCatch(
      
      harmonise_data(
        exposure_dat = exposure_dat,
        outcome_dat = out_dat,
        action = 2
      ),
      
      error = function(e) {
        
        warning(
          "Harmonization failed for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (is.null(dat) || nrow(dat) < 2) {
      
      warning(
        "Insufficient harmonized SNPs for ",
        idp_name,
        " -> ",
        trait
      )
      
      next
    }
    
    
    # ------------------------------
    # 8.3.6 Retain valid instruments
    # ------------------------------
    
    valid_dat <- dat[
      !is.na(dat$mr_keep) &
        dat$mr_keep == TRUE,
    ]
    
    if (nrow(valid_dat) < 2) {
      
      warning(
        "Fewer than 2 valid instruments after harmonization for ",
        idp_name,
        " -> ",
        trait
      )
      
      next
    }
    
    
    # ------------------------------
    # 8.3.7 Steiger directionality test
    # ------------------------------
    
    steiger_result <- tryCatch(
      
      directionality_test(
        valid_dat
      ),
      
      warning = function(w) {
        
        message(
          "Steiger warning for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(w)
        )
        
        suppressWarnings(
          directionality_test(valid_dat)
        )
      },
      
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
    
    
    # Do not exclude the association from the primary MR analysis
    # solely on the basis of the Steiger result. The Steiger test
    # is retained as a directionality and reverse-causality
    # assessment.
    
    if (!is.null(steiger_result)) {
      
      steiger_df <- data.frame(
        id.exposure = unique(valid_dat$id.exposure),
        id.outcome = unique(valid_dat$id.outcome),
        exposure = idp_name,
        outcome = trait,
        nsnp = nrow(valid_dat),
        snp_r2.exposure =
          steiger_result$snp_r2.exposure,
        snp_r2.outcome =
          steiger_result$snp_r2.outcome,
        correct_causal_direction =
          steiger_result$correct_causal_direction,
        steiger_pval =
          steiger_result$steiger_pval,
        stringsAsFactors = FALSE
      )
      
      directionality_results <- rbind(
        directionality_results,
        steiger_df
      )
    }
    
    
    # ------------------------------
    # 8.3.8 Save harmonized data
    # ------------------------------
    
    safe_idp_name <- gsub(
      "[^A-Za-z0-9._-]",
      "_",
      idp_name
    )
    
    safe_trait <- gsub(
      "[^A-Za-z0-9._-]",
      "_",
      trait
    )
    
    file_prefix <- paste0(
      safe_idp_name,
      "_on_",
      safe_trait
    )
    
    harmonized_file <- file.path(
      data_dir,
      paste0(
        file_prefix,
        "_harmonized.csv"
      )
    )
    
    fwrite(
      valid_dat,
      harmonized_file
    )
    
    
    # ------------------------------
    # 8.3.9 Primary MR analysis
    # ------------------------------
    
    mr_methods <- c(
      "mr_ivw",
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_egger_regression"
    )
    
    res <- tryCatch(
      
      mr(
        valid_dat,
        method_list = mr_methods
      ),
      
      error = function(e) {
        
        warning(
          "Primary MR failed for ",
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
    
    
    # ------------------------------
    # 8.3.10 MR-RAPS
    # ------------------------------
    
    mr_raps_result <- tryCatch(
      
      mr.raps(
        valid_dat$beta.exposure,
        valid_dat$beta.outcome,
        valid_dat$se.exposure,
        valid_dat$se.outcome
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
    
    
    if (!is.null(mr_raps_result)) {
      
      mr_raps_df <- data.frame(
        id.exposure =
          unique(valid_dat$id.exposure),
        id.outcome =
          unique(valid_dat$id.outcome),
        outcome = trait,
        exposure = idp_name,
        method = "MR-RAPS",
        nsnp = nrow(valid_dat),
        b = mr_raps_result$beta.hat,
        se = mr_raps_result$beta.se,
        pval = mr_raps_result$beta.p.value,
        stringsAsFactors = FALSE
      )
      
      res <- rbind(
        res,
        mr_raps_df
      )
    }
    
    
    # ------------------------------
    # 8.3.11 MR-PRESSO
    # ------------------------------
    
    if (nrow(valid_dat) >= 4) {
      
      presso_result <- tryCatch(
        
        run_mr_presso(
          dat = valid_dat
        ),
        
        error = function(e) {
          
          warning(
            "MR-PRESSO failed for ",
            idp_name,
            " -> ",
            trait,
            ": ",
            conditionMessage(e)
          )
          
          NULL
        }
      )
      
      
      if (!is.null(presso_result)) {
        
        presso_main <-
          presso_result[[1]]$`Main MR results`
        
        presso_global <-
          presso_result[[1]]$
          `MR-PRESSO results`$
          `Global Test`
        
        presso_outlier <-
          presso_result[[1]]$
          `MR-PRESSO results`$
          `Outlier Test`
        
        presso_distortion <-
          presso_result[[1]]$
          `MR-PRESSO results`$
          `Distortion Test`
        
        
        # Save MR-PRESSO main results
        
        presso_main_file <- file.path(
          presso_dir,
          paste0(
            file_prefix,
            "_main_results.csv"
          )
        )
        
        fwrite(
          as.data.frame(presso_main),
          presso_main_file
        )
        
        
        # Save global test
        
        presso_global_file <- file.path(
          presso_dir,
          paste0(
            file_prefix,
            "_global_test.csv"
          )
        )
        
        fwrite(
          as.data.frame(presso_global),
          presso_global_file
        )
        
        
        # Save outlier test
        
        if (!is.null(presso_outlier)) {
          
          presso_outlier_file <- file.path(
            presso_dir,
            paste0(
              file_prefix,
              "_outlier_test.csv"
            )
          )
          
          fwrite(
            as.data.frame(presso_outlier),
            presso_outlier_file
          )
        }
        
        
        # Save distortion test
        
        if (!is.null(presso_distortion)) {
          
          presso_distortion_file <- file.path(
            presso_dir,
            paste0(
              file_prefix,
              "_distortion_test.csv"
            )
          )
          
          fwrite(
            as.data.frame(presso_distortion),
            presso_distortion_file
          )
        }
        
        
        # Add outlier-corrected estimate to MR results
        # when available.
        
        if (!is.null(presso_main) &&
            "MR Analysis" %in% names(presso_main)) {
          
          outlier_corrected <-
            presso_main[
              presso_main$`MR Analysis` ==
                "Outlier-corrected",
              ,
              drop = FALSE
            ]
          
          if (nrow(outlier_corrected) > 0) {
            
            mr_presso_outlier_corrected <-
              data.frame(
                id.exposure =
                  unique(valid_dat$id.exposure),
                id.outcome =
                  unique(valid_dat$id.outcome),
                outcome = trait,
                exposure = idp_name,
                method =
                  "MR-PRESSO Outlier-corrected",
                nsnp = nrow(valid_dat),
                b =
                  outlier_corrected$
                  `Causal Estimate`[1],
                se =
                  outlier_corrected$Sd[1],
                pval =
                  outlier_corrected$
                  `P-value`[1],
                stringsAsFactors = FALSE
              )
            
            res <- rbind(
              res,
              mr_presso_outlier_corrected
            )
          }
        }
      }
    }
    
    
    # ------------------------------
    # 8.3.12 Generate odds ratios
    # ------------------------------
    
    or_df <- tryCatch(
      
      generate_odds_ratios(res),
      
      error = function(e) {
        
        warning(
          "OR calculation failed for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    
    if (!is.null(or_df)) {
      
      result <- merge(
        res,
        or_df,
        by = c(
          "id.exposure",
          "id.outcome",
          "outcome",
          "exposure",
          "method",
          "nsnp",
          "b",
          "se",
          "pval"
        ),
        all.x = TRUE
      )
      
    } else {
      
      result <- res
    }
    
    
    # ------------------------------
    # 8.3.13 Calculate overall R2 and F
    # ------------------------------
    
    if ("R2" %in% names(valid_dat)) {
      
      overall_R2 <- sum(
        valid_dat$R2,
        na.rm = TRUE
      )
      
      k <- nrow(valid_dat)
      
      overall_F_statistic <-
        (
          overall_R2 *
            (N - k - 1)
        ) /
        (
          k *
            (1 - overall_R2)
        )
      
      result$overall_R2 <-
        overall_R2
      
      result$overall_F_statistic <-
        overall_F_statistic
      
    } else {
      
      result$overall_R2 <- NA_real_
      result$overall_F_statistic <- NA_real_
    }
    
    
    # ------------------------------
    # 8.3.14 Add analysis metadata
    # ------------------------------
    
    result$exposure <- idp_name
    result$outcome <- trait
    result$exposure_sample_size <- N
    result$outcome_sample_size <-
      outcome_traits$Sample_size[j]
    
    
    # ------------------------------
    # 8.3.15 Save complete MR results
    # ------------------------------
    
    ivw_result <- result[
      result$method ==
        "Inverse variance weighted",
      ,
      drop = FALSE
    ]
    
    if (nrow(ivw_result) > 0) {
      
      ivw_pval <- ivw_result$pval[1]
      
    } else {
      
      ivw_pval <- NA_real_
    }
    
    
    result_file <- file.path(
      result_dir,
      paste0(
        file_prefix,
        "_MR_results.csv"
      )
    )
    
    fwrite(
      result,
      result_file
    )
    
    
    # ------------------------------
    # 8.3.16 Store significant results
    # ------------------------------
    
    if (!is.na(ivw_pval) &&
        ivw_pval < p_threshold) {
      
      significant_results <- rbind(
        significant_results,
        result
      )
    }
    
    
    # ------------------------------
    # 8.3.17 Heterogeneity analysis
    # ------------------------------
    
    het <- tryCatch(
      
      mr_heterogeneity(
        valid_dat
      ),
      
      error = function(e) {
        
        warning(
          "Heterogeneity analysis failed for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    
    if (!is.null(het)) {
      
      het$ivw_pval <- ivw_pval
      het$exposure <- idp_name
      het$outcome <- trait
      
      heterogeneity_results <- rbind(
        heterogeneity_results,
        het
      )
    }
    
    
    # ------------------------------
    # 8.3.18 Horizontal pleiotropy
    # ------------------------------
    
    pleio <- tryCatch(
      
      mr_pleiotropy_test(
        valid_dat
      ),
      
      error = function(e) {
        
        warning(
          "Pleiotropy analysis failed for ",
          idp_name,
          " -> ",
          trait,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    
    if (!is.null(pleio)) {
      
      pleio$ivw_pval <- ivw_pval
      pleio$exposure <- idp_name
      pleio$outcome <- trait
      
      pleiotropy_results <- rbind(
        pleiotropy_results,
        pleio
      )
    }
    
    
    # ------------------------------
    # 8.3.19 Scatter plot
    # ------------------------------
    
    scatter_plot <- tryCatch(
      
      mr_scatter_plot(
        res,
        valid_dat
      ),
      
      error = function(e) {
        
        warning(
          "Scatter plot failed for ",
          idp_name,
          " -> ",
          trait
        )
        
        NULL
      }
    )
    
    
    if (!is.null(scatter_plot)) {
      
      ggsave(
        filename = file.path(
          plot_dir,
          paste0(
            file_prefix,
            "_scatter.pdf"
          )
        ),
        plot = scatter_plot[[1]],
        width = 7,
        height = 7
      )
    }
    
    
    # ------------------------------
    # 8.3.20 Single-SNP and forest plot
    # ------------------------------
    
    res_single <- tryCatch(
      
      mr_singlesnp(
        valid_dat
      ),
      
      error = function(e) {
        
        warning(
          "Single-SNP analysis failed for ",
          idp_name,
          " -> ",
          trait
        )
        
        NULL
      }
    )
    
    
    if (!is.null(res_single)) {
      
      forest_plot <- tryCatch(
        mr_forest_plot(res_single),
        error = function(e) NULL
      )
      
      if (!is.null(forest_plot)) {
        
        ggsave(
          filename = file.path(
            plot_dir,
            paste0(
              file_prefix,
              "_forest.pdf"
            )
          ),
          plot = forest_plot[[1]],
          width = 7,
          height = 7
        )
      }
      
      
      # ------------------------------
      # 8.3.21 Funnel plot
      # ------------------------------
      
      funnel_plot <- tryCatch(
        mr_funnel_plot(res_single),
        error = function(e) NULL
      )
      
      if (!is.null(funnel_plot)) {
        
        ggsave(
          filename = file.path(
            plot_dir,
            paste0(
              file_prefix,
              "_funnel.pdf"
            )
          ),
          plot = funnel_plot[[1]],
          width = 7,
          height = 7
        )
      }
    }
    
    
    # ------------------------------
    # 8.3.22 Leave-one-out analysis
    # ------------------------------
    
    loo_result <- tryCatch(
      
      mr_leaveoneout(
        valid_dat
      ),
      
      error = function(e) {
        
        warning(
          "Leave-one-out analysis failed for ",
          idp_name,
          " -> ",
          trait
        )
        
        NULL
      }
    )
    
    
    if (!is.null(loo_result)) {
      
      loo_file <- file.path(
        result_dir,
        paste0(
          file_prefix,
          "_leave_one_out.csv"
        )
      )
      
      fwrite(
        loo_result,
        loo_file
      )
      
      
      loo_plot <- tryCatch(
        mr_leaveoneout_plot(loo_result),
        error = function(e) NULL
      )
      
      if (!is.null(loo_plot)) {
        
        ggsave(
          filename = file.path(
            plot_dir,
            paste0(
              file_prefix,
              "_leave_one_out.pdf"
            )
          ),
          plot = loo_plot[[1]],
          width = 7,
          height = 7
        )
      }
    }
    
    
    # ------------------------------
    # 8.3.23 Store MR results
    # ------------------------------
    
    result$key <- paste(
      idp_name,
      trait,
      sep = " -> "
    )
    
    all_mr_results[[length(all_mr_results) + 1]] <-
      result
    
    
    # ------------------------------
    # 8.3.24 Release memory
    # ------------------------------
    
    rm(
      VD_dat,
      out_dat,
      dat,
      valid_dat
    )
    
    gc()
    
    
    message(
      "Finished: ",
      idp_name,
      " -> ",
      trait
    )
  }
}


# ------------------------------
# 9. Combine and save results
# ------------------------------

if (length(all_mr_results) > 0) {
  
  all_mr_results_df <- bind_rows(
    all_mr_results
  )
  
  fwrite(
    all_mr_results_df,
    file.path(
      result_dir,
      "all_global_IDP_to_VaD_MR_results.csv"
    )
  )
}


# ------------------------------
# 10. Save significant results
# ------------------------------

if (nrow(significant_results) > 0) {
  
  fwrite(
    significant_results,
    file.path(
      result_dir,
      "significant_global_IDP_to_VaD_results.csv"
    )
  )
}


# ------------------------------
# 11. Save heterogeneity results
# ------------------------------

if (nrow(heterogeneity_results) > 0) {
  
  fwrite(
    heterogeneity_results,
    file.path(
      result_dir,
      "heterogeneity_results.csv"
    )
  )
}


# ------------------------------
# 12. Save pleiotropy results
# ------------------------------

if (nrow(pleiotropy_results) > 0) {
  
  fwrite(
    pleiotropy_results,
    file.path(
      result_dir,
      "pleiotropy_results.csv"
    )
  )
}


# ------------------------------
# 13. Save Steiger results
# ------------------------------

if (nrow(directionality_results) > 0) {
  
  fwrite(
    directionality_results,
    file.path(
      result_dir,
      "directionality_results.csv"
    )
  )
}


# ------------------------------
# 14. Completion message
# ------------------------------

message("\n========================================")
message("Global cortical IDP -> VaD MR completed.")
message("Results directory: ", result_dir)
message("========================================")
