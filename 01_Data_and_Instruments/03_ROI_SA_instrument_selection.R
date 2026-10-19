# ============================================================
# ROI cortical surface area instrument selection
#
# Purpose:
#   1. Convert MR-ready GWAS summary statistics to TwoSampleMR
#      exposure format
#   2. Select genome-wide significant instruments
#   3. Perform LD clumping
#   4. Calculate SNP-level R2 and F statistics
#   5. Remove weak instruments with F <= 10
#
# Input:
#   data/processed/IDP_ROI_SA_CAMIDE/
#
# Output:
#   data/processed/exp_dat_3_ROI_SA_CAMIDE/
#   data/processed/exp_dat_eaf_Ffilter_ROI_SA_CAMIDE/
# ============================================================

library(here)
library(vroom)
library(data.table)
library(TwoSampleMR)

# ------------------------------------------------------------
# 1. Load the formal R2/F calculation function
# ------------------------------------------------------------

source(
  here::here(
    "functions",
    "calculate_R2_F.R"
  )
)

# ------------------------------------------------------------
# 2. Define parameters and directories
# ------------------------------------------------------------

Ffilter <- 10

input_dir <- here::here(
  "data",
  "processed",
  "IDP_ROI_SA_CAMIDE"
)

clumped_output_dir <- here::here(
  "data",
  "processed",
  "exp_dat_3_ROI_SA_CAMIDE"
)

filtered_output_dir <- here::here(
  "data",
  "processed",
  "exp_dat_eaf_Ffilter_ROI_SA_CAMIDE"
)

dir.create(
  clumped_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  filtered_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 3. Identify QC files
# ------------------------------------------------------------

qc_files <- list.files(
  input_dir,
  pattern = "_rsID_QC\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(qc_files) == 0) {
  stop("No ROI QC files were found.")
}

qc_files <- sort(qc_files)

message(
  "Number of ROI exposure files: ",
  length(qc_files)
)

# ------------------------------------------------------------
# 4. Process each ROI
# ------------------------------------------------------------

for (file_path in qc_files) {
  
  file_name <- basename(file_path)
  
  idp_name <- sub(
    "_rsID_QC\\.tsv\\.gz$",
    "",
    file_name
  )
  
  message(
    "\nProcessing: ",
    idp_name
  )
  
  # ----------------------------------------------------------
  # 4.1 Read MR-ready summary statistics
  # ----------------------------------------------------------
  
  IDP <- vroom(
    file_path,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  # ----------------------------------------------------------
  # 4.2 Convert to TwoSampleMR exposure format
  # ----------------------------------------------------------
  
  exp_dat1 <- tryCatch(
    
    format_data(
      IDP,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "EAF",
      pval_col = "P",
      samplesize_col = "SAMPLESIZE"
    ),
    
    error = function(e) {
      
      warning(
        "format_data failed for ",
        idp_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(exp_dat1)) {
    next
  }
  
  # ----------------------------------------------------------
  # 4.3 Genome-wide significance threshold
  # ----------------------------------------------------------
  
  exp_dat2 <- exp_dat1 %>%
    filter(
      !is.na(pval.exposure),
      pval.exposure < 5e-8
    )
  
  if (nrow(exp_dat2) < 3) {
    
    message(
      "Fewer than 3 genome-wide significant SNPs: ",
      idp_name
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 4.4 LD clumping
  #
  #   P-value threshold: 5e-8
  #   LD r2 threshold:   0.001
  #   Window:            10,000 kb
  #   Population:        EUR
  # ----------------------------------------------------------
  
  exp_dat3 <- tryCatch(
    
    clump_data(
      exp_dat2,
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p1 = 5e-8,
      pop = "EUR"
    ),
    
    error = function(e) {
      
      warning(
        "LD clumping failed for ",
        idp_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(exp_dat3) || nrow(exp_dat3) < 2) {
    
    message(
      "Fewer than 2 independent instruments after clumping: ",
      idp_name
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 4.5 Save clumped instruments
  # ----------------------------------------------------------
  
  clumped_file <- file.path(
    clumped_output_dir,
    paste0(
      "exp_dat_",
      idp_name,
      ".csv"
    )
  )
  
  fwrite(
    exp_dat3,
    clumped_file
  )
  
  # ----------------------------------------------------------
  # 4.6 Calculate R2 and F statistics
  # ----------------------------------------------------------
  
  exp_dat3 <- tryCatch(
    
    calculate_R2_F(exp_dat3),
    
    error = function(e) {
      
      warning(
        "R2/F calculation failed for ",
        idp_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(exp_dat3)) {
    next
  }
  
  # ----------------------------------------------------------
  # 4.7 Remove weak instruments
  # ----------------------------------------------------------
  
  if (!"F" %in% names(exp_dat3)) {
    
    warning(
      "F-statistic column not found for ",
      idp_name
    )
    
    next
  }
  
  filtered_data <- exp_dat3 %>%
    filter(
      !is.na(F),
      F > Ffilter
    )
  
  if (nrow(filtered_data) < 2) {
    
    message(
      "F > ",
      Ffilter,
      " leaves fewer than 2 instruments: ",
      idp_name
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 4.8 Save F-filtered instruments
  # ----------------------------------------------------------
  
  filtered_file <- file.path(
    filtered_output_dir,
    paste0(
      "exp_dat_",
      idp_name,
      "_Ffiltered.csv"
    )
  )
  
  fwrite(
    filtered_data,
    filtered_file
  )
  
  message(
    "Finished: ",
    idp_name,
    " | instruments after clumping = ",
    nrow(exp_dat3),
    " | instruments after F > ",
    Ffilter,
    " = ",
    nrow(filtered_data)
  )
}

message(
  "\nROI surface-area instrument selection completed."
)