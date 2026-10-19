# ============================================================
# Cortical IDP genetic instrument selection
#
# Purpose:
#   1. Read processed global cortical IDP GWAS summary statistics
#   2. Select genome-wide significant SNPs
#   3. Perform LD clumping
#   4. Calculate SNP-level R2 and F statistics
#   5. Exclude weak instruments with F <= 10
#   6. Export MR-ready exposure datasets
#
# Input:
#   data/processed/IDP_Global_CAMIDE/*_rsID_QC.tsv.gz
#
# Output:
#   results/instruments/global_idp/
#       - All clumped instruments
#       - F-filtered instruments
#
# Instrument selection:
#   P < 5e-8
#   LD r2 < 0.001
#   Clumping window = 10,000 kb
#   F > 10
# ============================================================


# ------------------------------
# 1. Load packages
# ------------------------------

library(vroom)
library(data.table)
library(TwoSampleMR)


# ------------------------------
# 2. Load project-specific functions
# ------------------------------

source(
  file.path(
    "functions",
    "calculate_R2_F.R"
  )
)


# ------------------------------
# 3. Define analysis parameters
# ------------------------------

p_threshold <- 5e-8
clump_r2 <- 0.001
clump_kb <- 10000
f_threshold <- 10

population <- "EUR"


# ------------------------------
# 4. Define directories
# ------------------------------

input_dir <- file.path(
  "data",
  "processed",
  "IDP_Global_CAMIDE"
)

output_clumped_dir <- file.path(
  "results",
  "instruments",
  "global_idp",
  "clumped"
)

output_filtered_dir <- file.path(
  "results",
  "instruments",
  "global_idp",
  "F_filtered"
)

dir.create(
  output_clumped_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_filtered_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------
# 5. Check input directory
# ------------------------------

if (!dir.exists(input_dir)) {
  stop(
    "Input directory not found: ",
    input_dir
  )
}


# ------------------------------
# 6. Identify processed GWAS files
# ------------------------------

qc_files <- list.files(
  input_dir,
  pattern = "_rsID_QC\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(qc_files) == 0) {
  stop(
    "No processed cortical IDP files were found in: ",
    input_dir
  )
}

message(
  "Number of cortical IDP files identified: ",
  length(qc_files)
)


# ------------------------------
# 7. Process each cortical IDP
# ------------------------------

for (file_path in qc_files) {
  
  file_name <- basename(file_path)
  
  idp_name <- sub(
    "_meta_rsID_QC\\.tsv\\.gz$",
    "",
    file_name
  )
  
  message("\n========================================")
  message("Processing IDP: ", idp_name)
  message("========================================")
  
  
  # ------------------------------
  # 7.1 Read summary statistics
  # ------------------------------
  
  IDP <- vroom(
    file_path,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  
  # ------------------------------
  # 7.2 Check required columns
  # ------------------------------
  
  required_columns <- c(
    "SNP",
    "BETA",
    "SE",
    "A1",
    "A2",
    "EAF",
    "P",
    "SAMPLESIZE"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(IDP)
  )
  
  if (length(missing_columns) > 0) {
    
    warning(
      "Skipping ",
      idp_name,
      ". Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
    
    next
  }
  
  
  # ------------------------------
  # 7.3 Format exposure data
  # ------------------------------
  
  exp_dat1 <- format_data(
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
  )
  
  
  # ------------------------------
  # 7.4 Select genome-wide significant SNPs
  # ------------------------------
  
  exp_dat2 <- exp_dat1[
    !is.na(exp_dat1$pval.exposure) &
      exp_dat1$pval.exposure < p_threshold,
  ]
  
  message(
    "Genome-wide significant SNPs: ",
    nrow(exp_dat2)
  )
  
  if (nrow(exp_dat2) < 3) {
    
    warning(
      "Fewer than 3 genome-wide significant SNPs for ",
      idp_name,
      ". Skipping."
    )
    
    next
  }
  
  
  # ------------------------------
  # 7.5 LD clumping
  # ------------------------------
  
  # Explicitly match the LD parameters specified in the Methods.
  
  exp_dat3 <- tryCatch(
    
    clump_data(
      exp_dat2,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      pop = population
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
    
    warning(
      "Insufficient SNPs after LD clumping for ",
      idp_name,
      ". Skipping."
    )
    
    next
  }
  
  message(
    "SNPs after LD clumping: ",
    nrow(exp_dat3)
  )
  
  
  # ------------------------------
  # 7.6 Calculate R2 and F statistics
  # ------------------------------
  
  exp_dat3 <- calculate_R2_F(
    exp_dat3
  )
  
  
  # ------------------------------
  # 7.7 Export all clumped instruments
  # ------------------------------
  
  clumped_file <- file.path(
    output_clumped_dir,
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
  
  
  # ------------------------------
  # 7.8 Exclude weak instruments
  # ------------------------------
  
  if (!"F" %in% names(exp_dat3)) {
    
    stop(
      "Column 'F' was not generated by calculate_R2_F() for ",
      idp_name
    )
  }
  
  filtered_data <- exp_dat3[
    !is.na(exp_dat3$F) &
      exp_dat3$F > f_threshold,
  ]
  
  message(
    "SNPs with F > ",
    f_threshold,
    ": ",
    nrow(filtered_data)
  )
  
  
  # ------------------------------
  # 7.9 Export F-filtered instruments
  # ------------------------------
  
  if (nrow(filtered_data) < 2) {
    
    warning(
      "Fewer than 2 SNPs passed the F-statistic threshold for ",
      idp_name,
      ". No F-filtered dataset will be exported."
    )
    
    next
  }
  
  filtered_file <- file.path(
    output_filtered_dir,
    paste0(
      "exp_dat_",
      idp_name,
      "_Ffiltered.csv"
    )
  )
  
  # Important:
  # Export the F-filtered dataset rather than the unfiltered
  # clumped dataset.
  
  fwrite(
    filtered_data,
    filtered_file
  )
  
  message(
    "Saved F-filtered instruments to: ",
    filtered_file
  )
}


# ------------------------------
# 8. Completion message
# ------------------------------

message("\nCortical IDP instrument construction completed.")
message(
  "F-filtered instruments are stored in: ",
  output_filtered_dir
)
