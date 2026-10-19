# ============================================================
# Cortical IDP GWAS preprocessing
#
# Purpose:
#   1. Read raw CAMIDE global cortical IDP GWAS summary statistics
#   2. Add sample size according to the original GWAS coding
#   3. Match SNPs to allele-frequency information
#   4. Convert genomic variants to rsIDs using GRCh37/dbSNP144
#   5. Derive minor allele frequency (EAF)
#   6. Export MR-ready summary statistics
#
# Input:
#   data/external/IDP_Global_CAMIDE/*_meta.txt
#   data/external/IDP_Global_CAMIDE/Readme_SNP_afreq.txt
#
# Output:
#   data/processed/IDP_Global_CAMIDE/*_rsID_QC.tsv.gz
#
# Notes:
#   Raw GWAS summary statistics are not redistributed with the
#   repository. Users should obtain the original data from the
#   corresponding data source and place them in the specified
#   directory.
# ============================================================

# ------------------------------
# 1. Load packages
# ------------------------------

library(vroom)
library(dplyr)
library(MungeSumstats)
library(DelayedArray)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)


# ------------------------------
# 2. Define project directories
# ------------------------------

# All paths are project-relative.
# The repository should be opened as the project root.

input_dir <- file.path(
  "data",
  "external",
  "IDP_Global_CAMIDE"
)

output_dir <- file.path(
  "data",
  "processed",
  "IDP_Global_CAMIDE"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------
# 3. Check required input files
# ------------------------------

afreq_file <- file.path(
  input_dir,
  "Readme_SNP_afreq.txt"
)

if (!dir.exists(input_dir)) {
  stop(
    "Input directory not found: ",
    input_dir,
    "\nPlease place the CAMIDE GWAS files in the specified directory."
  )
}

if (!file.exists(afreq_file)) {
  stop(
    "Allele-frequency file not found: ",
    afreq_file
  )
}


# ------------------------------
# 4. Read allele-frequency data
# ------------------------------

message("Reading allele-frequency reference file...")

afreq <- vroom(
  afreq_file,
  col_names = TRUE,
  show_col_types = FALSE
)


# ------------------------------
# 5. Identify raw GWAS files
# ------------------------------

meta_files <- list.files(
  input_dir,
  pattern = "_meta\\.txt$",
  full.names = TRUE
)

if (length(meta_files) == 0) {
  stop(
    "No *_meta.txt files were found in: ",
    input_dir
  )
}

message(
  "Number of GWAS files identified: ",
  length(meta_files)
)


# ------------------------------
# 6. Process each global cortical IDP
# ------------------------------

for (file_path in meta_files) {
  
  message("\n----------------------------------------")
  message("Processing: ", basename(file_path))
  message("----------------------------------------")
  
  
  # ------------------------------
  # 6.1 Read GWAS summary statistics
  # ------------------------------
  
  sumstats <- vroom(
    file_path,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  
  # ------------------------------
  # 6.2 Add sample size
  # ------------------------------
  
  # The original CAMIDE GWAS encodes study sample size using
  # the variable N. Here, N is mapped to the corresponding
  # total sample size.
  
  if (!"N" %in% names(sumstats)) {
    stop(
      "Column 'N' was not found in: ",
      basename(file_path)
    )
  }
  
  sumstats <- sumstats %>%
    mutate(
      SAMPLESIZE = case_when(
        N == 1 ~ 31797,
        N == 2 ~ 36663,
        TRUE ~ NA_real_
      )
    )
  
  
  # ------------------------------
  # 6.3 Match allele-frequency information
  # ------------------------------
  
  required_sumstats_columns <- c(
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2"
  )
  
  missing_columns <- setdiff(
    required_sumstats_columns,
    names(sumstats)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "The following required columns are missing from ",
      basename(file_path),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  afreq_required_columns <- c(
    "SNP",
    "CHR",
    "POS",
    "A1",
    "A2"
  )
  
  missing_afreq_columns <- setdiff(
    afreq_required_columns,
    names(afreq)
  )
  
  if (length(missing_afreq_columns) > 0) {
    stop(
      "The following required columns are missing from the ",
      "allele-frequency reference file: ",
      paste(missing_afreq_columns, collapse = ", ")
    )
  }
  
  
  sumstats2 <- sumstats %>%
    left_join(
      afreq,
      by = c(
        "SNP",
        "CHR" = "CHR",
        "BP" = "POS",
        "A1" = "A1",
        "A2" = "A2"
      )
    )
  
  
  # ------------------------------
  # 6.4 Define output filename
  # ------------------------------
  
  base_name <- tools::file_path_sans_ext(
    basename(file_path)
  )
  
  out_file <- file.path(
    output_dir,
    paste0(
      base_name,
      "_rsID.tsv.gz"
    )
  )
  
  
  # ------------------------------
  # 6.5 Convert variants to rsIDs
  # ------------------------------
  
  message("Converting genomic variants to rsIDs...")
  
  format_sumstats(
    path = sumstats2,
    ref_genome = "GRCh37",
    dbSNP = 144,
    save_path = out_file
  )
  
  
  # ------------------------------
  # 6.6 Read converted summary statistics
  # ------------------------------
  
  sumstats_rsID <- vroom(
    out_file,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  
  # ------------------------------
  # 6.7 Derive effect-allele frequency
  # ------------------------------
  
  # FRQ is converted to the frequency of the minor allele.
  # This is used as EAF in downstream MR analyses.
  
  if (!"FRQ" %in% names(sumstats_rsID)) {
    stop(
      "Column 'FRQ' was not found after formatting: ",
      basename(out_file)
    )
  }
  
  sumstats_rsID <- sumstats_rsID %>%
    mutate(
      EAF = if_else(
        FRQ > 0.5,
        1 - FRQ,
        FRQ
      )
    )
  
  
  # ------------------------------
  # 6.8 Ensure sample size is available
  # ------------------------------
  
  # MungeSumstats may standardize column names during formatting.
  # Preserve a consistent SAMPLESIZE column for downstream
  # TwoSampleMR analyses.
  
  if (!"SAMPLESIZE" %in% names(sumstats_rsID)) {
    
    if ("N" %in% names(sumstats_rsID)) {
      
      sumstats_rsID <- sumstats_rsID %>%
        mutate(
          SAMPLESIZE = case_when(
            N == 1 ~ 31797,
            N == 2 ~ 36663,
            N > 2 ~ N,
            TRUE ~ NA_real_
          )
        )
      
    } else {
      
      # Use the original GWAS-specific sample size when the
      # standardized output no longer contains the coding variable.
      
      sample_size_value <- unique(
        sumstats$SAMPLESIZE[
          !is.na(sumstats$SAMPLESIZE)
        ]
      )
      
      if (length(sample_size_value) == 1) {
        
        sumstats_rsID$SAMPLESIZE <- sample_size_value
        
      } else {
        
        warning(
          "Unable to determine a unique sample size for ",
          basename(file_path),
          "."
        )
      }
    }
  }
  
  
  # ------------------------------
  # 6.9 Export QC-ready summary statistics
  # ------------------------------
  
  qc_file <- file.path(
    output_dir,
    paste0(
      base_name,
      "_rsID_QC.tsv.gz"
    )
  )
  
  vroom_write(
    sumstats_rsID,
    qc_file,
    delim = "\t"
  )
  
  message(
    "Saved processed summary statistics to: ",
    qc_file
  )
}


# ------------------------------
# 7. Completion message
# ------------------------------

message("\nAll cortical IDP GWAS files have been processed.")
message("Output directory: ", output_dir)
