# ============================================================
# ROI cortical surface area GWAS preprocessing
#
# Purpose:
#   1. Read CAM-CAN/IDE ROI surface-area GWAS summary statistics
#   2. Add sample-size information
#   3. Add SNP frequency information
#   4. Convert genomic coordinates to rsIDs using MungeSumstats
#   5. Generate MR-ready summary statistics
#
# Input:
#   data/raw/CAM_IDE_Suface_Area/
#
# Output:
#   data/processed/IDP_ROI_SA_CAMIDE/
#
# Notes:
#   - Raw GWAS files are not included in the public repository.
#   - All paths are project-root relative.
#   - EAF represents the effect-allele frequency.
# ============================================================

library(here)
library(vroom)
library(dplyr)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# ------------------------------------------------------------
# 1. Define directories
# ------------------------------------------------------------

input_dir <- here::here(
  "data",
  "raw",
  "CAM_IDE_Suface_Area"
)

output_dir <- here::here(
  "data",
  "processed",
  "IDP_ROI_SA_CAMIDE"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Read SNP frequency information
# ------------------------------------------------------------

afreq_file <- file.path(
  input_dir,
  "Readme_SNP_afreq.txt"
)

if (!file.exists(afreq_file)) {
  stop("SNP frequency file not found: ", afreq_file)
}

afreq <- vroom(
  afreq_file,
  col_names = TRUE,
  show_col_types = FALSE
)

required_afreq_cols <- c(
  "SNP",
  "CHR",
  "POS",
  "A1",
  "A2"
)

missing_afreq_cols <- setdiff(
  required_afreq_cols,
  names(afreq)
)

if (length(missing_afreq_cols) > 0) {
  stop(
    "Missing required columns in SNP frequency file: ",
    paste(missing_afreq_cols, collapse = ", ")
  )
}

# AF1 is preferred because it represents the frequency of A1.
if (!"AF1" %in% names(afreq)) {
  
  if ("FRQ" %in% names(afreq)) {
    afreq <- afreq %>%
      rename(AF1 = FRQ)
  } else {
    stop(
      "Neither AF1 nor FRQ was found in the SNP frequency file."
    )
  }
}

# Keep only the information needed for allele-frequency recovery.
afreq <- afreq %>%
  select(
    SNP,
    CHR,
    POS,
    A1,
    A2,
    AF1
  ) %>%
  distinct()

# ------------------------------------------------------------
# 3. Identify ROI GWAS files
# ------------------------------------------------------------

meta_files <- list.files(
  input_dir,
  pattern = "^plinkmeta[0-9]+\\.txt$",
  full.names = TRUE
)

if (length(meta_files) == 0) {
  stop(
    "No ROI GWAS files matching 'plinkmeta[0-9]+.txt' were found."
  )
}

meta_files <- sort(meta_files)

message(
  "Number of ROI GWAS files detected: ",
  length(meta_files)
)

# ------------------------------------------------------------
# 4. Process each ROI GWAS
# ------------------------------------------------------------

for (file_path in meta_files) {
  
  file_name <- basename(file_path)
  
  message("\nProcessing: ", file_name)
  
  # ----------------------------------------------------------
  # 4.1 Read GWAS summary statistics
  # ----------------------------------------------------------
  
  sumstats <- vroom(
    file_path,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  required_sumstats_cols <- c(
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "BETA",
    "SE",
    "P",
    "N"
  )
  
  missing_sumstats_cols <- setdiff(
    required_sumstats_cols,
    names(sumstats)
  )
  
  if (length(missing_sumstats_cols) > 0) {
    
    warning(
      "Skipping ",
      file_name,
      ". Missing columns: ",
      paste(missing_sumstats_cols, collapse = ", ")
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 4.2 Convert encoded N values to actual sample sizes
  # ----------------------------------------------------------
  
  sumstats <- sumstats %>%
    mutate(
      SAMPLESIZE = case_when(
        N == 1 ~ 31797,
        N == 2 ~ 36663,
        TRUE ~ NA_real_
      )
    )
  
  if (all(is.na(sumstats$SAMPLESIZE))) {
    
    warning(
      "No valid sample-size mapping was found for ",
      file_name
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # 4.3 Add allele-frequency information
  # ----------------------------------------------------------
  
  sumstats2 <- sumstats %>%
    left_join(
      afreq,
      by = c(
        "SNP",
        "CHR",
        "BP" = "POS",
        "A1",
        "A2"
      )
    )
  
  # ----------------------------------------------------------
  # 4.4 Write temporary input for MungeSumstats
  # ----------------------------------------------------------
  
  temp_input <- tempfile(
    pattern = "roi_sa_",
    fileext = ".tsv.gz"
  )
  
  vroom_write(
    sumstats2,
    temp_input,
    delim = "\t"
  )
  
  # ----------------------------------------------------------
  # 4.5 Convert genomic coordinates to rsIDs
  # ----------------------------------------------------------
  
  base_name <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(file_name)
  )
  
  rsid_file <- file.path(
    output_dir,
    paste0(base_name, "_rsID.tsv.gz")
  )
  
  format_sumstats(
    path = temp_input,
    ref_genome = "GRCh37",
    dbSNP = 144,
    save_path = rsid_file
  )
  
  unlink(temp_input)
  
  # ----------------------------------------------------------
  # 4.6 Read converted summary statistics
  # ----------------------------------------------------------
  
  sumstats_rsID <- vroom(
    rsid_file,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  # ----------------------------------------------------------
  # 4.7 Recover sample size
  # ----------------------------------------------------------
  
  sample_size_lookup <- sumstats %>%
    select(SNP, SAMPLESIZE) %>%
    distinct(SNP, .keep_all = TRUE)
  
  sumstats_rsID <- sumstats_rsID %>%
    select(-any_of("SAMPLESIZE")) %>%
    left_join(
      sample_size_lookup,
      by = "SNP"
    )
  
  # ----------------------------------------------------------
  # 4.8 Recover effect-allele frequency
  #
  # Important:
  #   EAF must correspond to the effect allele, not simply
  #   the minor allele frequency.
  # ----------------------------------------------------------
  
  allele_lookup <- sumstats %>%
    select(
      SNP,
      A1,
      A2,
      AF1
    ) %>%
    distinct()
  
  sumstats_rsID <- sumstats_rsID %>%
    left_join(
      allele_lookup %>%
        rename(
          source_A1 = A1,
          source_A2 = A2,
          source_AF1 = AF1
        ),
      by = "SNP"
    ) %>%
    mutate(
      EAF = case_when(
        A1 == source_A1 &
          A2 == source_A2 ~ source_AF1,
        
        A1 == source_A2 &
          A2 == source_A1 ~ 1 - source_AF1,
        
        TRUE ~ suppressWarnings(as.numeric(FRQ))
      )
    ) %>%
    select(
      -source_A1,
      -source_A2,
      -source_AF1
    )
  
  # ----------------------------------------------------------
  # 4.9 Basic QC checks
  # ----------------------------------------------------------
  
  if (!"EAF" %in% names(sumstats_rsID)) {
    stop(
      "EAF could not be generated for: ",
      file_name
    )
  }
  
  invalid_eaf <- sum(
    !is.na(sumstats_rsID$EAF) &
      (sumstats_rsID$EAF < 0 |
         sumstats_rsID$EAF > 1)
  )
  
  if (invalid_eaf > 0) {
    
    warning(
      invalid_eaf,
      " variants have invalid EAF values in ",
      file_name
    )
    
    sumstats_rsID$EAF[
      sumstats_rsID$EAF < 0 |
        sumstats_rsID$EAF > 1
    ] <- NA_real_
  }
  
  # ----------------------------------------------------------
  # 4.10 Save MR-ready QC file
  # ----------------------------------------------------------
  
  qc_file <- file.path(
    output_dir,
    paste0(base_name, "_rsID_QC.tsv.gz")
  )
  
  vroom_write(
    sumstats_rsID,
    qc_file,
    delim = "\t"
  )
  
  message(
    "Saved: ",
    qc_file
  )
}

message("\nROI surface-area preprocessing completed.")