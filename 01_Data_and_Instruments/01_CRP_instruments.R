# ============================================================
# CRP genetic instrument selection
#
# Purpose:
#   1. Read CRP GWAS summary statistics
#   2. Select genome-wide significant variants
#   3. Perform LD clumping
#   4. Calculate SNP-level R2 and F statistics
#   5. Remove weak instruments with F <= 10
#
# Input:
#   data/raw/CRP/CRP_with_EAF.tsv
#
# Output:
#   data/processed/CRP/
#       CRP_instruments_clumped.csv
#       CRP_instruments_Ffiltered.csv
# ============================================================

library(here)
library(vroom)
library(data.table)
library(dplyr)
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

p_threshold <- 5e-8
F_threshold <- 10

input_file <- here::here(
  "data",
  "raw",
  "CRP",
  "CRP_with_EAF.tsv"
)

output_dir <- here::here(
  "data",
  "processed",
  "CRP"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 3. Read CRP summary statistics
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "CRP summary-statistics file not found: ",
    input_file
  )
}

crp_dat <- vroom(
  input_file,
  col_names = TRUE,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 4. Basic input validation
# ------------------------------------------------------------

required_cols <- c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "eaf.exposure",
  "pval.exposure",
  "samplesize.exposure"
)

missing_cols <- setdiff(
  required_cols,
  names(crp_dat)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required CRP columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 5. Select genome-wide significant instruments
# ------------------------------------------------------------

crp_candidates <- crp_dat %>%
  filter(
    !is.na(pval.exposure),
    pval.exposure < p_threshold
  )

message(
  "Genome-wide significant variants: ",
  nrow(crp_candidates)
)

if (nrow(crp_candidates) < 3) {
  stop(
    "Fewer than 3 genome-wide significant CRP variants."
  )
}

# ------------------------------------------------------------
# 6. LD clumping
# ------------------------------------------------------------

crp_clumped <- tryCatch(
  
  clump_data(
    crp_candidates,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = p_threshold,
    pop = "EUR"
  ),
  
  error = function(e) {
    
    stop(
      "LD clumping failed: ",
      conditionMessage(e)
    )
  }
)

message(
  "Independent instruments after clumping: ",
  nrow(crp_clumped)
)

if (nrow(crp_clumped) < 2) {
  stop(
    "Fewer than 2 independent CRP instruments after clumping."
  )
}

# ------------------------------------------------------------
# 7. Calculate R2 and F statistics
# ------------------------------------------------------------

crp_clumped <- calculate_R2_F(
  crp_clumped
)

if (!"F" %in% names(crp_clumped)) {
  stop(
    "F-statistic was not generated."
  )
}

# ------------------------------------------------------------
# 8. Save clumped instruments before F filtering
# ------------------------------------------------------------

clumped_file <- file.path(
  output_dir,
  "CRP_instruments_clumped.csv"
)

fwrite(
  crp_clumped,
  clumped_file
)

# ------------------------------------------------------------
# 9. Remove weak instruments
# ------------------------------------------------------------

crp_final <- crp_clumped %>%
  filter(
    !is.na(F),
    F > F_threshold
  )

message(
  "Instruments with F > ",
  F_threshold,
  ": ",
  nrow(crp_final)
)

if (nrow(crp_final) < 2) {
  stop(
    "F > ",
    F_threshold,
    " leaves fewer than 2 CRP instruments."
  )
}

# ------------------------------------------------------------
# 10. Save final CRP instruments
# ------------------------------------------------------------

final_file <- file.path(
  output_dir,
  "CRP_instruments_Ffiltered.csv"
)

fwrite(
  crp_final,
  final_file
)

# ------------------------------------------------------------
# 11. Save analysis parameters
# ------------------------------------------------------------

parameters <- data.frame(
  parameter = c(
    "genome_wide_significance",
    "LD_r2",
    "LD_window_kb",
    "F_threshold",
    "instruments_before_clumping",
    "instruments_after_clumping",
    "instruments_after_F_filter"
  ),
  
  value = c(
    p_threshold,
    0.001,
    10000,
    F_threshold,
    nrow(crp_candidates),
    nrow(crp_clumped),
    nrow(crp_final)
  )
)

write.csv(
  parameters,
  file.path(
    output_dir,
    "CRP_instrument_selection_parameters.csv"
  ),
  row.names = FALSE
)

message(
  "\nCRP instrument selection completed."
)