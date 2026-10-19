# ============================================================
# Cortical IDP–VaD Genetic Correlation Using LDSC
# ============================================================
#
# Purpose:
#   Estimate genetic correlations between prespecified
#   cortical surface-area phenotypes and subcortical vascular
#   dementia (VaD) using LD score regression.
#
# Candidate IDPs:
#   - surface_area-SA_plinkmeta101
#   - surface_area-SA_plinkmeta105
#
# These IDPs were prespecified based on the preceding MR
# localization and mediation workflow.
#
# Notes:
#   - Cortical IDPs are continuous traits.
#   - Subcortical VaD is a binary trait and is converted to the
#     liability scale using sample and population prevalence.
#   - HapMap3 SNPs and European LD scores are used.
#   - Cross-trait LDSC intercepts are retained for QC.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(GenomicSEM)
library(dplyr)
library(vroom)
library(R.utils)
library(here)


# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

# Subcortical VaD sample information
VD_CASES <- 1031
VD_CONTROLS <- 474702

VD_TOTAL_N <- VD_CASES + VD_CONTROLS

VD_SAMPLE_PREVALENCE <- VD_CASES / VD_TOTAL_N

# Prespecified population prevalence assumption.
#
# IMPORTANT:
# This value should be justified in the Methods or
# Supplementary Methods.
VD_POPULATION_PREVALENCE <- 0.002


# ------------------------------------------------------------
# 3. Project paths
# ------------------------------------------------------------

# Raw FinnGen summary statistics
vd_file <- here::here(
  "data",
  "raw",
  "finngen_R12_VD_SC.gz"
)

# Munged cortical IDP summary statistics
munged_dir <- here::here(
  "results",
  "LDSC",
  "munged"
)

# Munged VaD summary statistics
vd_sumstats <- file.path(
  munged_dir,
  "VD_SC.sumstats.gz"
)

# European LD scores
ld_ref_dir <- here::here(
  "data",
  "reference",
  "eur_w_ld_chr"
)

hm3_snplist <- file.path(
  ld_ref_dir,
  "w_hm3.snplist"
)

# LDSC output directory
output_dir <- here::here(
  "results",
  "LDSC"
)

dir.create(
  munged_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Prespecified candidate IDPs
# ------------------------------------------------------------

candidate_idps <- data.frame(
  
  IDP = c(
    "surface_area-SA_plinkmeta101",
    "surface_area-SA_plinkmeta105"
  ),
  
  file = c(
    file.path(
      munged_dir,
      "surface_area-SA_plinkmeta101.sumstats.gz"
    ),
    file.path(
      munged_dir,
      "surface_area-SA_plinkmeta105.sumstats.gz"
    )
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Input validation
# ------------------------------------------------------------

required_files <- c(
  vd_file,
  vd_sumstats,
  hm3_snplist
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "The following required files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

missing_idp_files <- candidate_idps$file[
  !file.exists(candidate_idps$file)
]

if (length(missing_idp_files) > 0) {
  
  stop(
    "The following candidate IDP LDSC files were not found:\n",
    paste(
      missing_idp_files,
      collapse = "\n"
    ),
    "\nPlease run the CRP–IDP LDSC munging workflow first."
  )
}

if (!dir.exists(ld_ref_dir)) {
  
  stop(
    "LD score reference directory was not found: ",
    ld_ref_dir
  )
}


# ------------------------------------------------------------
# 6. Basic parameter checks
# ------------------------------------------------------------

if (
  VD_CASES <= 0 ||
  VD_CONTROLS <= 0
) {
  
  stop(
    "VD_CASES and VD_CONTROLS must both be positive."
  )
}

if (
  VD_POPULATION_PREVALENCE <= 0 ||
  VD_POPULATION_PREVALENCE >= 1
) {
  
  stop(
    "VD_POPULATION_PREVALENCE must be between 0 and 1."
  )
}

message(
  "Subcortical VaD sample prevalence: ",
  round(
    VD_SAMPLE_PREVALENCE,
    6
  )
)

message(
  "Subcortical VaD population prevalence: ",
  VD_POPULATION_PREVALENCE
)


# ------------------------------------------------------------
# 7. Prepare and munge subcortical VaD
# ------------------------------------------------------------

prep_and_munge_vdsc <- function(
    vd_path,
    hm3,
    output_dir,
    out_prefix = "VD_SC",
    sample_size
) {
  
  message(
    "Preparing ",
    out_prefix,
    " summary statistics..."
  )
  
  vd <- vroom(
    vd_path,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "rsids",
    "#chrom",
    "pos",
    "alt",
    "ref",
    "beta",
    "sebeta",
    "pval"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(vd)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Missing VaD columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  vd <- vd %>%
    transmute(
      
      SNP = rsids,
      
      CHR = `#chrom`,
      
      BP = pos,
      
      A1 = toupper(alt),
      
      A2 = toupper(ref),
      
      Z = beta / sebeta,
      
      P = pval,
      
      N = sample_size
      
    ) %>%
    
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(A1),
      !is.na(A2),
      A1 %in% c(
        "A",
        "C",
        "G",
        "T"
      ),
      A2 %in% c(
        "A",
        "C",
        "G",
        "T"
      ),
      A1 != A2,
      !is.na(Z),
      is.finite(Z),
      !is.na(P),
      is.finite(P),
      P > 0,
      P <= 1,
      !(
        (A1 == "A" & A2 == "T") |
          (A1 == "T" & A2 == "A") |
          (A1 == "C" & A2 == "G") |
          (A1 == "G" & A2 == "C")
      )
    ) %>%
    
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  tmpfile <- file.path(
    output_dir,
    paste0(
      out_prefix,
      ".tmp"
    )
  )
  
  write.table(
    vd,
    file = tmpfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  gzip(
    tmpfile,
    overwrite = TRUE
  )
  
  GenomicSEM::munge(
    files = paste0(
      tmpfile,
      ".gz"
    ),
    hm3 = hm3,
    trait.names = out_prefix,
    N = NULL
  )
  
  unlink(
    c(
      tmpfile,
      paste0(
        tmpfile,
        ".gz"
      )
    ),
    force = TRUE
  )
  
  return(
    paste0(
      out_prefix,
      ".sumstats.gz"
    )
  )
}


# ------------------------------------------------------------
# 8. Generate VaD LDSC summary statistics
# ------------------------------------------------------------

if (!file.exists(vd_sumstats)) {
  
  prep_and_munge_vdsc(
    vd_path = vd_file,
    hm3 = hm3_snplist,
    output_dir = munged_dir,
    out_prefix = "VD_SC",
    sample_size = VD_TOTAL_N
  )
}


# ------------------------------------------------------------
# 9. Pairwise LDSC function
# ------------------------------------------------------------

run_rg <- function(
    trait1,
    trait1_name,
    trait2,
    trait2_name
) {
  
  message(
    "Running LDSC: ",
    trait1_name,
    " <-> ",
    trait2_name
  )
  
  res <- GenomicSEM::ldsc(
    
    traits = c(
      trait1,
      trait2
    ),
    
    sample.prev = c(
      NA,
      VD_SAMPLE_PREVALENCE
    ),
    
    population.prev = c(
      NA,
      VD_POPULATION_PREVALENCE
    ),
    
    ld = ld_ref_dir,
    
    wld = ld_ref_dir,
    
    trait.names = c(
      trait1_name,
      trait2_name
    ),
    
    stand = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Genetic correlation
  # ----------------------------------------------------------
  
  rg_est <- res$S_Stand[1, 2]
  
  rg_se <- sqrt(
    res$V_Stand[2, 2]
  )
  
  rg_z <- rg_est / rg_se
  
  rg_p <- 2 * pnorm(
    -abs(rg_z)
  )
  
  
  # ----------------------------------------------------------
  # Cross-trait LDSC intercept
  # ----------------------------------------------------------
  
  intercept_est <- res$I[1, 2]
  
  intercept_se <- sqrt(
    res$V[2, 2]
  )
  
  
  # ----------------------------------------------------------
  # Number of SNPs
  # ----------------------------------------------------------
  
  n_snp <- res$m
  
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  
  data.frame(
    
    Trait1 = trait1_name,
    
    Trait2 = trait2_name,
    
    rg = rg_est,
    
    SE = rg_se,
    
    Z = rg_z,
    
    P = rg_p,
    
    cross_trait_intercept =
      intercept_est,
    
    cross_trait_intercept_SE =
      intercept_se,
    
    n_SNP = n_snp,
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 10. Run LDSC for prespecified candidate IDPs
# ------------------------------------------------------------

results_list <- vector(
  "list",
  nrow(candidate_idps)
)

for (i in seq_len(nrow(candidate_idps))) {
  
  results_list[[i]] <- run_rg(
    
    trait1 =
      candidate_idps$file[i],
    
    trait1_name =
      candidate_idps$IDP[i],
    
    trait2 =
      vd_sumstats,
    
    trait2_name =
      "VD_SC"
    
  )
}


# ------------------------------------------------------------
# 11. Combine results
# ------------------------------------------------------------

results <- bind_rows(
  results_list
)


# ------------------------------------------------------------
# 12. Save results
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "IDP101_105_VD_SC_LDSC.csv"
)

write.csv(
  results,
  output_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Print results
# ------------------------------------------------------------

print(results)

message(
  "IDP–subcortical VaD LDSC analysis completed."
)

message(
  "Results saved to: ",
  output_file
)