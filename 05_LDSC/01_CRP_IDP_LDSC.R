# ============================================================
# CRP–Cortical IDP Genetic Correlation Using LDSC
# ============================================================
#
# Purpose:
#   Estimate pairwise genetic correlations between genetically
#   predicted C-reactive protein (CRP) and prespecified cortical
#   imaging-derived phenotypes (IDPs) using LD score regression.
#
# Analysis:
#   CRP <-> cortical IDP
#
# Notes:
#   - Summary statistics are munged using GenomicSEM::munge().
#   - HapMap3 SNPs are used for LDSC.
#   - European LD scores are used.
#   - Cross-trait LDSC intercepts are retained for QC.
#   - Raw GWAS files should not be committed to the repository.
# ============================================================

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(GenomicSEM)
library(dplyr)
library(stringr)
library(vroom)
library(R.utils)
library(here)

# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

CRP_N <- 575531

# Project-relative paths
idp_dir <- here::here(
  "data",
  "raw",
  "CAMIDE_GWAS_for_VaD"
)

crp_file <- here::here(
  "data",
  "raw",
  "CRP_GCST90029070_buildGRCh37.tsv.gz"
)

ld_ref_dir <- here::here(
  "data",
  "reference",
  "eur_w_ld_chr"
)

hm3_snplist <- file.path(
  ld_ref_dir,
  "w_hm3.snplist"
)

# Output directories
ldsc_dir <- here::here(
  "results",
  "LDSC"
)

munged_dir <- file.path(
  ldsc_dir,
  "munged"
)

dir.create(
  ldsc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  munged_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 3. Identify cortical IDP summary-statistic files
# ------------------------------------------------------------

idp_files <- list.files(
  idp_dir,
  pattern = "_rsID_QC\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(idp_files) == 0) {
  stop(
    "No cortical IDP summary-statistic files were found in: ",
    idp_dir
  )
}

message(
  "Number of cortical IDP files identified: ",
  length(idp_files)
)

# ------------------------------------------------------------
# 4. Input validation
# ------------------------------------------------------------

required_files <- c(
  crp_file,
  hm3_snplist
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (!dir.exists(ld_ref_dir)) {
  stop(
    "LD score reference directory was not found: ",
    ld_ref_dir
  )
}

# ------------------------------------------------------------
# 5. Input-level QC function
# ------------------------------------------------------------

filter_for_munge <- function(df) {
  
  df %>%
    mutate(
      SNP = as.character(SNP),
      A1 = toupper(as.character(A1)),
      A2 = toupper(as.character(A2))
    ) %>%
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(A1),
      !is.na(A2),
      A1 %in% c("A", "C", "G", "T"),
      A2 %in% c("A", "C", "G", "T"),
      A1 != A2,
      !is.na(Z),
      is.finite(Z),
      !is.na(P),
      is.finite(P),
      P > 0,
      P <= 1,
      !is.na(N),
      is.finite(N),
      N > 0
    ) %>%
    filter(
      !(
        (A1 == "A" & A2 == "T") |
          (A1 == "T" & A2 == "A") |
          (A1 == "C" & A2 == "G") |
          (A1 == "G" & A2 == "C")
      )
    ) %>%
    distinct(SNP, .keep_all = TRUE)
}

# ------------------------------------------------------------
# 6. Prepare and munge CRP summary statistics
# ------------------------------------------------------------

prep_and_munge_crp <- function(
    crp_path,
    hm3,
    output_dir,
    out_prefix = "CRP",
    sample_size
) {
  
  message("Preparing CRP summary statistics...")
  
  crp <- vroom(
    crp_path,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "variant_id",
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "standard_error",
    "p_value"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(crp)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing CRP columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  crp <- crp %>%
    transmute(
      SNP = variant_id,
      CHR = chromosome,
      BP = base_pair_location,
      A1 = toupper(effect_allele),
      A2 = toupper(other_allele),
      BETA = beta,
      SE = standard_error,
      P = p_value
    ) %>%
    mutate(
      Z = BETA / SE,
      N = sample_size
    ) %>%
    filter_for_munge() %>%
    select(
      SNP,
      CHR,
      BP,
      A1,
      A2,
      Z,
      N,
      P
    )
  
  tmpfile <- file.path(
    output_dir,
    paste0(out_prefix, ".tmp")
  )
  
  write.table(
    crp,
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
    files = paste0(tmpfile, ".gz"),
    hm3 = hm3,
    trait.names = out_prefix,
    N = NULL
  )
  
  unlink(
    c(
      tmpfile,
      paste0(tmpfile, ".gz")
    ),
    force = TRUE
  )
  
  return(
    paste0(out_prefix, ".sumstats.gz")
  )
}

# ------------------------------------------------------------
# 7. Prepare and munge cortical IDP summary statistics
# ------------------------------------------------------------

prep_and_munge_idp <- function(
    idp_path,
    hm3,
    output_dir,
    out_prefix
) {
  
  message(
    "Preparing IDP: ",
    out_prefix
  )
  
  idp <- vroom(
    idp_path,
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "BETA",
    "SE",
    "P"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(idp)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing IDP columns for ",
      out_prefix,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  if (!"SAMPLESIZE" %in% names(idp)) {
    stop(
      "SAMPLESIZE column is missing for IDP: ",
      out_prefix
    )
  }
  
  idp <- idp %>%
    transmute(
      SNP,
      CHR,
      BP,
      A1 = toupper(A1),
      A2 = toupper(A2),
      Z = BETA / SE,
      N = SAMPLESIZE,
      P
    ) %>%
    filter_for_munge()
  
  tmpfile <- file.path(
    output_dir,
    paste0(out_prefix, ".tmp")
  )
  
  write.table(
    idp,
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
    files = paste0(tmpfile, ".gz"),
    hm3 = hm3,
    trait.names = out_prefix,
    N = NULL
  )
  
  unlink(
    c(
      tmpfile,
      paste0(tmpfile, ".gz")
    ),
    force = TRUE
  )
  
  return(
    paste0(out_prefix, ".sumstats.gz")
  )
}

# ------------------------------------------------------------
# 8. Munging CRP
# ------------------------------------------------------------

crp_sumstats <- prep_and_munge_crp(
  crp_path = crp_file,
  hm3 = hm3_snplist,
  output_dir = munged_dir,
  out_prefix = "CRP",
  sample_size = CRP_N
)

# ------------------------------------------------------------
# 9. Initialize results table
# ------------------------------------------------------------

ldsc_results <- data.frame(
  IDP = character(),
  rg = numeric(),
  SE = numeric(),
  Z = numeric(),
  P = numeric(),
  cross_trait_intercept = numeric(),
  cross_trait_intercept_SE = numeric(),
  n_SNP = integer(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. Pairwise CRP–IDP LDSC
# ------------------------------------------------------------

for (f in idp_files) {
  
  idp_name <- basename(f)
  
  idp_name <- str_remove(
    idp_name,
    "_rsID_QC\\.tsv\\.gz$"
  )
  
  message(
    "--------------------------------------------"
  )
  
  message(
    "Processing: ",
    idp_name
  )
  
  tryCatch({
    
    idp_sumstats <- prep_and_munge_idp(
      idp_path = f,
      hm3 = hm3_snplist,
      output_dir = munged_dir,
      out_prefix = idp_name
    )
    
    res <- GenomicSEM::ldsc(
      traits = c(
        idp_sumstats,
        crp_sumstats
      ),
      sample.prev = c(
        NA,
        NA
      ),
      population.prev = c(
        NA,
        NA
      ),
      ld = ld_ref_dir,
      wld = ld_ref_dir,
      trait.names = c(
        idp_name,
        "CRP"
      ),
      stand = TRUE
    )
    
    # --------------------------------------------------------
    # Genetic correlation
    # --------------------------------------------------------
    
    rg_est <- res$S_Stand[1, 2]
    
    rg_se <- sqrt(
      res$V_Stand[2, 2]
    )
    
    rg_z <- rg_est / rg_se
    
    rg_p <- 2 * pnorm(
      -abs(rg_z)
    )
    
    # --------------------------------------------------------
    # Cross-trait LDSC intercept
    # --------------------------------------------------------
    
    intercept_est <- res$I[1, 2]
    
    # The corresponding sampling variance of the
    # cross-trait intercept is stored in the LDSC
    # covariance matrix V.
    #
    # For two traits, the covariance parameter is the
    # second element in the lower-triangular ordering.
    
    intercept_se <- sqrt(
      res$V[2, 2]
    )
    
    # --------------------------------------------------------
    # Number of SNPs retained for LDSC
    # --------------------------------------------------------
    
    n_snp <- res$m
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    ldsc_results <- rbind(
      ldsc_results,
      data.frame(
        IDP = idp_name,
        rg = rg_est,
        SE = rg_se,
        Z = rg_z,
        P = rg_p,
        cross_trait_intercept = intercept_est,
        cross_trait_intercept_SE = intercept_se,
        n_SNP = n_snp,
        stringsAsFactors = FALSE
      )
    )
    
    message(
      sprintf(
        "%s: rg = %.4f, SE = %.4f, P = %.4g",
        idp_name,
        rg_est,
        rg_se,
        rg_p
      )
    )
    
  }, error = function(e) {
    
    warning(
      "LDSC failed for ",
      idp_name,
      ": ",
      conditionMessage(e)
    )
    
  })
  
  gc()
}

# ------------------------------------------------------------
# 11. Save results
# ------------------------------------------------------------

output_file <- file.path(
  ldsc_dir,
  "CRP_IDP_LDSC_results.csv"
)

write.csv(
  ldsc_results,
  output_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Print results
# ------------------------------------------------------------

print(ldsc_results)

message(
  "CRP–IDP LDSC analysis completed."
)

message(
  "Results saved to: ",
  output_file
)