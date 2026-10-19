# Genetic Associations Between C-Reactive Protein, Cortical Surface Area, and Vascular Dementia

## Summary

This repository contains the R code used for a two-sample Mendelian randomization (MR) study investigating the genetic associations between C-reactive protein (CRP), cortical imaging phenotypes, and vascular dementia (VaD).

The analytical framework includes univariable MR, two-step MR mediation, reverse MR, multivariable MR, linkage disequilibrium score regression (LDSC), and sensitivity analyses.

All analyses use publicly available GWAS summary statistics.

## File List

```text
01_Data_and_Instruments/
├── 01_CRP_instruments.R
├── 02_Cortical_Global_IDP_preprocessing.R
├── 02_Cortical_Global_IDP_instruments.R
├── 02_ROI_SA_sumstats_preprocessing.R
└── 03_ROI_SA_instrument_selection.R

02_Univariable_MR/
├── 01_CRP_to_VaD.R
├── 02_Global_IDP_to_VaD.R
├── 03_ROI_SA_to_VaD.R
└── 04_CRP_to_Candidate_IDP.R

03_Mediation/
├── 01_Two_step_MR_mediation.R
└── 02_Reverse_MR.R

04_Multivariable_MR/
└── 01_CRP_Covariates_on_VD_MVMR.R

05_LDSC/
├── 01_CRP_IDP_LDSC.R
├── 02_CRP_VaD_LDSC.R
└── 03_IDP_VaD_LDSC.R

06_Sensitivity_Analysis/
├── 01_MR_Robustness.R
├── 02_MR_PRESSO.R
├── 03_Heterogeneity.R
├── 04_Pleiotropy.R
├── 05_Leave_one_out.R
└── 06_Steiger.R
```

## Analysis Workflow

### Step 1. Genetic Instrument Preparation

GWAS summary statistics were processed and genetic instruments were selected using genome-wide significance and LD clumping criteria.

* P < 5 × 10^-8
* LD r² < 0.001
* Clumping window = 10,000 kb
* European ancestry

### Step 2. Univariable MR

The following associations were evaluated:

* CRP → VaD
* Global cortical phenotypes → VaD
* Regional cortical phenotypes → VaD
* CRP → candidate cortical phenotypes

Random-effects inverse-variance weighted (IVW) MR was used as the primary analysis.

Sensitivity analyses included MR-Egger, weighted median, weighted mode, MR-RAPS, and MR-PRESSO.

### Step 3. Two-Step MR Mediation

Potential mediation was assessed using:

CRP → Cortical Surface Area → VaD


The indirect effect was estimated using the product of the two MR estimates, with the standard error calculated using the delta method.

### Step 4. Reverse MR

Reverse-direction MR analyses were performed where sufficient genetic instruments were available to assess potential reverse causality.

### Step 5. Multivariable MR

Separate two-exposure MVMR analyses were performed to evaluate whether the CRP–VaD association was independent of prespecified metabolic, vascular, and lifestyle-related traits.

### Step 6. LDSC

Pairwise genetic correlations were assessed among:

* CRP and cortical phenotypes
* CRP and VaD
* Candidate cortical phenotypes and VaD

Cross-trait LDSC intercepts were also evaluated.

### Step 7. Sensitivity Analyses

| File                 | Analysis                                          |
| -------------------- | ------------------------------------------------- |
| `01_MR_Robustness.R` | MR-Egger, weighted median, weighted mode, MR-RAPS |
| `02_MR_PRESSO.R`     | MR-PRESSO                                         |
| `03_Heterogeneity.R` | Cochran's Q                                       |
| `04_Pleiotropy.R`    | MR-Egger intercept                                |
| `05_Leave_one_out.R` | Leave-one-out analysis                            |
| `06_Steiger.R`       | Steiger directionality assessment                 |

These analyses use the same harmonized datasets as the primary MR analyses and do not independently re-select or re-harmonize genetic instruments.

## Software

Analyses were performed using R 4.5.2.

Major packages include:

* `TwoSampleMR` 0.6.29
* `mr.raps` 0.4.3
* `MRPRESSO` 1.0
* `MendelianRandomization` 0.10.0
* `MVMR` 0.4.2
* `GenomicSEM` 0.0.5

The computational environment is recorded in `renv.lock`.

## Data Availability

The GWAS summary statistics used in this study were obtained from publicly available resources.

Raw GWAS summary statistics are not included in this repository. Users should obtain the corresponding datasets from the original data sources and comply with their respective data-use requirements.

## Reproducibility

The scripts are organized according to the analytical workflow described above.

Project-relative paths are used to facilitate reproducibility across computing environments. Before running the analyses, users should obtain the required GWAS summary statistics and place them in the input locations specified in the corresponding scripts.

The R package environment can be restored using `renv.lock`.


## Citation

If you use this code, please cite the associated publication and this repository.
