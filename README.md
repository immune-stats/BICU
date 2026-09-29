# FeLV and FIV clinical decision support models

Data and R code accompanying *Sequential Reassessment of FeLV and FIV Clinical Decision Support Models Using an Extended Veterinary Isolation Unit Dataset*.

## Repository contents

- `dataset.csv` — de-identified admission-level data used in the analyses.
- `data_dictionary.csv` — variable definitions, permitted categories and missing-value meanings.
- `analysis.R` — single script reproducing the numerical tables and Figure 1.
- `results/` — numerical CSV outputs, Figure 1, recorded warnings and R session information.

## Reproduce the analysis

Use R with the `pROC` and `ResourceSelection` packages. Set this folder as the working directory and run:

```r
install.packages(c("pROC", "ResourceSelection"))  # Only needed once
source("analysis.R")
```

The equivalent terminal command, run from this folder, is `Rscript analysis.R`. The script reads `dataset.csv` and writes its outputs to `results/`, replacing files with the same names. The recorded analysis used R 4.5.1, pROC 1.19.0.1 and ResourceSelection 0.3.6; see `results/sessionInfo.txt` and `results/package_versions.csv`. The CSVs contain numerical results at full precision; the article tables apply publication formatting.

## Understand the data

The dataset contains 1,304 records and 14 variables. One row represents one admission record. There is no identifier linking records from the same animal, so rows with identical categories must not be deduplicated. A blank `felv` or `fiv` value means that the record is ineligible for that outcome analysis; it does **not** mean a negative test. Retrovirus-negative non-cases may contribute to both analyses, and records positive for both retroviruses may contribute to both positive groups. The data dictionary specifies the remaining missing values and category codes.

The development period is 2013–2022 and the later testing period is 2023–2025. Development samples require complete candidate variables. Later FeLV testing requires complete original FeLV predictors; later FIV testing also requires complete sex and neuter information. The combined models use complete cases for their own predictors. This gives FeLV development n = 630 (126 positive), testing n = 336 (58 positive), and combined n = 966 before additional model requirements; FIV development n = 638 (134 positive), testing n = 344 (67 positive), and combined n = 982. The updated FeLV model uses n = 965.

The original model coefficients and decision thresholds are specified in `analysis.R` and applied unchanged to later admissions. The combined-period models are fitted separately. The FIV comparison using joint sex and neuter categories is provided in supplementary Table S1.

## Find the article results

| Article item | File in `results/` |
| --- | --- |
| Tables 1 and 2 — population comparisons | `Table_1_Population.csv`, `Table_2_Population.csv` |
| Table 3 — fixed models across periods | `Table_3_Fixed_models.csv` |
| Tables 4 and 5 — regression results | `Table_4_Regression.csv`, `Table_5_Regression.csv` |
| Table S1 — FIV model comparison | `Table_S1_Coefficients.csv`, `Table_S1_Performance.csv` |
| Table S2 — original FeLV model in combined data | `Table_S2_Coefficients.csv` |
| Table S3 — candidate screening | `Table_S3_Global_screening.csv`, `Table_S3_Univariable_coefficients.csv` |
| Table S4 — updated model coefficients | `Table_S4_Coefficients.csv` |
| Figure 1 — ROC curves | `Figure_1_ROC.png`, `Figure_1_ROC.pdf` |

`Cohort_counts.csv`, `Variable_availability.csv` and `Combined_model_performance.csv` provide further numerical details. `Model_and_test_warnings.txt` records statistical warnings. AUC intervals use DeLong's method; coefficient intervals and fixed-threshold sensitivity, specificity and balanced-accuracy intervals use Wald estimates. Population comparisons use Pearson chi-squared tests without continuity correction. Some comparisons have expected cell counts below five, so their asymptotic p-values require caution.
