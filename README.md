### Reproducible Statistical Analysis README

Instructions for running code associated with "Synchronized seasonal excretion of multiple coronaviruses in Australian Pteropus spp is associated with co-infections in juvenile and sub-adult bats." Figures from the manuscript can be recreated using the "CovOZ_Figures_Submission_Clean.R" script.


### Coinfection Analysis

Runs chi-squared tests on coinfections of beta 2d.iv and beta 2d.v. Generates summary statistics, test statistics, and p-values from manuscript.

- input files: individual_variant_covariates.csv
- script file: coinfection_final.R

### Individual Level Dynamics of Infection: Dynamic Binary Regression

Runs individual level dynamic binary regression models. Produces output file that can recreate figures.

- input files: individual_variant_covariates.csv
- script files:
  - logistic_curves_final.R
  - GP_regression.stan
- output files: logistic_curve_out.RData
  
  
### Dynamics of Circulation at the Population Level

Runs combined (individual and pooled data) dynamic models. Produces output file that can recreate figures.

- input files: combined_out_variant.csv
- script files:
  - cluster_curves_final.R
  - GP_withLL.stan
- output files: cluster_curves.csv
