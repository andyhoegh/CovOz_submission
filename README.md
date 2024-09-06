### Reproducible Statistical Analysis README

Instructions for running code associated with "Synchronized seasonal excretion of multiple coronaviruses in Australian Pteropus spp is associated with co-infections in juvenile and sub-adult bats"


### Coinfection Analysis

- input files: individual_variant_covariates.csv
  - needed: species + age + infection status
- script file: coinfection_final.R

### Individual Level Dynamics of Infection: Dynamic Binary Regression

- input files: individual_variant_covariates.csv
  - needed: species + sampling_date + infection status
- script files:
  - logistic_curves_final.R
  - GP_regression.stan
  
  
### Dynamics of Circulation at the Population Level