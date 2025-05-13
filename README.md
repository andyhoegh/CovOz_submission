## Synchronized seasonal excretion of multiple coronaviruses in Australian Pteropus spp is associated with co-infections in juvenile and sub-adult bats.

### General Information

This repo contains instructions and source code for reproducing the statistical analyses in the manuscript.
Pteropus spp is associated with co-infections in juvenile and subadult bats"

Author/Principal Investigator Information
Name: Plowright, Raina K.
ORCID: https://orcid.org/0000-0002-3338-6590
Institution: Cornell University
Email: rkp57@cornell.edu

Author/Corresponding Authors
Name: Peel, Alison J.
ORCID: https://orcid.org/0000-0003-3538-3550
Institution: University of Sydney
Email: alison.peel@sydney.edu.au

Name: Hoegh, Andrew
ORCID: https://orcid.org/0000-0003-1176-4965
Institution: Montana State University
Email: andrew.hoegh@montana.edu

Information about funding sources that supported the collection of the data: 
The project was supported by the National Science Foundation (DEB1716698, EF2133763, EF-2231624), and the
DARPA PREEMPT program Cooperative Agreement # D18AC00031. The content of the information does not
necessarily reflect the position or the policy of the U.S. government, and no official endorsement should be inferred.
AJP was supported by an ARC DECRA fellowship (DE190100710) and a University of Sydney Horizons Fellowship.

### Sharing / Access Information

Licenses/restrictions placed on the data: 
This dataset is shared under a Creative Commons 1.0 Universal Public Domain Dedication (https://creativecommons.org/publicdomain/zero/1.0/). The material can be copied, modified and used without permission, but attribution to the original authors is always appreciated. 

Recommended citation for this dataset: 
Alison J. Peel, Manuel Ruiz-Aravena, Karan Kim, Braden Scherting, Caylee A. Falvo, Daniel E. Crowley, Vincent J. Munster, Edward J. Annand, Karren Plain, Devin N. Jones, Tamika J. Lunn, Adrienne S. Dale, Andrew Hoegh, John-Sebastian Eden, Raina K. Plowright. (2024) Data from: Synchronized seasonal excretion of multiple coronaviruses in Australian Pteropus spp is associated with co-infections in juvenile and subadult bats. [dataset] Cornell University Library eCommons Repository. https://doi.org/10.7298/w7sw-6161

### Repo Contents

- scripts: contains the source .R and .stan files to reproduce the anaysis. Each file is detailed below in the specific sections corresponding to the statistical analyses.
- data: contains the raw source data and model generated output.
- figures: contains the final output figures from the manuscript. These can be recreated with the CovOZ_Figures_Submission_Clean.R script.

### 1. System Requirements

#### Hardware Requirements

Our source code requires only a standard computer. Much of the Markov chain Monte Carlo code is run in parallel so a computer with ample memory and multiple cores can be advantageous. The runtimes below are generated using a macbook with the recommended specs (64 GB RAM, 8 cores at 2.7 GHz). The code will also work on linux or windows computer.

#### Software Requirements

Reproducing the statistical analyses requires a current version of R and stan. We use version 4.4.1 of R and version 2.32.2 of stan.

#### Package dependencies and versions

Users will need the following packages install the following packages to execute the code. Our versions are effective October 1, 2024

```
tidyverse 2.0.0
lubridate 1.9.3
stringr 1.5.1
rstan 2.32.6
cowplot 1.1.3
ggtext 0.1.2
jpeg 0.1-10
scales 1.3.0
tictoc 1.2.1
```


### 2. Installation Guide

Running the analysis requires:

- installing R. Depending on wifi speeds, installing R usually takes a few minutes.
- installing stan. Depending on wifi speeds, installing stan usually takes a few minutes.
- installing the necessary R packages (listed above). Depending on wifi speeds, installing packages usually takes about 30 seconds per package.




### 3. Demo

This source code is not an R package with a formal demo, but rather source code is included for the various analyses in section 4.


### 4. Instructions for Use

#### 4.1 Coinfection Analysis

Runs chi-squared tests on coinfections of beta 2d.iv and beta 2d.v. Generates summary statistics, test statistics, and p-values from manuscript.

- input files: individual_variant_covariates.csv
- script file: coinfection_final.R
- run time: approximately 1 second

#### 4.2 Individual Level Dynamics of Infection: Dynamic Binary Regression

Runs individual level dynamic binary regression models. Produces output file that can recreate figures.

- input files: individual_variant_covariates.csv
- script files:
  - logistic_curves_final.R
  - GP_regression.stan
- output files: logistic_curve_out.RData
- run time: approximately 66 minutes
  
  
#### 4.3 Dynamics of Circulation at the Population Level

Runs combined (individual and pooled data) dynamic models. Produces output file that can recreate figures.

- input files: combined_out_variant.csv
- script files:
  - cluster_curves_final.R
  - GP_withLL.stan
- output files: cluster_curves.csv
- run time: approximately 25 minutes

#### 4.4 Manuscript Figures

Combined script that uses output files created by previous scripts to recreate all figures in the manuscript.


- input files: 
  - model_output/cluster_curves.csv
  - combined_out_variant.csv
  - individual_variant_covariates.csv
  - model_output/logistic_curve_out.RData

- script files:
  - CovOZ_Figures_Submission_Clean.R

- output files: 
  - Figure2_final.png
  - Figure3_final.png
  - Figure4_A-D_final.png
  - Figure6_AP.png
  - Figure7.png
  - SIFigure8.png
  - SIFigure9.png

- run time: approximately 16 seconds

#### 4.5 Model Comparison Integrated

Compares LOOIC values for sets of model frameworks.

- input files: combined_out_variant.csv
- script files:
  - Pred_Comparisons.R
  - GP_withLL.stan
- output files: preds.RData
- run time: approximately 2 hours

#### 4.6 Model Comparison Individual

Compares LOOIC values for sets of model frameworks.

- input files: combined_out_variant.csv
- script files:
  - logistic_curves_loo.R
  - GP_regression.stan
  - GP_regression_add.stan
  - GP_regression_interact.stan
- output files: 
  - logistic_curve_loo_age.RData
  - logistic_curve_loo_age_add_sex.RData
  - logistic_curve_loo_age_interact_sex.RData
- run time: approximately 6:45 hours
