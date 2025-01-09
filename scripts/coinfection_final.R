############################################################################################################
############################################################################################################
# R code for "Recaptured and co-infected individuals" Section
############################################################################################################
############################################################################################################

library(tidyverse)
library(lubridate)
library(tictoc)
tic('coinfection')

individual_variant_covariates <- read_csv("data/individual_variant_covariates.csv")

####################################
### Calculate Summary Statistics
####################################
individual_variants_wide <- individual_variant_covariates |> 
  pivot_wider(names_from = type, values_from = variant_positive)

bff_wide <- individual_variants_wide |>
  filter(bat_species == 'bff') 

bff_coinfections <- bff_wide |>
  mutate(infections = beta_2d_i + beta_2d_ii + beta_2d_iii + beta_2d_iv + beta_2d_v + beta_2d_vi) |>
  filter(infections > 1) |>
  tally()

ghff_coinfections <- individual_variants_wide |>
  filter(bat_species == 'ghff') |>
  mutate(infections = beta_2d_i + beta_2d_ii + beta_2d_iii + beta_2d_iv + beta_2d_v + beta_2d_vi) |>
  filter(infections > 1) |>
  tally()

bff_beta2d_iv <- bff_wide |>
  filter(beta_2d_iv) |>
  tally()
  
bff_beta2d_v <- bff_wide |>
  filter(beta_2d_v) |>
  tally()

bff_beta2d_iv_v_coinfection <- bff_wide |>
  filter(beta_2d_iv & beta_2d_v) |>
  tally()

####################################
### 2- way Chi-Squared Test
####################################

table(bff_wide$beta_2d_iv, bff_wide$beta_2d_v)
Xsq <- chisq.test(bff_wide$beta_2d_iv, bff_wide$beta_2d_v, simulate.p.value = TRUE)

####################################
### 3- way Chi-Squared Test with age
####################################

bats_iv <- individual_variant_covariates %>% 
  filter(type == 'beta_2d_iv') %>% 
  rename(positive_iv = variant_positive) %>% 
  mutate(key = paste(accession_update, '-', sample_id, sep = '')) |>
           select(key, positive_iv, bat_age) %>%
  filter(bat_age %in% c('adult','juve','sub_adult'))

bats_v <- individual_variant_covariates %>% 
  filter(type == 'beta_2d_v') %>% 
  rename(positive_v = variant_positive) %>% 
  mutate(key = paste(accession_update, '-', sample_id, sep = '')) |>
  select(key, positive_v, bat_age) %>%
  filter(bat_age %in% c('adult','juve','sub_adult'))

comb <- bats_iv %>% left_join(bats_v, by = 'key')
comb_adult <- comb %>% filter(bat_age.x == 'adult') 
comb_subadult <- comb %>% filter(bat_age.x == 'sub_adult') 
comb_juve <- comb %>% filter(bat_age.x == 'juve') 

iv_rate <- bats_iv %>% group_by(bat_age) %>% summarise(rate = mean(positive_iv)) 
iv_adult <- iv_rate %>% filter(bat_age == 'adult') %>% select(rate) %>% pull()
iv_juve <- iv_rate %>% filter(bat_age == 'juve') %>% select(rate) %>% pull()
iv_subadult <- iv_rate %>% filter(bat_age == 'sub_adult') %>% select(rate) %>% pull()

v_rate <- bats_v %>% group_by(bat_age) %>% summarise(rate = mean(positive_v)) 
v_adult <- v_rate %>% filter(bat_age == 'adult') %>% select(rate) %>% pull()
v_juve <- v_rate %>% filter(bat_age == 'juve') %>% select(rate) %>% pull()
v_subadult <- v_rate %>% filter(bat_age == 'sub_adult') %>% select(rate) %>% pull()

n_adult <- nrow(bats_iv %>% filter(bat_age == 'adult'))
n_subadult <- nrow(bats_iv %>% filter(bat_age == 'sub_adult'))
n_juve <- nrow(bats_iv %>% filter(bat_age == 'juve'))

## Adult
expected_matrix_adult <- n_adult * matrix(c((1- iv_adult) * (1- v_adult), 
                                iv_adult * (1- v_adult),
                                (1-iv_adult) * v_adult,
                                iv_adult * v_adult), 2, 2)
observed_matrix_adult <- matrix(c(comb_adult %>% filter(!positive_iv, !positive_v) %>% tally() %>% pull(),
                            comb_adult %>% filter(positive_iv, !positive_v) %>% tally() %>% pull(),
                            comb_adult %>% filter(!positive_iv, positive_v) %>% tally() %>% pull(),
                            comb_adult %>% filter(positive_iv, positive_v) %>% tally() %>% pull()),
                          2,2)


chisq_adult <- sum(((observed_matrix_adult - expected_matrix_adult) / sqrt(expected_matrix_adult))^2)
1- pchisq(chisq_adult,1)


## Subadult

expected_matrix_subadult <- n_subadult * matrix(c((1- iv_subadult) * (1- v_subadult), 
                                            iv_subadult * (1- v_subadult),
                                            (1-iv_subadult) * v_subadult,
                                            iv_subadult * v_subadult), 2, 2)
observed_matrix_subadult <- matrix(c(comb_subadult %>% filter(!positive_iv, !positive_v) %>% tally() %>% pull(),
                                  comb_subadult %>% filter(positive_iv, !positive_v) %>% tally() %>% pull(),
                                  comb_subadult %>% filter(!positive_iv, positive_v) %>% tally() %>% pull(),
                                  comb_subadult %>% filter(positive_iv, positive_v) %>% tally() %>% pull()),
                                2,2)


chisq_subadult <- sum(((observed_matrix_subadult - expected_matrix_subadult) / sqrt(expected_matrix_subadult))^2)
1- pchisq(chisq_subadult,1)

## Juveniles

expected_matrix_juve <- n_juve * matrix(c((1- iv_juve) * (1- v_juve), 
                                                  iv_juve * (1- v_juve),
                                                  (1-iv_juve) * v_juve,
                                                  iv_juve * v_juve), 2, 2)
observed_matrix_juve <- matrix(c(comb_juve %>% filter(!positive_iv, !positive_v) %>% tally() %>% pull(),
                                     comb_juve %>% filter(positive_iv, !positive_v) %>% tally() %>% pull(),
                                     comb_juve %>% filter(!positive_iv, positive_v) %>% tally() %>% pull(),
                                     comb_juve %>% filter(positive_iv, positive_v) %>% tally() %>% pull()),
                                   2,2)


chisq_juve <- sum(((observed_matrix_juve - expected_matrix_juve) / sqrt(expected_matrix_juve))^2)
1- pchisq(chisq_juve,1)
toc()