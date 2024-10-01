############################################################################################################
############################################################################################################
# R code for "Individual Level Dynamics of Infection: Dynamic Binary Regression" Section
############################################################################################################
############################################################################################################

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(lubridate)
library(tictoc)

tic('logistic_curves_final')

bat_list <- read_csv('data/individual_variant_covariates.csv') %>%
  filter(!is.na(bat_age), bat_age != 'BLANK') |>
  group_by(bat_species, type) %>%
  tally()  %>% ungroup()


########################################################################
### Run dynamic binary regression for species - clade combinations
### Note: run time can take several hours
########################################################################

out_fitted_curves <- NULL
out_mu <- NULL
out_beta <- NULL

for (iter in 1:nrow(bat_list)){
  print(iter) 
  print(bat_list[iter,])
  bats <- read_csv('data/individual_variant_covariates.csv') %>%
    filter(bat_age %in% c('adult','juve','sub_adult'), type == bat_list[iter,'type'] %>% pull(), bat_species == bat_list[iter, 'bat_species'] %>% pull() )  %>%
    arrange(sampling_date) %>% mutate(n = 1, y = as.numeric(variant_positive)) %>%
    mutate(t = substr(sampling_date, 1, 10) %>% ymd() %>% as.integer())
  
  tnew <- seq(min(bats$t-30), max(bats$t + 30), by=10)
  
  tmp <- bats %>%
    group_by(sampling_date, t, bat_age) %>%
    summarize(n = n(), y = sum(variant_positive), .groups = 'drop')
  
  X <- as.matrix(model.matrix(y~bat_age -1 , data = tmp))
  tmp <- tmp %>% bind_cols(X)
  
  fit <- stan("scripts/GP_regression.stan",
              data = list(N1 = nrow(tmp),
                          N2 = length(tnew),
                          time_int = c(tmp$t, tnew),
                          n = tmp$n,
                          y = tmp$y,
                          x1 = tmp$bat_ageadult,
                          x2 = tmp$bat_agejuve,
                          x3 = tmp$bat_agesub_adult,
                          ig_alpha = 18.51,
                          ig_beta = 1198.18),
              iter = 2000, chains=4, seed=1222021)
  
  z_vals <- rstan::extract(fit)$z
  beta1_vals <-  rstan::extract(fit)$beta1
  beta2_vals <-  rstan::extract(fit)$beta2
  beta3_vals <-  rstan::extract(fit)$beta3
  beta1_2 <- beta1_vals - beta2_vals
  beta1_3 <- beta1_vals - beta3_vals
  beta2_3 <- beta2_vals - beta3_vals
  
  num_pts <- length(c(tmp$t,tnew))
  out_curves <- tibble(mean = c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean))),
                       lower =  c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025))),
                       upper =  c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975))),
                       date = rep(as_date(c(tmp$t,tnew)), 3),
                       type = rep(bat_list[iter,'type'] %>% pull(),num_pts *3) ,
                       age  = c(rep('adult',num_pts) ,
                                rep('juve',num_pts),
                                rep('sub_adult', num_pts)),
                       species = rep(bat_list[iter, 'bat_species']%>% pull(), num_pts*3))
  
  out_fitted_curves <- out_fitted_curves %>% bind_rows(out_curves)
  
  out_vals <- tibble(mean = c(mean(beta1_vals),
                              mean(beta2_vals),
                              mean(beta3_vals)),
                     lower = c(quantile(beta1_vals, probs = .025),
                               quantile(beta2_vals, probs = .025),
                               quantile(beta3_vals, probs = .025)),
                     upper = c(quantile(beta1_vals, probs = .975),
                               quantile(beta2_vals, probs = .975),
                               quantile(beta3_vals, probs = .975)),
                     type = rep(bat_list[iter,'type'] %>% pull(),3) ,
                     age  = c('adult','juve','sub_adult'),
                     species = rep(bat_list[iter, 'bat_species']%>% pull(), 3))
  out_mu <- out_mu %>% bind_rows(out_vals)
  
  beta_diff <- tibble(vals = c(beta1_2, beta1_3, beta2_3),
                        type = rep(c('beta1 - beta2', 'beta1 - beta3', 'beta2 - beta3'), each = 4000),
                                                clade = rep(bat_list[iter,'type'] |> pull(), 4000 * 3),
                      species = rep(bat_list[iter, 'bat_species']%>% pull(), 4000 * 3))
  out_beta <- out_beta %>% bind_rows(beta_diff)
}

out_fitted_curves <- out_fitted_curves |>
  mutate(key = paste(as.character(date), type, age, species)) |>
  filter(!duplicated(key)) |>
  select(-key)

save(out_mu, out_fitted_curves,  out_beta, file = 'data/model_output/logistic_curve_out.RData')

toc()
