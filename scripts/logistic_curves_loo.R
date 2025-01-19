library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(lubridate)
library(tictoc)
tic('start')
bat_list <- read_csv('data/individual_variant_covariates.csv') %>%
  filter(!is.na(bat_age), bat_age != 'BLANK') |>
  group_by(bat_species, type) %>%
  tally()  %>% ungroup()

####### Age Model


out_fitted_curves <- NULL
out_mu <- NULL
looic <- rep(0, nrow(bat_list))

for (iter in 1:nrow(bat_list)){
  print(iter)
  print(bat_list[iter,])
  bats <- read_csv('data/individual_variant_covariates.csv') %>%
    filter(bat_age %in% c('adult','juve','sub_adult'), 
           type == bat_list[iter,'type'] %>% pull(), 
           bat_species == bat_list[iter, 'bat_species'] %>% pull() )  %>%
    arrange(sampling_date) %>% 
    mutate(n = 1, y = as.numeric(variant_positive)) %>%
    mutate(t = substr(sampling_date, 1, 10) %>% ymd() %>% as.integer())
  
  # change the `by` argument for smoother visuals or shorter computation
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
  
  looic[iter] <- loo::loo(fit)$estimates['looic','Estimate']
  z_vals <- extract(fit)$z
  beta1_vals <- extract(fit)$beta1
  beta2_vals <- extract(fit)$beta2
  beta3_vals <- extract(fit)$beta3
  
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
  
  out_curves %>%
    filter(species == 'bff') %>%
    ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = mean )) +
    geom_ribbon(alpha = .3) + facet_grid(type~factor(age, levels = c('juve','sub_adult','adult'))) +
    geom_line(linetype = 3) +
    theme_bw() +
    xlab('') + ylab('Prevalence') +
    theme(legend.position = 'none') +
    # geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff') ) +
    #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
    xlim(min(out_curves$date), max(out_curves$date)) +
    ggtitle("Black Flying Fox Individual Data")
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
}

save(out_mu, out_fitted_curves, looic, file = 'data/model_output/logistic_curve_loo_age.RData')
load(file = 'data/model_output/logistic_curve_loo_age.RData')
sum(looic)
####### Age + Sex Model 

# bat_list_AP <- read_csv('Data/processed_data/individual_variant_covariates_AP.csv') %>%
#   mutate(bat_age = case_when(
#     bat_age == 'adult' ~ 'adult',
#     bat_age == 'juve' ~ 'juve',
#     bat_age == 'sub_adult' ~ 'sub_adult',
#     bat_age == 'dep_pup'~ 'juve'
#   )) |>
#   filter(!is.na(bat_age), bat_sex != 'BLANK') |>
#   group_by(bat_species, type, bat_age, bat_sex) %>%
#   tally()  %>% ungroup()
#  
out_fitted_curves <- NULL
out_mu <- NULL
looic <- rep(0, nrow(bat_list))

for (iter in 1:nrow(bat_list)){
  print(iter)
  print(bat_list[iter,])
  
  bats <- read_csv('data/individual_variant_covariates.csv') %>%
    filter(bat_age %in% c('adult','juve','sub_adult'), 
           type == bat_list[iter,'type'] %>% pull(), 
           bat_species == bat_list[iter, 'bat_species'] %>% pull() )  %>%
    arrange(sampling_date) %>% 
    mutate(n = 1, y = as.numeric(variant_positive)) %>%
    mutate(t = substr(sampling_date, 1, 10) %>% ymd() %>% as.integer())


  # change the `by` argument for smoother visuals or shorter computation
  tnew <- seq(min(bats$t-30), max(bats$t + 30), by=10)

  tmp <- bats %>%
    group_by(sampling_date, t, bat_age, bat_sex) %>%
    summarize(n = n(), y = sum(variant_positive), .groups = 'drop')

  X <- as.matrix(model.matrix(y~bat_sex + bat_age  , data = tmp))
  tmp <- tmp %>% bind_cols(X)

  fit <- stan("scripts/GP_regression_add.stan",
              data = list(N1 = nrow(tmp),
                          N2 = length(tnew),
                          time_int = c(tmp$t, tnew),
                          n = tmp$n,
                          y = tmp$y,
                          x1 = tmp$`(Intercept)`,
                          x2 = tmp$bat_sexmale,
                          x3 = tmp$bat_agejuve,
                          x4 = tmp$bat_agesub_adult,
                          ig_alpha = 18.51,
                          ig_beta = 1198.18),
              iter = 2000, chains=4, seed=1222021)
  looic[iter] <- loo::loo(fit)$estimates['looic','Estimate']
  #system('say done done')
  z_vals <- extract(fit)$z
  beta1_vals <- extract(fit)$beta1
  beta2_vals <- extract(fit)$beta2
  beta3_vals <- extract(fit)$beta3
  beta4_vals <- extract(fit)$beta4

  num_pts <- length(c(tmp$t,tnew))
  out_curves <- tibble(mean = c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) +  matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean))),
                       lower =  c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) +  matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025))),
                       upper =  c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) +  matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975))),
                       date = rep(as_date(c(tmp$t,tnew)), 6),
                       type = rep(bat_list[iter,'type'] %>% pull(),num_pts *6) ,
                       age  = c(rep('female adult', num_pts) ,
                                rep('female juve', num_pts),
                                rep('female sub_adult', num_pts),
                                rep('male adult', num_pts),
                                rep('male juve', num_pts),
                                rep('male sub_adult', num_pts)),
                       species = rep(bat_list[iter, 'bat_species']%>% pull(), num_pts*6))


  out_curves %>%
    filter(species == 'bff') %>%
    ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = mean )) +
    geom_ribbon(alpha = .3) + facet_grid(type~factor(age, levels = c('female adult','male adult','female sub_adult',
                                                                     'male sub_adult', 'female juve','male juve'))) +
    geom_line(linetype = 3) +
    theme_bw() +
    xlab('') + ylab('Prevalence') +
    theme(legend.position = 'none') +
   # geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff') ) +
  #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
    xlim(min(out_curves$date), max(out_curves$date)) +
    ggtitle("Black Flying Fox Individual Data")
  
  out_fitted_curves <- out_fitted_curves %>% bind_rows(out_curves)

  out_vals <- tibble(mean = c(mean(beta1_vals),
                  mean(beta2_vals),
                  mean(beta3_vals),
                  mean(beta4_vals)),
                lower = c(quantile(beta1_vals, probs = .025),
                          quantile(beta2_vals, probs = .025),
                          quantile(beta3_vals, probs = .025),
                          quantile(beta4_vals, probs = .025)),
                upper = c(quantile(beta1_vals, probs = .975),
                          quantile(beta2_vals, probs = .975),
                          quantile(beta3_vals, probs = .975),
                          quantile(beta4_vals, probs = .975)),
                type = rep(bat_list[iter,'type'] %>% pull(),4) ,
                age  = c('(Intercept)','bat_sexmale','bat_agejuve', 'bat_agesub_adult'),
                species = rep(bat_list[iter, 'bat_species']%>% pull(), 4))
  
  out_mu <- out_mu %>% bind_rows(out_vals)
}

save(out_mu, out_fitted_curves, looic, file = 'data/model_output/logistic_curve_loo_age_add_sex.RData')
load(file = 'data/model_output/logistic_curve_loo_age_add_sex.RData')
sum(looic) # 654
####### Age x Sex Model 


out_fitted_curves <- NULL
out_mu <- NULL
looic <- rep(0, nrow(bat_list))

for (iter in 1:nrow(bat_list)){
  print(iter)
  print(bat_list[iter,])
  bats <- read_csv('data/individual_variant_covariates.csv') %>%
    filter(bat_age %in% c('adult','juve','sub_adult'), 
           type == bat_list[iter,'type'] %>% pull(), 
           bat_species == bat_list[iter, 'bat_species'] %>% pull() )  %>%
    arrange(sampling_date) %>% 
    mutate(n = 1, y = as.numeric(variant_positive)) %>%
    mutate(t = substr(sampling_date, 1, 10) %>% ymd() %>% as.integer())
  
  
  # change the `by` argument for smoother visuals or shorter computation
  tnew <- seq(min(bats$t-30), max(bats$t + 30), by=10)
  
  tmp <- bats %>%
    group_by(sampling_date, t, bat_age, bat_sex) %>%
    summarize(n = n(), y = sum(variant_positive), .groups = 'drop')
  
  X <- as.matrix(model.matrix(y~bat_sex * bat_age  , data = tmp))
  tmp <- tmp %>% bind_cols(X)
  
  fit <- stan("scripts/GP_regression_interact.stan",
              data = list(N1 = nrow(tmp),
                          N2 = length(tnew),
                          time_int = c(tmp$t, tnew),
                          n = tmp$n,
                          y = tmp$y,
                          x1 = tmp$`(Intercept)`,
                          x2 = tmp$bat_sexmale,
                          x3 = tmp$bat_agejuve,
                          x4 = tmp$bat_agesub_adult,
                          x5 = tmp$`bat_sexmale:bat_agejuve`,
                          x6 = tmp$`bat_sexmale:bat_agesub_adult`,
                          ig_alpha = 18.51,
                          ig_beta = 1198.18),
              iter = 2000, chains=4, seed=1222021)
  
  looic[iter] <- loo::loo(fit)$estimates['looic','Estimate']
  z_vals <- extract(fit)$z
  beta1_vals <- extract(fit)$beta1
  beta2_vals <- extract(fit)$beta2
  beta3_vals <- extract(fit)$beta3
  beta4_vals <- extract(fit)$beta4
  beta5_vals <- extract(fit)$beta5
  beta6_vals <- extract(fit)$beta6
  
  num_pts <- length(c(tmp$t,tnew))
  out_curves <- tibble(mean = c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean)),
                                pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) +  matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, mean))),
                       lower =  c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) +  matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .025))),
                       upper =  c(pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta3_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975)),
                                  pnorm(apply(z_vals + matrix(beta1_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) + matrix(beta2_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)) +  matrix(beta4_vals, nrow = nrow(z_vals), ncol = ncol(z_vals)), 2, quantile, prob = .975))),
                       date = rep(as_date(c(tmp$t,tnew)), 6),
                       type = rep(bat_list[iter,'type'] %>% pull(),num_pts *6) ,
                       age  = c(rep('female adult', num_pts) ,
                                rep('female juve', num_pts),
                                rep('female sub_adult', num_pts),
                                rep('male adult', num_pts),
                                rep('male juve', num_pts),
                                rep('male sub_adult', num_pts)),
                       species = rep(bat_list[iter, 'bat_species']%>% pull(), num_pts*6))
  
  
  out_curves %>%
    filter(species == 'bff') %>%
    ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = mean )) +
    geom_ribbon(alpha = .3) + facet_grid(type~factor(age, levels = c('female adult','male adult','female sub_adult',
                                                                     'male sub_adult', 'female juve','male juve'))) +
    geom_line(linetype = 3) +
    theme_bw() +
    xlab('') + ylab('Prevalence') +
    theme(legend.position = 'none') +
    # geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff') ) +
    #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
    xlim(min(out_curves$date), max(out_curves$date)) +
    ggtitle("Black Flying Fox Individual Data")
  
  out_fitted_curves <- out_fitted_curves %>% bind_rows(out_curves)
  
  out_vals <- tibble(mean = c(mean(beta1_vals),
                              mean(beta2_vals),
                              mean(beta3_vals),
                              mean(beta4_vals)),
                     lower = c(quantile(beta1_vals, probs = .025),
                               quantile(beta2_vals, probs = .025),
                               quantile(beta3_vals, probs = .025),
                               quantile(beta4_vals, probs = .025)),
                     upper = c(quantile(beta1_vals, probs = .975),
                               quantile(beta2_vals, probs = .975),
                               quantile(beta3_vals, probs = .975),
                               quantile(beta4_vals, probs = .975)),
                     type = rep(bat_list[iter,'type'] %>% pull(),4) ,
                     age  = c('(Intercept)','bat_sexmale','bat_agejuve', 'bat_agesub_adult'),
                     species = rep(bat_list[iter, 'bat_species']%>% pull(), 4))
  
  out_mu <- out_mu %>% bind_rows(out_vals)
}

save(out_mu, out_fitted_curves, looic, file = 'data/model_output/logistic_curve_loo_age_interact_sex.RData')
load( file = 'data/model_output/logistic_curve_loo_age_interact_sex.RData')
sum(looic) # 643
toc() #6:45
# ######
# 
# load('Data/processed_data/logistic_curve_out_AP.RData')
# bats <- read_csv('Data/processed_data/individual_variant_covariates.csv')
# bats %>% filter(bat_age %in% c('adult','juve','sub_adult') ) %>%
#   group_by(type, bat_age, bat_species) %>% 
#   summarize(n = n(), n_pos = sum(variant_positive), mean_pos = mean(variant_positive))
# 
# bat_prop  <- bats %>% filter(bat_age %in% c('adult','juve','sub_adult') ) %>%
#   group_by(type, bat_age, bat_species, date) %>% 
#   summarize(n = n(), n_pos = sum(variant_positive), mean_pos = mean(variant_positive),
#             .groups = 'drop')
# 
# 
# out_fitted_curves %>%
#   filter(species == 'bff') %>%
#   ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = mean )) + 
#   geom_ribbon(alpha = .3) +
#   geom_point(inherit.aes = F, aes(x = date, y = mean_pos, size = n), data = bat_prop %>% filter(bat_species == 'bff') %>% mutate(age = bat_age), alpha = .1) +
#   facet_grid(type~factor(age, levels = c('juve','sub_adult','adult'))) + 
#   geom_line(linetype = 3) +
#   theme_bw() + 
#   xlab('') + ylab('Prevalence') +
#  # theme(legend.position = 'none') +
#   # geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff') ) +
#   #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
#   xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
#   ylim(0, 1) +
#   ggtitle("Black Flying Fox Individual Data") 
# 
# 
# bats <- read_csv('Data/processed_data/individual_variant_covariates.csv') %>%
#   filter(bat_species == 'bff', type %in% c('beta 2d.iv', 'beta 2d.v'), variant_positive)
# 
# dup_keys <- bats %>% group_by(key) %>% tally() %>% filter(n == 2)
# bat_co <- bats %>% filter(key %in% dup_keys$key) %>% arrange(sampling_date) %>%
#   mutate(age = factor(bat_age, levels = c('juve','sub_adult','adult')))
# 
# out_fitted_curves %>%
#   filter(species == 'bff', type %in% c('beta 2d.v', 'beta 2d.iv')) %>%
#   ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = mean )) + 
#   geom_ribbon(alpha = .3) +
#   geom_point(inherit.aes = F, aes(x = date, y = mean_pos, size = n), 
#              data = bat_prop %>% filter(bat_species == 'bff', type %in% c('beta 2d.v', 'beta 2d.iv')) %>% mutate(age = bat_age),
#              alpha = .1) +
#   facet_grid(type~factor(age, levels = c('juve','sub_adult','adult'))) + 
#   geom_line(linetype = 3) +
#   theme_bw() + 
#   xlab('') + ylab('Prevalence') +
#   # theme(legend.position = 'none') +
#   geom_rug(inherit.aes = F, aes(x = date), data = bat_co, color = 'red' ) +
#   #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
#   xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
#   ylim(0, 1) +
#   ggtitle("Black Flying Fox Individual Data") 
# 
# 
# 
# out_fitted_curves %>%
#   filter(species == 'ghff') %>%
#   ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = mean )) + 
#   geom_ribbon(alpha = .3) +
#   geom_point(inherit.aes = F, aes(x = date, y = mean_pos, size = n), data = bat_prop %>% filter(bat_species == 'ghff') %>% mutate(age = bat_age), alpha = .1) +
#   facet_grid(type~factor(age, levels = c('juve','sub_adult','adult'))) + 
#   geom_line(linetype = 3) +
#   theme_bw() + 
#   xlab('') + ylab('Prevalence') +
#   #theme(legend.position = 'none') +
#   # geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff') ) +
#   #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
#   xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
#   ylim(0, 1) +
#   ggtitle("Grey Headed Flying Fox Individual Data") 
# 
# 
# out_mu %>% mutate(age = factor(age, levels = c('juve','sub_adult','adult'))) %>%
#   ggplot(aes(x = lower, xend = upper, y = type, yend = type)) +
#   geom_segment() + facet_grid(species~age) + theme_bw() +
#   geom_point(aes(x = mean)) + theme(legend.position = 'none') +
#   ylab("Clade") + xlab(expression(beta)) + 
#   ggtitle('Dynamic Binary Regression Coefficients') +
#   labs(caption = '95% credibility intervals')
# 
# # out_curves %>%
# #   filter(species == 'bff') %>%
# #   ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, y = median )) +
# #   geom_ribbon(alpha = .3) + facet_grid(type~factor(age, levels = c('juve','sub_adult','adult'))) +
# #   geom_line(linetype = 3) +
# #   theme_bw() +
# #   xlab('') + ylab('Prevalence') +
# #   theme(legend.position = 'none') +
# #   # geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff') ) +
# #   #  geom_rug(inherit.aes = F, aes(x = date), data = bat_samples %>% filter(bat_species == 'bff', variant_positive == 1), color = 'red', sides = 'b') +
# #   xlim(min(out_curves$date), max(out_curves$date)) +
# #   ggtitle("Black Flying Fox Individual Data")
