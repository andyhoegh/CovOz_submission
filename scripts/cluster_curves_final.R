############################################################################################################
############################################################################################################
# R code for "Dynamics of Circulation at the Population Level" Section
############################################################################################################
############################################################################################################
library(tidyverse)
library(lubridate)
library(stringr)
library(rstan)
options(mc.cores = parallel::detectCores())

combined_out_variant <- read_csv('data/combined_out_variant.csv')

all_clusters <- c('beta 2d.i', "beta 2d.ii", 'beta 2d.iii', 'beta 2d.iv', 'beta 2d.v', 'beta 2d.vi')
out_fitted_curves <- NULL

####### Model for cluster combination #############

for (clst in all_clusters){
  print(clst)
  bats <- read_csv('data/combined_out_variant.csv') %>%
    mutate(positive = ifelse(variant_positive=="FALSE", 0, 1)) %>%
    filter(type == clst) %>% arrange(date_end)
  
  t_all <- substr(bats$date_end, 1, 10) %>% ymd() %>% as.integer()
  t <- unique(t_all)
  tnew <- seq(min(t-30), max(t + 30), by=10)
  
  N <- length(t)
  k <- table(t_all)
  max_k <- max(k)
  m <- matrix(rep(0, N * max(k)), nrow=N)
  
  m[1,] <- c(bats$n[1:k[1]], rep(0, max_k-k[1]))
  for (i in 2:N){
    m[i,1:k[i]] <- bats$n[(sum(k[1:(i-1)])+1):sum(k[1:i])]
  }
  
  y <- matrix(rep(0, N * max(k)), nrow=N)
  
  y[1,] <- c(bats$positive[1:k[1]], rep(0, max_k-k[1]))
  for (i in 2:N){
    y[i,1:k[i]] <- bats$positive[(sum(k[1:(i-1)])+1):sum(k[1:i])]
  }
  row_index <- matrix(rep(1:nrow(m), ncol(m)), nrow = nrow(m), ncol = ncol(m))
  col_index <- matrix(rep(1:ncol(m), nrow(m)), nrow = nrow(m), ncol = ncol(m), byrow = T)
  
  fit <- stan("scripts/GP_withLL.stan",
              data = list(N1 = N,
                          N2 = length(tnew),
                          N_tot = nrow(bats),
                          time_int = c(t, tnew),
                          max_k = max_k,
                          k = k,
                          m = m,
                          y = y,
                          row_index = row_index[m>0],
                          col_index = col_index[m>0],
                          ig_alpha = 18.51,
                          ig_beta = 1198.18),
              iter = 2000, chains=4, seed=1222021)
  z_vals3 <- extract(fit)$z
  out_curves <- tibble(median = pnorm(apply(z_vals3[1:4000,], 2, median)),
                       lower =  pnorm(apply(z_vals3[1:4000,], 2, quantile, prob = .025)),
                       upper =  pnorm(apply(z_vals3[1:4000,], 2, quantile, prob = .975)),
                       date = as_date(c(t,tnew)),
                       cluster = clst)
  out_fitted_curves <- out_fitted_curves %>% bind_rows(out_curves)
  
}

write_csv(out_fitted_curves, file ="data/model_output/cluster_curves.csv")
