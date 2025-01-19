library(tidyverse)
library(rstan)
library(knitr)
library(lubridate)
library(stringr)
tic('LOO')
options(mc.cores = 4)


site_list <- c('Redcliffe', 'Sunnybank','Burleigh','Clunes','Toowoomba')
all_clusters <- c('beta 2d.i', "beta 2d.ii", 'beta 2d.iii',
                  'beta 2d.iv', 'beta 2d.v', 'beta 2d.vi')
out_loo_clust_site <- NULL
out_loo_clust <- NULL
out_loo_site <- NULL


####### Model for site-cluster combination #############

for (site_loc in site_list){
  for (clst in all_clusters){
    print(clst)
    print(site_loc)  
    bats <- read_csv('data/combined_out_variant.csv') %>%
      mutate(positive = ifelse(variant_positive=="FALSE", 0, 1)) %>%
      filter(type == clst, site == site_loc) %>% arrange(date_end)
    
    t_all <- substr(bats$date_end, 1, 10) %>% ymd() %>% as.integer()
    t <- unique(t_all)
    
    # change the `by` argument for smoother visuals or shorter computation
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
    
  
    loo_val <- tibble(site = site_loc,
                      cluster = clst,
                      loo = loo(fit)$looic)
    out_loo_clust_site <- out_loo_clust_site %>% bind_rows(loo_val)
  }
}

####### Model for site #############

for (site_loc in site_list){
  print(site_loc)  
  bats <- read_csv('data/combined_out_variant.csv') %>%
    mutate(positive = ifelse(variant_positive=="FALSE", 0, 1)) %>%
    filter(site == site_loc) %>% arrange(date_end)
    
  t_all <- substr(bats$date_end, 1, 10) %>% ymd() %>% as.integer()
  t <- unique(t_all)
    
  # change the `by` argument for smoother visuals or shorter computation
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
    
  loo_val <- tibble(site = site_loc,
                    loo = loo(fit)$looic)
  out_loo_site <- out_loo_site %>% bind_rows(loo_val)
}


####### Model for cluster combination #############

for (clst in all_clusters){
  print(clst)
  bats <- read_csv('data/combined_out_variant.csv') %>%
    mutate(positive = ifelse(variant_positive=="FALSE", 0, 1)) %>%
    filter(type == clst) %>% arrange(date_end)
    
  t_all <- substr(bats$date_end, 1, 10) %>% ymd() %>% as.integer()
  t <- unique(t_all)
  
  # change the `by` argument for smoother visuals or shorter computation
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
    
  loo_val <- tibble(cluster = clst,
                    loo = loo(fit)$looic)
  out_loo_clust <- out_loo_clust %>% bind_rows(loo_val)
}

####### A Single Combined Model #############
bats <- read_csv('data/combined_out_variant.csv') %>%
      mutate(positive = ifelse(variant_positive=="FALSE", 0, 1)) %>% arrange(date_end)
    
t_all <- substr(bats$date_end, 1, 10) %>% ymd() %>% as.integer()
t <- unique(t_all)
    
# change the `by` argument for smoother visuals or shorter computation
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
    
out_loo_all <- loo(fit)$looic

out_loo_all / -2
out_loo_clust %>% summarize(sum(loo)) / -2
out_loo_site %>% summarize(sum(loo)) / -2
out_loo_clust_site %>% summarize(sum(loo)) / -2

save(out_loo_all, out_loo_clust, out_loo_site, out_loo_clust_site, file = "data/model_output/preds.RData")
toc()