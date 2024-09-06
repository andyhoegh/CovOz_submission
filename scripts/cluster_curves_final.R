############################################################################################################
############################################################################################################
# R code for "Dynamics of Circulation at the Population Level" Section
############################################################################################################
############################################################################################################

library(tidyverse)
library(lubridate)
library(stringr)

options(mc.cores = 4)

all_clusters <- c('beta 2d.i', "beta 2d.ii", 'beta 2d.iii', 'beta 2d.iv', 'beta 2d.v', 'beta 2d.vi')
out_fitted_curves <- NULL

####### Model for cluster combination #############

for (clst in all_clusters){
  print(clst)
  bats <- read_csv('Data/processed_data/combined_results_variant.csv') %>%
    mutate(positive = ifelse(positive=="FALSE", 0, 1)) %>%
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
  
  fit <- stan("Scripts/stan/GP_withLL.stan",
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


write_csv(out_fitted_curves, file ="Data/results/cluster_curves.csv")
out_fitted_curves <- read_csv('Data/results/cluster_curves.csv')
#
bats <- read_csv('Data/processed_data/combined_results_variant.csv') %>%
  mutate(positive = ifelse(positive=="FALSE", 0, 1))
#
png("ClusterCurves5_v2.png", width = 16, height = 9, units = 'in', res = 1200)
out_fitted_curves %>%
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = cluster, y = median )) +
  geom_ribbon(alpha = .3) + facet_grid(rows = vars(cluster)) +
  geom_line(linetype = 3) +
  theme_bw() +
  xlab('') + ylab('Prevalence') +
  theme(legend.position = 'none')+
#dev.off() +
  geom_rug(inherit.aes = F, aes(x = date), data = bats %>% mutate(date = as_date(date_end))) +
  geom_rug(inherit.aes = F, aes(x = date), data = bats %>% filter(positive == 1) %>% mutate(date = as_date(date_end), cluster = type), color = 'red', sides = 'b') +
  xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
  scale_y_continuous(breaks=c(0,.10, .20)) +
  annotate("rect", xmin = as_date('2018-06-21'), xmax = as_date('2018-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = as_date('2019-06-21'), xmax = as_date('2019-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = as_date('2020-06-21'), xmax = as_date('2020-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = as_date('2017-07-15'), xmax = as_date('2017-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red")
dev.off()

pdf("ClusterCurves5_v2.pdf", width = 16, height = 9)
out_fitted_curves %>%
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = cluster, y = median )) +
  geom_ribbon(alpha = .3) + facet_grid(rows = vars(cluster)) +
  geom_line(linetype = 3) +
  theme_bw() +
  xlab('') + ylab('Prevalence') +
  theme(legend.position = 'none')+
  #dev.off() +
  geom_rug(inherit.aes = F, aes(x = date), data = bats %>% mutate(date = as_date(date_end))) +
  geom_rug(inherit.aes = F, aes(x = date), data = bats %>% filter(positive == 1) %>% mutate(date = as_date(date_end), cluster = type), color = 'red', sides = 'b') +
  xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
  scale_y_continuous(breaks=c(0,.10, .20)) +
  annotate("rect", xmin = as_date('2018-06-21'), xmax = as_date('2018-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = as_date('2019-06-21'), xmax = as_date('2019-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = as_date('2020-06-21'), xmax = as_date('2020-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red") +
  annotate("rect", xmin = as_date('2017-07-15'), xmax = as_date('2017-09-21'), ymin = 0, ymax = .2,
           alpha = .1,fill = "red")
dev.off()

# pdf("ClusterCurves6.pdf", width = 16, height = 9)
# out_fitted_curves %>%
#   ggplot(aes(x = date, ymin = lower, ymax = upper, fill = cluster, y = median, color = cluster )) +
#   geom_ribbon(alpha = .3) +
#   geom_line(linetype = 3) +
#   theme_bw() +
#   xlab('') + ylab('Prevalence') +
#   theme(legend.position = 'none')+
#   #dev.off() +
#   geom_rug(inherit.aes = F, aes(x = date), data = bats %>% mutate(date = as_date(date_end))) +
#   geom_rug(inherit.aes = F, aes(x = date), data = bats %>% filter(positive == 1) %>% mutate(date = as_date(date_end), cluster = type), color = 'red', sides = 'b') +
#   xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
#   scale_y_continuous(breaks=c(0,.10, .20))
# dev.off()
