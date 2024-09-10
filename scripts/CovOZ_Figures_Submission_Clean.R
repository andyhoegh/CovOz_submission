library(tidyverse)
library(lubridate)
library(cowplot)
library(ggtext)
library(jpeg)
library(scales)

#### Load colours and names ####
colours_clades6 <- c("#191919",viridis_pal(option = "B")(6)[2:5], "#fff207") 
colours_clades4 <- colours_clades6[2:5]
red_viridis <- "#DD513AFF"
colours_repro <- c(rep("#696969",4))
alpha_repro <- c(0.02, 0.08, 0.14, 0.22)

### ### ### ### ### ###
#### Figure 2 ####
### ### ### ### ### ###
out_fitted_curves <- read_csv('data/model_output/cluster_curves.csv')
bats <- read_csv('data/combined_out_variant.csv') %>%
  mutate(positive = ifelse(variant_positive=="FALSE", 0, 1))

png("Figure2_final.png", width = 5, height = 6, units = 'in', res = 1500) 
out_fitted_curves %>%
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = cluster, colour = cluster,y = median )) +
  # REPRODUCTIVE SHADING
  # mating
  annotate("rect", xmin = as_date('2018-03-01'), xmax = as_date('2018-03-31'), ymin = 0, ymax = .25, alpha = alpha_repro[4],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-03-01'), xmax = as_date('2019-03-31'), ymin = 0, ymax = .25, alpha = alpha_repro[4],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-03-01'), xmax = as_date('2020-03-31'), ymin = 0, ymax = .25, alpha = alpha_repro[4],fill = "#696969") +
  geom_text(x = as_date('2018-03-15'), y = .23, label = 'M', color = 'grey70', size =2.3)+
  # pregnancy shading
  annotate("rect", xmin = as_date('2018-04-01'), xmax = as_date('2018-09-30'), ymin = 0, ymax = .25, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-04-01'), xmax = as_date('2019-09-30'), ymin = 0, ymax = .25, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-04-01'), xmax = as_date('2020-09-30'), ymin = 0, ymax = .25, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-07-01'), xmax = as_date('2017-09-30'), ymin = 0, ymax = .25, alpha = alpha_repro[1],fill = "#696969") +
  geom_text(x = as_date('2017-08-01'), y = .23, label = 'P', color = 'grey70', size =2.3)+
  # birth shading
  annotate("rect", xmin = as_date('2018-10-01'), xmax = as_date('2018-10-31'), ymin = 0, ymax = .25, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-10-01'), xmax = as_date('2019-10-31'), ymin = 0, ymax = .25, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-10-01'), xmax = as_date('2020-10-31'), ymin = 0, ymax = .25, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-10-01'), xmax = as_date('2017-10-31'), ymin = 0, ymax = .25, alpha = alpha_repro[2],fill = "#696969") +
  geom_text(x = as_date('2017-10-15'), y = .23, label = 'B', color = 'grey70', size =2.3)+
  # lactation shading
  annotate("rect", xmin = as_date('2018-11-01'), xmax = as_date('2019-03-01'), ymin = 0, ymax = .25, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-11-01'), xmax = as_date('2020-03-01'), ymin = 0, ymax = .25, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-11-01'), xmax = as_date('2021-02-01'), ymin = 0, ymax = .25, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-11-01'), xmax = as_date('2018-03-01'), ymin = 0, ymax = .25, alpha = alpha_repro[3],fill = "#696969") +
  geom_text(x = as_date('2017-12-31'), y = .23, label = 'L', color = 'grey70', size =2.3)+
  facet_grid(rows = vars(cluster)) +
  geom_ribbon(alpha = .4,linewidth=0) +
  geom_line(lty = 1,linewidth = 1.2) +
  scale_fill_manual(values = colours_clades6) +
  scale_colour_manual(values = colours_clades6) +
  theme_cowplot() + 
  theme(legend.position = "none") +
  xlab('') + ylab('Prevalence') +
  geom_rug(inherit.aes = F, aes(x = date), data = bats %>% mutate(date = as_date(date_end))) +
  geom_rug(inherit.aes = F, aes(x = date), data = bats %>% filter(positive == 1) %>% mutate(date = as_date(date_end), cluster = type), color = 'red', sides = 'b') +
  xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
  scale_y_continuous(breaks=c(0,.10, .20)) +
  scale_x_date(breaks = seq(as.Date("2017-09-01"), as.Date("2020-09-01"), by = "6 months"),
               date_labels = "%b\n%Y",limits = c(as.Date("2017-07-01"), as.Date("2021-02-01"))) 

dev.off()

### ### ### ### ### ###
#### Figure 3 ####
### ### ### ### ### ###
out_fitted_curves <- read_csv('data/model_output/cluster_curves.csv') %>% 
  group_by(cluster) %>%
  filter(!duplicated(date)) %>%
  ungroup()

tmp <- out_fitted_curves %>% 
  bind_rows(tibble(median = 0, lower = 0, upper = 0, date = mdy("01-01-2020"), cluster = 'beta 2d.iii')) %>% # workaround to get January axis to show
  mutate(new_date = substr(decimal_date(date),5,8), 
         time = paste('2020', new_date, sep = ''), 
         date2 = date_decimal(as.numeric(time)),
         date3 = date(date2),
         year = as.character(year(date))) 

png("Figure3_final.png", width = 9, height = 4, units = 'in', res = 800)

tmp %>% 
  filter(!cluster %in% c('beta 2d.i', 'beta 2d.vi') ) %>%
  filter(year %in% c(2018,2019,2020)) %>%
  ggplot(aes(x = date3, y = median, ymax = upper, ymin = lower, color = cluster, fill = cluster)) +
  # REPRODUCTIVE SHADING
  # mating shading
  annotate("rect", xmin = as_date('2020-03-01'), xmax = as_date('2020-03-31'), ymin = 0, ymax = .3, alpha = alpha_repro[4],fill = "#696969") +
  geom_text(x = as_date('2020-03-15'), y = .29, label = 'M', color = 'grey70', size =2)+
  # pregnancy shading
  annotate("rect", xmin = as_date('2020-04-01'), xmax = as_date('2020-09-30'), ymin = 0, ymax = .3, alpha = alpha_repro[1],fill = "#696969") +
  geom_text(x = as_date('2020-07-01'), y = .29, label = 'P', color = 'grey70', size =2)+
  # birth shading
  annotate("rect", xmin = as_date('2020-10-01'), xmax = as_date('2020-10-31'), ymin = 0, ymax = .3, alpha = alpha_repro[2],fill = "#696969") +
  geom_text(x = as_date('2020-10-15'), y = .29, label = 'B', color = 'grey70', size =2)+
  # lactation shading
  annotate("rect", xmin = as_date('2020-11-01'), xmax = as_date('2020-12-31'), ymin = 0, ymax = .3, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-01-01'), xmax = as_date('2020-03-01'), ymin = 0, ymax = .3, alpha = alpha_repro[3],fill = "#696969") +
  geom_text(x = as_date('2020-12-07'), y = .29, label = 'L', color = 'grey70', size =2)+
  geom_text(x = as_date('2020-01-31'), y = .29, label = 'L', color = 'grey70', size =2)+
  facet_grid(cols = vars(cluster)) + 
  geom_ribbon(aes(lty = year),alpha = 0.2,linewidth=0) +
  geom_line(aes(lty = year),linewidth = 1) +
  scale_linetype_manual(values = c('solid', 'dotted', 'dashed'))+
  theme_set(theme_cowplot())+
  xlab('') + ylab('Prevalence') +
  theme(legend.position = 'right')+
  labs(linetype = "Year") + 
  guides(color = "none", fill = "none") +
  scale_fill_manual(values = colours_clades4)+
  scale_colour_manual(values = colours_clades4)+
  scale_y_continuous(limits=c(-.01,.3)) +
  scale_x_date(date_labels = "%b", breaks = "3 month", date_minor_breaks = "1 month") 

dev.off()


### ### ### ### ### ###
### Figure 4 A-B (individual bar plot prevalence)  ####
### ### ### ### ### ###

# read data
df <- read.csv('data/individual_variant_covariates.csv')

### create a barplot with prevalence ###
# get positive results 
df1 <- df %>%
  group_by(bat_age, bat_species,type, variant_positive) %>%
  filter(variant_positive == 'TRUE') %>%
  summarize(positive = n()) %>% 
  select('bat_age', 'bat_species', 'type', 'positive') %>%
  as.data.frame()   %>%
  ungroup() 
# calculate total tested 
df2 <- df %>% 
  group_by(bat_age, bat_species,type) %>%
  summarize(total_tested = n()) %>%
  group_by(bat_age, bat_species) %>%
  slice(1) %>%
  select(1,2,4) %>% 
  filter(bat_age == 'adult' | bat_age == 'sub_adult' | bat_age == 'juve')%>%
  ungroup() 

# create a full grid of combinations 
df_plotNew <- expand.grid(bat_species = c('bff', 'ghff'), 
                          bat_age = c('juve', 'sub_adult', 'adult'),
                          type = c('beta 2d.i','beta 2d.ii','beta 2d.iii', 'beta 2d.iv', 'beta 2d.v', 'beta 2d.vi'))

###
df_plot1 <- merge(df_plotNew, df1, by = c('bat_species', 'bat_age','type'), all.x=TRUE) %>%
  mutate(positive = ifelse(is.na(positive), 0,positive)) %>%
  right_join(., df2, by = c('bat_species', 'bat_age')) %>%
  mutate(bat_age= factor(bat_age, levels = c('juve','sub_adult','adult'), labels=c('Juvenile','Subadult','Adult'))) %>% 
  mutate(prevalence = (signif(positive/total_tested, digits = 2)*100)) %>% 
  rowwise() %>% 
  mutate(lci = binom.test(x = positive,n = total_tested)$conf.int[[1]]*100,
         uci = binom.test(x = positive,n = total_tested)$conf.int[[2]]*100)

levels(df_plot1$type) <- gsub("beta ", "", levels(df_plot1$type), fixed=TRUE)

df_plot1$bat_species <- ifelse(df_plot1$bat_species == 'bff', 'Black Flying Fox',
                               'Grey-Headed Flying Fox')

bff_prev_plot <- ggplot(df_plot1 %>% filter(bat_species=='Black Flying Fox')) +  
  aes(y = fct_rev(type), x = prevalence, fill = type) +
  geom_col(alpha = 0.6) +
  geom_errorbar(aes(xmin = lci, xmax = uci),alpha = 0.3,width = 0.15, position = position_dodge(0.5))+
  geom_vline(xintercept=0)+
  geom_segment(aes(y=0.41,yend=0.41, x=0, xend=20))+
  geom_text(aes(label= paste(sprintf("%.1f", round(prevalence, 1)), '%', sep = ''), x = -1.2), size = 3.5,
            position=position_dodge(width=0.5),
            vjust=0.4) +
  scale_fill_manual(values = colours_clades6)+
  labs(y = "", x = "Prevalence (%)",fill="Clade") +
  theme_set(theme_cowplot())+
  theme(legend.position = "none",
        panel.spacing.x = unit(0.5, "lines"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(face = "bold"),
        plot.margin = margin(5.5, 5.5, 0, 5.5)) +
  scale_x_continuous(breaks = seq(0,25, 5), limits = c(-1.2,20)) +
  facet_grid(bat_age~bat_species, scales = 'free')

ghff_prev_plot <- ggplot(df_plot1 %>% filter(bat_species=='Grey-Headed Flying Fox')) +  
  aes(y = fct_rev(type), x = prevalence, fill = type) +
  geom_col(alpha = 0.6) +
  geom_errorbar(aes(xmin = lci, xmax = uci),alpha = 0.3,width = 0.15, position = position_dodge(0.5))+
  geom_vline(xintercept=0)+
  geom_segment(aes(y=0.41,yend=0.41, x=0, xend=100))+
  geom_text(aes(label= paste(sprintf("%.1f", round(prevalence, 0)), '%', sep = ''), x = -6), size = 3.5,
            position=position_dodge(width=0.5),
            vjust=0.4) +
  scale_fill_manual(values = colours_clades6)+
  labs(y = "", x = "Prevalence (%)",fill="Clade") +
  theme_set(theme_cowplot())+
  theme(#legend.position = "none",
    panel.spacing.x = unit(0.5, "lines"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)) +
  scale_x_continuous(breaks = seq(0,100, 25), limits = c(-6,100)) +
  facet_grid(bat_age~bat_species, scales = 'free')


### ### ### ### ### ###
### Figures 4 C-D BFF and GHFF fitted curves  ####
### ### ### ### ### ###
load('data/logistic_curve_out.RData') #loads out_fitted_curves
bats <- read_csv('data/individual_variant_covariates.csv')

# GROUP BY FIRST DATE OF EACH SAMPLING SESSION
bat_prop  <- bats %>% filter(bat_age %in% c('adult','juve','sub_adult') ) %>%
  group_by(type, bat_age, bat_species, date_end) %>%
  summarize(n = n(), n_pos = sum(variant_positive), mean_pos = mean(variant_positive),
            .groups = 'drop') %>% 
  rename(date = date_end)

bff_curves <- out_fitted_curves %>%
  filter(species == 'bff', 
         type %in% c('beta 2d.ii', 'beta 2d.iv','beta 2d.v', 'beta 2d.vi')) %>% 
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, colour = type, y = mean )) +
  geom_point(inherit.aes = F, aes(x = date, y = mean_pos, size = n), shape = 1, 
             data = bat_prop %>% filter(bat_species == 'bff',
                                        type %in% c('beta 2d.ii', 'beta 2d.iv','beta 2d.v', 'beta 2d.vi')) %>% mutate(age = bat_age), alpha = .25) +
  geom_ribbon(alpha = .4,linewidth=0) +
  geom_line(lty = 1,linewidth = 1.2) +
  facet_grid(factor(age, levels = c('juve','sub_adult','adult'), labels=c('Juvenile','Subadult','Adult'))~type) +
  scale_fill_manual(values = colours_clades4)+
  scale_colour_manual(values = colours_clades4)+
  scale_size_continuous(range = c(1/3,25/3), breaks = seq(0, 25, by = 5),'# samples')+
  scale_x_date(breaks = seq(as.Date("2018-03-01"), as.Date("2020-03-01"), by = "12 months"), date_labels = "%b\n%Y",limits = c(as.Date("2017-07-01"), as.Date("2021-02-01"))) +
  theme_set(theme_cowplot())+
  xlab('') + ylab('Prevalence') +
  labs(subtitle = "Black flying foxes\n")+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  ylim(0, 1) +
  guides(color = "none", fill = "none") 

ghff_curves <- out_fitted_curves %>%
  filter(species == 'ghff', 
         type %in% c('beta 2d.iii')) %>% 
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, colour = type, y = mean )) +
  geom_point(inherit.aes = F, aes(x = date, y = mean_pos, size = n),  shape = 1,
             data = bat_prop %>% filter(bat_species == 'ghff',
                                        type %in% c('beta 2d.iii')) %>% mutate(age = bat_age), alpha = .25) +
  geom_ribbon(alpha = .4,linewidth=0) +
  geom_line(lty = 1,linewidth = 1.2) +
  facet_grid(factor(age, levels = c('juve','sub_adult','adult'), labels=c('Juvenile','Subadult','Adult'))~type) +
  scale_fill_manual(values = colours_clades4)+
  scale_colour_manual(values = colours_clades4)+
  scale_size_continuous(range = c(1/5,25/5), breaks = seq(0, 25, by = 5),'# samples')+
  scale_x_date(breaks = seq(as.Date("2018-07-01"), as.Date("2020-07-01"), by = "12 months"), date_labels = "%b\n%Y",limits = c(as.Date("2017-11-01"), as.Date("2019-06-01"))) +
  theme_set(theme_cowplot())+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  xlab('') + ylab('') +
  labs(subtitle = "Grey-headed\nflying foxes")+
  ylim(0, 1) +
  guides(color = "none", fill = "none")

### ### ### ### ### ###
### Figure 4 - all combined ####
### ### ### ### ### ###
png("Figure4_A-D_final.png", width = 12, height = 11.5, units = 'in', res = 1200)

a_b <- plot_grid(bff_prev_plot,ghff_prev_plot, ncol=2, rel_widths = c(0.95,1), labels = c("A","B"))
# Create a blank spacer plot
spacer_plot <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 2, 0))  # Adjust bottom margin 
c_d <- plot_grid(bff_curves +theme(legend.position = "none"),ghff_curves+theme(legend.position = "none"),get_legend(bff_curves), ncol=3, rel_widths = c(16,3.2,2),labels = c("C","D",""))
plot_grid(a_b,spacer_plot, c_d,nrow = 3, rel_heights = c(4.5,0.4,6))
dev.off()

### ### ### ### ### ###
### Figure 5 - Andy (or Dan?) has final code. Do we need to update this with new logistic curve output?
### ### ### ### ### ###

load('data/logistic_curve_out_417.RData') # should this change?

out_beta <-   out_beta %>%
  bind_cols(tibble(clade_spec= rep(c('beta 2d.ii', 'beta 2d.iv', 'beta 2d.v'), each = 12000))) %>%
  mutate(vals2 = case_when(
    type == 'beta1 - beta2' ~ -1 * vals,
    type == 'beta1 - beta3' ~ -1 * vals,
    type == 'beta2 - beta3' ~ vals
  )) %>%
  mutate(contrast = case_when(
    type == 'beta1 - beta2' ~ 'Juvenile - Adult',
    type == 'beta1 - beta3' ~ 'Subadult - Adult',
    type == 'beta2 - beta3' ~ 'Juvenile - Subadult',
  ))

probs <- out_beta %>% group_by(clade_spec, type) %>%
  summarize(prob = mean(vals2 > 0), .groups = 'drop' ) %>%
  mutate(contrast = case_when(
    type == 'beta1 - beta2' ~ 'Juvenile - Adult',
    type == 'beta1 - beta3' ~ 'Subadult - Adult',
    type == 'beta2 - beta3' ~ 'Juvenile - Subadult',
  ) )

png("Figures/Figure6_AP.png", width = 8, height = 6, units = 'in', res = 1200)
out_beta %>% ggplot(aes( x= vals2)) + geom_histogram(aes(y = after_stat(density)), binwidth = 0.1) +
  geom_vline(xintercept = 0, color = red_viridis) +
  geom_text(aes(label = paste('Prob = ',round(prob,3), sep = ''), y = 1.5, x = -.75), data = probs,size = 3) +
  facet_grid(contrast~ clade_spec, scales = "free_y") + theme_bw() +
  geom_segment(aes(x=-.75, y=1.45, xend=.2, yend=.1), arrow = arrow(length=unit(0.25, 'cm'))) +
  theme_cowplot()+
  xlab('')+
  ylab('Density')

dev.off()


### ### ### ### ### ###
### Supp Figure 8 
### ### ### ### ### ###

out_fitted_curves <- read_csv('data/model_output/cluster_curves.csv')

bats <- read_csv('data/combined_out_variant.csv') %>%
  mutate(variant_positive = ifelse(variant_positive=="FALSE", 0, 1))

bat_prop  <- bats %>% 
  group_by(type, site, date_end,session_type) %>%
  mutate(session_type = factor(session_type, levels = c("C", "R"), labels = c("Ind", "UR"))) %>% 
  summarize(n = n(), n_pos = sum(variant_positive), mean_pos = mean(variant_positive),
            .groups = 'drop') %>% 
  rename(date=date_end,
         cluster=type)  

alpha_reproSI = c(0.01, 0.05, 0.10, 0.15) # lighter shading to allow for points

png("SIFigure8.png", width = 6, height = 7, units = 'in', res = 1500) 

out_fitted_curves %>% 
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = cluster, colour = cluster,y = median )) +
  annotate("rect", xmin = as_date('2018-03-01'), xmax = as_date('2018-03-31'), ymin = 0, ymax = .6, alpha = alpha_repro[4],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-03-01'), xmax = as_date('2019-03-31'), ymin = 0, ymax = .6, alpha = alpha_repro[4],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-03-01'), xmax = as_date('2020-03-31'), ymin = 0, ymax = .6, alpha = alpha_repro[4],fill = "#696969") +
  geom_text(x = as_date('2018-03-15'), y = .56, label = 'M', color = 'grey70', size =2.3)+
  # pregnancy shading
  annotate("rect", xmin = as_date('2018-04-01'), xmax = as_date('2018-09-30'), ymin = 0, ymax = .6, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-04-01'), xmax = as_date('2019-09-30'), ymin = 0, ymax = .6, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-04-01'), xmax = as_date('2020-09-30'), ymin = 0, ymax = .6, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-07-01'), xmax = as_date('2017-09-30'), ymin = 0, ymax = .6, alpha = alpha_repro[1],fill = "#696969") +
  geom_text(x = as_date('2017-08-01'), y = .56, label = 'P', color = 'grey70', size =2.3)+
  # birth shading
  annotate("rect", xmin = as_date('2018-10-01'), xmax = as_date('2018-10-31'), ymin = 0, ymax = .6, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-10-01'), xmax = as_date('2019-10-31'), ymin = 0, ymax = .6, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-10-01'), xmax = as_date('2020-10-31'), ymin = 0, ymax = .6, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-10-01'), xmax = as_date('2017-10-31'), ymin = 0, ymax = .6, alpha = alpha_repro[2],fill = "#696969") +
  geom_text(x = as_date('2017-10-15'), y = .56, label = 'B', color = 'grey70', size =2.3)+
  # lactation shading
  annotate("rect", xmin = as_date('2018-11-01'), xmax = as_date('2019-03-01'), ymin = 0, ymax = .6, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-11-01'), xmax = as_date('2020-03-01'), ymin = 0, ymax = .6, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-11-01'), xmax = as_date('2021-02-01'), ymin = 0, ymax = .6, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-11-01'), xmax = as_date('2018-03-01'), ymin = 0, ymax = .6, alpha = alpha_repro[3],fill = "#696969") +
  geom_text(x = as_date('2017-12-31'), y = .56, label = 'L', color = 'grey70', size =2.3)+  
  facet_grid(rows = vars(cluster)) +
  geom_ribbon(alpha = .4, linewidth=0, show.legend = FALSE) +
  geom_line(lty = 1,linewidth = 1.2, show.legend = FALSE) +
  scale_fill_manual(values = colours_clades6) +
  scale_colour_manual(values = colours_clades6) +
  theme_set(theme_cowplot())+
  xlab('') + ylab('Prevalence') + labs(size = "Sample\n  size", shape = "Session\n  type")+
  xlim(min(out_fitted_curves$date), max(out_fitted_curves$date)) +
  scale_y_continuous(breaks=c(0,.2,.4,0.6)) +
  geom_point(inherit.aes = F, aes(x = as.Date(date), y = mean_pos, size = n,shape = session_type), data = bat_prop , alpha = .3) +#
  scale_shape_manual(values = c(1, 19))+ 
  scale_size_continuous(range = c(0.1, 3)) +  # Adjust the range to reduce overall point size
  scale_x_date(breaks = seq(as.Date("2017-09-01"), as.Date("2020-09-01"), by = "6 months"),
               date_labels = "%b\n%Y",limits = c(as.Date("2017-07-01"), as.Date("2021-02-01"))) 


dev.off()

### ### ### ### ### ###
### Supp Figure 9 
### ### ### ### ### ###
load('data/logistic_curve_out.RData') 
bats <- read_csv('data/individual_variant_covariates.csv')

# GROUP BY FIRST DATE OF EACH SAMPLING SESSION
bat_prop  <- bats %>% filter(bat_age %in% c('adult','juve','sub_adult') ) %>% 
  group_by(type, bat_age, bat_species, date_end) %>%
  summarize(n = n(), n_pos = sum(variant_positive), mean_pos = mean(variant_positive),
            .groups = 'drop') %>% 
  rename(date = date_end)

png("SIFigure9.png", width = 16, height = 7, units = 'in', res = 300)

bff_curves <- out_fitted_curves %>% 
  filter(species == 'bff', 
         type %in% c('beta 2d.ii', 'beta 2d.iv','beta 2d.v', 'beta 2d.vi')) %>% 
  ggplot(aes(x = date, ymin = lower, ymax = upper, fill = type, colour = type, y = mean )) +
  # REPRODUCTIVE SHADING
  # mating shading
  annotate("rect", xmin = as_date('2018-02-01'), xmax = as_date('2018-03-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[4],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-02-01'), xmax = as_date('2019-03-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[4],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-02-01'), xmax = as_date('2020-03-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[4],fill = "#696969") +
  geom_text(x = as_date('2018-03-01'), y = .6, label = 'M', color = 'grey70', size =2.3)+
  # pregnancy shading
  annotate("rect", xmin = as_date('2018-04-01'), xmax = as_date('2018-09-30'), ymin = 0, ymax = 0.65, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-04-01'), xmax = as_date('2019-09-30'), ymin = 0, ymax = 0.65, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-04-01'), xmax = as_date('2020-09-30'), ymin = 0, ymax = 0.65, alpha = alpha_repro[1],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-07-01'), xmax = as_date('2017-09-30'), ymin = 0, ymax = 0.65, alpha = alpha_repro[1],fill = "#696969") +
  geom_text(x = as_date('2017-07-01'), y = .6, label = 'P', color = 'grey70', size =2.3)+
  # birth shading
  annotate("rect", xmin = as_date('2018-10-01'), xmax = as_date('2018-10-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-10-01'), xmax = as_date('2019-10-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-10-01'), xmax = as_date('2020-10-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[2],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-10-01'), xmax = as_date('2017-10-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[2],fill = "#696969") +
  geom_text(x = as_date('2017-10-15'), y = .6, label = 'B', color = 'grey70', size =2.3)+
  # lactation shading
  annotate("rect", xmin = as_date('2018-11-01'), xmax = as_date('2019-01-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2019-11-01'), xmax = as_date('2020-01-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2020-11-01'), xmax = as_date('2021-01-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[3],fill = "#696969") +
  annotate("rect", xmin = as_date('2017-11-01'), xmax = as_date('2018-01-31'), ymin = 0, ymax = 0.65, alpha = alpha_repro[3],fill = "#696969") +
  geom_text(x = as_date('2017-12-15'), y = .6, label = 'L', color = 'grey70', size =2.3)+
  geom_line(lty = 1,linewidth = 1.2) +
  geom_jitter(inherit.aes = F, aes(x = date, y = mean_pos, size = n, colour = type),shape = 1, width = 4,height = 0,
              data = bat_prop %>% filter(bat_species == 'bff',
                                         type %in% c('beta 2d.ii', 'beta 2d.iv','beta 2d.v', 'beta 2d.vi')) %>% mutate(age = bat_age)) + 
  facet_grid(factor(age, levels = c('juve','sub_adult','adult'), labels=c('Juvenile','Subadult','Adult'))~type) +
  scale_fill_manual(values = colours_clades6[c(2,4,5,6)])+
  scale_colour_manual(values = colours_clades6[c(2,4,5,6)])+
  scale_size_continuous(range = c(1/3,25/3), breaks = seq(0, 25, by = 5),'# samples')+
  scale_x_date(breaks = seq(as.Date("2017-09-01"), as.Date("2020-09-01"), by = "6 months"), date_labels = "%b\n%Y",limits = c(as.Date("2017-07-01"), as.Date("2021-02-01"))) +
  theme_set(theme_cowplot())+
  xlab('') + ylab('Prevalence') +
  labs(colour = "Clade")+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  ylim(0, 0.7) +
  guides(fill = "none") 


# BFF ONLY
plot_grid(bff_curves, ncol=1) 
dev.off()
 

