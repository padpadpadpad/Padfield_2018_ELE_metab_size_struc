#---------------------------------------------------------------#
# Analysis of total abundance and estimated P per individual ####
#---------------------------------------------------------------#

# This script analyses the relationship between total abundance and estimated individual metabolic rate across all communities
# Recreates Figure 4 of Padfield et al. (2018) Ecology Letters

# clear workspace
mise::mise(vars = TRUE, pkgs = TRUE, console = TRUE, figs = TRUE)

# load in packages
library(smatr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(purrr)

# load in raw count data and process ####
d_counts <- readRDS('raw_data/counts.rds')

# output of mle GPP analysis
# ln_r0_eco_GP = -3.46 95% CI: -6.27 - -1.02
# a_GP = 0.88 95% CI: 0.57 - 1.17
# E_GP = 0.61 eV 95% CI: 0.11 - 1.12

# create size_corrected biomass using value of alpha
d_count_sum <- separate(d_counts, id, c('pond', 'ancestral', 'treatment'), sep = "_") %>%
  select(., pond, treatment, ancestral, carbon, df) %>%
  group_by(., pond, treatment, ancestral, df) %>%
  dplyr::summarise(., Mtot_cor = sum(carbon^0.88),
            count = n()) %>%
  data.frame() %>%
  mutate(., Mtot_cor = Mtot_cor * df,
         count = count * df)

# add column for estimated GPP, estimated P per individual and log count
d_count_sum <-mutate(d_count_sum, 
                      log_ntot = log(count),
                      temp = ifelse(treatment == 'A', 16, 20),
                      ln_Mtot_cor_GP = log(Mtot_cor*exp(-3.46)*exp(0.61*(1/(8.62e-5*(18+273.15)) - (1/(8.62e-5*(temp+273.15)))))),
                      log_Mind = log(exp(ln_Mtot_cor_GP)/count))

# quickplot
ggplot(d_count_sum) +
  geom_point(aes(log_Mind, log_ntot, col = treatment, shape = ancestral)) +
  scale_color_manual(values = c('black', 'red')) +
  stat_smooth(aes(log_Mind, log_ntot), method = 'lm', se = FALSE, col = 'black') +
  xlab(expression(log[10]~M[ind]^alpha)) +
  ylab(expression(log[10]~count)) +
  theme_bw(base_size = 14, base_family = 'Helvetica')

# use sma to look at log ntot vs log Mind
fit_sma <- sma(log_ntot ~ log_Mind, d_count_sum, slope.test = -1, robust = TRUE)

# boostrap predictions for confidence intervals of sma plot ####

# get high resolution of prediction variable
boot_preds <- data.frame(expand.grid(log_Mind = seq(min(d_count_sum$log_Mind, na.rm = T), max(d_count_sum$log_Mind, na.rm = T), length.out = 200), stringsAsFactors = FALSE))

# bootstrap data and fits
fit_boots <- d_count_sum %>%
  modelr::bootstrap(n = 1000, id = 'boot_num') %>%
  group_by(boot_num) %>%
  mutate(., fit = map(strap, ~ sma(log_ntot ~ log_Mind, data.frame(.), robust = TRUE))) %>%
  ungroup() %>%
  mutate(., intercept = map_dbl(fit, ~coef(.x)[1]),
         slope = map_dbl(fit, ~coef(.x)[2])) %>%
  select(., -fit) %>%
  group_by(boot_num) %>%
  do(data.frame(fitted = .$intercept + .$slope*boot_preds$log_Mind, 
                log_Mind = boot_preds$log_Mind)) %>%
  ungroup() %>%
  group_by(., log_Mind) %>%
  dplyr::summarise(., conf_low = quantile(fitted, 0.025),
            conf_high = quantile(fitted, 0.975)) %>%
  ungroup()

# get predictions from the bootstrapped models
preds <- data.frame(expand.grid(log_Mind = seq(min(d_count_sum$log_Mind, na.rm = T), max(d_count_sum$log_Mind, na.rm = T), length.out = 40), stringsAsFactors = FALSE)) %>%
  mutate(., preds = coef(fit_sma)[1] + coef(fit_sma)[2]*log_Mind)

# make final plot #####
ggplot() +
  geom_point(aes(log_Mind, log_ntot, col = ancestral), d_count_sum) +
  geom_line(aes(log_Mind, preds), preds) +
  geom_ribbon(aes(x = log_Mind, ymin = conf_low, ymax = conf_high), alpha = 0.1, fit_boots) +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  xlab(expression(ln~Average~individual~gross~photosynthesis(µmol~O[2]~L^-1~hr^-1))) +
  ylab(expression(ln~Total~community~abundance~(count~L^-1))) +
  xlim(-13.6, -10.7) +
  ylim(14.5, 20.5)+
  theme(legend.position = 'none')

# save final plot ####
ggsave('plots/Figure_4.pdf', plot = last_plot(), width = 6, height = 5)