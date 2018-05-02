#------------------------------------------#
# analysis of MST model fits from rStan ####
#------------------------------------------#

# This script analyses the MST model fits by looking at the relationship between temperature-corrected rate and mass-corrected biomass using the models fitted in a Bayesian framework using rStan
# Recreates Figure 2, Figure 3 and Figure 4 of Padfield et al. (2018) Ecology Letters from the Bayesian models (not included in the paper)

# clear workspace
mise::mise(vars = TRUE, pkgs = TRUE, console = TRUE, figs = TRUE)

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)
library(loo)
library(rethinking)
library(purrr)
library(readr)

# custom functions for manipulating rStan model objects ####

# create dplyr functions to return credible intervals of stan predictions
CI_stan_preds <- function(data, conf_level = 0.95){
  data %>% 
    gather(., 'Nnew', 'samps') %>%
    mutate(., Nnew = readr::parse_number(Nnew)) %>%
    group_by(., Nnew) %>%
    summarise(., mu = quantile(samps, probs = 0.5),
              lwr_CI = quantile(samps, probs = 0 + (0.5 - conf_level/2)),
              upr_CI = quantile(samps, probs = 1 - (0.5 - conf_level/2))) %>%
    ungroup()
}

# function to take extract parameter samples into dataframe
get_params <- function(param, samples){
  temp <- samples[[param]]
  t <- data.frame(param = param, mu = temp, stringsAsFactors = FALSE)
  return(t)
}

# load in data raw rate data ####
d_rates <- readRDS('raw_data/rates.rds') %>%
  mutate(., ancest_code = ifelse(ancestral == 'amb', 0, 1),
         ikt = (1/8.62e-05/(273.15 + 18)) - (1/8.62e-05/(temp + 273.15))) %>%
  arrange(., id) %>%
  filter(., id != '20_amb_A_16')

# GP model selection ####

# load in models
GP_mods <- readRDS('processed_data/GP_stan_fits.rds')

# extract the log likelihood of each model
loglik1 <- extract_log_lik(GP_mods$model1)
loglik2 <- extract_log_lik(GP_mods$model2)
loglik3 <- extract_log_lik(GP_mods$model3)
loglik4 <- extract_log_lik(GP_mods$model4)

# compare models using loo, according to Gelman
loo::compare(loo(loglik1),
        loo(loglik2))
loo::compare(loo(loglik2),
        loo(loglik3))
loo::compare(loo(loglik3),
        loo(loglik4))
# best model is model 4

GP_comp_mods <- rethinking::compare(GP_mods$model1,
                    GP_mods$model2,
                    GP_mods$model3,
                    GP_mods$model4)
plot(GP_comp_mods)
# best model is model 4

# CR model selection ####

# load in models
R_mods <- readRDS('processed_data/stan_output_R.rds')

loglik1 <- extract_log_lik(R_mods$model1)
loglik2 <- extract_log_lik(R_mods$model2)
loglik3 <- extract_log_lik(R_mods$model3)
loglik4 <- extract_log_lik(R_mods$model4)

# compare models
loo::compare(loo(loglik1),
        loo(loglik2))
loo::compare(loo(loglik2),
        loo(loglik3))
loo::compare(loo(loglik3),
        loo(loglik4))

# compare models
R_comp_mods <- rethinking::compare(R_mods$model1,
                    R_mods$model2,
                    R_mods$model3,
                    R_mods$model4,
                    func = WAIC)
plot(R_comp_mods)

R_comp_mods2 <- rethinking::compare(R_mods$model1,
                                   R_mods$model2,
                                   R_mods$model3,
                                   R_mods$model4,
                                   func = LOO)
plot(R_comp_mods2)

# best model changes from model 3 (loo) to model 4 (waic)
# best advice says we cannot trust the loo estimate for models 1-3

# plot the parameters of model 3
plot(R_mods$model3, pars = c('ln_r0', 'ln_r0_heat', 'bT', 'a'))
# 95% CI of ln_r0_heat still encompass 0 so going to delete it from the model - does not improve predictive accuracy

# best model is model 4 - no interactions of long-term warming on alpha, the activation energy or metabolic normalisation constant

#-----------------------------------#
# Process favoured model for GPP ####
#-----------------------------------#

# extract samples, compute mean, 95% credible intervals for estimated parameters
GP_samples <- rstan::extract(GP_mods$model4)

# get names of all parameters
names(GP_samples)

# extract parameter estimates
params <- map_df(c('ln_r0', 'bT', 'a'), get_params, GP_samples)

# quick plot
ggplot(params) +
  geom_density(aes(mu), fill = 'green4', alpha = 0.5) +
  facet_wrap(~ param, scales = 'free')

# create 95% CI for each parameter
GP_param_CI <- group_by(params, param) %>%
  summarise(., lwr_CI = quantile(mu, 0.025),
            upr_CI = quantile(mu, 0.975),
            mu = quantile(mu, 0.5))

# mle GPP output 
# ln_r0_eco_GP = -3.46 95% CI: -6.27 - -1.02
# a_GP = 0.88 95% CI: 0.57 - 1.17
# E_GP = 0.61 eV 95% CI: 0.11 - 1.12

# extract temperature corrected rate of each community
temp_cor_rate <- data.frame(GP_samples[['temp_cor_mu']]) %>%
  CI_stan_preds(.) %>%
  rename(., mu_tempcor = mu, lwr_CI_tempcor = lwr_CI, upr_CI_tempcor = upr_CI) %>%
  mutate(., flux = 'GP',
         ancestral = d_rates$ancestral)
  
# extract mass corrected biomass of each community
size_cor_biomass <- data.frame(GP_samples[['size_cor_biom_mu']]) %>%
  CI_stan_preds(.) %>%
  rename(., mu_sizecor = mu, lwr_CI_sizecor = lwr_CI, upr_CI_sizecor = upr_CI) %>%
  mutate(., flux = 'GP',
         ancestral = d_rates$ancestral)

# combine these together
GP_preds <- cbind(select(temp_cor_rate, mu_tempcor, lwr_CI_tempcor, upr_CI_tempcor, flux, ancestral),
                  select(size_cor_biomass, mu_sizecor, lwr_CI_sizecor, upr_CI_sizecor)) %>%
  mutate(ln_rate = d_rates$ln_GP,
         treatment = d_rates$treatment)


#----------------------------------#
# Process favoured model for CR ####
#----------------------------------#

# extract samples, compute mean, 95% credible intervals of estimated parameters
R_samples <- rstan::extract(R_mods$model4)

# get names of all parameters
names(R_samples)

# extract parameter estimates
params <- map_df(c('ln_r0', 'bT', 'a'), get_params, R_samples)

# quick plot
ggplot(params) +
  geom_density(aes(mu), fill = 'red', alpha = 0.5) +
  facet_wrap(~ param, scales = 'free')

# create 95% CI for each parameter
R_param_CI <- group_by(params, param) %>%
  summarise(., lwr_CI = quantile(mu, 0.025),
            upr_CI = quantile(mu, 0.975),
            mu = quantile(mu, 0.5))

# extract temperature corrected rate of each community
temp_cor_rate <- data.frame(R_samples[['temp_cor_mu']]) %>%
  CI_stan_preds(.) %>%
  rename(., mu_tempcor = mu, lwr_CI_tempcor = lwr_CI, upr_CI_tempcor = upr_CI) %>%
  mutate(., flux = 'R',
         ancestral = d_rates$ancestral)

# extract mass corrected biomass of each community
size_cor_biomass <- data.frame(R_samples[['size_cor_biom_mu']]) %>%
  CI_stan_preds(.) %>%
  rename(., mu_sizecor = mu, lwr_CI_sizecor = lwr_CI, upr_CI_sizecor = upr_CI) %>%
  mutate(., flux = 'R',
         ancestral = d_rates$ancestral)

# create predictions dataframe
R_preds <- cbind(select(temp_cor_rate, mu_tempcor, lwr_CI_tempcor, upr_CI_tempcor, flux, ancestral),
                  select(size_cor_biomass, mu_sizecor, lwr_CI_sizecor, upr_CI_sizecor)) %>%
  mutate(., ln_rate = d_rates$ln_R,
         treatment = d_rates$treatment)

#---------------------------------------------------------------------#
# analysis of temperature-corrected rate vs mass-corrected biomass ####
#---------------------------------------------------------------------#

# GPP analysis 

# fit with SMA for GPP
fit_GP_sma <- smatr::sma(mu_tempcor ~ mu_sizecor, GP_preds, na.action = na.fail, slope.test = 1, robust = TRUE)

# predictions for GPP
GP_preds <- mutate(GP_preds, preds = coef(fit_GP_sma)[1] + coef(fit_GP_sma)[2]*mu_sizecor)

# GP plot
plot_GP <- ggplot(GP_preds, aes(mu_sizecor, mu_tempcor, col = ancestral)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr_CI_tempcor, ymax = upr_CI_tempcor)) +
  geom_errorbarh(aes(xmin = lwr_CI_sizecor, xmax = upr_CI_sizecor)) +
  geom_line(aes(mu_sizecor, preds), col = 'black') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  ylim(3, 7.25) +
  xlim(3, 7.25) +
  geom_segment(aes(x = 4, y = 4, xend = 7, yend = 7), linetype = 2, col = 'black') +
  xlab(expression(Mass~corrected~biomass~(µg~C^-1))) +
  ylab(expression(atop(Temperature~corrected, community~flux~(µmol~O[2]~L^-1~hr^-1)))) +
  ggtitle(expression((a)~Gross~primary~production))

# CR analysis

# fit with SMA for R
fit_R_sma <- smatr::sma(mu_tempcor ~ mu_sizecor, R_preds, na.action = na.fail, slope.test = 1, robust = TRUE)

# predictions for R
R_preds <- mutate(R_preds, preds = coef(fit_R_sma)[1] + coef(fit_R_sma)[2]*mu_sizecor)

# R plot
plot_R <- ggplot(R_preds, aes(mu_sizecor, mu_tempcor, col = ancestral)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr_CI_tempcor, ymax = upr_CI_tempcor)) +
  geom_errorbarh(aes(xmin = lwr_CI_sizecor, xmax = upr_CI_sizecor)) +
  geom_line(aes(mu_sizecor, preds), col = 'black') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  ylim(2.2,5.9) +
  xlim(2.2,5.9) +
  geom_segment(aes(x = 2.65, y = 2.65, xend = 5.5, yend = 5.5), linetype = 2, col = 'black') +
  xlab(expression(Mass~corrected~biomass~(µg~C^-1))) +
  ylab(expression(atop(Temperature~corrected, community~flux~(µmol~O[2]~L^-1~hr^-1)))) +
  ggtitle(expression((b)~Community~respiration))

# arrange final plot, equivalent of Figure 3
p1 <- gridExtra::grid.arrange(plot_GP + theme(legend.position = 'none'), plot_R + theme(legend.position = 'none'), ncol = 2)
ggsave('plots/Figure_3_stan.pdf', plot = p1, width = 10, height = 5)

#----------------------------------------------------------------------------#
# analysis of total abundance versus estimated individual metabolic rate  ####
#----------------------------------------------------------------------------#

# get individual photosynthesis out and save
ind_rate_cor <- data.frame(GP_samples[['ind_mu']]) %>%
  CI_stan_preds(.) %>%
  mutate(., flux = 'GP',
         id = d_rates$id,
         pond = d_rates$pond,
         df = d_rates$df) %>%
  separate(., id, c('pond', 'ancestral', 'treatment', 'temp'), sep = "_") %>%
  unite(., id, c(pond, ancestral, treatment), sep = '_')

# load in raw count data and process ####
d_counts <- readRDS('raw_data/counts.rds') %>%
  separate(., id, c('pond', 'ancestral', 'treatment'), sep = "_", remove = FALSE) %>%
  select(., id, pond, treatment, ancestral, carbon, df) %>%
  group_by(., pond, treatment, ancestral, df, id) %>%
  dplyr::summarise(., count = n()) %>%
  data.frame() %>%
  mutate(., count = count * df,
         log_ntot = log(count))

# merge predictions with log_ntot
ind_rate_cor <-  merge(ind_rate_cor, select(d_counts, id, log_ntot, ancestral), by = 'id')

# use sma to look at log ntot vs log Mind
fit_sma <- smatr::sma(log_ntot ~ mu, ind_rate_cor, slope.test = -1, robust = TRUE)

ind_rate_cor <- mutate(ind_rate_cor, preds = coef(fit_sma)[1] + coef(fit_sma)[2]*mu)

# plot (equivalent of Figure 4)
ggplot(ind_rate_cor, aes(mu, log_ntot)) +
  geom_point(aes(col = ancestral)) +
  geom_line(aes(mu, preds)) + 
  geom_errorbarh(aes(xmin = lwr_CI, xmax = upr_CI, col = ancestral)) +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  xlab(expression(ln~Average~individual~gross~photosynthesis(µmol~O[2]~L^-1~hr^-1))) +
  ylab(expression(ln~Total~community~abundance~(count~L^-1))) +
  xlim(-14.1, -10.7) +
  ylim(14.5, 20) +
  theme(legend.position = 'none')

ggsave('plots/Figure_4_stan.pdf', plot = last_plot(), width = 7, height = 5)



