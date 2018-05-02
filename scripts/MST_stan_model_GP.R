#--------------------------------------------------------------------------------------#
# fitting the MST model to GPP, temperature and size distribution data using rStan #####
#--------------------------------------------------------------------------------------#

# This script fits the MST model to the metabolism, temperature and size distribution data using a Bayesian framework. The code is essentially duplicated for GPP and CR
# Uses the model described in the theory section of Padfield et al. (2018) Ecology Letters

# The Bayesian analysis does not appear in the paper, but I feel it is a useful addition to the GitHub repository as it allows the propagation of uncertainty onto the predictions. This is especially useful when using the values of alpha and the activation energy to calculate mass-corrected biomass and temperature-corrected rate

# I am relatively new to Bayesian statistics and writing the models used here was challenging. Consequently, I cannot guarantee that I have adhered to best practice or that the code and models are completely correct. However, they do reassuringly give idential mean estimates to the maximum likelihood approach. And any deviations from best practice or errors were not through a lack of trying.

# clear workspace
mise::mise(vars = TRUE, pkgs = TRUE, console = TRUE, figs = TRUE)

# load in packages
library(rstan)
library(ggplot2)
library(tidyr)
library(dplyr)

# load in data and process ####
d_rates <- readRDS('raw_data/rates.rds')

# process rate data
d_rates <- mutate(d_rates, ancest_code = ifelse(ancestral == 'amb', 0, 1),
             ikt = (1/8.62e-05/(273.15 + 18)) - (1/8.62e-05/(K))) %>%
  arrange(., id)

# load in raw count data
d_count <- readRDS('raw_data/counts.rds')

# process count data
d_count <- separate(d_count, id, c('pond', 'ancestral', 'treatment'), sep = '_') %>%
  select(., pond, treatment, ancestral, carbon) %>%
  mutate(., temp = ifelse(treatment == 'W', 20, 16)) %>%
  unite(., id, c(pond, ancestral, treatment, temp), sep = '_', remove = FALSE) %>%
  group_by(., id) %>%
  mutate(., row = 1:n(),
         ancest_code = ifelse(ancestral == 'amb', 0, 1)) %>%
  data.frame()

# check which biomass sample was missed
d_rates$id[! d_rates$id %in% unique(cd$id)]

# drop 20_amb_A_16
d_rates <- filter(d_rates, id != '20_amb_A_16')

# create vector of when mass changes
# this is needed for the stan model
mass_n <- group_by(d_count, id) %>%
  summarise(n = n()) %>%
  data.frame()
df <- d_rates$df

# create vector position of the first mass for each sample
mass_n = mass_n$n
first_mass = c(1, head(cumsum(mass_n), 38) + 1)

# they are in the correct order
unique(d_count$id) == unique(d_rates$id)
# YIPPEE

# analysis of GPP ####

# create data
stan_datalist = list(N = nrow(d_rates),
                     temp = d_rates$ikt, 
                     ln_rate = d_rates$ln_GP,
                     K = nrow(d_count),
                     mass = log(d_count$carbon),
                     mass_n = mass_n,
                     first_mass = first_mass,
                     df = d_rates$df,
                     ancest = d_rates$ancest_code)

# create function for start values
stan_inits <- function() list(a = 0.75, bT = 0.7, ln_r0 = -3, tau_R = 1)

# run models

# 1. model with all interactions
model_stan1 <- rstan::stan(file = 'stan_models/MST_model_all_interactions.stan',
                          data = stan_datalist, 
                          init = stan_inits, 
                          iter = 5e3, 
                          chains = 3)

# 2. model with E*long-term warming removed
model_stan2 <- rstan::stan(file = 'stan_models/MST_model_bTc_alpha_interaction.stan',
                           data = stan_datalist, 
                           init = stan_inits, 
                           iter = 5e3, 
                           chains = 3)

# 3. model with a*long-term warming removed
model_stan3 <- rstan::stan(file = 'stan_models/MST_model_bTc_interaction.stan',
                           data = stan_datalist, 
                           init = stan_inits, 
                           iter = 5e3, 
                           chains = 3)

# 4. model with all interactions removed
model_stan4 <- rstan::stan(file = 'stan_models/MST_model_no_interactions.stan',
                           data = stan_datalist, 
                           init = stan_inits, 
                           iter = 5e3, 
                           chains = 3)

# save stan fits ####
stan_fits <- list(model1 = model_stan1,
                  model2 = model_stan2,
                  model3 = model_stan3,
                  model4 = model_stan4)

saveRDS(stan_fits, 'processed_data/stan_output_GP.rds')
