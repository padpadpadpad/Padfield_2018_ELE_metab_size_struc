#---------------------------------------------------------------------------------------------------------#
# fitting the MST model to metabolism, temperature and size distribution data using maximum likelihood ####
#---------------------------------------------------------------------------------------------------------#

# This script fits the MST model to the metabolism, temperature and size distribution data using maximum likelihood
# Uses the model described in the theory section of Padfield et al. (2018) Ecology Letters

# clear workspace
mise::mise(vars = TRUE, pkgs = TRUE, console = TRUE, figs = TRUE)

# load in packages ####
library(bbmle)
library(ggplot2)
library(tidyr)
library(dplyr)

# set.seed
set.seed(42)

# load in data and process ready for modelling ####

# load in rate data
pd <- readRDS('raw_data/rates.rds')

# create numeric code for ancestral treatment
pd <- mutate(pd, ancest_code = ifelse(ancestral == 'amb', 0, 1)) %>%
  # arrange by id
  arrange(., id) %>%
  # change temperature to inverse temperature, centred temperature is 18 ºC
  mutate(., ikt = (1/8.62e-05/(273.15 + 18)) - (1/8.62e-05/(K)))

# load in raw count data
cd <- readRDS('raw_data/counts.rds')

# process count data
cd <- cd %>%
  separate(., id, c('pond', 'ancestral', 'treatment'), sep = '_') %>%
  select(., pond, treatment, ancestral, carbon) %>%
  # add column for temperature
  mutate(., temp = ifelse(treatment == 'W', 20, 16)) %>%
  # add column for id
  unite(., id, c(pond, ancestral, treatment, temp), sep = '_') %>%
  # add column for row tally for each group
  group_by(., id) %>%
  mutate(., row = 1:n()) %>%
  data.frame()

# check which biomass sample was missed in the flow cytometry
pd$id[! pd$id %in% unique(cd$id)]

# drop 20_amb_A_16
pd <- filter(pd, id != '20_amb_A_16')

# create vector of first values to check everything has been done properly later on
start_vals <- group_by(cd, id) %>%
  summarise(first_val = carbon[1]) %>%
  data.frame()

# spread and transpose to get into correct format for maximum likelihood model
cd <- spread(cd, id, carbon) %>%
  select(., -row) %>%
  data.frame()
cd <- as.data.frame(t(cd))  

# check is all correct, should all be TRUE
cd[,1] == start_vals$first_val
# If these are all TRUE - rejoice

#----------------------------------#
# Gross Primary Production fits ####
#----------------------------------#

# global model with fits for GP

# full model includes:
# size distribution
# treatment temperature effect
# interactions with ancestry on those things
# ancest code - amb = 0

# bind together GP data and count data
f_data_GP <- cbind(select(pd, ln_GP, df, ikt, ancest_code), cd)

# delete rows with NAs in
f_data_GP <- f_data_GP[!is.na(f_data_GP[,4]),]

# set column positions for all the different data types
f_col = 1
df_col = 2
ikt_col = 3
ancest_col = 4
mass_col = 5:ncol(f_data_GP)

# filter the data for ancestral treatments
f_data_GP_amb <- f_data_GP[f_data_GP[,4] == 0,]
f_data_GP_warm <- f_data_GP[f_data_GP[,4] == 1,]

# write likelihood function
flux_predict <- function(data, ln_r0_eco, a, E, bAncest){
  
  # create size corrected biomass
  size_cor_dist <- log(data[,df_col] * apply(data[,mass_col]^a, 1, sum, na.rm = TRUE))
  
  # return model with effects
  return(ln_r0_eco + size_cor_dist + E*data[,ikt_col])
}

# write function to get deviances from actual values
# delta values are for differences between photosynthesis and respiration
ll_global <- function(ln_r0_eco, a, E, delta_ln_r0_eco, delta_a, delta_E){
  
  # actual values vs predicted values for ambient ponds
  devs_amb <- f_data_GP_amb[, f_col] - flux_predict(f_data_GP_amb, 
                                                   ln_r0_eco, 
                                                   a, 
                                                   E)
  
  # actual values vs predicted values for warm ponds
  devs_warm <- f_data_GP_warm[, f_col] - flux_predict(f_data_GP_warm, 
                                                     ln_r0_eco + delta_ln_r0_eco,
                                                     a + delta_a, 
                                                     E + delta_E)
  # deviations combined
  devs <- c(devs_amb, devs_warm)
  # get standard deviations 
  sigma <- sqrt(var(devs))
  
  return(-sum(dnorm(devs, mean = 0, sd = sigma, log = TRUE)))
  return(sum(probs))
}

# fit global model with all interactions
ll.fit_global <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                              a = 0.75,
                                              E = 1,
                                              delta_ln_r0_eco = 2, 
                                              delta_a = 0,
                                              delta_E = 0),
                      control = list(maxit = 10000),
                      method = 'Nelder-Mead')

# look at model
summary(ll.fit_global)

# remove delta_E - fix it as 0
ll.fit_global2 <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                               a = 0.75,
                                               E = 1,
                                               delta_ln_r0_eco = 2, 
                                               delta_a = 0),
                       control = list(maxit = 10000),
                       fixed = list(delta_E = 0),
                       method = 'Nelder-Mead')

summary(ll.fit_global2)

# remove delta a - fix it as 0
ll.fit_global3 <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                               a = 0.75,
                                               E = 1,
                                               delta_ln_r0_eco = 2),
                       control = list(maxit = 10000),
                       fixed = list(delta_a = 0,
                                    delta_E = 0),
                       method = 'Nelder-Mead')

summary(ll.fit_global3)

# remove delta_ln_r0_eco - fix it as 0
ll.fit_global4 <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                               a = 0.75,
                                               E = 1),
                       control = list(maxit = 10000),
                       fixed = list(delta_ln_r0_eco = 0,
                                    delta_a = 0,
                                    delta_E = 0),
                       method = 'Nelder-Mead')

# check these models
AIC(ll.fit_global,
    ll.fit_global2,
    ll.fit_global3,
    ll.fit_global4)

anova(ll.fit_global,
      ll.fit_global2,
      ll.fit_global3,
      ll.fit_global4)

# none of these are significantly better than the null model with no long-term warming effect
summary(ll.fit_global4)

# fit log likelihood for this best model - ll.fit_global4 - no ancestral effects
ll.fit_final_prof <- profile(ll.fit_global4)

# save models
to_save <- list(models = list(ll.fit_global,
                              ll.fit_global2,
                              ll.fit_global3,
                              ll.fit_global4),
                best_profile = ll.fit_final_prof)

saveRDS(to_save, 'processed_data/mle_output_GP.rds')

#-------------------------------#
# Community Respiration Fits ####
#-------------------------------#

# global model with fits for CR

# full model includes:
# size distribution
# treatment temperature effect
# interactions with ancestry on those things
# ancest code - amb = 0

f_data_R <- cbind(select(pd, ln_R, df, ikt, ancest_code), cd)

# delete rows with NAs in
f_data_R <- f_data_R[!is.na(f_data_R[,4]),]

# set column positions
f_col = 1
df_col = 2
ikt_col = 3
ancest_col = 4
mass_col = 5:ncol(f_data_R)

# filter the data for ancestral treatments
f_data_R_amb <- f_data_R[f_data_R[,4] == 0,]
f_data_R_warm <- f_data_R[f_data_R[,4] == 1,]

# write likelihood function ####
flux_predict <- function(data, ln_r0_eco, a, E, bAncest){
  
  # create size corrected biomass
  size_cor_dist <- log(data[,df_col] * apply(data[,mass_col]^a, 1, sum, na.rm = TRUE))
  
  # return model with effects
  return(ln_r0_eco + size_cor_dist + E*data[,ikt_col])
}

# write function to get deviances from actual values
# delta values are for differences between ambient and warm mesocosms
ll_global <- function(ln_r0_eco, a, E, delta_ln_r0_eco, delta_a, delta_E){
  
  # actual values vs predicted values for respiration
  devs_amb <- f_data_R_amb[, f_col] - flux_predict(f_data_R_amb, 
                                                    ln_r0_eco, 
                                                    a, 
                                                    E)
  # actual values vs predicted values for photosynthesis
  devs_warm <- f_data_R_warm[, f_col] - flux_predict(f_data_R_warm, 
                                                      ln_r0_eco + delta_ln_r0_eco,
                                                      a + delta_a, 
                                                      E + delta_E)
  # deviations combined
  devs <- c(devs_amb, devs_warm)
  # get standard deviations 
  sigma <- sqrt(var(devs))
  
  return(-sum(dnorm(devs, mean = 0, sd = sigma, log = TRUE)))
  return(sum(probs))
}

# fit global model
ll.fit_global <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                              a = 0.75,
                                              E = 1,
                                              delta_ln_r0_eco = 2, 
                                              delta_a = 0,
                                              delta_E = 0),
                      control = list(maxit = 10000),
                      method = 'Nelder-Mead')

# look at model - fully expect some terms to drop out
summary(ll.fit_global)

# remove delta_E - fixed to 0
ll.fit_global2 <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                               a = 0.75,
                                               E = 1,
                                               delta_ln_r0_eco = 2, 
                                               delta_a = 0),
                       control = list(maxit = 10000),
                       fixed = list(delta_E = 0),
                       method = 'Nelder-Mead')

summary(ll.fit_global2)

# remove delta_a - fixed to 0
ll.fit_global3 <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                               a = 0.75,
                                               E = 1,
                                               delta_ln_r0_eco = 2),
                       control = list(maxit = 10000),
                       fixed = list(delta_a = 0,
                                    delta_E = 0),
                       method = 'Nelder-Mead')

summary(ll.fit_global3)

# remove delta_ln_r0_eco - fixed to 0
ll.fit_global4 <- mle2(ll_global, start = list(ln_r0_eco = -4, 
                                               a = 0.75,
                                               E = 1),
                       control = list(maxit = 10000),
                       fixed = list(delta_ln_r0_eco = 0,
                                    delta_a = 0,
                                    delta_E = 0),
                       method = 'Nelder-Mead')

# check these models
AIC(ll.fit_global,
    ll.fit_global2,
    ll.fit_global3,
    ll.fit_global4)

anova(ll.fit_global,
      ll.fit_global2,
      ll.fit_global3,
      ll.fit_global4)

# Best model is ll.fit_global4
ll.fit_final_prof <- profile(ll.fit_global4)

# to save
to_save <- list(models = list(ll.fit_global,
                              ll.fit_global2,
                              ll.fit_global3,
                              ll.fit_global4),
                best_profile = ll.fit_final_prof)

saveRDS(to_save, 'processed_data/mle_output_R.rds')


