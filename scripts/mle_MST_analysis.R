#-------------------------------------------------------#
# analysis of raw metabolism data and MST model fits ####
#-------------------------------------------------------#

# This script analyses the raw metabolism data and also analyses the MST model fits by looking at the relationship between temperature-corrected rate and mass-corrected biomass
# Recreates Figure 2 and Figure 3 of Padfield et al. (2018) Ecology Letters

# clear workspace
mise::mise(vars = TRUE, pkgs = TRUE, console = TRUE, figs = TRUE)

# load in packages ####
library(smatr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(ggplot2)
library(bbmle)
library(purrr)

# -----------------------------------#
# analysis of raw metabolism data ####
#------------------------------------#

# load data and process
d_rates <- readRDS('raw_data/rates.rds') 

# linear mixed model for GPP
GP_fit <- lmer(ln_GP ~ treatment*ancestral + (1|pond), d_rates)
GP_fit2 <- lmer(ln_GP ~ treatment + ancestral + (1|pond), d_rates)
GP_fit3 <- lmer(ln_GP ~ ancestral + (1|pond), d_rates)
GP_fit4 <- lmer(ln_GP ~ treatment + (1|pond), d_rates)
GP_fit5 <- lmer(ln_GP ~ 1 + (1|pond), d_rates)

# compare models
MuMIn::AICc(GP_fit,
            GP_fit2,
            GP_fit3,
            GP_fit4,
            GP_fit5)

anova(GP_fit3, GP_fit5)

# linear mixed model for CR
R_fit <- lmer(ln_R ~ treatment*ancestral + (1|pond), d_rates)
R_fit2 <- lmer(ln_R ~ treatment + ancestral + (1|pond), d_rates)
R_fit3 <- lmer(ln_R ~ treatment + (1|pond), d_rates)
R_fit4 <- lmer(ln_R ~ ancestral + (1|pond), d_rates)
R_fit5 <- lmer(ln_R ~ 1 + (1|pond), d_rates)

# compare models
anova(R_fit3, R_fit5)

MuMIn::AICc(R_fit,
            R_fit2,
            R_fit3,
            R_fit4,
            R_fit5)

# make plot of raw metabolism data ####
P_plot_raw <- ggplot(d_rates, aes(treatment, ln_GP, fill = ancestral, col = ancestral)) +
  geom_boxplot(aes(fill = ancestral, col = ancestral), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.6)) +
  stat_summary(position = position_dodge(width = 0.6), geom = 'crossbar', fatten = 0, color = 'white', width = 0.45, fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(treatment, ln_GP, col = ancestral), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.2), size = 1.5) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  ggtitle('(a) Gross primary production') +
  ylab(expression(ln~Community~flux~(µmol~O[2]~L^-1~hr^-1))) +
  xlab(' ') +
  scale_x_discrete(labels=c(expression(atop(Ambient,Incubator)), expression(atop(Warmed,Incubator)))) +
  theme(axis.text.x = element_text(size = 12))

R_plot_raw <- ggplot(d_rates, aes(treatment, ln_R, fill = ancestral, col = ancestral)) +
  geom_boxplot(aes(fill = ancestral, col = ancestral), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.6)) +
  stat_summary(position = position_dodge(width = 0.6), geom = 'crossbar', fatten = 0, color = 'white', width = 0.45, fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(treatment, ln_R, col = ancestral), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.2), size = 1.5) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  ggtitle('(b) Community respiration') +
  xlab('') +
  ylab(expression(ln~Community~flux~(µmol~O[2]~L^-1~hr^-1))) +
  scale_x_discrete(labels=c(expression(atop(Ambient,Incubator)), expression(atop(Warmed,Incubator)))) +
  theme(axis.text.x = element_text(size = 12))

raw_plot <- gridExtra::grid.arrange(P_plot_raw + theme(legend.position = 'none') + ylim(3, 7.25), R_plot_raw + theme(legend.position = 'none') + ylim(3, 7.25), ncol = 2)

ggsave('plots/Figure_2.pdf', plot = raw_plot, width = 10, height = 5)

#------------------------------------------------------------#
# analysis of output of MST maximum likelihood model fits ####
#------------------------------------------------------------#

# load in GP model output
GP_mle <- readRDS('processed_data/mle_output_GP.rds')

# unload output
GPmodel1 <- GP_mle$models[[1]]
GPmodel2 <- GP_mle$models[[2]]
GPmodel3 <- GP_mle$models[[3]]
GPmodel4 <- GP_mle$models[[4]]

GP_ll <- GP_mle$best_profile

# model comparisons
AIC(GPmodel1, GPmodel2, GPmodel3, GPmodel4)
anova(GPmodel1, GPmodel2, GPmodel3, GPmodel4)

# best model is model 4
summary(GPmodel4)
confint(GP_ll)

# ln_r0_eco_GP = -3.46 95% CI: -6.27 - -1.02
# a_GP = 0.88 95% CI: 0.57 - 1.17
# E_GP = 0.61 eV 95% CI: 0.11 - 1.12

# load in R model output
R_mle <- readRDS('processed_data/mle_output_R.rds')

# unload output
Rmodel1 <- R_mle$models[[1]]
Rmodel2 <- R_mle$models[[2]]
Rmodel3 <- R_mle$models[[3]]
Rmodel4 <- R_mle$models[[4]]

R_ll <- R_mle$best_profile

# model comparisons
AIC(Rmodel1, Rmodel2, Rmodel3, Rmodel4)
anova(Rmodel1, Rmodel2, Rmodel3, Rmodel4)
anova(Rmodel2, Rmodel3)

# best model is model 4 - was refit with model 5
summary(Rmodel4)
confint(R_ll)

# ln_r0_eco_R = - 5.50 95% CI: -10.1 - -2.22
# a_R = 0.8 95% CI: 0.31 - 1.18
# E_R = 1.27 eV 95% CI: 0.70 - 1.83

# create mass corrected biomass estimates ####

# load in raw count data and process
d_count <- readRDS('raw_data/counts.rds') %>%
  separate(., id, c('pond', 'ancestral', 'treatment'), sep = "_") %>%
  select(., pond, treatment, ancestral, carbon) %>%
  mutate(., temp = ifelse(treatment == 'W', 20, 16)) %>%
  unite(., id, c(pond, ancestral, treatment, temp), sep = '_') %>%
  group_by(., id) %>%
  # create mass corrected biomass by raising all the carbon values by alpha for GPP (0.88)
  dplyr::summarise(., Mtot_cor_GP = sum(carbon^0.88),
                   Mtot_cor_R = sum(carbon^0.8)) %>%
  data.frame()

# create temperature-corrected rate ####

# load in rate data 
d_rates <- readRDS('raw_data/rates.rds')

# merge dataframe
d_rates <- merge(d_rates, d_count, by = 'id')

# process dataframe and calculate temperature corrected rate
# uses parameter estimates from the maximum likelihood approach - see values above
# ln_Mtot_cor = mass corrected biomass
# ln_GP_cor = temperature corrected GP
# ln_R_cor = temperature corrected R

d_rates <- mutate(d_rates,
                Mtot_cor_GP = Mtot_cor_GP*df,
                Mtot_cor_R = Mtot_cor_R*df,
                ln_Mtot_cor_GP = log(Mtot_cor_GP*exp(-3.46)),
                ln_Mtot_cor_R = log(Mtot_cor_R*exp(-5.5)),
                GP_pred = log(exp(0.61*(1/(8.62e-5*(18+273.15)) - (1/(8.62e-5*(temp+273.15)))))),
                ln_GP_cor = ln_GP - GP_pred,
                R_pred = log(exp(1.27*(1/(8.62e-5*(18+273.15)) - (1/(8.62e-5*(temp+273.15)))))),
                ln_R_cor = ln_R - R_pred,
                ikt = (1/8.62e-05/(273.15 + 18)) - (1/8.62e-05/(temp + 273.15)))

# GP analysis ####
d_GP <- select(d_rates, pond, ancestral,ln_GP, ln_GP_cor, ikt, ln_Mtot_cor_GP, treatment)
d_GP_clean <- d_GP[complete.cases(d_GP),]

# fit GP sma
fit_GP_sma <- sma(ln_GP_cor ~ ln_Mtot_cor_GP, d_GP, na.action = na.fail, slope.test = 1, robust = TRUE)

# predictions ####

# bootstrap confidence intervals
# get high resolution of prediction variable
boot_preds <- data.frame(expand.grid(ln_Mtot_cor_GP = seq(min(d_GP$ln_Mtot_cor_GP, na.rm = T), max(d_GP$ln_Mtot_cor_GP, na.rm = T), length.out = 200), stringsAsFactors = FALSE))

# bootstrap data and fits
fit_boots <- d_GP %>%
  modelr::bootstrap(n = 1000, id = 'boot_num') %>%
  group_by(boot_num) %>%
  mutate(., fit = map(strap, ~ sma(ln_GP_cor ~ ln_Mtot_cor_GP, data.frame(.), robust = TRUE))) %>%
  ungroup() %>%
  mutate(., intercept = map_dbl(fit, ~coef(.x)[1]),
         slope = map_dbl(fit, ~coef(.x)[2])) %>%
  select(., -fit) %>%
  group_by(boot_num) %>%
  do(data.frame(fitted = .$intercept + .$slope*boot_preds$ln_Mtot_cor_GP, 
                ln_Mtot_cor_GP = boot_preds$ln_Mtot_cor_GP)) %>%
  ungroup() %>%
  group_by(., ln_Mtot_cor_GP) %>%
  dplyr::summarise(., conf_low = quantile(fitted, 0.025),
                   conf_high = quantile(fitted, 0.975)) %>%
  ungroup()

# line of best fit
GP_preds <- data.frame(expand.grid(ln_Mtot_cor_GP = seq(min(d_GP$ln_Mtot_cor_GP, na.rm = T), max(d_GP_clean$ln_Mtot_cor_GP, na.rm = T), length.out = 40), ancestral = unique(d_GP$ancestral), stringsAsFactors = FALSE)) %>%
  mutate(., preds = coef(fit_GP_sma)[1] + coef(fit_GP_sma)[2]*ln_Mtot_cor_GP)

# final GP plot
GP_plot <- ggplot(d_GP) +
  geom_segment(aes(x = 3.87, y = 3.87, xend = 6.84, yend = 6.84), linetype = 2) +
  geom_point(aes(ln_Mtot_cor_GP, ln_GP_cor, col = ancestral)) +
  geom_line(aes(ln_Mtot_cor_GP, preds), GP_preds) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high, x = ln_Mtot_cor_GP), alpha = 0.1, fit_boots) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  ylim(3.25,7.25) +
  xlim(3.25,7.25) +
  xlab(expression(ln~Mass~corrected~biomass~(µg~C^-1))) +
  ylab(expression(atop(ln~Temperature~corrected, community~flux~(µmol~O[2]~L^-1~hr^-1)))) +
  ggtitle(expression((a)~Gross~primary~production))

# R analysis ####

# select columns
d_R <- select(d_rates, pond, ancestral,ln_R, ln_R_cor, ikt, ln_Mtot_cor_R, treatment)
d_R_clean <- d_R[complete.cases(d_R),]

# fit sma
fit_R_sma <- sma(ln_R_cor ~ ln_Mtot_cor_R, d_R, na.action = na.fail, slope.test = 1, robust = TRUE)

# predictions for R ####

# get high resolution of prediction variable
boot_preds <- data.frame(expand.grid(ln_Mtot_cor_R = seq(min(d_R$ln_Mtot_cor_R, na.rm = T), max(d_R$ln_Mtot_cor_R, na.rm = T), length.out = 200), stringsAsFactors = FALSE))

# bootstrap data and fits
fit_boots <- d_R %>%
  modelr::bootstrap(n = 1000, id = 'boot_num') %>%
  group_by(boot_num) %>%
  mutate(., fit = map(strap, ~ sma(ln_R_cor ~ ln_Mtot_cor_R, data.frame(.), robust = TRUE))) %>%
  ungroup() %>%
  mutate(., intercept = map_dbl(fit, ~coef(.x)[1]),
         slope = map_dbl(fit, ~coef(.x)[2])) %>%
  select(., -fit) %>%
  group_by(boot_num) %>%
  do(data.frame(fitted = .$intercept + .$slope*boot_preds$ln_Mtot_cor_R, 
                ln_Mtot_cor_R = boot_preds$ln_Mtot_cor_R)) %>%
  ungroup() %>%
  group_by(., ln_Mtot_cor_R) %>%
  dplyr::summarise(., conf_low = quantile(fitted, 0.025),
                   conf_high = quantile(fitted, 0.975)) %>%
  ungroup()

R_preds <- data.frame(expand.grid(ln_Mtot_cor_R = seq(min(d_R$ln_Mtot_cor_R, na.rm = T), max(d_R$ln_Mtot_cor_R, na.rm = T), length.out = 40), ancestral = unique(d_R$ancestral), stringsAsFactors = FALSE)) %>%
  mutate(., preds = coef(fit_R_sma)[1] + coef(fit_R_sma)[2]*ln_Mtot_cor_R)

# final R plot
R_plot <- ggplot(d_R) +
  geom_segment(aes(x = 2.51, y = 2.51, xend = 5.47, yend = 5.47), linetype = 2) +
  geom_point(aes(ln_Mtot_cor_R, ln_R_cor, col = ancestral)) +
  geom_line(aes(ln_Mtot_cor_R, preds), R_preds) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high, x = ln_Mtot_cor_R), alpha = 0.1, fit_boots) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  ylim(2.15,5.75) +
  xlim(2.15,5.75) +
  xlab(expression(ln~Mass~corrected~biomass~(µg~C^-1))) +
  ylab(expression(atop(ln~Temperature~corrected, community~flux~(µmol~O[2]~L^-1~hr^-1)))) +
  ggtitle(expression((b)~Community~respiration))

# arrange final plot ####
p1 <- gridExtra::grid.arrange(GP_plot + theme(legend.position = 'none'), R_plot + theme(legend.position = 'none'), ncol = 2)
ggsave('plots/Figure_3.pdf', plot = p1, width = 10, height = 5)