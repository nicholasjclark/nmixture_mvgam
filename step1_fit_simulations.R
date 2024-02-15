#### Simulate imperfectly-detected multispecies abundance data ####
# Required libraries:
# insight
# reshape2
# mvnfast
# scales
# extraDistr
# mvgam
# spAbundance
# ggplot2
library(mvgam)
library(ggplot2)

# Source the simulation functions
source('sim_nmix.R')

# Simulate with: 
# 8 years (i.e. timepoints)
# 8 sample sites
# 2 species with correlated environmental responses
# base detection probability = 50%
# base mean latent abundance = 40
# No observations missing
# 4 replicates per year (at medium range detection probability and with base
# lambdas around 30 - 40, we generally may need at least 4 - 6 visits per year
# for estimating abundance: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14054)
simdat <- simulate_nmix(n_sites = 8,
                        n_timepoints = 6,
                        n_species = 2,
                        n_replicates = 4,
                        prop_missing = 0,
                        base_detprob = 0.20,
                        base_lambda = 50)

# Inspect the observed and true count distributions
plot_counts(simdat)

# Inspect the returned modeling data.frame
model_df <- simdat$model_df
dplyr::glimpse(model_df)
length(which(is.na(model_df$obs_count)))

# The simulated rainfall effects on detection probability
plot_rainfuncs(simdat)

# The simulated linear effects of temperature
simdat$abund_params$temps

ggplot(model_df, aes(x = temperature, 
                     y = true_count, 
                     col = species,
                     fill = species)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~site)

# Inspect the trend_map data.frame
trend_map <- simdat$trend_map
dplyr::glimpse(trend_map)

# Number of unique species * site combinations, which form the 
# basis of the latent (possibly dynamic) factors
length(unique(trend_map$trend))

# Add a 'cap' variable to the data.frame to be used for setting the
# upper limit for marginalising over the latent abundance states;
# this needs to be large enough so that we can get a reasonable idea
# of support for different detection / latent abundance combinations
model_df %>%
  dplyr::group_by(species, site) %>%
  dplyr::mutate(max_obs = max(obs_count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cap = max(150, 
                          as.integer(floor(max_obs * 2.5)))) -> model_df

# A quick plot to show the replicate observations over time for 
# some series, with the 'true' counts also shown somehow
# Set default colours to look nicer, and use white outlines
ggplot(model_df %>%
         dplyr::filter(site %in% c('site_1', 'site_2',
                                   'site_3', 'site_4',
                                   'site_5', 'site_6')), 
       aes(x = time, y = true_count, fill = species)) +
  geom_line(colour = 'black', size = 1) +
  geom_point(shape = 21, colour = 'white', fill = 'black',
             size = 2.5) +
  geom_line(aes(x = time, y = cap), colour = 'gray',
            size = 1) +
  geom_jitter(aes(x = time, y = obs_count),
              width = 0.08,
              shape = 21, colour = 'white', size = 2.5) +
  facet_grid(rows = vars(species),
             cols = vars(site), 
             scales = 'free') +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(y = 'Count', x = 'Year')
  
# Model the data using mvgam()
mod <- mvgam(formula = obs_count ~
               # Each species' detection prob can have a different nonlinear 
               # association with rainfall; estimate with a hierarchical
               # GAM to allow for partial pooling of functions
               s(rainfall, k = 4) +
               s(rainfall, species, bs = 'fs', k = 4),
             
             trend_formula = ~ 
               # Hierarchical intercepts for each species*site combination
               s(trend, bs = 're') +
               
               # Effects of temperature for each species
               # (could use s(temperature, by = species, bs = 're') if
               # more than 2 species are simulated)
               species * temperature,
             
             # Independent dynamic processes for each species*site
             # combination to capture unmodelled autocorrelation
             trend_model = AR(),
             
             # Regularize key model parameters using informative priors
             priors = c(prior(normal(-0.5, 0.5), class = Intercept),
                        prior(std_normal(), class = b),
                        prior(normal(20, 10), class = lambda),
                        prior(exponential(1), class = sigma_raw_trend),
                        prior(normal(4, 1.5), class = Intercept_trend),
                        prior(exponential(1), class = sigma)),
             trend_map = trend_map,
             family = nmix(),
             algorithm = 'meanfield',
             return_model_data = TRUE,
             samples = 1000,
             data = model_df)

# Parameter summaries (ignoring spline coefficients)
summary(mod)
mcmc_plot(mod, variable = 'ar1', regex = TRUE)
mcmc_plot(mod, variable = 'sigma', regex = TRUE)


# Modeled nonlinear functions vs the true simulated functions
# (ignore differences in scale here; plot_predictions considers the 
# intercepts while plot_tempfuncs does not)
plot_predictions(mod, condition = list('species'),
                 type = 'link',
                 conf_level = 0.2) +
  ylab('Expected latent abundance (N)')

plot_predictions(mod, condition = list('temperature',
                                    'species'),
                 type = 'link',
                 conf_level = 0.2) +
  ylab('Expected latent abundance (N)')
mcmc_plot(mod, variable = 'temperature', regex = TRUE)
simdat$abund_params$temps

plot_predictions(mod, condition = c('rainfall',
                                    'species'),
                 type = 'detection') +
  ylab('Pr(detection)') +
  ylim(c(0, 1))
plot_rainfuncs(simdat)

# Plot some of the latent N predictions vs the simulated truths
hc <- hindcast(mod, type = 'latent_N')
plot_latentN(hc, model_df, 
             species = 'sp_1',
             site = 'site_1')
plot_latentN(hc, model_df, 
             species = 'sp_1',
             site = 'site_2')
plot_latentN(hc, model_df, 
             species = 'sp_1',
             site = 'site_3')
plot_latentN(hc, model_df, 
             species = 'sp_1',
             site = 'site_4')

plot_latentN(hc, model_df, 
             species = 'sp_2',
             site = 'site_1')
plot_latentN(hc, model_df, 
             species = 'sp_2',
             site = 'site_2')
plot_latentN(hc, model_df, 
             species = 'sp_2',
             site = 'site_3')
plot_latentN(hc, model_df, 
             species = 'sp_2',
             site = 'site_4')

# Proposed simulation strategy
# 20 iterations at each of:
# 8 sites
# 2 species
# Normal(50, 5)[25, 100] base lambda per site*species combo
# 2, 3, 4 replicates per time period
# 6, 8, 10 total time periods 
# 0.2, 0.4, 0.6 base p
# mod1 (spline of rainfall), mod2 (rainfall + rainfall^2)
expand.grid(reps = 2:4,
            times = c(6,8,10),
            p = c(0.2,0.4,0.6),
            detcov = c('linear', 'nonlinear'),
            iters = 1:20) %>%
  NROW()
# Want to calculate:
#1. Correlation between each latent_N temporal sample and the true 
#   temporal path (how well does it capture the true trajectory?)

#2. DRPS of latent_N estimates (how well does it estimate true abundance?)

#3. Variogram score of estimated function 1st derivatives for each 
# function sample (z-scored) compared to the true derivative
# (how well does it estimate the overall shape of nonlinear effects?)


