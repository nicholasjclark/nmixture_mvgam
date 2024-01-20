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

# A custom ggplot2 theme
theme_set(theme_classic(base_size = 12, base_family = 'serif') +
            theme(axis.line.x.bottom = element_line(colour = "black",
                                                    size = 1),
                  axis.line.y.left = element_line(colour = "black",
                                                  size = 1)))
options(ggplot2.discrete.colour = c("#A25050",
                                    "#00008b",
                                    'darkred',
                                    "#010048"),
        ggplot2.discrete.fill = c("#A25050",
                                  "#00008b",
                                  'darkred',
                                  "#010048"))

# Source the simulation functions
source('sim_nmix.R')

# Simulate with: 
# 6 years
# 10 sites
# 2 species
# 3 replicates per year
# base detection probability = 0.6
# base mean latent abundance = 30
# 10% of observations missing
set.seed(1)
simdat <- simulate_nmix(n_sites = 10,
                        n_timepoints = 6,
                        n_species = 2,
                        n_replicates = 3,
                        prop_missing = 0.1,
                        base_lambda = 30)
plot_tempfuncs(simdat, 1)
plot_tempfuncs(simdat, 2)
simdat$abund_params

plot_rainfuncs(simdat, 1)
plot_rainfuncs(simdat, 2)
simdat$detection_params

# Inspect the returned modeling data.frame
model_df <- simdat$model_df
dplyr::glimpse(model_df)
length(which(is.na(model_df$obs_count)))
hist(model_df$obs_count)
hist(model_df$true_count)

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
  dplyr::mutate(cap = as.integer(floor(max_obs * 3))) -> model_df

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
  
# Have our true temp associations come through in
# the simulated data?
mod_temp <- gam(true_count ~ interaction(species, site) +
              s(temperature, by = species, k = 10) +
              s(time, by = interaction(species, site), k = 4),
            family = poisson(),
            data = model_df)
plot_predictions(mod_temp, condition = c('temperature', 
                                         'species',
                                         'species'),
                 type = 'response', points = 0.5) +
  ylab('True abundance') +
  theme_bw() + theme(legend.position = 'none')

# Have true detection associations come through?
mod_rain <- gam(true_detect ~ 
              species + 
             s(rainfall, by = species),
           family = betar(),
           data = model_df)
plot_predictions(mod_rain, condition = c('rainfall', 
                                         'species',
                                         'species'),
                 type = 'response', points = 0.5) +
  theme_bw() + theme(legend.position = 'none') +
  ylim(c(0, 1))

# Try a model
mod <- mvgam(formula = obs_count ~
               # rainfall can impact detection for each 
               # species in nonlinear ways
               species + 
               s(rainfall, by = species, k = 6) - 1,
             
             trend_formula = ~ 
               
               # each species' abundance can have a different nonlinear 
               # association with temperature
               species + 
               s(temperature, by = species, k = 6) + 
               s(time, by = trend, k = 4),
             trend_map = trend_map,
             family = nmix(),
             data = model_df,
             priors = c(prior(std_normal(), class = b),
                        prior(normal(2, 2), class = Intercept_trend)),
             return_model_data = TRUE,
             algorithm = 'meanfield',
             samples = 500)

summary(mod, include_betas = FALSE)
plot(mod, type = 'smooths')
plot(mod, type = 'smooths', trend_effects = TRUE)

# Modeled nonlinear functions vs the true simulated functions
plot_predictions(mod, condition = c('temperature',
                                    'species',
                                    'species'),
                 type = 'link') +
  ylab('Expected abundance') +
  theme_bw() +
  theme(legend.position = 'none')

plot_predictions(mod_temp, condition = c('temperature',
                                         'species',
                                         'species'),
                 type = 'response', points = 0.5) +
  ylab('True abundance') +
  theme_bw() + theme(legend.position = 'none')

plot_predictions(mod, condition = c('rainfall',
                                    'species',
                                    'species'),
                 type = 'detection') +
  ylab('Pr(detection)') +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = 'none')

plot_predictions(mod_rain, condition = c('rainfall', 
                                         'species',
                                         'species'),
                 type = 'response', points = 0.5) +
  theme_bw() + 
  ylim(c(0, 1)) +
  theme(legend.position = 'none')


plot_latentN(model_df,
             trend_map,
             mod,
             trend = 1)
plot_latentN(model_df,
             trend_map,
             mod,
             trend = 2)
plot_latentN(model_df,
             trend_map,
             mod,
             trend = 3)
plot_latentN(model_df,
             trend_map,
             mod,
             trend = 7)
plot_latentN(model_df,
             trend_map,
             mod,
             trend = 8)

