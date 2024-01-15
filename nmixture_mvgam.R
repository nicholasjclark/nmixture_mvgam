# Simulate data from a Poisson-Binomial N-Mixture model
library(mvgam)
myplot = function(...){
  plot(..., bty = 'l', pch = 16,
       col = 'gray30', cex = 0.8)
  box(bty = 'l', lwd = 2)
}

# True abundance is predicted by a single nonlinear function of temperature
# as well as a nonlinear long-term Gaussian Process trend
set.seed(999)
gamdat <- gamSim(n = 80); N <- NROW(gamdat)
abund_linpred <- gamdat$y; temperature <- gamdat$x2
myplot(abund_linpred ~ temperature)
trend <-  mvgam:::sim_gp(rnorm(3, 0, 0.5),
                        alpha_gp = 5,
                        rho_gp = 14, h = N)
true_abund <- floor(12 + abund_linpred + trend)

# Detection probability changes non-linearly with rainfall, which
# is a highly seasonal variable
rainfall <- sin(2 * pi * (1:N) / 12) +
  0.5 * cos(2 * pi * (1:N) / 12) +
  mvgam:::sim_gp(rnorm(3, 0, 0.1),
                 alpha_gp = 0.1,
                 rho_gp = 25, h = N)
detect_linpred <- 0.7 + (1.25 * rainfall) +
  (-0.90 * rainfall^2) + (-0.60 * rainfall^3)
detect_prob <- plogis(detect_linpred)
myplot(detect_prob ~ rainfall,
       ylim = c(0, 1))

# Simulate observed counts
obs_abund <- rbinom(N, size = true_abund, prob = detect_prob)

# Plot true and observed time series
plot(true_abund,
     type = 'l',
     ylab = 'Abundance',
     xlab = 'Time',
     ylim = c(0, max(true_abund)),
     bty = 'l',
     lwd = 2)
box(bty = 'l', lwd = 2)
lines(obs_abund, col = 'darkred', lwd = 2)
title(main = 'True = black; observed = red')

# Gather data into a dataframe suitable for mvgam modelling;
# This will require a 'cap' variable specifying the maximum K to marginalise
# over when estimating latent abundance (it does NOT have to be a fixed value)
model_dat <- data.frame(obs_abund,
                        temperature,
                        rainfall,
                        cap = max(obs_abund) + 20,
                        time = 1:N,
                        series = as.factor('series1'))

# Training and testing folds
data_train <- model_dat %>% dplyr::filter(time <= 75)
data_test <- model_dat %>% dplyr::filter(time > 75)

# Fit a model with informative priors on the two intercept parameters
# and on the length scale of the GP temporal trend parameter
# Note that the 'trend_formula' applies to the latent count process
# (a Poisson process with log-link), while the 'formula' applies to the
# detection probability (logit link)
mod <- mvgam(formula = obs_abund ~ s(rainfall, k = 6),
             trend_formula = ~ s(temperature, k = 6) +
               gp(time, k = 12, c = 5/4, scale = FALSE),
             family = nmix(),
             data = data_train,
             newdata = data_test,
             priors = c(prior(std_normal(), class = '(Intercept)'),
                        prior(normal(2, 2), class = '(Intercept)_trend'),
                        prior(normal(12, 3), class = 'rho_gp_trend(time)')))

# Model summary and diagnostics
summary(mod)
plot(mod, type = 'residuals')

# Intercept parameters
mcmc_plot(mod,
          variable = "Intercept",
          regex = TRUE,
          type = 'hist')

# Hindcasts and forecasts of latent abundance (with truth overlain)
fc <- forecast(mod, type = 'latent_N')
plot(fc); points(true_abund, pch = 16, cex = 0.8)

# Latent abundance predictions are restricted based on 'cap'
max(model_dat$cap); range(fc$forecasts[[1]])

# Hindcasts and forecasts of detection probability (with truth overlain)
fc <- forecast(mod, type = 'detection')
plot(fc); points(detect_prob, pch = 16, cex = 0.8)

# Hindcasts and forecasts of observations
# (after accounting for detection error)
fc <- forecast(mod, type = 'response')
plot(fc)

# Hindcasts and forecasts of response expectations
# (with truth overlain)
fc <- forecast(object = mod, type = 'expected')
plot(fc); points(detect_prob * true_abund, pch = 16, cex = 0.8)

# Plot conditional effects
library(ggplot2)

# Effects on true abundance can be visualised using type = 'link'
abund_plots <- plot(conditional_effects(mod,
                    type = 'link',
                    effects = c('temperature', 'time')),
                    plot = FALSE)

# Effect of temperature on abundance
abund_plots[[1]] +
  ylab('Latent abundance')
myplot(true_abund ~ temperature)

# Long-term trend in abundance
abund_plots[[2]] +
  ylab('Latent abundance')
myplot(trend)

# Effect of rainfall on detection probability
det_plot <- plot(conditional_effects(mod,
                    type = 'detection',
                    effects = 'rainfall'),
                    plot = FALSE)
det_plot[[1]] +
  ylab('Pr(detection)')
myplot(detect_prob ~ rainfall)

# More targeted plots can use marginaleffects capabilities;
# Here visualise how response predictions might change
# if we considered different possible 'cap' limits on latent
# abundance and different environmental measurements
plot_predictions(mod, condition = list('temperature',
                                       cap = c(15, 50),
                                       rainfall = c(-1.25, 0, 1.25)),
                 type = 'response',
                 conf_level = 0.5) +
  ylab('Observed abundance') +
  theme_classic()

#### A simulation with multiple replicates at each timepoint ####
# 4 timepoints, 3 reps per timepoint, 2 species with different
# temporal trends and detection probabilities
library(dplyr)
set.seed(123)
data.frame(site = 1,
           trend = 1,
           replicate = rep(1:3, 4),
           time = sort(rep(1:4, 3)),
           species = 'sp_1',
           truth = c(rep(28, 3),
                     rep(26, 3),
                     rep(23, 3),
                     rep(19, 3)),
           obs = c(rbinom(3, 28, 0.5),
                   rbinom(3, 26, 0.5),
                   rbinom(3, 23, 0.5),
                   rbinom(3, 19, 0.5))) %>%
  dplyr::mutate(series = paste0('site_', site,
                                '_', species,
                                '_rep_', replicate),
                time = as.numeric(time),
                cap = 35) %>%
  dplyr::select(- replicate) -> testdat

# Add another species
testdat = testdat %>%
  dplyr::bind_rows(data.frame(site = 1,
                              trend = 2,
                              replicate = rep(1:3, 4),
                              time = sort(rep(1:4, 3)),
                              species = 'sp_2',
                              truth = c(rep(4, 3),
                                        rep(7, 3),
                                        rep(15, 3),
                                        rep(8, 3)),
                              obs = c(rbinom(3, 4, 0.8),
                                      rbinom(3, 7, 0.8),
                                      rbinom(3, 15, 0.8),
                                      rbinom(3, 8, 0.8))) %>%
                     dplyr::mutate(series = paste0('site_', site,
                                                   '_', species,
                                                   '_rep_', replicate),
                                   time = as.numeric(time),
                                   cap = 18) %>%
                     dplyr::select(- replicate))

# Factors for species and series IDs
testdat$species <- as.factor(testdat$species)
testdat$series <- factor(testdat$series, 
                         levels = unique(testdat$series))

# The trend_map operates by species here
trend_map <- testdat %>%
  dplyr::select(trend, series) %>%
  dplyr::distinct()
trend_map

# Model allows species to have different detection probabilities
# and different temporal trends
mod <- mvgam(obs ~ species,
             trend_formula = ~ s(time, by = species, k = 4) +
               species,
             trend_model = 'None',
             trend_map = trend_map,
             family = nmix(),
             data = testdat,
             priors = c(prior(std_normal(), class = b),
                        prior(std_normal(), class = Intercept),
                        prior(normal(1, 1.5), class = Intercept_trend)),
             return_model_data = TRUE)

# Summary and other diagnostics work fine
summary(mod)
loo(mod)

# Smooths of time
plot(mod, type = 'smooths', trend_effects = TRUE)

# Detection probs per species
plot_predictions(mod, condition = list('species',
                                       cap = 10),
                 type = 'detection') +
  ylab('Pr(detection)') +
  theme_classic() +
  theme(legend.position = 'none')

# Latent N predictions per species (ignoring the observations)
plot_predictions(mod, 
                 condition = list(time = 1:4, 
                                  species = 'sp_1',
                                  cap = 35),
                 type = 'latent_N')

plot_predictions(mod, 
                 condition = list(time = 1:4, 
                                  species = 'sp_2',
                                  cap = 18),
                 type = 'latent_N')

# Predictions conditional on the observations
ypreds <- mvgam:::mcmc_chains(mod$model_output, 'latent_ypred')
dim(ypreds)

# Pull out the time series of each replicate and plot
plot(1,
     type = 'n',
     xlim = c(1, 4),
     ylim = range(ypreds[,which(testdat$trend == 1)]))
for(i in 1:3){
  for(j in 1:300){
    points(x = jitter(1:4, 0.2),
           ypreds[j,which(as.numeric(testdat$series) == i)],
           col = "#BEBEBE4C", pch = 16, cex = 0.7)
  }
}
points(x = 1:4,
       y = c(28, 26, 23, 19),
       pch = 16, cex = 1.1, col = 'white')
points(x = 1:4,
       y = c(28, 26, 23, 19),
       pch = 16, cex = 0.9,col = 'darkred')


#### Reproducing Jeff Doser's example: ####
# https://www.jeffdoser.com/files/spabundance-web/articles/nmixturemodels
library(mvgam)
load(url('https://github.com/doserjef/spAbundance/raw/main/data/dataNMixSim.rda'))
data.one.sp <- dataNMixSim
data.one.sp$y <- data.one.sp$y[1, , ]
abund.cov <- dataNMixSim$abund.covs[, 1]
abund.factor <- as.factor(dataNMixSim$abund.covs[, 2])
det.cov <- dataNMixSim$det.covs$det.cov.1[,]
det.cov[is.na(det.cov)] <- rnorm(length(which(is.na(det.cov))))
det.cov2 <- dataNMixSim$det.covs$det.cov.2
det.cov2[is.na(det.cov2)] <- rnorm(length(which(is.na(det.cov2))))

mod_data <- do.call(rbind,
                    lapply(1:NROW(data.one.sp$y), function(x){
                      data.frame(y = data.one.sp$y[x,],
                                 abund_cov = abund.cov[x],
                                 abund_fac = abund.factor[x],
                                 det_cov = det.cov[x,],
                                 det_cov2 = det.cov2[x,],
                                 time = 1:NCOL(data.one.sp$y),
                                 site = paste0('site', x))
                    })) %>%
  dplyr::mutate(series = as.factor(paste0(site, '_', time)),
                site = factor(site, levels = unique(site))) %>%
  dplyr::mutate(time = 1,
                cap = max(data.one.sp$y, na.rm = TRUE) + 8)
dplyr::glimpse(mod_data)
head(mod_data)
levels(mod_data$site)

trend_map <- mod_data %>%
  dplyr::mutate(trend = as.numeric(site)) %>%
  dplyr::select(series, trend)
dplyr::glimpse(trend_map)
head(trend_map)

mod <- mvgam(
  # overall effects of covariates on detection probability
  formula = y ~ det_cov + det_cov2,
  # overall effects of the covariates on latent abundance
  trend_formula = ~ abund_cov +
    s(abund_fac, bs = 're'),
  # linking multiple observations to each site
  trend_map = trend_map,
  family = nmix(),
  trend_model = 'None',
  data = mod_data,
  # standard normal priors on key regression params
  priors = c(prior(std_normal(), class = 'b'),
             prior(std_normal(), class = 'Intercept'),
             prior(std_normal(), class = 'Intercept_trend')),
  run_model = TRUE,
  burnin = 400,
  samples = 300)
code(mod)
summary(mod)

# Overall detection probability
avg_predictions(mod, type = 'detection')

library(ggplot2)
# Effect of covariates on latent abundance
abund_plots <- plot(conditional_effects(mod,
                                        type = 'link',
                                        effects = c('abund_cov',
                                                    'abund_fac')),
                    plot = FALSE)
abund_plots[[1]] +
  ylab('Latent abundance')
abund_plots[[2]] +
  ylab('Latent abundance')

det_plots <- plot(conditional_effects(mod,
                                      type = 'detection',
                                      effects = c('det_cov',
                                                  'det_cov2')),
                  plot = FALSE)
det_plots[[1]] +
  ylab('Pr(detection)')
det_plots[[2]] +
  ylab('Pr(detection)')
