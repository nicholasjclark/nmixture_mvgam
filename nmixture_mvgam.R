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
plot(trend)
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
