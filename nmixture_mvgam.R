# Simulate data from a Poisson-Binomial N-Mixture model
library(mvgam)

# True abundance is predicted by a single nonlinear function of temperature
# as well as a moderate autoregressive process
set.seed(123)
gamdat <- gamSim(n = 80)
N <- NROW(gamdat)
abund_linpred <- as.vector(scale(gamdat$y))
temperature <- gamdat$x2
true_abund <- rpois(N,
                    exp(1.5 +
                          abund_linpred +
                          mvgam:::sim_ar3(ar1 = 0.7, tau = 100, h = N)))
plot(true_abund ~ temperature)

# Detection probability increases linearly with decreasing rainfall
rainfall <- rnorm(N)
detect_linpred <- 0.5 + -0.95 * rainfall
detect_prob <- plogis(detect_linpred)

# Simulate observed counts
obs_abund <- rbinom(N, size = true_abund, prob = detect_prob)
plot(obs_abund ~ temperature)
plot(obs_abund ~ rainfall)
plot(obs_abund ~ true_abund)

# Gather data into a dataframe suitable for mvgam modelling;
# This will require a 'cap' variable specifying the maximum K to marginalise
# over when estimating latent abundance
model_dat <- data.frame(obs_abund,
                        temperature,
                        rainfall,
                        cap = max(obs_abund) + 20,
                        time = 1:N,
                        series = as.factor('series1'))

# Fit a model with a more useful prior on the rainfall coefficient
mod <- mvgam(formula = obs_abund ~ rainfall,
             trend_formula = ~ s(temperature, k = 7),
             trend_model = AR(),
             family = nmix(),
             data = model_dat,
             priors = prior(std_normal(), class = 'rainfall'))
code(mod)
summary(mod)

plot(mod, type = 'smooths', trend_effects = TRUE)
plot(mod, type = 'residuals')
plot(mod, type = 'trend')
plot(mod, type = 'forecast')

hc <- hindcast(mod, type = 'expected')
plot(hc)

loo(mod)

conditional_effects(mod, type = 'link')
mcmc_plot(mod, variable = 'rainfall', type = 'hist')

# Need to allow different types in marginaleffects
plot_predictions(mod, condition = list('rainfall',
                                       cap = 30),
                 type = 'detection')

newdata <- data.frame(rainfall = runif(100, -2, 2),
                      cap = 30,
                      temperature = mean(model_dat$temperature))
preds <- predict(mod, newdata = newdata, type = 'detection')
plot(preds[1,] ~ newdata$rainfall)
