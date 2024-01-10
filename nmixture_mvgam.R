# Simulate data from a Poisson-Binomial N-Mixture model
library(mvgam)

# True abundance is predicted by a single nonlinear function of temperature
# as well as a nonlinear long-term trend
set.seed(123)
gamdat <- gamSim(n = 80)
N <- NROW(gamdat)
abund_linpred <- gamdat$y
plot(abund_linpred, type = 'l')
temperature <- gamdat$x2
trend <- mvgam:::sim_gp(rnorm(3, 0, 0.1),
                        alpha_gp = 3,
                        rho_gp = 16, h = N)
plot(trend, type = 'l')
true_abund <- floor(10 + abund_linpred + trend)
plot(true_abund ~ temperature)
plot(true_abund, type = 'l')

# Detection probability increases linearly with decreasing rainfall
rainfall <- rnorm(N)
detect_linpred <- 0.4 + -0.55 * rainfall
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
                        cap = max(obs_abund) + 10,
                        time = 1:N,
                        series = as.factor('series1'))

# Add some missing observations
model_dat$obs_abund[c(3,8,12)] <- NA

# Training and testing
data_train <- model_dat %>%
  dplyr::filter(time <= 70)
data_test <- model_dat %>%
  dplyr::filter(time > 70)

# Ensure get_mvgam_priors() allows the trend intercept to be modified
get_mvgam_priors(formula = obs_abund ~ s(rainfall, k = 4),
                 trend_formula = ~ s(temperature, k = 5) +
                   gp(time, k = 10, c = 5/4),
                 trend_model = 'None',
                 family = nmix(),
                 data = data_train)

# Fit a model with a more useful prior on the rainfall coefficient
mod <- mvgam(formula = obs_abund ~ s(rainfall, k = 4),
             trend_formula = ~ s(temperature, k = 5) +
               gp(time, k = 10, c = 5/4),
             trend_model = 'None',
             family = nmix(),
             data = data_train,
             newdata = data_test,
             priors = c(prior(normal(-0.5, 1), class = '(Intercept)'),
                        prior(std_normal(), class = '(Intercept)_trend')),
             chains = 4,
             burnin = 500,
             samples = 350,
             run_model = TRUE)
code(mod)

# Ensure summary works (needs some work for some reason)
summary(mod)

# All standard plots should work
plot(mod, type = 'smooths', trend_effects = TRUE)
plot(mod, type = 'smooths')
plot(mod, type = 'residuals')
plot(mod, type = 'forecast')
plot(mod, type = 'trend')
mcmc_plot(mod,
          variable = "Intercept",
          regex = TRUE,
          type = 'hist')

# link hindcasts are the mean of the latent abundance process
hc <- hindcast(mod, type = 'link')
plot(hc)

# latent_N hindcast is the computed latent abundance
hc <- hindcast(mod, type = 'latent_N')
plot(hc)

# Detection probabilities
hc <- hindcast(mod, type = 'detection')
plot(hc)
plot(hc$hindcasts[[1]][1,] ~ rainfall[1:70])

# Forecasts of latent abundance should be restricted
# based on the 'cap' variable in the data
fc <- forecast(mod, type = 'latent_N')
plot(fc)
points(true_abund, pch = 16, cex = 0.8)
range(fc$forecasts[[1]], na.rm = TRUE)

# Forecasts of detection probability
fc <- forecast(mod, type = 'detection')
plot(fc)
points(detect_prob, pch = 16, cex = 0.8)

####                                                             ####
# TO DO: Need to fit the model without test data and ensure forecasts
# all work #
####                                                             ####

# loo should also work, though need to do some testing to ensure
# it is stable
loo(mod)

# This can be a bit confusing, as the detection probability effects
# are difficult to visualise; it will be better if we can use specific
# types in marginaleffects functions (type = 'detection' and type = 'latent_N')
conditional_effects(mod, type = 'link')
mcmc_plot(mod, variable = 'ar1', regex = TRUE, type = 'hist')

# Need to allow different types in marginaleffects
plot_predictions(mod, condition = list('rainfall',
                                       cap = 30),
                 type = 'detection')

# Easy enough to plot some detection probability predictions
newdata <- data.frame(rainfall = runif(100, min(rainfall), max(rainfall)),
                      cap = 30,
                      temperature = mean(model_dat$temperature),
                      series = 'series1',
                      time = 1)
preds <- predict(mod, newdata = newdata, type = 'detection')
plot(preds[1,] ~ newdata$rainfall, pch = 16, col = 'white',
     ylim = c(0, 1), ylab = 'Pr(detection)',
     xlab = 'rainfall')
for(i in 1:100){
  lines(preds[i,] ~ newdata$rainfall,
        col = 'grey70')
}
lines(y = data.frame(detect_prob, rainfall) %>%
        dplyr::arrange(rainfall) %>%
        dplyr::pull(detect_prob),
      x = data.frame(detect_prob, rainfall) %>%
        dplyr::arrange(rainfall) %>%
        dplyr::pull(rainfall), lwd = 3)
