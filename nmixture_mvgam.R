# Simulate data from a Poisson-Binomial N-Mixture model
library(mvgam)

# True abundance is predicted by a single nonlinear function of temperature
# as well as a moderate autoregressive process
set.seed(123)
gamdat <- gamSim(n = 100)
N <- NROW(gamdat)
abund_linpred <- as.vector(scale(gamdat$y))
temperature <- gamdat$x2
true_abund <- rpois(N,
                    exp(1.5 +
                          abund_linpred +
                          mvgam:::sim_ar3(ar1 = 0.7, tau = 100)))
plot(true_abund ~ temperature)

# Detection probability increases linearly with decreasing rainfall
rainfall <- rnorm(N)
detect_linpred <- 0.5 + -0.65 * rainfall
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

# Throw in some NAs just to be sure the model runs appropriately
model_dat$obs_abund[c(5, 7, 14)] <- NA
head(model_dat, 8)

# Generate skeleton model code and data that can be modified to accommodate
# the N-mixture process
mod_skeleton <- mvgam(formula = obs_abund ~ rainfall,
                      trend_formula = ~ s(temperature, k = 8),
                      trend_model = AR(),
                      family = poisson(),
                      data = model_dat,
                      run_model = FALSE,
                      autoformat = FALSE)

# Make additions for Poisson-Binomial N-mixture, using code from
# Martijn Bollen's examples
# (https://github.com/MartijnUH/RWsim_abundance_models/blob/main/Stan/NMM_multithreaded.stan)

# 1. Add functions
# 2. Add data lines for the cap, and also add to the model_data
mod_data <- mod_skeleton$model_data
mod_data$cap <- model_dat$cap
mod_data$cap <- as.vector(mod_data$cap)[which(as.vector(mod_data$y_observed) == 1)]

# 3. Add detection probability parameter and model to transformed parameters
# 4. Update likelihood function in model
# 5. Add generated quantities for latent predictions

# Autoformat and compile the model
mod_file <- mvgam:::.autoformat('nmixture_mvgam.stan', overwrite_file = FALSE)
mod_file <- readLines(textConnection(mod_file), n = -1)
cmd_mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(mod_file),
                                   stanc_options = list('O1'))

# Fit the model
fit1 <- cmd_mod$sample(data = mod_data,
                       chains = 4,
                       parallel_chains = 4,
                       refresh = 100,
                       iter_sampling = 500,
                       iter_warmup = 800)

# Convert to mvgam
out_gam_mod <- mvgam:::read_csv_as_stanfit(fit1$output_files())
out_gam_mod <- mvgam:::repair_stanfit(out_gam_mod)
mod_skeleton$model_output <- out_gam_mod
class(mod_skeleton) <- 'mvgam'
plot(mod_skeleton, type = 'smooths', trend_effects = TRUE)
plot(y = abund_linpred, x = temperature)

mcmc_plot(mod_skeleton, variable = 'ar1', regex = TRUE,
          type = 'hist') # true = 0.7
mcmc_plot(mod_skeleton, variable = 'sigma', regex = TRUE,
          type = 'hist') # true = 0.1
mcmc_plot(mod_skeleton, variable = 'rainfall',
          type = 'hist') # true = -0.65
pairs(mod_skeleton, variable = c('(Intercept)', 'rainfall'))
plot(mod_skeleton, type = 'forecast')

# Predictions of true (latent) abundance
Npreds <- as.matrix(mod_skeleton$model_output, 'latent_ypred')
plot(true_abund, pch = 16, col = 'white',
     ylim = range(Npreds))
for(i in 1:500){
  lines(Npreds[i,], col = 'grey80')
}
points(true_abund, pch = 16, col = 'black')  

hc <- hindcast(mod_skeleton, type = 'link')
plot(hc)

#### Notes ####
# Will need a new family function nmix()
# use family$family = 'nmix'
# use family$link = 'log' for brms default priors to work well

# Will need to allow non-dynamic process models (i.e. trend_model == 'None')
# as an option, but only for this particular family

# Will need a different Intercept prior for the observation model, as this
# should relate to average detection probability rather than to the observed counts

# Will have to think about how to calculate residuals 
# (using the Binomial observations?)

# Will need a predict function for new observations
poisbin_pred = function(cap, Xp, betas, p, 
                        truth = NULL, 
                        type = 'link',
                        density = FALSE){
  
  lambdas <- as.vector((matrix(Xp, ncol = NCOL(Xp)) %*%
                      betas) + attr(Xp, 'model.offset'))
  
  if(type ==  'link'){
    # 'link' predictions are expectations of the latent abundance
    out <- lambdas
    
    if(density){
      out <- vector(length = length(truth))
      for(i in seq_along(truth)){
        ks <- truth[i]:cap[i]
        lik_binom <- dbinom(truth, size = ks, p = p[i], log = TRUE)
        lik_poisson <- dpois(x = ks, lambda = lambdas[i], log = TRUE)
        loglik = lik_binom + lik_poisson
        out[i] <- mvgam:::log_sum_exp(loglik)
      }
    }
    
  } else if(type == 'response'){
    xpred <- extraDistr::rtpois(n = length(lambdas), 
                                lambda = lambdas, 
                                b = cap)
    out <- rbinom(length(lambdas), size = xpred, prob = p)
    
  } else if(type == 'variance'){
    xpred <- extraDistr::rtpois(n = length(lambdas), 
                                lambda = lambdas, 
                                b = cap)
    # Variance of a Binomial distribution
    out <- xpred * p * (1 - p)
  } else {
    # Expectations
    xpred <- extraDistr::rtpois(n = length(lambdas), 
                                lambda = lambdas, 
                                b = cap)
    out <- xpred * p
  }
  
  return(out)
}

Xp <- matrix(c(32, 38))
attr(Xp, 'model.offset') <- 0
betas <- 1
poisbin_pred(cap = rep(40, 2),
             Xp = Xp,
             betas = 1,
             p = c(0.7, 0.8),
             type = 'link',
             truth = c(41, 45),
             density = TRUE)

poisbin_pred(cap = rep(40, 2),
             Xp = Xp,
             betas = 1,
             p = c(0.7, 0.8),
             type = 'expected')

poisbin_pred(cap = rep(40, 2),
             Xp = Xp,
             betas = 1,
             p = c(0.7, 0.8),
             type = 'response')


# Will likely also need a predict function for latent 
# true abundance given an observed abundance, following here:
# https://discourse.mc-stan.org/t/calculating-abundance-from-marginalized-n-mixture-model/12751/7
p <- .75
cap <- 50
lambda <- 31
y <- 30
ks <- y:cap

lik_binom <- dbinom(x = y, size = ks, p = p, log = TRUE)
lik_poisson <- dpois(x = ks, lambda = lambda, log = TRUE)
loglik = lik_binom + lik_poisson

# Density for the given observation
mvgam:::log_sum_exp(loglik)

# Predictions
lik <- exp(loglik)
probs <- lik/sum(lik)
N <- sample(x = ks, size = 1, prob = probs)
plot(ks, probs)

