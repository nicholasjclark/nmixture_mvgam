# Function to plot latent abundance estimates vs truth
plot_latentN = function(model_df, trend_map,
                        model, trend = 1){
  ypreds <- mvgam:::mcmc_chains(model$model_output, 'latent_ypred')
  ypreds <- predict(model, newdata = model_df %>%
                      dplyr::arrange(time, series),
                    type = 'latent_N')
  
  truthdat <- model_df %>%
    dplyr::left_join(trend_map) %>%
    dplyr::arrange(time, series) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::filter(trend == !!trend) %>%
    dplyr::mutate(series = droplevels(series))
  
  ids = truthdat %>%
    dplyr::pull(id)
  
  truths <- truthdat %>%
    dplyr::select(time, true_count) %>%
    dplyr::distinct() %>%
    dplyr::arrange(time)
  tot_replicates <- length(unique(truthdat$series))
  
  # Pull out the time series of each replicate and plot
  plot(1,
       type = 'n',
       xlim = c(1, NROW(truths)),
       ylim = c(0, max(c(max(truths$true_count),
                         max(ypreds[,ids])))))
  for(i in 1:tot_replicates){
    for(j in 1:300){
      points(x = jitter(1:6, 0.2),
             ypreds[j,ids[which(as.numeric(truthdat$series) == i)]],
             col = "#BEBEBE4C", pch = 16, cex = 0.7)
    }
  }
  points(x = 1:NROW(truths),
         y = truths$true_count,
         pch = 16, cex = 1.1, col = 'white')
  points(x = 1:NROW(truths),
         y = truths$true_count,
         pch = 16, cex = 0.9,col = 'darkred')
  
}

# Function to construct species' nonlinear responses to
# temperature (which must be scaled to [0,1])
temperature_function = function(a, b, c, temperature){
  (exp(a * temperature +
         b * a * temperature ^ 2 +
         a / 4 * temperature ^ 4) - 
     a / 1.5) * c
}

plot_tempfuncs = function(simdat, species = 1){
  simtemp <- seq(min(simdat$model_df$temperature), 
                 max(simdat$model_df$temperature), 
                 length.out = 500)
  y <- temperature_function(simdat$abund_params$temps1[species],
                            simdat$abund_params$temps2[species],
                            simdat$abund_params$temps3[species],
                            simtemp)
  plot(x = simtemp, y = y,
       type = 'l', lwd = 2,
       xlab = 'temperature',
       ylab = 'F(temperature)',
       bty = 'l')
  box(bty = 'l', lwd = 2)
}

# Function to construct nonlinear responses of 
# detection probability to rainfall variation
rainfall_function = function(a, b, c, d, rainfall){
  plogis(a + (b * rainfall) +
           (c * rainfall ^ 2) + 
           (d * rainfall ^ 3))
}

plot_rainfuncs = function(simdat, species = 1){
  simrain <- seq(min(simdat$model_df$rainfall), 
                 max(simdat$model_df$rainfall), 
                 length.out = 500)
  y <- rainfall_function(simdat$detection_params$alphas[species],
                         simdat$detection_params$linrainfalls[species],
                         simdat$detection_params$quadrainfalls[species],
                         simdat$detection_params$cubrainfalls[species],
                         simrain)
  plot(x = simrain, y = y,
       type = 'l', lwd = 2,
       xlab = 'rainfall',
       ylab = 'F(rainfall)',
       bty = 'l', ylim = c(0, 1))
  box(bty = 'l', lwd = 2)
}

# Function to simulate from a squared exponential Gaussian Process
sim_gp = function(N, alpha, rho){
  Sigma <- alpha ^ 2 *
    exp(-0.5 * ((outer(1:N, 1:N, "-") / rho) ^ 2)) +
    diag(1e-9, N)
  mvnfast::rmvn(1,
                mu = rep(0, N),
                sigma = Sigma)[1,]
}

# Functions to melt arrays into the appropriate 'long'
# dataframe format
melt_sp_array = function(sp_array, value_name = 'obs_count'){
  reshape2::melt(sp_array,
                 varnames = c('site', 'species',
                              'time', 'replicate'),
                 value.name = value_name) %>%
    dplyr::mutate(species = paste0('sp_', species),
                  site = paste0('site_', site)) %>%
    dplyr::mutate(series = paste0(site,
                                  '_', species,
                                  '_rep_', replicate))
}

melt_env_array = function(env_array, value_name = 'temperature'){
  reshape2::melt(env_array,
                 varnames = c('site', 'time', 'replicate'),
                 value.name = value_name) %>%
    dplyr::mutate(site = paste0('site_', site))
}


# Simulation function
simulate_nmix = function(n_sites = 5,
                         n_timepoints = 8,
                         n_species = 2,
                         n_replicates = 4,
                         prop_missing = 0.1,
                         base_detprob = 0.6,
                         base_lambda = 15){

  if(base_lambda < 10L){
    stop('Not recommended to use base abundance < 10')
  }
  if(base_detprob < 0.2){
    stop('Not recommended to use base detection probability < 0.2')
  }
  
  # Required libraries (besides mvgam and spAbundance)
  requireNamespace('insight')
  insight::check_if_installed('reshape2')
  insight::check_if_installed('scales')
  insight::check_if_installed('mvnfast')
  insight::check_if_installed('extraDistr')
  
  # Simulate rainfall for every day during the potential sampling window
  N <- 365 * n_timepoints
  rainfall <- matrix(NA, 
                     nrow = n_sites,
                     ncol = N)
  for(i in 1:n_sites){
    trend <- mvgam:::sim_ar3(ar1 = 0.9,
                             h = N, 
                             tau = 10)
    
    rainfall[i, ] <- as.vector(scale(sin(2 * pi * (1:N) / 365) +
      runif(1, 0.8, 1.1) * cos(2 * pi * (1:N) / 365) + 
      trend)) 
  }

  # Species' base detection probabilities hierarchically
  sp_det_alphas <- extraDistr::rtnorm(n_species, 
                                      qlogis(base_detprob), 
                                      0.1, 
                                      a = 0.3, 
                                      b = 0.85)
  
  # All species will become more detectable as rainfall
  # reaches a 'happy medium'
  sp_rain_1s <- runif(n_species, 0.75, 1.4)
  sp_rain_2s <- runif(n_species, -0.55, -0.10)
  sp_rain_3s <- runif(n_species, -0.35, -0.10) 
  
  # Detection probabilities at each site and timepoint
  sp_det_probs <- array(NA , dim = c(n_sites, n_species, N))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      sp_det_probs[i, j, ] <- rainfall_function(sp_det_alphas[j],
                                                sp_rain_1s[j],
                                                sp_rain_2s[j],
                                                sp_rain_3s[j],
                                                rainfall[i, ])
    }
  }

  # Simulate a long-term moving average temperature pattern at each
  # site that revolves around a primary long-term trend
  set.seed(111)
  main_temp <- sim_gp(N = N, alpha = 1, rho = N / 3)
  temperature <- matrix(NA, 
                     nrow = n_sites,
                     ncol = N)
  for(i in 1:n_sites){
    trend <- main_temp + 0.15 * 
      scale(cumsum(rnorm(N)))
    
    temperature[i, ] <- (trend - min(trend)) /
      (max(trend) - min(trend))
  }

  # Simulate species base abundances hierarchically
  sp_N_alphas <- floor(extraDistr::rtnorm(n_species, 
                                          base_lambda, 
                                          0.5,
                                          a = 5, b = 100))
  
  # All species respond nonlinearly to long-term temperature change
  sp_temp_1s <- runif(n_species, 10, 18)
  sp_temp_2s <- scales::rescale(rnorm(n_species),
                                to = c(-1.4, -1.1))
  sp_temp_3s <- sample(c(-1, 1), n_species, replace = TRUE)
  
  # Latent abundances at each site and timepoint are made up of
  # nonlinear responses to the moving average temperature variable as 
  # well as a random walk trend component
  scale2 = function(x){
    scales::rescale(x, to = c(-2, 2))
  }
  
  sp_abundances <- array(NA , dim = c(n_sites, n_species, N))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      sp_abundances[i, j, ] <- pmax(0, 
                                    floor(sp_N_alphas[j] +
                                            temperature_function(sp_temp_1s[j],
                                                                 sp_temp_2s[j],
                                                                 sp_temp_3s[j],
                                                                 temperature[i, ]) +
                                            scale2(cumsum(rnorm(N)))))
    }
  }

  # Now simulate the replicate samples, spanning the 'breeding' season
  # each year (roughly between May and July)
  sample_times <- array(NA, dim = c(n_timepoints, n_replicates))
  for(t in 1:n_timepoints){
    sample_times[t, ] <- sample(seq.int(118, 180),
                                n_replicates,
                                replace = FALSE) + 365 * (t - 1)
  }

  sp_observations <- array(NA , dim = c(n_sites, n_species, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      for(t in 1:n_timepoints){
        sp_observations[i, j, t, ] <- rbinom(n_replicates, 
                                             size = floor(mean(sp_abundances[i, j, sample_times[t,]])), 
                                             prob = sp_det_probs[i, j, sample_times[t,]])
      }
    }
  }

  # Set missing observations
  if(prop_missing > 0L){
    missing_obs <- sample(1:length(sp_observations),
                          floor(prop_missing * length(sp_observations)),
                          replace = FALSE)
    sp_observations[missing_obs] <- NA
  }

  sp_truths <- array(NA , dim = c(n_sites, n_species, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      for(t in 1:n_timepoints){
        sp_truths[i, j, t, ] <- floor(mean(sp_abundances[i, j, sample_times[t,]]))
      }
    }
  }
  
  #### Return necessary objects ####
  # The count observations as a dataframe
  obs_df <- melt_sp_array(sp_observations,
                          'obs_count')
  truth_df <- melt_sp_array(sp_truths,
                            'true_count')
  model_df <- obs_df %>%
    dplyr::left_join(truth_df)
  
  # The environmental measurements as dataframes
  temp_observations <- array(NA , dim = c(n_sites, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(t in 1:n_timepoints){
      temp_observations[i, t, ] <- mean(temperature[i, sample_times[t,]])
    }
  }
  temp_obs_df <- melt_env_array(temp_observations, 'temperature')
  
  rainfall_observations <- array(NA , dim = c(n_sites, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(t in 1:n_timepoints){
      rainfall_observations[i, t, ] <- rainfall[i, sample_times[t,]]
    }
  }
  rainfall_obs_df <- melt_env_array(env_array = rainfall_observations, 
                                    value_name = 'rainfall')
  
  # Joining environmental measurements to the counts
  model_df %>%
    dplyr::left_join(temp_obs_df) %>%
    dplyr::left_join(rainfall_obs_df) %>%
    dplyr::mutate(year = as.factor(time)) %>%
    dplyr::select(site, 
                  species, 
                  replicate,
                  series,
                  year,
                  time,
                  temperature,
                  rainfall,
                  obs_count,
                  true_count)%>%
    dplyr::mutate(series = as.factor(series),
                  site = as.factor(site),
                  species = as.factor(species)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(true_detect = rainfall_function(sp_det_alphas[as.numeric(species)],
                                    sp_rain_1s[as.numeric(species)],
                                    sp_rain_2s[as.numeric(species)],
                                    sp_rain_3s[as.numeric(species)],
                                    rainfall)) %>%
    dplyr::ungroup() -> model_df
  
  # The trend_map for mvgam modelling
  # Each species * site combination will be modeled as a 
  # unique latent factor to allow for multiple replicate observations
  # in each year
  trend_map <- model_df %>%
    dplyr::select(series, species, site) %>%
    dplyr::mutate(trend = as.factor(paste0(species,
                                           '_',
                                           site))) %>%
    dplyr::select(-species, -site) %>%
    dplyr::mutate(trend = as.numeric(trend)) %>%
    dplyr::distinct()

  # Return as a list
  return(list(model_df = model_df,
              trend_map = trend_map,
              rainfall_observations = rainfall_observations,
              temp_observations = temp_observations,
              count_observations = sp_observations,
              detection_params = list(alphas = sp_det_alphas,
                                      linrainfalls = sp_rain_1s,
                                      quadrainfalls = sp_rain_2s,
                                      cubrainfalls = sp_rain_3s),
              abund_params = list(temps1 = sp_temp_1s,
                                  temps2 = sp_temp_2s,
                                  temps3 = sp_temp_3s)
              
              ))
}

