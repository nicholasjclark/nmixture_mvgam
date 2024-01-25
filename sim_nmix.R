# Function to plot latent abundance estimates vs truth
plot_latentN = function(hindcasts, data, species = 'sp_1',
                        site = 'site_1'){
  all_series <- unique(data %>%
                         dplyr::filter(species == !!species,
                                       site == !!site) %>%
                         dplyr::pull(series))
  
  # Grab the first replicate that represents this series
  # so we can get the true simulated values
  series <- as.numeric(all_series[1])
  truths <- data %>%
    dplyr::arrange(time, series) %>%
    dplyr::filter(series == !!levels(data$series)[series]) %>%
    dplyr::pull(true_count)
  
  # In case some replicates have missing observations,
  # pull out predictions for ALL replicates and average over them
  hcs <- do.call(rbind, lapply(all_series, function(x){
    ind <- which(names(hindcasts$hindcasts) %in% as.character(x))
    hindcasts$hindcasts[[ind]]
  }))
  
  # Calculate posterior empirical quantiles of predictions
  pred_quantiles <- data.frame(t(apply(hcs, 2, function(x) 
    quantile(x, probs = c(0.05, 0.2, 0.3, 0.4, 
                          0.5, 0.6, 0.7, 0.8, 0.95)))))
  pred_quantiles$time <- 1:NROW(pred_quantiles)
  pred_quantiles$truth <- truths
  
  # Grab observations
  data %>%
    dplyr::filter(series %in% all_series) %>%
    dplyr::select(time, obs_count) -> observations
  ggplot(pred_quantiles, aes(x = time, group = 1)) +
    geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "#DCBCBC") + 
    geom_ribbon(aes(ymin = X20., ymax = X80.), fill = "#C79999") +
    geom_ribbon(aes(ymin = X30., ymax = X70.), fill = "#B97C7C") +
    geom_ribbon(aes(ymin = X40., ymax = X60.), fill = "#A25050") +
    geom_line(aes(x = time, y = truth),
               colour = 'black', size = 1) +
    geom_point(aes(x = time, y = truth),
               shape = 21, colour = 'white', fill = 'black',
               size = 2.5) +
    geom_jitter(data = observations, aes(x = time, y = obs_count),
               width = 0.06, 
               shape = 21, fill = 'darkred', colour = 'white', size = 2.5) +
    labs(y = 'Latent abundance (N)',
         x = 'Time',
         title = paste0(species, '; ', site))
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
  sp_N_alphas <- matrix(floor(extraDistr::rtnorm(n_species*n_sites, 
                                          base_lambda, 
                                          base_lambda / 5,
                                          a = 5, b = 100)),
                        ncol = n_species,
                        nrow = n_sites)
  
  # All species respond nonlinearly to long-term temperature change
  sp_temp_1s <- runif(n_species, 10, 18)
  sp_temp_2s <- scales::rescale(rnorm(n_species),
                                to = c(-1.4, -1.1))
  sp_temp_3s <- sample(c(-1, 1), n_species, replace = TRUE)
  
  # Latent abundances at each site and timepoint are made up of
  # nonlinear responses to the moving average temperature variable as 
  # site*species level intercepts and a random walk trend component
  scale2 = function(x, base_lambda){
    
    scales::rescale(x, to = c(-(base_lambda / 10), 
                              (base_lambda / 10)))
  }
  
  sp_abundances <- array(NA , dim = c(n_sites, n_species, N))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      sp_abundances[i, j, ] <- pmax(0, 
                                    floor(sp_N_alphas[i, j] +
                                            temperature_function(sp_temp_1s[j],
                                                                 sp_temp_2s[j],
                                                                 sp_temp_3s[j],
                                                                 temperature[i, ]) +
                                            scale2(cumsum(rnorm(N)), base_lambda)))
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

