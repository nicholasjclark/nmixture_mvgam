simulate_nmix = function(n_sites = 5,
                         n_timepoints = 8,
                         n_species = 2,
                         n_replicates = 4,
                         prop_missing = 0.1,
                         base_detprob = 0.6,
                         base_lambda = 8){

  # Required libraries (besides mvgam and spAbund)
  requireNamespace('insight')
  insight::check_if_installed('reshape2')
  insight::check_if_installed('mvnfast')
  insight::check_if_installed('extraDistr')
  
  # Function to simulate from a squared exponential Gaussian Process
  sim_gp = function(N, alpha, rho){
    Sigma <- alpha ^ 2 *
      exp(-0.5 * ((outer(1:N, 1:N, "-") / rho) ^ 2)) +
      diag(1e-9, N)
    mvnfast::rmvn(1,
                  mu = rep(0, N),
                  sigma = Sigma)[1,]
  }
  
  # Function to melt arrays into the appropriate 'long'
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
  
  # Simulate rainfall for every day during the potential sampling window
  N <- 365 * n_timepoints
  rainfall <- matrix(NA, 
                     nrow = n_sites,
                     ncol = N)
  for(i in 1:n_sites){
    trend <- sim_gp(N = N, alpha = 0.5, rho = 15) + 
      mvgam:::sim_ar3(ar1 = 0.5, h = N, tau = 10)
    
    rainfall[i, ] <- as.vector(scale(sin(2 * pi * (1:N) / 365) +
      runif(1, 0.8, 1.1) * cos(2 * pi * (1:N) / 365) + 
      trend)) 
  }

  # Species' base detection probabilities
  sp_det_alphas <- extraDistr::rtnorm(n_species, 
                                      qlogis(base_detprob), 
                                      0.1, 
                                      a = 0.3, 
                                      b = 0.85)
  
  # All species will become more detectable as rainfall
  # reaches a 'happy medium'
  sp_rain_1s <- runif(n_species, 0.90, 1.15)
  sp_rain_2s <- runif(n_species, -0.85, -0.50)
  sp_rain_3s <- runif(n_species, -0.55, -0.25) 
  
  # Detection probabilities at each site and timepoint
  sp_det_probs <- array(NA , dim = c(n_sites, n_species, N))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      sp_det_probs[i, j, ] <-  plogis(sp_det_alphas[j] + 
                                        (sp_rain_1s[j] * rainfall[i, ]) +
                                        (sp_rain_2s[j] * rainfall[i, ] ^ 2) + 
                                        (sp_cubcov_effects[j] * rainfall[i, ] ^ 3))
    }
  }

  # Simulate a long-term moving average temperature pattern at each
  # site that revolves around a primary long-term trend
  main_temp <- sim_gp(N = N, alpha = 0.5, rho = N / 6)
  temperature <- matrix(NA, 
                     nrow = n_sites,
                     ncol = N)
  for(i in 1:n_sites){
    trend <- main_temp + 0.1 * 
      mvgam:::sim_ar3(ar1 = 1, h = N, tau = 100)
    
    temperature[i, ] <- (trend-min(trend))/(max(trend)-min(trend))
  }

  # Simulate species base abundances
  sp_N_alphas <- floor(extraDistr::rtnorm(n_species, 
                                          base_lambda, 
                                          0.5,
                                          a = 3, b = 100))
  
  # All species respond nonlinearly to long-term temperature change
  sp_temp_1s <- runif(n_species, -0.6, 0.6)
  sp_temp_2s <- runif(n_species, 0.2, 0.4)
  sp_temp_4s <- runif(n_species, -0.05, 0.05)

  # Latent abundances at each site and timepoint are made up of
  # nonlinear responses to the moving average temperature variable as 
  # well as a random walk trend component
  sp_abundances <- array(NA , dim = c(n_sites, n_species, N))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      sp_abundances[i, j, ] <- pmax(0, floor(sp_N_alphas[j] +
                                       sp_temp_1s * temperature[i, ] +
                                       sp_temp_2s * temperature[i, ] ^ 2 +
                                       sp_temp_4s * temperature[i, ] ^ 4 +
                                       cumsum(rnorm(N, 0, 0.1))))
    }
  }

  # Now simulate the replicate samples, taken during the 'breeding' season
  # each year (between May and June)
  sample_times <- which(rep(1:365, n_timepoints) %in% c(121, 130, 145, 159))
  sp_observations <- array(NA , dim = c(n_sites, n_species, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      for(t in 1:n_timepoints){
        sp_observations[i, j, t, ] <- rbinom(n_replicates, 
                                             size = sp_abundances[i, j, sample_times[t]], 
                                             prob = sp_det_probs[i, j, sample_times[t]])
      }
    }
  }

  sp_truths <- array(NA , dim = c(n_sites, n_species, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      for(t in 1:n_timepoints){
        sp_truths[i, j, t, ] <- sp_abundances[i, j, sample_times[t]]
      }
    }
  }
  
  # Return necessary objects
  obs_df <- melt_sp_array(sp_observations,
                          'obs_count')
  truth_df <- melt_sp_array(sp_truths,
                            'true_count')
  model_df <- obs_df %>%
    dplyr::left_join(truth_df)
  
  temp_observations <- array(NA , dim = c(n_sites, n_timepoints))
  for(i in 1:n_sites){
    for(t in 1:n_timepoints){
      temp_observations[i, t] <- temperature[i, sample_times[t]]
    }
  }
  
  temp_obs_df <- do.call(rbind, lapply(seq_len(n_sites), function(i){
    data.frame(site = paste0('site_', i),
               time = 1:n_timepoints,
               temperature = temp_observations[i, ])
  }))
  
  rainfall_observations <- array(NA , dim = c(n_sites, n_timepoints))
  for(i in 1:n_sites){
    for(t in 1:n_timepoints){
      rainfall_observations[i, t] <- rainfall[i, sample_times[t]]
    }
  }
  
  rainfall_obs_df <- do.call(rbind, lapply(seq_len(n_sites), function(i){
    data.frame(site = paste0('site_', i),
               time = 1:n_timepoints,
               rainfall = rainfall_observations[i, ])
  }))
  
  model_df %>%
    dplyr::left_join(temp_obs_df) %>%
    dplyr::left_join(rainfall_obs_df) %>%
    dplyr::select(site, 
                  species, 
                  replicate,
                  series,
                  time,
                  temperature,
                  rainfall,
                  obs_count,
                  true_count)-> model_df
  
  return(list(model_df = model_df,
              rainfall_observations = rainfall_observations,
              temp_observations = temp_observations,
              count_observations = sp_observations,
              detection_params = list(alphas = sp_det_alphas,
                                      linrainfalls = sp_rain_1s,
                                      quadrainfalls = sp_rain_2s,
                                      cubrainfalls = sp_rain_3s),
              abund_params = list(lintemps = sp_temp_1s,
                                  quadtemps = sp_temp_2s,
                                  fourthtemps = sp_temp_4s)
              
              ))
}

