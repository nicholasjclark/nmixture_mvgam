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

# Function to plot observed counts vs truth
plot_counts = function(simdat){
  ggplot(rbind(simdat$model_df %>%
                 dplyr::mutate(count = obs_count,
                               type = 'Observed'),
               simdat$model_df %>%
                 dplyr::mutate(count = true_count,
                               type = 'Truth')), 
         aes(x = as.numeric(count),
             fill = type)) +
    geom_histogram(col = 'white') + labs(x = 'Population count',
                                         y = 'Frequency') +
    facet_wrap(~type, nrow = 2,
               scales = 'free_y') + theme(legend.position = 'none')
}

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
  temperature = temperature - 0.5
  c + (a * temperature +
     b * temperature ^ 2)
}

plot_tempfuncs = function(simdat){
  n_species <- length(unique(simdat$model_df$species))
  simtemp <- seq(min(simdat$model_df$temperature), 
                 max(simdat$model_df$temperature), 
                 length.out = 500)
  sim_funcs <- do.call(rbind, lapply(seq_len(n_species), function(species){
    data.frame(temp = simtemp,
               y = temperature_function(simdat$abund_params$temps1[species],
                                        simdat$abund_params$temps2[species],
                                        simdat$abund_params$temps3[species],
                                        simtemp),
               species = paste0('species_', species))
  }))

  ggplot(sim_funcs, aes(x = temp,
                        y = y,
                        col = species)) + 
    geom_line(size = 1) + 
    labs(x = 'Temperature (normalized)', y = 'F(temperature) on latent N')
}

# Function to construct nonlinear responses of 
# detection probability to rainfall variation
rainfall_function = function(a, b, c, d, rainfall){
  plogis(a + (b * rainfall) +
           (c * rainfall ^ 2) + 
           (d * rainfall ^ 3))
}

plot_rainfuncs = function(simdat){
  n_species <- length(unique(simdat$model_df$species))
  simrain <- seq(min(simdat$model_df$rainfall), 
                 max(simdat$model_df$rainfall), 
                 length.out = 500)
  sim_funcs <- do.call(rbind, lapply(seq_len(n_species), function(species){
    data.frame(rainfall = simrain,
               y = rainfall_function(simdat$detection_params$alphas[species],
                                     simdat$detection_params$linrainfalls[species],
                                     simdat$detection_params$quadrainfalls[species],
                                     simdat$detection_params$cubrainfalls[species],
                                     simrain),
               species = paste0('species_', species))
  }))
  
  ggplot(sim_funcs, aes(x = rainfall,
                        y = y,
                        col = species)) + 
    geom_line(size = 1) + 
    ylim(c(0, 1)) +
    labs(x = 'Rainfall (z-scored)', y = 'Pr(detection)')
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
  if(base_detprob < 0.1){
    stop('Not recommended to use base detection probability < 0.1')
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
                             tau = 5)
    
    rainfall[i, ] <- as.vector(scale(sin(2 * pi * (1:N) / 365) +
      runif(1, 0.8, 1.1) * cos(2 * pi * (1:N) / 365) + 
      trend)) 
  }

  # Species' base detection probabilities
  sp_det_alphas <- rep(qlogis(base_detprob), n_species)
  
  # All species will become more detectable as rainfall
  # reaches a 'happy medium'
  sp_rain_1s <- runif(n_species, 0.5, 1.1)
  sp_rain_2s <- runif(n_species, -0.55, -0.15)
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

  # Now simulate the replicate samples, spanning the 'breeding' season
  # each year (roughly between May and July)
  sample_times <- array(NA, dim = c(n_timepoints, n_replicates))
  for(t in 1:n_timepoints){
    sample_times[t, ] <- sample(seq.int(118, 180),
                                n_replicates,
                                replace = FALSE) + 365 * (t - 1)
  }
  
  # Simulate species base abundances hierarchically
  sp_N_alphas <- matrix(floor(extraDistr::rtnorm(n_species*n_sites, 
                                          base_lambda, 
                                          base_lambda / 10,
                                          a = base_lambda / 2, 
                                          b = base_lambda * 2)),
                        ncol = n_species,
                        nrow = n_sites)
  
  # Simulate moving average temperatures per site, which revolve
  # around a central 'global' temperature trend
  main_temp <- sim_gp(N = n_timepoints, alpha = 1, 
                      rho = n_timepoints / 3)
  
  closed_temp = function(main_temp, n_timepoints, n_replicates){
    trend <- main_temp + 0.65 * 
      scale(cumsum(rnorm(n_timepoints)))
    trend <- scales::rescale(trend, c(-1.85, 1.85))
    out <- do.call(cbind, lapply(seq_len(n_replicates), function(x){ 
      trend
    }))
    return(out)
  }
  
  temp_observations <- array(NA , dim = c(n_sites, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    temp_observations[i, , ] <- closed_temp(main_temp,
                                            n_timepoints,
                                            n_replicates)
  }
  
  # All species respond linearly to long-term temperature change
  sp_temps <- runif(n_species, 0.15, 0.25) *
    sample(c(-1, 1), n_species, TRUE)

  # Latent abundances at each site and timepoint are made up of
  # responses to the moving average temperature variable and
  # site*species level intercept, plus some year to year autocorrelation
  closed_abund = function(alpha, temps, sp_temp, 
                          n_timepoints, n_replicates){
    
    # Species abundances evolve as AR1s plus temperature effects
    ar1 <- runif(1, 0.25, 0.65)
    linpreds <- log(alpha) + temps * sp_temp
    errors <- rnorm(n_timepoints, 0, 0.25)
    trend <- vector(length = n_timepoints)
    trend[1] <- linpreds[1]
    for(t in 2:n_timepoints){
      trend[t] <- ar1 * (errors[t - 1]) +
        linpreds[t] +
        errors[t]
    }
    trend <- ceiling(exp(trend))
    out <- do.call(cbind, lapply(seq_len(n_replicates), function(x){ 
      trend
    }))
    return(out)
  }
  
  sp_truths <- array(NA , dim = c(n_sites, n_species, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      sp_truths[i, j, , ] <- closed_abund(alpha = sp_N_alphas[i, j],
                                          temps = temp_observations[i, , 1],
                                          sp_temp = sp_temps[j],
                                          n_timepoints = n_timepoints,
                                          n_replicates = n_replicates)
    }
  }
  
  # Take imperfect observations, which depend on rainfall
  sp_observations <- array(NA , dim = c(n_sites, n_species, n_timepoints, n_replicates))
  for(i in 1:n_sites){
    for(j in 1:n_species){
      for(t in 1:n_timepoints){
        sp_observations[i, j, t, ] <- rbinom(n_replicates, 
                                             size = sp_truths[i, j, t,], 
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
  
  #### Return necessary objects ####
  # The count observations as a dataframe
  obs_df <- melt_sp_array(sp_observations,
                          'obs_count')
  truth_df <- melt_sp_array(sp_truths,
                            'true_count')
  model_df <- obs_df %>%
    dplyr::left_join(truth_df)
  
  # The environmental measurements as dataframes
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
              abund_params = list(alphas = sp_N_alphas,
                                  temps = sp_temps)
              ))
}

