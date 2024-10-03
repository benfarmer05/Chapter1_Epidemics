
  # .rs.restartR(clean = TRUE)
  rm(list = ls())
  
  # Load required libraries
  library(here)
  library(tidyverse)
  library(DEoptim)
  library(deSolve)
  
  # Import workspace from FLKEYS_data_processing.R
  load(here("output", "data_processing_workspace.RData"))
  
  # Function to refactor site names
  refactor_site_names <- function(data) {
    data %>%
      mutate(
        Site = case_when(
          Site == "Offshore" ~ "off",
          Site == "Midchannel" ~ "mid",
          Site == "Nearshore" ~ "near",
          TRUE ~ Site
        )
      )
  }
  
  # Define the SIR model function
  SIR_model <- function(t, y, params) {
    S = y[1]
    I = y[2]
    R = y[3]
    
    with(as.list(params), {
      transmission_rate = b * (1 / (1 + exp(-l * (C))) + offset)
      
      dS.dt = -transmission_rate * S * I / N
      dI.dt = transmission_rate * S * I / N - g * I
      dR.dt = g * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  # Function to optimize SIR parameters
  optimize_parameters <- function(initial_state, days, coral_data, lambda_modifier) {
    objective_function <- function(params) {
      betas = params[1]
      gammas = params[2]
      lambdas = params[3]
      
      SIR_out <- tryCatch({
        data.frame(ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                       days, SIR_model, c(b = betas, g = gammas, N = initial_state[4], l = lambdas, C = initial_state[5])))
      }, error = function(e) {
        print("Error in ODE:")
        print(e)
        return(NA)
      })
      
      sim_inf = SIR_out[which(SIR_out$time %in% days), "I"]
      sim_rem = SIR_out[which(SIR_out$time %in% days), "R"]
      
      obs_inf = unlist(coral_data[1])
      obs_rem = unlist(coral_data[2])
      
      # Rescale observed removed values
      max_obs_inf = max(obs_inf)
      max_obs_rem = max(obs_rem)
      rem_inf_ratio = max_obs_rem / max_obs_inf
      obs_rem = obs_rem / rem_inf_ratio
      sim_rem = sim_rem / rem_inf_ratio
      
      # Calculate differences and minimize
      diff_inf = sim_inf - obs_inf
      diff_rem = sim_rem - obs_rem
      
      sum_squared_diff = sum(abs(diff_inf)) + sum(abs(diff_rem))
      return(sum_squared_diff)
    }
    
    # Define bounds and control settings for DEoptim
    lower_bounds = c(0, 0, lambda_modifier)
    upper_bounds = c(4, 4, lambda_modifier)
    control = list(itermax = 100)
    
    # Run the optimization
    result = DEoptim(fn = objective_function, lower = lower_bounds, upper = upper_bounds,
                     data = coral_data, time = days, initial_state = initial_state,
                     control = control)
    
    return(result$optim$bestmem)
  }
  
  # Main execution loop for each site
  obs.model <- refactor_site_names(obs)
  
  lambda.modifier = 1.0
  
  
  sites = unique(summary$site)
  my.SIRS = vector('list', length(sites))
  params = vector('list', length(sites))
  
  for (i in seq_along(sites)) {
    site.loop = sites[i]
    
    days = summary %>%
      filter(site == site.loop) %>%
      pull(days.survey) %>%
      na.omit()  # Ensure no NAs in days
    
    first_valid_idx = which(!is.na(days))[1]
    days = days[first_valid_idx:length(days)]
    
    N.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>%
      pull(N.site)
    cover.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>%
      pull(cover.site)
    
    inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
      pull(tissue)    
    
    remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
      pull(tissue)
    
    inftiss = inftiss[first_valid_idx:length(inftiss)]
    remtiss = remtiss[first_valid_idx:length(remtiss)]
    
    # Initial conditions
    I.tiss = inftiss[1]
    S.tiss = N.site - I.tiss
    R.tiss = 0
    
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    # Optimize parameters for the current site
    optimized_params = optimize_parameters(initial_state.tiss, days, coraldata.tiss, lambda.modifier)
    
    # Store optimized parameters in my.SIRS and params
    my.SIRS[[i]] = optimized_params
    params[[i]] = list(beta = optimized_params[1], gamma = optimized_params[2], lambda = optimized_params[3])
    
    # Print optimized parameters for the current site
    print(paste("Site:", site.loop, "Optimized parameters:", paste(optimized_params, collapse = ", ")))
  }
