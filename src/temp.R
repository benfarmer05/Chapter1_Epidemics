  
  # Optimization function
  optimize_SIR <- function(data, time, initial_state, lambda.modifier, offset, lower_bounds, upper_bounds, control){
    # Objective function is left unchanged
    objective_function <- function(params, data, time, initial_state){
      list_params <- setNames(params, c("beta", "gamma", "lambda"))
      
      # SIR model output
      SIR.out <- ode(y = c(S = initial_state$S, I = initial_state$I, R = initial_state$R), 
                     times = time, func = SIR, parms = list(
                       b = list_params$beta, 
                       g = list_params$gamma, 
                       N = initial_state$N, 
                       l = list_params$lambda, 
                       C = initial_state$cover
                     ))
      SIR.out <- as.data.frame(SIR.out)
      
      # Simulated infected and removed
      sim.inf <- SIR.out[SIR.out$time %in% data$days, "I"]
      sim.rem <- SIR.out[SIR.out$time %in% data$days, "R"]
      
      # Observed infected and removed
      obs.inf <- data$inftiss
      obs.rem <- data$remtiss
      
      # Scale and rescale removed
      rem.inf.ratio <- max(obs.rem) / max(obs.inf)
      obs.rem <- obs.rem / rem.inf.ratio
      sim.rem <- sim.rem / rem.inf.ratio
      
      # Difference calculation
      diff.inf <- abs(sim.inf - obs.inf)
      diff.rem <- abs(sim.rem - obs.rem)
      
      return(sum(diff.rem))
    }
    
    # Run the optimization
    result <- DEoptim(fn = objective_function, lower = lower_bounds, upper = upper_bounds, 
                      data = data, time = time, initial_state = initial_state, control = control)
    
    # Extract optimized parameters
    optimized_params <- setNames(result$optim$bestmem, c("beta", "gamma", "lambda"))
    return(optimized_params)
  }
  
  # Simulation Function
  simulate_SIR <- function(initial_state, optimized_params, time, lambda.modifier, offset){
    transmission_rate <- optimized_params$beta * (1 / (1 + exp(-optimized_params$lambda * initial_state$cover)) + offset)
    SIR.out <- ode(y = c(S = initial_state$S, I = initial_state$I, R = initial_state$R), 
                   times = time, func = SIR, parms = list(
                     b = transmission_rate, 
                     g = optimized_params$gamma, 
                     N = initial_state$N, 
                     l = optimized_params$lambda, 
                     C = initial_state$cover
                   ))
    return(as.data.frame(SIR.out))
  }
  
  # Main workflow - example for one site
  optimize_and_simulate_site <- function(site_data, lambda.modifier, offset, lower_bounds, upper_bounds, control){
    # Define initial states
    initial_state <- list(
      S = site_data$N - site_data$polyp_SA,
      I = site_data$polyp_SA,
      R = 0,
      N = site_data$N,
      cover = site_data$cover
    )
    
    # Time and observation data
    data <- list(
      days = site_data$days,
      inftiss = site_data$inftiss,
      remtiss = site_data$remtiss
    )
    
    # Optimize parameters
    optimized_params <- optimize_SIR(data = data, time = site_data$time, initial_state = initial_state, 
                                     lambda.modifier = lambda.modifier, offset = offset, 
                                     lower_bounds = lower_bounds, upper_bounds = upper_bounds, control = control)
    
    # Simulate the SIR model with optimized parameters
    simulation_result <- simulate_SIR(initial_state = initial_state, optimized_params = optimized_params, 
                                      time = site_data$time, lambda.modifier = lambda.modifier, offset = offset)
    
    # Calculate R0 and adjusted beta
    beta_adj <- optimized_params$beta * (1 / (1 + exp(-lambda.modifier * initial_state$cover)) + offset)
    R0 <- beta_adj / optimized_params$gamma
    
    return(list(optimized_params = optimized_params, simulation = simulation_result, R0 = R0, beta_adj = beta_adj))
  }
  
  # Example call for Nearshore
  site_data <- list(
    days = summary %>% filter(site == "Nearshore") %>% pull(days.survey),
    polyp_SA = polyp_SA.nearshore[[1]],
    N = summary %>% filter(site == "Nearshore", TP == "T11") %>% pull(tot.sustiss),
    cover = (summary %>% filter(site == "Nearshore", TP == "T11") %>% pull(tot.sustiss)) / 200 * 0.3,
    time = time_list[[1]],
    inftiss = obs %>% filter(Site == "Nearshore", Category == "LS", Compartment == "Infected") %>% pull(tissue) + 
      obs %>% filter(Site == "Nearshore", Category == "MS", Compartment == "Infected") %>% pull(tissue) + 
      obs %>% filter(Site == "Nearshore", Category == "HS", Compartment == "Infected") %>% pull(tissue),
    remtiss = obs %>% filter(Site == "Nearshore", Category == "LS", Compartment == "Dead") %>% pull(tissue) + 
      obs %>% filter(Site == "Nearshore", Category == "MS", Compartment == "Dead") %>% pull(tissue) + 
      obs %>% filter(Site == "Nearshore", Category == "HS", Compartment == "Dead") %>% pull(tissue)
  )
  
  # Run the optimization and simulation for Nearshore
  result_nearshore <- optimize_and_simulate_site(site_data, lambda.modifier = 1.0, offset = 1 - 1 / (1 + exp(-1.0)),
                                                 lower_bounds = c(0, 0, 1.0), upper_bounds = c(4, 4, 1.0), 
                                                 control = list(itermax = 100))
  
  # Access the results
  result_nearshore$optimized_params
  result_nearshore$R0
