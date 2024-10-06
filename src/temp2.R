  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(DEoptim)
  library(deSolve)
  library(purrr)  # for functional programming
  
  # Import workspace from FLKEYS_data_processing.R
  load(here("output", "plots_obs_workspace.RData"))
  
  # Refactor names in 'obs' to match 'summary'
  obs.model <- obs %>%
    mutate(
      Site = case_when(
        Site == "Offshore" ~ "off",
        Site == "Midchannel" ~ "mid",
        Site == "Nearshore" ~ "near",
        TRUE ~ Site
      )
    )
  
  # ONLY if intending to exclude degree-heating weeks (DHW or DHWs)
  DHW.onset.date = as.Date("2019-07-01")
  DHW.modifier = summary_unique %>%
    filter(date > DHW.onset.date) %>%
    nrow()
  
  # # Scenario 1 [maximum transmission modifier of 1.0, with 100% coral cover]
  #lambda of 3: R0 is extremely low (0.8 or something), no outbreak in Midchannel
  #lambda of 1.6: R0 is solidly below 1.0., at 0.98 for Midchannel
  #lambda of 1.3: R0 is juuust over 1.0, get very small and sustained infection, but essentially no outbreak in Midchannel. R0 < 1 in Offshore but still some infections ??
  #lambda of 1: R0 is low-ish (1.02), get strange and late outbreak in Midchannel
  #lambda of 0.8: R0 is moderate (1.04), again a late and too-strong outbreak in Midchannel. Offshore is late but removal is really close
  #lambda of 0.7: R0 is high (1.05), a late and too-strong outbreak in Midchannel still. offshore looks GREAT, though
  #lambda of 0.5: R0 is high (1.06), a somewhat late and too-strong outbreak in Midchannel
  
  lambda.modifier = 1.0
  offset = 1 - 1 / (1 + exp(-lambda.modifier))
  
  SIR = function(t, y, p){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p), {
      
      # Introduce a non-linear effect of cover on the transmission rate
      # transmission_rate = b * (1 + l * sqrt(C)) #setting lambda to zero nullifies the effect of cover and reverts the model to a basic SIR
      # transmission_rate = b * (l * sqrt(C))
      # transmission_rate = b * (l * (1-exp(-130*(C)))) #20
      transmission_rate = b * (1 / (1 + exp(-l * (C))) + offset)
      
      dS.dt = -transmission_rate * S * I / N
      dI.dt = transmission_rate * S * I / N - g * I
      dR.dt = g * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  sites = unique(summary$site)
  # list(time_list <- vector("list", length = length(sites)))
  
  # Define the function to apply for each site
  process_site <- function(site.loop) {
    
    days <- summary %>%
      # drop_na(days.survey) %>% # NOTE - area to return to after fixing backtracking with patient zero corals. we do want to drop this
      filter(site == site.loop) %>%
      pull(days.inf.site)
    
    # Adjust days to chop off observations after onset of Degree Heating Weeks
    days <- days[1:(length(days) - DHW.modifier)]
    
    # Find the first non-NA index in 'days'
    first_valid_idx <- which(!is.na(days))[1]
    
    # Trim 'days' starting from the first non-NA value
    days.obs = days[first_valid_idx:length(days)]
    days.model = seq(from = min(days.obs), to = max(days.obs), by = 1)
    
    N.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #all values for N.site are the same between categories, so slice first row
      pull(N.site)
    
    cover.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #same as above
      pull(cover.site)
    
    # Sequence of infected & removed SA's for each SusCat within site
    inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    
    remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)
    
    # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
    inftiss <- inftiss[first_valid_idx:length(inftiss)]
    remtiss <- remtiss[first_valid_idx:length(remtiss)]
    
    # Initial conditions
    I.tiss = inftiss[1] #first non-NA & non-zero infection entry
    # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
    S.tiss = N.site - I.tiss
    R.tiss = 0
    
    ############################## optimize parameters ##############################
    # Set up the data and initial conditions
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    objective_function = function(params, data, time, initial_state) {
      betas = params[1]
      gammas = params[2]
      lambdas = params[3]
      
      SIR.out = tryCatch({
        data.frame(ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                       time, SIR, c(b = betas, g = gammas, N = initial_state[4], l = lambdas, C = initial_state[5])))
      }, error = function(e) {
        return(NA)
      })
      
      if (is.na(SIR.out)) return(Inf)
      
      #extract simulated values at time points matching observations
      sim.inf = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I')]
      sim.rem = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
      
      #extract observed values [repetitive code, but works]
      obs.inf = unlist(data[1])
      obs.rem = unlist(data[2])
      
      #use ratio of maximum removed value to maximum infected value to rescale removed. fits the 'I' curve better this way
      max_obs_inf = max(obs.inf)
      max_obs_rem = max(obs.rem)
      rem.inf.ratio = max_obs_rem / max_obs_inf
      obs.rem = obs.rem / rem.inf.ratio
      sim.rem = sim.rem / rem.inf.ratio
      
      # #test to separately rescale removed simulated curve
      # max_sim_inf = max(sim.inf)
      # max_sim_rem = max(sim.rem)
      # rem.inf.ratio = max_sim_rem/max_sim_inf
      # sim.rem = sim.rem/rem.inf.ratio
      
      #minimize using sum of absolute differences
      sum_squared_diff_I = sum(sum(abs(diff.inf))) #can multiply this by 2 or similar to weight it extra
      sum_squared_diff_R = sum(sum(abs(diff.rem)))
      # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
      sum_squared_diff = sum_squared_diff_R
      
      # #minimize using sum of squared differences
      # sum_squared_diff_I = sum(sum(diff.inf^2))
      # sum_squared_diff_R = sum(sum(diff.rem^2))
      # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
      
      return(sum_squared_diff)
    }
    
    # uniform
    lower_bounds.tiss = c(0, 0, lambda.modifier)  #lower bounds for beta, gamma and lambda
    # lower_bounds.tiss = c(0, 0.6, 1)  #lower bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(1, 1, 0.3)  #upper bounds for beta, gamma and lambda
    upper_bounds.tiss = c(4, 4, lambda.modifier)  #upper bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(5, 5, 1)  #upper bounds for beta, gamma and lambda

    control = list(itermax = 100)
    
    # Run the optimization
    result.tiss = DEoptim(fn = objective_function, lower = lower_bounds.tiss, upper = upper_bounds.tiss,
                          data = coraldata.tiss, time = days.model, initial_state = initial_state.tiss,
                          control = control)
    
    # Extract the optimized parameters
    optimized_params.tiss = result.tiss$optim$bestmem
    
    # Print the optimized parameters
    min.beta.tiss = as.numeric(optimized_params.tiss[1])
    min.gamma.tiss = as.numeric(optimized_params.tiss[2])
    min.lambda.tiss = as.numeric(optimized_params.tiss[3])
    min.beta.tiss.adj = min.beta.tiss * (1 / (1 + exp(-lambda.modifier * (cover.site))) + offset)
    # min.beta.tiss.adj = min.beta.tiss * (min.lambda.tiss * (1-exp(-130*(cover.site))))
    # min.beta.tiss.adj = (min.beta.tiss * (1 + min.lambda.tiss * sqrt(cover.site)))
    # min.beta.tiss.adj = (min.beta.tiss * (min.lambda.tiss * sqrt(cover.site)))
    R0 = min.beta.tiss.adj / min.gamma.tiss
    cat("Optimized Tissue Model Parameters for", site.loop, " site:\n")
    cat("Beta:", min.beta.tiss, "\n")
    cat("Cover-adjusted beta:", min.beta.tiss.adj, "\n")
    cat("Gamma:", min.gamma.tiss, "\n")
    cat("Lambda:", min.lambda.tiss, "\n")
    cat("R0:", R0, '\n')
    
    return(c(min.beta.tiss, min.beta.tiss.adj, min.gamma.tiss, min.lambda.tiss, R0, cover.site))
  }
  
  # Apply the process_site function to each site using map
  params <- map(sites, process_site)
  
  # Convert list to dataframe if needed
  params_df <- do.call(rbind, params)
  
  params_df