  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(DEoptim)
  library(deSolve)
  
  ################################## Set-up ##################################
  
  #import workspace from upstream script
  load(here("output", "plots_obs_workspace.RData"))
  
  #all modeling output in this script is for a single-host SIR model
  curr.host = 'Single-host'
  
  #refactor names in 'obs' to match 'summary'
  obs.model <- obs %>%
    mutate(
      Site = case_when(
        Site == "Offshore" ~ "off",
        Site == "Midchannel" ~ "mid",
        Site == "Nearshore" ~ "near",
        TRUE ~ Site
      )
    )
  
  #parameters for adjusting modeled outbreak based on Degree Heating Weeks (DHWs)
  # NOTE - these do nothing and should be deprecated in favor of the now-updated approach to thermal stress
  DHW.onset.date = as.Date("2019-07-01")
  # DHW.modifier = summary_unique %>%
  #   filter(date > DHW.onset.date) %>%
  #   nrow()
  DHW.modifier = 0 #null model where DHWs do not matter
  
  #NOTE - 10 dec 2024
  #   - should maybe try a DHW threshold of 8.0? I think a problem right now is the phase lag in the second 
  #       infection peak. but this may also be a fundamental flaw in how the second peak is occurring AS
  #       bleaching is occurring. not sure how I would explicitly encourage the modeled infections to 
  #       kick in again after heat stress tapers out....but maybe what I have is actually somewhat realistic too?
  #       lick maybe R0 is really higher as temperature kicks in, causing an uptick in infections...but that 
  #       doesn't seem supported by the observations. not sure. the answer might in part be, next study,
  #       fully integrating temperature in to affect R0 from start to end of outbreak
  #
  #   - something else to think about is if this hints at, since we're only fitting removal, that "infections"
  #       must in some way be kicking up prior to the observed bump in removal? perhaps suggesting that 
  #       the vector itself (maybe symbionts) is surging, precipitating losses...seems a bit contrived. either
  #       way I think explicitly modeling the vector pool some day may make sense
  #
  #   - and another alternative takeaway may be that the system could just be simpler - thermal stress wears off,
  #       disease kicks back into gear. what we're missing is a smart way to test that which informs what is
  #       happening with R0 (possibly fuel for future work)
  
  #NOTE - ideally, all preparation of days and days.model would happen outside of the loops below, making adjustment of fit to include/exclude
  #         days above SST/DHW thresholds easier. for now, this works
  # NOTE - also, date_threshold should be coded dynamically based on the threshold value being used
  SST_threshold_value = 30.5 #the SST on the date that patient-zero SCTLD was backtracked to [could also try 30.5C, thermal stress threshold in corals]
  # DHW_threshold_value = 4.0 #4 is a threshold for coral bleaching in the literature; could try 3 (Whitaker 2024) or 2 (Gierz 2020)
  DHW_threshold_value = NA
  
  #SST approach
  # Find the date when the minimum SST occurs, to split apart 2018 and 2019's excessive heating events
  min_sst_date <- DHW.CRW %>%
    filter(SST.90th_HS == min(SST.90th_HS, na.rm = TRUE)) %>%
    arrange(date) %>% 
    slice(1) %>%  # In case there are multiple minimum values, pick the earliest one
    pull(date)
  
  # Find the date (following the minimum SST date) just before SST exceeds the threshold
  date_threshold <- DHW.CRW %>%
    filter(date > min_sst_date, SST.90th_HS >= SST_threshold_value) %>%
    arrange(date) %>%
    slice(1) %>%
    pull(date) - 1
  
  # Update Error_eval with metrics and thresholds
  error_eval <- error_eval %>%
    mutate(
      SST_threshold = if_else(host == 'Single-host', SST_threshold_value, SST_threshold),
      date_thresh = if_else(host == 'Single-host', date_threshold, date_thresh)
    )
  
  # # Scenario 1 [maximum transmission modifier of 1.0, with 100% coral cover]
  #lambda of 3: R0 is extremely low (0.8 or something), no outbreak in Midchannel
  #lambda of 1.6: R0 is solidly below 1.0., at 0.98 for Midchannel
  #lambda of 1.3: R0 is juuust over 1.0, get very small and sustained infection, but essentially no outbreak in Midchannel. R0 < 1 in Offshore but still some infections ??
  #lambda of 1: R0 is low-ish (1.02), get strange and late outbreak in Midchannel
  #lambda of 0.8: R0 is moderate (1.04), again a late and too-strong outbreak in Midchannel. Offshore is late but removal is really close
  #lambda of 0.7: R0 is high (1.05), a late and too-strong outbreak in Midchannel still. offshore looks GREAT, though
  #lambda of 0.5: R0 is high (1.06), a somewhat late and too-strong outbreak in Midchannel
  
  lambda.modifier = 1.0 # NOTE - 7 feb 2025 - used this for a long time. revert to it if there are issues
  # lambda.modifier = 1.15
  offset = 1 - 1 / (1 + exp(-lambda.modifier))
  
  sites = unique(summary$site)
  
  # Prepare 'days.obs' and 'days.model' based on 'summary' data
  days_sites <- summary %>%
    group_by(site) %>%
    summarize(
      days = list(days.inf.site[1:(length(days.inf.site))]),
      # Extract days.obs and days.model for each site as lists to maintain lengths
      days.obs = list({
        days_cleaned <- na.omit(days.inf.site[1:(length(days.inf.site))])
        days_cleaned[which(!is.na(days_cleaned))[1]:length(days_cleaned)]
      }),
      days.model = list(seq(from = min(na.omit(days.inf.site)), to = max(na.omit(days.inf.site)), by = 1))
    )
  
  # pull SST from summary table
  SST_sites <- summary %>%
    group_by(site) %>%
    filter(date >= first(first.infdate.site)) %>%
    # # Calculate the time variable starting from day 0
    # mutate(time = as.integer(difftime(date, first(first.infdate.site), units = "days"))) %>%
    ungroup() %>%
    select(date, site, days.inf.site, SST.90th_HS) # Adjust the columns as needed
    # rename(SST = SST.90th_HS, time = days.inf.site)
  
  # Create a sequence of all dates between the first and last date for each site
  all_dates <- SST_sites %>%
    group_by(site) %>%
    summarise(start_date = min(date), end_date = max(date)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(date_sequence = list(seq.Date(as.Date(start_date), as.Date(end_date), by = "day"))) %>%
    unnest(date_sequence) %>%
    select(site, date_sequence) %>%
    rename(date = date_sequence)
  
  # use date sequence to pull all SST's between start and end dates of each site's infection observations
  SST_sites <- all_dates %>%
    left_join(SST_sites %>% select(site, date, days.inf.site), by = c("site", "date")) %>%
    left_join(DHW.CRW %>% select(date, SST.90th_HS), by = "date") %>%
    rename(SST = SST.90th_HS) %>%
    group_by(site) %>%
    mutate(time = row_number() - 1) %>%
    ungroup() %>%
    select(site, date, days.inf.site, time, SST)
  
  #for DHW
  DHW_sites <- summary %>%
    group_by(site) %>%
    filter(date >= first(first.infdate.site)) %>%
    # # Calculate the time variable starting from day 0
    # mutate(time = as.integer(difftime(date, first(first.infdate.site), units = "days"))) %>%
    ungroup() %>%
    select(date, site, days.inf.site, DHW_from_90th_HS.1) # Adjust the columns as needed
  # rename(DHW = DHW_from_90th_HS.1, time = days.inf.site)
  
  # Create a sequence of all dates between the first and last date for each site
  all_dates <- DHW_sites %>%
    group_by(site) %>%
    summarise(start_date = min(date), end_date = max(date)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(date_sequence = list(seq.Date(as.Date(start_date), as.Date(end_date), by = "day"))) %>%
    unnest(date_sequence) %>%
    select(site, date_sequence) %>%
    rename(date = date_sequence)
  
  # use date sequence to pull all DHW's between start and end dates of each site's infection observations
  DHW_sites <- all_dates %>%
    left_join(DHW_sites %>% select(site, date, days.inf.site), by = c("site", "date")) %>%
    left_join(DHW.CRW %>% select(date, DHW_from_90th_HS.1), by = "date") %>%
    rename(DHW = DHW_from_90th_HS.1) %>%
    group_by(site) %>%
    mutate(time = row_number() - 1) %>%
    ungroup() %>%
    select(site, date, days.inf.site, time, DHW)
  
  ################################## Model: single-host ##################################
  SIR = function(t,y,p){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p),{
      
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
  
  my.SIRS.basic = vector('list', length(sites))
  params.basic = vector('list', length(sites))
  curr.type = 'Fitted' #the below is for a basic fitting model for single-host transmission (no DHW or projection)
  
  for(i in 1:length(sites)){
    
    site.loop = sites[i]
  
    # site.loop = "mid" #for testing purposes
    # i = 1
    # site.loop = "near" #for testing purposes
    # i = 2
    # site.loop = "off" #for testing purposes
    # i = 3
    
    days <- days_sites %>% # NOTE - make sure this is working right with backtracked patient zero corals
      filter(site == site.loop) %>%
      pull(days) %>%
      unlist()
    
    days.obs <- days_sites %>%
      filter(site == site.loop) %>%
      pull(days.obs) %>%
      unlist() 
    
    days.model <- days_sites %>%
      filter(site == site.loop) %>%
      pull(days.model) %>%
      unlist()
    
    N.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #all values for N.site are the same between categories, so slice first row
      pull(N.site)
    
    cover.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #same as above
      pull(cover.site)
    
    #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
    inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)    
    
    remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)
    
    # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
    first_valid_idx <- which(!is.na(days))[1] #find the first non-NA index
    inftiss <- inftiss[first_valid_idx:length(inftiss)]
    remtiss <- remtiss[first_valid_idx:length(remtiss)]
    
    #initial conditions
    I.tiss = inftiss[1] #first non-NA & non-zero infection entry
    # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
    S.tiss = N.site - I.tiss
    R.tiss = 0
    
    ################################## Optimize single-host model ##################################
    # Set up the data and initial conditions
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    objective_function = function(params, data, time, initial_state){
      
      #testing
      betas = 4 #betas = 3.91
      gammas = 3.12 #gammas = 3.01
      lambdas = 1.0
      initial_state = initial_state.tiss
      time = days.model
      data = coraldata.tiss

      betas = params[1]
      gammas = params[2]
      lambdas = params[3]
      
      SIR.out = tryCatch({
        data.frame(ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                       time, SIR, c(b = betas, g = gammas,
                                    N = initial_state[4],
                                    l = lambdas,
                                    C = initial_state[5])))
      }, error = function(e) {
        print("Error in ODE:")
        print(e)
        return(NA)
      })
      
      #SST approach
      time_cutoff <- SST_sites %>%
        filter(site == site.loop, date >= date_threshold, SST >= SST_threshold_value) %>%
        summarize(first_exceed_time = min(time, na.rm = TRUE)) %>%
        pull(first_exceed_time)
      
      # #DHW approach
      # time_cutoff <- DHW_sites %>%
      #   filter(site == site.loop, date >= date_threshold, DHW >= DHW_threshold_value) %>%
      #   summarize(first_exceed_time = min(time, na.rm = TRUE)) %>%
      #   pull(first_exceed_time)
      
      # Trim days.obs and SIR.out based on the time_cutoff
      days.obs_trimmed <- days.obs[days.obs < time_cutoff]
      
      # Extract simulated values
      sim.inf <- SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'I')]
      sim.rem <- SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'R')]
      sim.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I')]
      sim.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
      
      # Extract observed values
      obs.inf <- unlist(data[[1]])[days.obs < time_cutoff]  # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.rem <- unlist(data[[2]])[days.obs < time_cutoff]  # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.inf.total = unlist(data[[1]])
      obs.rem.total = unlist(data[[2]])
      
      # # NOTE - this was a version where fit was not constrained to the first epidemic wave
      # #extract simulated values at time points matching observations
      # sim.inf = SIR.out[which(SIR.out$time %in% head(days.obs, -1)), which(colnames(SIR.out) %in% 'I')] # NOTE - remove final timepoint for infecteds (because of backtracking, see processing script)
      # sim.rem = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
      # 
      # #extract observed values [repetitive code, but works]
      # obs.inf = unlist(data[[1]])[-length(data[[1]])] # NOTE - remove final timepoint for infecteds (because of backtracking, see processing script)
      # obs.rem = unlist(data[2])
      
      # # NOTE - this as a version where rescaling was required, because infections were part of the fitting process (not *just* removal)
      # #use ratio of maximum removed value to maximum infected value to rescale removed. fits the 'I' curve better this way
      # max_obs_inf = max(obs.inf)
      # max_obs_rem = max(obs.rem)
      # rem.inf.ratio = max_obs_rem/max_obs_inf
      # obs.rem = obs.rem/rem.inf.ratio
      # sim.rem = sim.rem/rem.inf.ratio
      # 
      # # #test to separately rescale removed simulated curve
      # # max_sim_inf = max(sim.inf)
      # # max_sim_rem = max(sim.rem)
      # # rem.inf.ratio = max_sim_rem/max_sim_inf
      # # sim.rem = sim.rem/rem.inf.ratio
      
      # Calculate residuals
      # diff.inf = (sim.inf - obs.inf)
      diff.rem = (sim.rem - obs.rem)
      diff.rem.total = (sim.rem.total - obs.rem.total)
      
      # #aggregate sum of absolute residuals
      # # NOTE - this is not sum of squares and should be clearly stated/defended in the manuscript if used
      # sum_abs_diff_I = sum(sum(abs(diff.inf))) #can multiply this by 2 or similar to weight it extra
      # sum_abs_diff_R = sum(sum(abs(diff.rem)))
      # # sum_diff = sum_abs_diff_I + sum_abs_diff_R #this is the version where infections AND removal are fitted, not *just* removal
      # sum_diff = sum_abs_diff_R
      
      # #decided to just focus on removal for fit
      # #minimize using sum of squared residuals
      # sum_squared_diff_I = sum(diff.inf^2)
      # sum_squared_diff_R = sum(diff.rem^2)
      # # sum_diff = sum_squared_diff_I + sum_squared_diff_R #this is the version where infections AND removal are fitted, not *just* removal
      # sum_diff = sum_squared_diff_R
      
      #minimize using sum of squared residuals
      sum_diff = sum(diff.rem^2)
      sum_diff.total = sum(diff.rem.total^2)
      
      # NOTE - see Kalizhanova et al. 2024 (TB SIR) for other error assessments - including mean absolute error (MAE)
      # Total Sum of Squares (TSS) for removal only
      mean_obs_rem = mean(obs.rem, na.rm = TRUE)  # Compute mean of observed removals
      tss_rem = sum((obs.rem - mean_obs_rem)^2)   # Sum of squared differences from mean
      mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
      tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
      
      #R-squared
      r_squared_rem = 1 - (sum_diff / tss_rem)
      r_squared_rem.total = 1 - (sum_diff.total / tss_rem.total)
      
      error_eval <<- error_eval %>%
        mutate(
          SSR = if_else(site == site.loop & host == curr.host & type == curr.type, 
                        if_else(wave == 'Pre-heat', sum_diff, if_else(wave == 'Both', sum_diff.total, SSR)), SSR),
          TSS = if_else(site == site.loop & host == curr.host & type == curr.type, 
                        if_else(wave == 'Pre-heat', tss_rem, if_else(wave == 'Both', tss_rem.total, TSS)), TSS),
          R_squared = if_else(site == site.loop & host == curr.host & type == curr.type, 
                              if_else(wave == 'Pre-heat', r_squared_rem, if_else(wave == 'Both', r_squared_rem.total, R_squared)), R_squared),
          
          # Update list-columns with vectors
          sim_inf = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Pre-heat', list(sim.inf), if_else(wave == 'Both', list(sim.inf.total), sim_inf)), sim_inf),
          sim_rem = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Pre-heat', list(sim.rem), if_else(wave == 'Both', list(sim.rem.total), sim_rem)), sim_rem),
          obs_inf = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Pre-heat', list(obs.inf), if_else(wave == 'Both', list(obs.inf.total), obs_inf)), obs_inf),
          obs_rem = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Pre-heat', list(obs.rem), if_else(wave == 'Both', list(obs.rem.total), obs_rem)), obs_rem)
        )
      
      return(sum_diff) #return only the residual metric for the epidemic wave being fit to
    }
    
    # uniform
    lower_bounds.tiss = c(0, 0, lambda.modifier)  #lower bounds for beta, gamma and lambda
    # lower_bounds.tiss = c(0, 0.6, 1)  #lower bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(1, 1, 0.3)  #upper bounds for beta, gamma and lambda
    upper_bounds.tiss = c(4, 4, lambda.modifier)  #upper bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(5, 5, 1)  #upper bounds for beta, gamma and lambda
    
    control = list(itermax = 100)  # Maximum number of iterations. 200 is default
    
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
    
    params.basic[[i]] = c(min.beta.tiss, min.beta.tiss.adj, min.gamma.tiss, min.lambda.tiss, R0, cover.site)
    
    #simulation using initial state variables and best-fit beta/gamma parameters
    SIR.out.tiss = data.frame(ode(c(S = initial_state.tiss[1], I = initial_state.tiss[2], R = initial_state.tiss[3]),
                                  days.model, SIR, c(b = min.beta.tiss, g = min.gamma.tiss,
                                               N = initial_state.tiss[4],
                                               l = min.lambda.tiss,
                                               C = initial_state.tiss[5])))
    my.SIRS.basic[[i]] = SIR.out.tiss
  }
  
  ################################## Model: thermal stress w/ single-host ##################################
  SIR.DHW = function(t,y,p,SST,DHW){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p),{
      
      # Find the closest index for t in the SST data and fetch the associated SST and DHW
      closest_index <- which.min(abs(SST$time - t))
      if (closest_index < 1 || closest_index > nrow(SST)) {
        stop("Closest index is out of bounds for SST.")
      }
      current_SST = SST$SST[closest_index]
      current_DHW = DHW$DHW[closest_index]
      
      # Introduce a non-linear effect of cover on the transmission rate [base transmission rate with cover factor]
      # transmission_rate = b * (1 + l * sqrt(C)) #setting lambda to zero nullifies the effect of cover and reverts the model to a basic SIR
      # transmission_rate = b * (l * sqrt(C))
      # transmission_rate = b * (l * (1-exp(-130*(C)))) #20
      transmission_rate = b * (1 / (1 + exp(-l * (C))) + offset)
      removal_rate = g
      
      # SST approach
      # Only apply SST and DHW effects in a way that can simulate a second wave
      thermal_onset_time = min(
        SST$time[SST$SST >= SST_threshold & SST$date >= date_threshold],
        na.rm = TRUE
      )

      if (t >= thermal_onset_time) {

        # SST-based logic for increasing or decreasing rates
        if (current_SST < SST_threshold) {
          # SST below threshold - increase rates
          # transmission_rate <- transmission_rate * (1 + z * (SST_threshold - current_SST))
          # removal_rate <- removal_rate * (1 + e * (SST_threshold - current_SST))
          transmission_rate <- transmission_rate
          removal_rate <- removal_rate
        } else {
          # SST above threshold - decrease rates
          transmission_rate <- transmission_rate * (1 / (1 + exp(-z * (current_SST - SST_threshold))))
          removal_rate <- removal_rate * (1 / (1 + exp(-e * (current_SST - SST_threshold))))
        }
      }
      
      # # DHW approach
      # # Only apply SST and DHW effects in a way that can simulate a second wave
      # thermal_onset_time = min(
      #     DHW$time[DHW$DHW >= DHW_threshold & DHW$date >= date_threshold],
      #     na.rm = TRUE
      #   )
      # 
      # if (t >= thermal_onset_time) {
      # 
      #   # SST-based logic for increasing or decreasing rates
      #   if (current_DHW < DHW_threshold) {
      #     # SST below threshold - increase rates
      #     # transmission_rate <- transmission_rate * (1 + z * (SST_threshold - current_SST))
      #     # removal_rate <- removal_rate * (1 + e * (SST_threshold - current_SST))
      #     transmission_rate <- transmission_rate
      #     removal_rate <- removal_rate
      #   } else {
      #     # SST above threshold - decrease rates
      #     transmission_rate <- transmission_rate * (1 / (1 + exp(-z * (current_DHW - DHW_threshold))))
      #     removal_rate <- removal_rate * (1 / (1 + exp(-e * (current_DHW - DHW_threshold))))
      #   }
      # }
      
      # Ensure rates are non-negative
      transmission_rate <- max(transmission_rate, 0)
      removal_rate <- max(removal_rate, 0)
      
      dS.dt = -transmission_rate * S * I / N 
      dI.dt = transmission_rate * S * I / N - removal_rate * I
      dR.dt = removal_rate * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  my.SIRS.basic.DHW = vector('list', length(sites))
  params.basic.DHW = vector('list', length(sites))
  curr.type = 'DHW' #the below is for a model of single-host transmission, with consideration of thermal stress
  
  for(i in 1:length(sites)){
    
    site.loop = sites[i]
    
    # site.loop = "mid" #for testing purposes
    # i = 1
    # site.loop = "near" #for testing purposes
    # i = 2
    # site.loop = "off" #for testing purposes
    # i = 3
    
    days <- days_sites %>% # NOTE - make sure this is working right with backtracked patient zero corals
      filter(site == site.loop) %>%
      pull(days) %>%
      unlist()
    
    days.obs <- days_sites %>%
      filter(site == site.loop) %>%
      pull(days.obs) %>%
      unlist() 
    
    days.model <- days_sites %>%
      filter(site == site.loop) %>%
      pull(days.model) %>%
      unlist()
    
    SST_df = SST_sites %>%
      filter(site == site.loop) %>%
      select(date, time, SST)
    
    DHW_df <- DHW_sites %>%
      filter(site == site.loop) %>%
      select(date, time, DHW)
    
    SST_values = SST_df %>%
      pull(SST)
    
    DHW_values = DHW_df %>%
      pull(DHW)
    
    # Ensure that the lengths of SST_values and DHW_values match the length of days.model
    #   NOTE - this error checker really should be verifying dates, not the sequence - would likely require comparing to summary
    if (length(SST_values) != length(days.model)) {
      stop("Length of SST_values does not match length of days.model.")
    }
    if (length(DHW_values) != length(days.model)) {
      stop("Length of DHW_values does not match length of days.model.")
    }
    
    N.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #all values for N.site are the same between categories, so slice first row
      pull(N.site)
    
    cover.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #same as above
      pull(cover.site)
    
    #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
    inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)    
    
    remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)
    
    # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
    first_valid_idx <- which(!is.na(days))[1] #find the first non-NA index
    inftiss <- inftiss[first_valid_idx:length(inftiss)]
    remtiss <- remtiss[first_valid_idx:length(remtiss)]
    
    #initial conditions
    I.tiss = inftiss[1] #first non-NA & non-zero infection entry
    # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
    S.tiss = N.site - I.tiss
    R.tiss = 0
    
    ################################## Optimize thermal-stress single-host model ##################################
    # Set up the data and initial conditions
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    pre.fitted.params = params.basic[[i]] # NOTE - relies on hard-coding params.basic & could be revised. works fine now, but check here if bugs
    betas = pre.fitted.params[1]
    gammas = pre.fitted.params[3]
    lambdas = pre.fitted.params[4]
    
    objective_function = function(params, data, time, initial_state){
      
      # #testing
      # betas = 3.91
      # gammas = 3.01
      # lambdas = 1.0
      # zetas = 0.89
      # etas = 0.70
      # initial_state = initial_state.tiss
      # time = days.model
      # data = coraldata.tiss
      
      zetas = params[1]
      etas = params[2]
      # SST_threshold_value = params[3]
      
      SIR.out = tryCatch({
        data.frame(ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                       time, SIR.DHW, c(b = betas, g = gammas,
                                    N = initial_state[4],
                                    z = zetas,
                                    e = etas,
                                    l = lambdas,
                                    C = initial_state[5],
                                    SST_threshold = SST_threshold_value,
                                    DHW_threshold = DHW_threshold_value), # NOTE - hard-coded
                       SST = SST_df,
                       DHW = DHW_df)) # NOTE - this is also hard-coded right now, ideally would not be
      }, error = function(e) {
        print("Error in ODE:")
        print(e)
        return(NA)
      })
      
      #SST approach
      time_cutoff <- SST_sites %>%
        filter(site == site.loop, date >= date_threshold, SST >= SST_threshold_value) %>%
        summarize(first_exceed_time = min(time, na.rm = TRUE)) %>%
        pull(first_exceed_time)
      
      # #DHW approach
      # time_cutoff <- DHW_sites %>%
      #   filter(site == site.loop, date >= date_threshold, DHW >= DHW_threshold_value) %>%
      #   summarize(first_exceed_time = min(time, na.rm = TRUE)) %>%
      #   pull(first_exceed_time)
      
      # Trim days.obs and SIR.out based on the time_cutoff
      #   - this is the reverse of the cutoff-based approach used for the initial fit: this time, only observational data *including and
      #       after* the first date in which the SST (or DHW) threshold is breached, is being fitted to
      days.obs_trimmed <- days.obs[days.obs >= time_cutoff]
      
      # Extract simulated values at matching time points in SIR.out for trimmed days.obs
      sim.inf <- SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'I')]
      sim.rem <- SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'R')]
      sim.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I')]
      sim.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]

      # Extract observed values for the trimmed days.obs
      obs.inf <- unlist(data[[1]])[days.obs >= time_cutoff]  # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.rem <- unlist(data[[2]])[days.obs >= time_cutoff]  # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.inf.total = unlist(data[[1]])
      obs.rem.total = unlist(data[[2]])
      
      # # NOTE - this was a version where fit was not constrained to the first epidemic wave
      # #extract simulated values at time points matching observations
      # sim.inf = SIR.out[which(SIR.out$time %in% head(days.obs, -1)), which(colnames(SIR.out) %in% 'I')] # NOTE - remove final timepoint for infecteds (because of backtracking, see processing script)
      # sim.rem = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
      # 
      # #extract observed values [repetitive code, but works]
      # obs.inf = unlist(data[[1]])[-length(data[[1]])] # NOTE - remove final timepoint for infecteds (because of backtracking, see processing script)
      # obs.rem = unlist(data[2])
      
      # # NOTE - this as a version where rescaling was required, because infections were part of the fitting process (not *just* removal)
      # #use ratio of maximum removed value to maximum infected value to rescale removed. fits the 'I' curve better this way
      # max_obs_inf = max(obs.inf)
      # max_obs_rem = max(obs.rem)
      # rem.inf.ratio = max_obs_rem/max_obs_inf
      # obs.rem = obs.rem/rem.inf.ratio
      # sim.rem = sim.rem/rem.inf.ratio
      # 
      # # #test to separately rescale removed simulated curve
      # # max_sim_inf = max(sim.inf)
      # # max_sim_rem = max(sim.rem)
      # # rem.inf.ratio = max_sim_rem/max_sim_inf
      # # sim.rem = sim.rem/rem.inf.ratio
      
      # Calculate residuals
      # diff.inf = (sim.inf - obs.inf)
      diff.rem = (sim.rem - obs.rem)
      diff.rem.total = (sim.rem.total - obs.rem.total)
      
      # #aggregate sum of absolute differences
      # # NOTE - this is not sum of squares and should be clearly stated/defended in the manuscript
      # sum_abs_diff_I = sum(sum(abs(diff.inf))) #can multiply this by 2 or similar to weight it extra
      # sum_abs_diff_R = sum(sum(abs(diff.rem)))
      # # sum_diff = sum_abs_diff_I + sum_abs_diff_R #this is the version where infections AND removal are fitted, not *just* removal
      # sum_diff = sum_abs_diff_R
      
      # # # decided to just focus on removal for fit
      # #minimize using sum of squared differences
      # sum_squared_diff_I = sum(sum(diff.inf^2))
      # sum_squared_diff_R = sum(sum(diff.rem^2))
      # sum_diff = sum_squared_diff_I + sum_squared_diff_R #this is the version where infections AND removal are fitted, not *just* removal
      # sum_diff = sum_squared_diff_R
      
      #minimize using sum of squared residuals
      sum_diff = sum(diff.rem^2)
      sum_diff.total = sum(diff.rem.total^2)
      
      # NOTE - see Kalizhanova et al. 2024 (TB SIR) for other error assessments - including mean absolute error (MAE)
      # Total Sum of Squares (TSS) for removal only
      mean_obs_rem = mean(obs.rem, na.rm = TRUE)  # Compute mean of observed removals
      tss_rem = sum((obs.rem - mean_obs_rem)^2)   # Sum of squared differences from mean
      mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)  # Compute mean of observed removals
      tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)   # Sum of squared differences from mean
      
      #R-squared
      r_squared_rem = 1 - (sum_diff / tss_rem)
      r_squared_rem.total = 1 - (sum_diff.total / tss_rem.total)
      
      # Update the Error_eval dataframe with current settings using dplyr
      error_eval <<- error_eval %>%
        mutate(
          SSR = if_else(site == site.loop & host == curr.host & type == curr.type, 
                        if_else(wave == 'Post-heat', sum_diff, if_else(wave == 'Both', sum_diff.total, SSR)), SSR),
          TSS = if_else(site == site.loop & host == curr.host & type == curr.type, 
                        if_else(wave == 'Post-heat', tss_rem, if_else(wave == 'Both', tss_rem.total, TSS)), TSS),
          R_squared = if_else(site == site.loop & host == curr.host & type == curr.type, 
                              if_else(wave == 'Post-heat', r_squared_rem, if_else(wave == 'Both', r_squared_rem.total, R_squared)), R_squared),
          
          # Update list-columns with vectors
          sim_inf = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Post-heat', list(sim.inf), if_else(wave == 'Both', list(sim.inf.total), sim_inf)), sim_inf),
          sim_rem = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Post-heat', list(sim.rem), if_else(wave == 'Both', list(sim.rem.total), sim_rem)), sim_rem),
          obs_inf = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Post-heat', list(obs.inf), if_else(wave == 'Both', list(obs.inf.total), obs_inf)), obs_inf),
          obs_rem = if_else(site == site.loop & host == curr.host & type == curr.type, 
                            if_else(wave == 'Post-heat', list(obs.rem), if_else(wave == 'Both', list(obs.rem.total), obs_rem)), obs_rem)
        )
      
      return(sum_diff)
    }

    # lower_bounds.tiss = c(0.01, 0.01)  # Lower bounds for zeta and eta
    # upper_bounds.tiss = c(1.0, 1.0)          # Upper bounds for zeta and eta
    # lower_bounds.tiss = c(0.01, 0.01)  # Lower bounds for zeta and eta
    # upper_bounds.tiss = c(3.0, 3.0)          # Upper bounds for zeta and eta
    lower_bounds.tiss = c(0.0001, 0.0001)  # Lower bounds for zeta and eta
    upper_bounds.tiss = c(3.0, 3.0)          # Upper bounds for zeta and eta
    
    control = list(itermax = 100)  # Maximum number of iterations. 200 is default
    
    # Run the optimization
    result.tiss = DEoptim(fn = objective_function, lower = lower_bounds.tiss, upper = upper_bounds.tiss,
                          data = coraldata.tiss, time = days.model, initial_state = initial_state.tiss,
                          control = control)
    
    # Extract the optimized parameters
    optimized_params.tiss = result.tiss$optim$bestmem
    
    # Print the optimized parameters
    min.zeta.tiss = as.numeric(optimized_params.tiss[1])
    min.eta.tiss = as.numeric(optimized_params.tiss[2])
    
    params.basic.DHW[[i]] = c(min.zeta.tiss, min.eta.tiss)
    
    #simulation using initial state variables and best-fit beta/gamma (and now zeta) parameters
    SIR.out.tiss = data.frame(ode(c(S = initial_state.tiss[1], I = initial_state.tiss[2], R = initial_state.tiss[3]),
                                  days.model, SIR.DHW, c(b = betas, g = gammas,
                                                     N = initial_state.tiss[4],
                                                     z = min.zeta.tiss,
                                                     e = min.eta.tiss,
                                                     l = lambdas,
                                                     C = initial_state.tiss[5],
                                                     SST_threshold = SST_threshold_value,
                                                     DHW_threshold = DHW_threshold_value),
                                  SST = SST_df,
                                  DHW = DHW_df))
    my.SIRS.basic.DHW[[i]] = SIR.out.tiss
  }
  
  ################################## Model: full outbreak ##################################
  my.SIRS.basic.full = vector('list', length(sites))
  params.basic.full = vector('list', length(sites))
  curr.type = 'Fitted' #the below is for a basic fitting model for single-host transmission (no DHW or projection)
  
  for(i in 1:length(sites)){
    
    site.loop = sites[i]
    
    # site.loop = "mid" #for testing purposes
    # i = 1
    # site.loop = "near" #for testing purposes
    # i = 2
    # site.loop = "off" #for testing purposes
    # i = 3
    
    days <- days_sites %>% # NOTE - make sure this is working right with backtracked patient zero corals
      filter(site == site.loop) %>%
      pull(days) %>%
      unlist()
    
    days.obs <- days_sites %>%
      filter(site == site.loop) %>%
      pull(days.obs) %>%
      unlist() 
    
    days.model <- days_sites %>%
      filter(site == site.loop) %>%
      pull(days.model) %>%
      unlist()
    
    N.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #all values for N.site are the same between categories, so slice first row
      pull(N.site)
    
    cover.site = susceptible_ref %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #same as above
      pull(cover.site)
    
    #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
    inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)    
    
    remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)
    
    # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
    first_valid_idx <- which(!is.na(days))[1] #find the first non-NA index
    inftiss <- inftiss[first_valid_idx:length(inftiss)]
    remtiss <- remtiss[first_valid_idx:length(remtiss)]
    
    #initial conditions
    I.tiss = inftiss[1] #first non-NA & non-zero infection entry
    # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
    S.tiss = N.site - I.tiss
    R.tiss = 0
    
    ################################## Optimize full outbreak ##################################
    # Set up the data and initial conditions
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    objective_function = function(params, data, time, initial_state){
      
      # #testing
      # betas = 4 #betas = 3.91
      # gammas = 3.12 #gammas = 3.01
      # lambdas = 1.0
      # initial_state = initial_state.tiss
      # time = days.model
      # data = coraldata.tiss
      
      betas = params[1]
      gammas = params[2]
      lambdas = params[3]
      
      SIR.out = tryCatch({
        data.frame(ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                       time, SIR, c(b = betas, g = gammas,
                                    N = initial_state[4],
                                    l = lambdas,
                                    C = initial_state[5])))
      }, error = function(e) {
        print("Error in ODE:")
        print(e)
        return(NA)
      })
      
      # Extract simulated values
      sim.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I')]
      sim.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
      
      # Extract observed values
      obs.inf.total = unlist(data[[1]])
      obs.rem.total = unlist(data[[2]])
      
      # Calculate residuals
      # diff.inf = (sim.inf - obs.inf)
      diff.rem.total = (sim.rem.total - obs.rem.total)
      
      #minimize using sum of squared residuals
      sum_diff.total = sum(diff.rem.total^2)
      
      # NOTE - see Kalizhanova et al. 2024 (TB SIR) for other error assessments - including mean absolute error (MAE)
      # Total Sum of Squares (TSS) for removal only
      mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
      tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
      
      #R-squared
      r_squared_rem.total = 1 - (sum_diff.total / tss_rem.total)
      
      error_eval <<- error_eval %>%
        mutate(
          SSR = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', sum_diff.total, SSR),
          TSS = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', tss_rem.total, TSS),
          R_squared = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', r_squared_rem.total, R_squared),
          
          # Update list-columns with vectors
          sim_inf = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(sim.inf.total), sim_inf),
          sim_rem = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(sim.rem.total), sim_rem),
          obs_inf = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(obs.inf.total), obs_inf),
          obs_rem = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(obs.rem.total), obs_rem)
        )
      
      return(sum_diff.total) #return only the residual metric for the epidemic wave being fit to
    }
    
    # uniform
    lower_bounds.tiss = c(0, 0, lambda.modifier)  #lower bounds for beta, gamma and lambda
    # lower_bounds.tiss = c(0, 0.6, 1)  #lower bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(1, 1, 0.3)  #upper bounds for beta, gamma and lambda
    upper_bounds.tiss = c(4, 4, lambda.modifier)  #upper bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(5, 5, 1)  #upper bounds for beta, gamma and lambda
    
    control = list(itermax = 100)  # Maximum number of iterations. 200 is default
    
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
    
    params.basic.full[[i]] = c(min.beta.tiss, min.beta.tiss.adj, min.gamma.tiss, min.lambda.tiss, R0, cover.site)
    
    #simulation using initial state variables and best-fit beta/gamma parameters
    SIR.out.tiss = data.frame(ode(c(S = initial_state.tiss[1], I = initial_state.tiss[2], R = initial_state.tiss[3]),
                                  days.model, SIR, c(b = min.beta.tiss, g = min.gamma.tiss,
                                                     N = initial_state.tiss[4],
                                                     l = min.lambda.tiss,
                                                     C = initial_state.tiss[5])))
    my.SIRS.basic.full[[i]] = SIR.out.tiss
  }
  
  ################################## Save output ##################################
  
  # STOPPING POINT - 12 NOV 2024
  # 1.) Still need to try this with DHWs and see if I can produce a less wiggly output
  #       - DONE, looks okay
  # 2.) Try only fitting pre-threshold-break for the initial fit and see what happens
  #       - DONE, looks great
  # 3.) Figure out why Nearshore is still not looking great
  #       - DONE, still fits strangely. may just be because of its late timing of outbreak. or differences in transmission due to ocean currents, connectivity, reservoirs nearby, etc.
  # 4.) Test different thresholds, including 28.5C (SST when first wave began offshore)
  #       - DONE, interesting output. largely similar, but the fitting is MUCH slower and still poorly constrained at Nearshore
  # 5.) Try this on multi-group SIR?
  #       - Not done yet!
  # 6.) Compare to pixel-level Coral Reef Watch data
  #       - Not done yet!
  # 7.) Try smoothing the SST data to isolate the trend and eliminate some noise
  #       - Not done yet!
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "basic_SIR_workspace.RData"))
  