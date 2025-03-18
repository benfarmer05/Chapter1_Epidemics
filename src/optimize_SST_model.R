  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(ggpmisc)
  library(patchwork)
  library(deSolve)
  library(dplyr)
  library(DEoptim)
  
  #import workspace from upstream script
  load(here("output/plots_basic_workspace.RData"))
  
  # ################################## More mature sigmoid & logistic behaviors of temperature ##################################
  # # AKA - zeta / eta parameter plotting sandbox
  # 
  # # Load necessary libraries
  # library(ggplot2)
  # library(viridis)
  # library(patchwork)  # For arranging multiple plots
  # 
  # # Define temperature range from 23 to 33
  # T <- seq(23, 33, length.out = 100)
  # 
  # #APPROACH 1: Inflection point at 30.5C (original)
  # temp_effect_inflection <- function(temp, T_inflection = 30.5, zeta_bounded = 0.5) {
  #   
  #   # Convert bounded zeta to effective zeta
  #   # zeta_effective <- 1 * zeta_bounded/(1-1*zeta_bounded)
  #   zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
  #   
  #   # Logistic function with inflection point at T_inflection
  #   modifier <- 1 / (1 + exp(-zeta_effective * (temp - T_inflection)))
  #   
  #   return(modifier)
  # }
  # 
  # #APPROACH 2: Original model clarified
  # temp_effect_original <- function(temp, T_min = 23, zeta_bounded = 0.5) {
  #   
  #   # Convert bounded zeta to effective zeta
  #   
  #   # zeta_effective <- 5 * zeta_bounded/(1-0.7*zeta_bounded)
  #   zeta_effective <- 1 * zeta_bounded/(1 - 1*zeta_bounded)
  #   
  #   # Original formula clarified
  #   # This is a non-centered sigmoid that:
  #   # - Starts near 0 at T_min
  #   # - Increases asymptotically toward 1 as temp increases
  #   # - Has no inflection point assumption
  #   modifier <- (1 - exp(-zeta_effective * (temp - T_min))) / (1 + exp(-zeta_effective * (temp - T_min)))
  #   
  #   return(modifier)
  # }
  # 
  # #APPROACH 3: Decreasing function with temperature
  # temp_effect_decreasing <- function(temp, T_max = 33, T_min = 23, zeta_bounded = 0.5) {
  #   
  #   # Convert bounded zeta to effective zeta
  #   # zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
  #   zeta_effective <- 1 * zeta_bounded/(1-1*zeta_bounded)
  #   
  #   # Calculate midpoint of temperature range for centering the sigmoid
  #   T_mid <- (T_min + T_max) / 2
  #   
  #   # Decreasing sigmoid that starts near 1 (at low temps) and approaches 0 as temp increases
  #   # Note the positive sign before zeta_effective which makes the function decrease with temperature
  #   modifier <- 1 / (1 + exp(zeta_effective * (temp - T_mid)))
  #   
  #   return(modifier)
  # }
  # 
  # # APPROACH 4: Non-centered sigmoid decreasing with temperature (FIXED)
  # temp_effect_noncentered_decreasing <- function(temp, T_min = 23, T_max = 33, zeta_bounded = 0.5) {
  # 
  #   # Convert bounded zeta to effective zeta
  #   zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded) # apply more tuning to behavior
  #   # zeta_effective <- zeta_bounded/(1-zeta_bounded)
  # 
  #   # Modified non-centered sigmoid that:
  #   # - Equals 1 at T_min (23°C)
  #   # - Decreases toward 0 as temperature increases
  #   # - Maintains the non-centered sigmoid shape but inverted
  # 
  #   # Calculate the relative position in the temperature range
  #   rel_temp <- (temp - T_min) / (T_max - T_min)
  # 
  #   # Apply the decreasing non-centered sigmoid
  #   modifier <- 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))
  #   # modifier <- 1 - ((1 - exp(-zeta_bounded * rel_temp)) / (1 + exp(-zeta_bounded * rel_temp)))
  # 
  #   return(modifier)
  # }
  # 
  # # temp_effect_noncentered_decreasing <- function(temp, T_min = 23, T_max = 33, zeta_bounded = 0.5) {
  # #   # Convert bounded zeta to effective zeta
  # #   zeta_effective <- zeta_bounded/(1-zeta_bounded)
  # # 
  # #   # Calculate the relative position in the temperature range
  # #   rel_temp <- (temp - T_min) / (T_max - T_min)
  # # 
  # #   # Most simplified form - standard logistic function
  # #   modifier <- 1 / (1 + exp(zeta_effective * rel_temp))
  # # 
  # #   return(modifier)
  # # }
  # 
  # # Define bounded zeta values from 0 to 1 (avoiding 1 which would create infinity)
  # zeta_values <- c(0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
  # 
  # # Create data for all approaches
  # plot_data_inflection <- data.frame()
  # plot_data_original <- data.frame()
  # plot_data_decreasing <- data.frame()
  # plot_data_noncentered_decreasing <- data.frame()
  # 
  # for (zeta in zeta_values) {
  #   # Data for approach 1
  #   modifiers_inflection <- sapply(T, function(t) temp_effect_inflection(t, zeta_bounded = zeta))
  #   temp_data_inflection <- data.frame(
  #     Temperature = T,
  #     Modifier = modifiers_inflection,
  #     Zeta = as.factor(zeta)
  #   )
  #   plot_data_inflection <- rbind(plot_data_inflection, temp_data_inflection)
  #   
  #   # Data for approach 2
  #   modifiers_original <- sapply(T, function(t) temp_effect_original(t, zeta_bounded = zeta))
  #   temp_data_original <- data.frame(
  #     Temperature = T,
  #     Modifier = modifiers_original,
  #     Zeta = as.factor(zeta)
  #   )
  #   plot_data_original <- rbind(plot_data_original, temp_data_original)
  #   
  #   # Data for approach 3
  #   modifiers_decreasing <- sapply(T, function(t) temp_effect_decreasing(t, zeta_bounded = zeta))
  #   temp_data_decreasing <- data.frame(
  #     Temperature = T,
  #     Modifier = modifiers_decreasing,
  #     Zeta = as.factor(zeta)
  #   )
  #   plot_data_decreasing <- rbind(plot_data_decreasing, temp_data_decreasing)
  #   
  #   # Data for approach a4
  #   modifiers_noncentered_decreasing <- sapply(T, function(t) temp_effect_noncentered_decreasing(t, zeta_bounded = zeta))
  #   temp_data_noncentered_decreasing <- data.frame(
  #     Temperature = T,
  #     Modifier = modifiers_noncentered_decreasing,
  #     Zeta = as.factor(zeta)
  #   )
  #   plot_data_noncentered_decreasing <- rbind(plot_data_noncentered_decreasing, temp_data_noncentered_decreasing)
  # }
  # 
  # # Plot 1: Original inflection point model
  # p1 <- ggplot(plot_data_inflection, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
  #   geom_line(size = 1) +
  #   labs(
  #     title = "Temperature Threshold with Classic Logistic Function",
  #     x = "Temperature (°C)",
  #     y = "Transmission Modifier"
  #   ) +
  #   theme_minimal(base_size = 12) +
  #   theme(
  #     legend.position = "right",
  #     panel.grid.minor = element_blank(),
  #     plot.title = element_text(face = "bold")
  #   ) +
  #   scale_color_viridis_d(option = "inferno") +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  #   scale_x_continuous(breaks = seq(23, 33, 2)) +
  #   geom_vline(xintercept = 30.5, linetype = "dashed", alpha = 0.5)
  # 
  # # Plot 2: Original increasing model
  # p2 <- ggplot(plot_data_original, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
  #   geom_line(size = 1) +
  #   labs(
  #     title = "Temperature Effect with non-centered Sigmoid Function",
  #     x = "Temperature (°C)",
  #     y = "Transmission Modifier"
  #   ) +
  #   theme_minimal(base_size = 12) +
  #   theme(
  #     legend.position = "right",
  #     panel.grid.minor = element_blank(),
  #     plot.title = element_text(face = "bold")
  #   ) +
  #   scale_color_viridis_d(option = "inferno") +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  #   scale_x_continuous(breaks = seq(23, 33, 2)) +
  #   geom_vline(xintercept = 23, linetype = "dashed", alpha = 0.5)
  # 
  # # Plot 3: Classic sigmoid decreasing model
  # p3 <- ggplot(plot_data_decreasing, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
  #   geom_line(size = 1) +
  #   labs(
  #     title = "Temperature Effect with Decreasing Sigmoid Function",
  #     x = "Temperature (°C)",
  #     y = "Transmission Modifier"
  #   ) +
  #   theme_minimal(base_size = 12) +
  #   theme(
  #     legend.position = "right",
  #     panel.grid.minor = element_blank(),
  #     plot.title = element_text(face = "bold")
  #   ) +
  #   scale_color_viridis_d(option = "inferno") +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  #   scale_x_continuous(breaks = seq(23, 33, 2)) +
  #   geom_vline(xintercept = 28, linetype = "dashed", alpha = 0.5)
  # 
  # # Plot 4: Fixed non-centered sigmoid decreasing model
  # p4 <- ggplot(plot_data_noncentered_decreasing, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
  #   geom_line(size = 1) +
  #   labs(
  #     title = "Temperature Effect with Inverted Non-centered Sigmoid",
  #     x = "Temperature (°C)",
  #     y = "Transmission Modifier"
  #   ) +
  #   theme_minimal(base_size = 12) +
  #   theme(
  #     legend.position = "right",
  #     panel.grid.minor = element_blank(),
  #     plot.title = element_text(face = "bold")
  #   ) +
  #   scale_color_viridis_d(option = "inferno") +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  #   scale_x_continuous(breaks = seq(23, 33, 2)) +
  #   geom_vline(xintercept = 23, linetype = "dashed", alpha = 0.5)
  # 
  # # Combine plots - using a 2x2 grid for better viewing
  # (p1 + p2) / (p3 + p4)
  # 
  # # # Let's focus only on the fourth plot
  # # p4
  # 
  # # Function explanations for mathematical clarity
  # cat("Approach 1: Inflection Point at 30.5°C\n")
  # cat("Formula: modifier = 1 / (1 + exp(-zeta_effective * (temp - T_inflection)))\n")
  # cat("- Classic logistic function centered at T_inflection\n")
  # cat("- At T_inflection (30.5°C), modifier = 0.5\n")
  # cat("- Assumes symmetric behavior around inflection point\n\n")
  # 
  # cat("Approach 2: Original Model Clarified\n")
  # cat("Formula: modifier = (1 - exp(-zeta_effective * (temp - T_min))) / (1 + exp(-zeta_effective * (temp - T_min)))\n")
  # cat("- At T_min (23°C), modifier ≈ 0\n")
  # cat("- No assumption of inflection point\n")
  # cat("- Asymmetric behavior that rises more steeply at first\n\n")
  # 
  # cat("Approach 3: Decreasing Function with Temperature\n")
  # cat("Formula: modifier = 1 / (1 + exp(zeta_effective * (temp - T_mid)))\n")
  # cat("- At low temperatures (23°C), modifier ≈ 1\n")
  # cat("- At high temperatures (33°C), modifier approaches 0\n")
  # cat("- Higher zeta values create steeper drop in transmission\n")
  # cat("- Centered at T_mid (28°C) where modifier = 0.5\n\n")
  # 
  # cat("Approach 4: Inverted Non-centered Sigmoid with Temperature\n")
  # cat("Formula: modifier = 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))\n")
  # cat("- At T_min (23°C), modifier = 1 exactly for all zeta values\n")
  # cat("- Decreases toward 0 as temperature increases\n")
  # cat("- Higher zeta values create steeper drop in transmission\n")
  # cat("- Preserves the asymmetric shape but inverted from approach 2\n")
  # 
  ################################## SST optimization ##################################
  
  #making choice to pull boundaries of SST from entire dataset (1985 - 2024). could constrain
  T_min = DHW.CRW.full %>%
    pull(SST.90th_HS) %>%
    min
  T_max = DHW.CRW.full %>%
    pull(SST.90th_HS) %>%
    max
  
  SIR_project = function(t,y,p,SST,DHW){ # 'p' is parameters or params
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
      
      #host density-null conditions
      transmission_modifier = 1
      removal_rate = g
      transmission_rate = b * transmission_modifier
      
      # # Convert bounded zeta to effective zeta (scales function behavior from 0 to 1)
      # zeta_effective = z / (1 - z) # zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
      # eta_effective = e / (1 - e)
      # 
      # # Calculate the relative position in the temperature range (sets domain window of function to realistic temperatures)
      # rel_temp <- (current_SST - T_min) / (T_max - T_min)
      # 
      # #modify epidemic rates by sea surface temperature
      # beta_scaling_factor = 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))
      # gamma_scaling_factor = 1 - ((1 - exp(-eta_effective * rel_temp)) / (1 + exp(-eta_effective * rel_temp)))
      # # beta_scaling_factor = 1
      # # gamma_scaling_factor = 1
      # transmission_rate <- transmission_rate * beta_scaling_factor
      # removal_rate <- removal_rate * gamma_scaling_factor
      
      

        zeta_effective <- z / (1 - z)
        modifier <- 1 / (1 + exp(zeta_effective * (current_SST - 30.5)))
        
        # NOTE - still working on this....can't quite seem to nail something down
        

      dS.dt = -transmission_rate * S * I / N
      dI.dt = transmission_rate * S * I / N - removal_rate * I
      dR.dt = removal_rate * I
  
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  my.SIRS.SST = vector('list', length(sites))
  params.SST = vector('list', length(sites))
  curr.type = 'Fitted' #the below is for a fitting model for multi-host transmission (no DHW or projection)
  susceptible_ref_SST = susceptible_ref %>%
    mutate(Site = case_when(
      Site == "Offshore" ~ "off",
      Site == 'Midchannel' ~ 'mid',
      Site == 'Nearshore' ~ 'near'
    ))
  
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
    
    N.site = susceptible_ref_SST %>%
      filter(Site == site.loop) %>%
      slice(1) %>% #all values for N.site are the same between categories, so slice first row
      pull(N.site)
    
    cover.site = susceptible_ref_SST %>%
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
    
    ################################## Optimize single-host SST model ##################################
    # Set up the data and initial conditions
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    # Define the objective function for optimization
    objective_function = function(params, data, time, initial_state){
      
      # # Print parameters to track what values are being tested
      # cat("Testing params:", params, "\n")
      
      # #testing
      # betas = 0.76
      # gammas = 0.56
      # lambdas = 1.0
      # zetas = 0.89
      # etas = 0.70
      # initial_state = initial_state.tiss
      # time = days.model
      # data = coraldata.tiss
      
      betas = params[1]
      gammas = params[2]
      lambdas = params[3]
      zetas = params[4]
      etas = params[5]
      
      # Check if SST_df and DHW_df are valid
      if(any(is.na(SST_df$time)) || any(is.na(SST_df$SST)) || 
         any(is.na(DHW_df$time)) || any(is.na(DHW_df$DHW))) {
        cat("Warning: NAs in SST or DHW data\n")
      }
      
      # Expand tryCatch to provide more information
      SIR.out = tryCatch({
        result <- ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                      time, SIR_project, c(b = betas, g = gammas,
                                           N = initial_state[4],
                                           z = zetas,
                                           e = etas,
                                           l = lambdas,
                                           C = initial_state[5]),
                      SST = SST_df,
                      DHW = DHW_df)
        # Check if ODE output contains NAs
        if(any(is.na(result))) {
          cat("Warning: ODE solution contains NAs\n")
        }
        data.frame(result)
      }, error = function(e) {
        cat("Error in ODE:", conditionMessage(e), "\n")
        return(NA)
      })
      
      # Check if ODE solving failed
      if(all(is.na(SIR.out))) {
        cat("ODE solving failed completely\n")
        return(Inf)  # Return Inf instead of NaN to avoid optimizer error
      }
      
      # Extract simulated values
      sim.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I')]
      sim.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
      
      # Extract observed values
      obs.inf.total = unlist(data[[1]])
      obs.rem.total = unlist(data[[2]])
      
      # Calculate residuals
      diff.inf.total = (sim.inf.total - obs.inf.total)
      diff.rem.total = (sim.rem.total - obs.rem.total)
      
      #aggregate sum of absolute residuals
      # NOTE - this is not sum of squares and should be clearly stated/defended in the manuscript if used
      sum_abs_diff_I.total = sum(sum(abs(diff.inf.total)))
      sum_abs_diff_R.total = sum(sum(abs(diff.rem.total)))
      sum_diff.abs.total = sum_abs_diff_R.total
      
      #minimize using sum of squared residuals
      sum_diff.total = sum(diff.rem.total^2)
      
      return(sum_diff.abs.total) #return only the residual metric for the epidemic wave being fit to
    }
    
    ############################## OPTIMIZE PARAMETERS ############################################################
    
    #boundary conditions
    lower_bounds.tiss = c(0, 0, lambda.modifier, 0, 0)  # Lower bounds for betas, gammas, lambdas, zetas, etas
    upper_bounds.tiss = c(4, 4, lambda.modifier, 1, 1)  # Upper bounds for betas, gammas, lambdas, zetas, etas
    
    control = list(itermax = 50)  # number of optimizer iterations. 200 is default
    
    # Run the optimization
    result.tiss = DEoptim(fn = objective_function, lower = lower_bounds.tiss, upper = upper_bounds.tiss,
                          data = coraldata.tiss, time = days.model, initial_state = initial_state.tiss,
                          control = control)
    
    # Extract the optimized parameters
    optimized_params.tiss = result.tiss$optim$bestmem
    
    # Print the optimized parameters
    min.beta.tiss = as.numeric(optimized_params.tiss[1])
    min.gamma.tiss = as.numeric(optimized_params.tiss[2])
    min.zeta.tiss = as.numeric(optimized_params.tiss[4])
    min.eta.tiss = as.numeric(optimized_params.tiss[5])
    
    params.SST[[i]] = c(min.beta.tiss, min.gamma.tiss, min.zeta.tiss, min.eta.tiss)
    
    #simulation using initial state variables and best-fit beta/gamma/zeta/eta parameters
    SIR.out.tiss = data.frame(ode(c(S = initial_state.tiss[1], I = initial_state.tiss[2], R = initial_state.tiss[3]),
                                  days.model, SIR_project, c(b = min.beta.tiss, g = min.gamma.tiss,
                                                         N = initial_state.tiss[4],
                                                         z = min.zeta.tiss,
                                                         e = min.eta.tiss,
                                                         l = lambdas,
                                                         C = initial_state.tiss[5]),
                                  SST = SST_df,
                                  DHW = DHW_df))
    
    my.SIRS.SST[[i]] = SIR.out.tiss
  }
  
  ################################## optimized plots ##################################
  
  #Nearshore
  site.loop = 'Nearshore'
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 2
  
  output.basic.nearshore.SST = my.SIRS.SST[[order]]
  params.basic.nearshore.SST = params.SST[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.nearshore = params.basic.nearshore.SST[1]
  eta.nearshore = params.basic.nearshore.SST[2]
  
  tab.nearshore = tibble(round(beta.nearshore, 2), round(beta.nearshore.adj, 2), round(gamma.nearshore, 2),
                         round(R0.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.nearshore.SST[which(output.basic.nearshore.SST$time %in% days.obs), which(colnames(output.basic.nearshore.SST) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore.SST = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.nearshore.SST,
                              R_squared))
  
  output.basic.nearshore.SST = pivot_longer(output.basic.nearshore.SST, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.nearshore.SST)[1] = 'days.model'
  colnames(output.basic.nearshore.SST)[3] = 'tissue'
  
  p.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.nearshore.SST$days.model), y = min(output.basic.nearshore.SST$tissue)*0.7, label = list(tab.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #Midchannel
  site.loop = 'Midchannel'
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 1
  
  output.basic.midchannel.SST = my.SIRS.SST[[order]]
  params.basic.midchannel.SST = params.SST[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.midchannel = params.basic.midchannel.SST[1]
  eta.midchannel = params.basic.midchannel.SST[2]
  
  tab.midchannel = tibble(round(beta.midchannel, 2), round(beta.midchannel.adj, 2), round(gamma.midchannel, 2),
                          round(R0.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.midchannel.SST[which(output.basic.midchannel.SST$time %in% days.obs), which(colnames(output.basic.midchannel.SST) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.midchannel.SST = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.midchannel.SST,
                              R_squared))
  
  output.basic.midchannel.SST = pivot_longer(output.basic.midchannel.SST, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.midchannel.SST)[1] = 'days.model'
  colnames(output.basic.midchannel.SST)[3] = 'tissue'
  
  p.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.midchannel.SST$days.model), y = min(output.basic.midchannel.SST$tissue)*0.7, label = list(tab.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_line(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #Offshore
  site.loop = 'Offshore'
  curr.site = 'off'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 3
  
  output.basic.offshore.SST = my.SIRS.SST[[order]]
  params.basic.offshore.SST = params.SST[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.offshore = params.basic.offshore.SST[1]
  eta.offshore = params.basic.offshore.SST[2]
  
  tab.offshore = tibble(round(beta.offshore, 2), round(beta.offshore.adj, 2), round(gamma.offshore, 2),
                        round(R0.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.offshore.SST[which(output.basic.offshore.SST$time %in% days.obs), which(colnames(output.basic.offshore.SST) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.offshore.SST = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.offshore.SST,
                              R_squared))
  
  output.basic.offshore.SST = pivot_longer(output.basic.offshore.SST, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.offshore.SST)[1] = 'days.model'
  colnames(output.basic.offshore.SST)[3] = 'tissue'
  
  p.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.offshore.SST$days.model), y = min(output.basic.offshore.SST$tissue)*0.7, label = list(tab.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic.SST
  p.D.fit.nearshore.basic.SST
  
  p.I.fit.midchannel.basic.SST
  p.D.fit.midchannel.basic.SST
  
  p.I.fit.offshore.basic.SST
  p.D.fit.offshore.basic.SST
  
  p.fit.nearshore.basic.SST
  p.fit.midchannel.basic.SST
  p.fit.offshore.basic.SST
  
  
  
  
  
  
  # ############################## sandbox for manually tweaking zeta/eta ##################################
  # 
  # #sandbox conditions
  # beta.nearshore.sand.SST = 0.65137 # 0.76 #0.64
  # gamma.nearshore.sand.SST = 0.5622153 # 0.56
  # zeta_val = 0.5 #0.07107107 # 0.03 # 0.07107107 # 0
  # eta_val = 0.7 #3
  # lambda_val = 1.0
  # offset_val = 1 - 1 / (1 + exp(-lambda * 1.0))
  # 
  # SIR.sand.no_cover.SST = function(t,y,p,SST,DHW){ # 'p' is parameters or params
  #   {
  #     S = y[1]
  #     I = y[2]
  #     R = y[3]
  #   }
  #   with(as.list(p),{
  #     
  #     # Find the closest index for t in the SST data and fetch the associated SST and DHW
  #     closest_index <- which.min(abs(SST$time - t))
  #     if (closest_index < 1 || closest_index > nrow(SST)) {
  #       stop("Closest index is out of bounds for SST.")
  #     }
  #     current_SST = SST$SST[closest_index]
  #     current_DHW = DHW$DHW[closest_index]
  #     
  #     #host density-null conditions
  #     transmission_modifier = 1
  #     removal_rate = g
  #     transmission_rate = b * transmission_modifier
  #     
  #     # Convert bounded zeta to effective zeta (scales function behavior from 0 to 1)
  #     zeta_effective = zeta_val / (1 - zeta_val) # zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
  #     eta_effective = eta_val / (1 - eta_val)
  #     
  #     # Calculate the relative position in the temperature range (sets domain window of function to realistic temperatures)
  #     rel_temp <- (current_SST - T_min) / (T_max - T_min)
  #     
  #     
  #     # transmission_modifier = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C)) / (1 - exp(-k_val)))
  #     beta_scaling_factor = 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))
  #     gamma_scaling_factor = 1 - ((1 - exp(-eta_effective * rel_temp)) / (1 + exp(-eta_effective * rel_temp)))
  #     transmission_rate <- transmission_rate * beta_scaling_factor
  #     removal_rate <- removal_rate * gamma_scaling_factor
  #     
  #     
  #     dS.dt = -transmission_rate * S * I / N
  #     dI.dt = transmission_rate * S * I / N - removal_rate * I
  #     dR.dt = removal_rate * I
  #     
  #     return(list(c(dS.dt, dI.dt, dR.dt)))
  #   })
  # }
  # 
  # #run nearshore model
  # curr.site = 'near'
  # days.obs <- days_sites %>%
  #   filter(site == curr.site) %>%
  #   pull(days.obs) %>%
  #   unlist()
  # 
  # output.basic.nearshore.sand.SST = data.frame(ode(c(S = S.nearshore, I = I.nearshore, R = R.nearshore),
  #                                              days.model.nearshore, SIR.sand.no_cover.SST, c(b = beta.nearshore.sand.SST, g = gamma.nearshore.sand.SST,
  #                                                                                         N = N.nearshore,
  #                                                                                         z = zeta_val,
  #                                                                                         e = eta_val,
  #                                                                                         l = lambda.nearshore,
  #                                                                                         C = cover.nearshore),
  #                                              SST = SST_df,
  #                                              DHW = DHW_df))
  # 
  # sim.rem.total = output.basic.nearshore.sand.SST[which(output.basic.nearshore.sand.SST$time %in% days.obs), which(colnames(output.basic.nearshore.sand.SST) %in% 'R')]
  # obs.rem.total = obs.model %>%
  #   filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
  #   slice(head(row_number(), n()-DHW.modifier)) %>%
  #   pull(tissue)
  # if (length(obs.rem.total) > length(sim.rem.total)) {
  #   obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  # }
  # 
  # diff.rem.total = (sim.rem.total - obs.rem.total)
  # sum_diff.total = sum(diff.rem.total^2)
  # mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  # tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  # r_squared.basic.nearshore.sand.no_cover.SST = 1 - (sum_diff.total / tss_rem.total)
  # 
  # 
  # 
  # output.basic.nearshore.sand.SST = pivot_longer(output.basic.nearshore.sand.SST, cols = -1, names_to = c("Compartment")) %>%
  #   mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
  #   mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
  #                               ifelse(Compartment == 'I', 'Infected',
  #                                      ifelse(Compartment == 'R', 'Dead', Compartment))))
  # 
  # colnames(output.basic.nearshore.sand.SST)[1] = 'days.model'
  # colnames(output.basic.nearshore.sand.SST)[3] = 'tissue'
  # 
  # p.fit.nearshore.basic.DHW.sand.SST = ggplot(data = output.basic.nearshore.sand.SST, aes(days.model, tissue, colour = Compartment)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of tissue (m2)") +
  #   ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   annotate(geom = "table", x = min(output.basic.nearshore.sand.SST$days.model), y = min(output.basic.nearshore.sand.SST$tissue)*0.7, label = list(tab.nearshore),
  #            vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
  #   theme_classic(base_family = 'Georgia') +
  #   theme(panel.background = element_rect(fill = "gray90"))
  # 
  # p.S.fit.nearshore.basic.DHW.sand.SST = ggplot(data = output.basic.nearshore.sand.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of susceptible tissue") +
  #   ggtitle(paste(c("", site.loop), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.I.fit.nearshore.basic.DHW.sand.SST = ggplot(data = output.basic.nearshore.sand.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of infected tissue") +
  #   ggtitle(paste(c("", site.loop), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.D.fit.nearshore.basic.DHW.sand.SST = ggplot(data = output.basic.nearshore.sand.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of dead tissue") +
  #   ggtitle(paste(c("", site.loop), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.fit.nearshore.basic.DHW.sand.SST
  # p.I.fit.nearshore.basic.DHW.sand.SST
