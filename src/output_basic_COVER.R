  
  
  # some thoughts:
  #   - this seemed to work better when I was working with numbers really close to, but below, 1.0. struggling
  #       when working with values close to, but ABOVE, 1.0.
  #   - why might that be?
  #   - I was getting R0's that made more sense on January 26th, and also earlier in the afternoon Jan 29 (check Box version control)
  #
  #   - I tested something below where 1.0 is the limit for a coral cover of 100%...but not the same as what worked well before.
  #       - with a lambda of 1.0, you end up with a very large effect of coral cover on transmission. this means that wonky things can
  #           happen
  #   - theoretically, should there not be some very large lambda value for which coral cover explains so much of transmission that
  #       a cross-site transfer of beta from Nearshore to Offshore / Midchannel simply results in zero outbreak?
  #         - that's what I was seeing before (like with the sqrt function directly), where past certain thresholds, there is either
  #           zero outbreak after the transfer, or an outbreak that reduces to exactly the Nearshore outbreak, proportionally
  #             - I should probably get this back to that place. and maybe hit a stopping point and move on to the multi-group model again
  #   - okay! I sort of hit that point again: a lambda of 3.0 does indeed result in no outbreak for Midchannel & Offshore
  #
  #   - one takeway I've got is, it seems like it might be really important to have values close to 1.0 so that the infection
  #       can take off quickly early in the outbreak. the total removal can get close with current methods below when transfering
  #       between sites, but there's usually a phase difference (temporal mismatch). consider looking at the other functions again, like
  #       a simple square root or the exponential one Dan did. these might provide clues too
  #         - part of the problem is that when I change lambda, it's changing the relationship of ALL sites to cover, but also
  #             their relationship with each other
  
  rm(list=ls())
  
  load("FLKEYS_workspace.RData")
  
  library(tidyverse)
  library(DEoptim)
  library(deSolve)

  # #ONLY if you want to exclude degree-heating weeks (DHW or DHWs)
  # DHW.modifier = 8 #10
  # tissue.summary = tissue.summary %>%
  #   group_by(Site_type) %>%
  #   slice(head(row_number(), n()-DHW.modifier)) #group_size()
  # obs = obs %>%
  #   group_by(Site, Compartment, Category) %>%
  #   slice(head(row_number(), n()-DHW.modifier))
  # obs.basic = obs.basic %>%
  #   group_by(Site, Compartment) %>%
  #   slice(head(row_number(), n()-DHW.modifier))
  
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
  
  my.SIRS = vector('list', length(sites))
  params = vector('list', length(sites))
  
  for(i in 1:length(sites)){
    site = sites[i]
  
    # site = "Midchannel" #for testing purposes
    # i = 1
    # site = "Nearshore" #for testing purposes
    # i = 2
    # site = "Offshore" #for testing purposes
    # i = 3
    
    days = tissue.summary %>% na.omit() %>% filter(Site_type == site) %>% pull(Day)
    
    ###tissue SA & counts of susceptible (healthy) hosts at the start of monitoring
    # also, selecting the correct starting polyp_SA for later in the loop
    prev.timepoint = ''
    if(site == 'Nearshore'){
      prev.timepoint = 'T11' #timepoint before first documented infection timepoint [site specific]
      polyp_SA = polyp_SA.nearshore[[1]]
      area.site = 200 #meters squared, 2D. this could be changed to reflect available 3D substrate within the site, but we don't have that
    } else if(site == "Midchannel"){
      prev.timepoint = 'T7'
      polyp_SA = polyp_SA.midchannel[[1]]
      area.site = 200
    } else{ #offshore
      prev.timepoint = 'T5'
      polyp_SA = polyp_SA.offshore[[1]]
      area.site = 200
    }
    
    time = time_list[[i]]
    
    # polyp_SA = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
    
    S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
    N.site = S0.snapshot$tot.sustiss
    cover.site = N.site / area.site * 0.3 #0.2840909 (~1:3X) is the ratio of our tissue density to CPCe-based cover (Williams et al. 2021)
    
    #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
    LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
    MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
    HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
    inftiss = LS.inftiss+MS.inftiss+HS.inftiss
    
    # polyp_SA = min(inftiss[1:5])/5
    # polyp_SA = inftiss[1]/1.5
    polyp_SA = min(inftiss[1:5])/1.5
    
    LS.remtiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Dead") %>% pull(prop)
    MS.remtiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Dead") %>% pull(prop)
    HS.remtiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Dead") %>% pull(prop)
    remtiss = LS.remtiss+MS.remtiss+HS.remtiss
    
    #initial conditions
    I.tiss = polyp_SA
    S.tiss = N.site - I.tiss
    R.tiss = 0
    
    ############################## optimize parameters ##############################
    # Set up the data and initial conditions
    coraldata.tiss = list(inftiss, remtiss)
    initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
    
    objective_function = function(params, data, time, initial_state){
      betas = params[1]
      gammas = params[2]
      lambdas = params[3]
      
      SIR.out = data.frame(ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
                               time, SIR, c(b = betas, g = gammas,
                                            N = initial_state[4],
                                            l = lambdas,
                                            C = initial_state[5])))
      
      #extract simulated values at time points matching observations
      sim.inf = SIR.out[which(SIR.out$time %in% days), which(colnames(SIR.out) %in% 'I')]
      sim.rem = SIR.out[which(SIR.out$time %in% days), which(colnames(SIR.out) %in% 'R')]
      
      #extract observed values [repetitive code, but works]
      obs.inf = unlist(data[1])
      obs.rem = unlist(data[2])
      
      #use ratio of maximum removed value to maximum infected value to rescale removed. fits the 'I' curve better this way
      max_obs_inf = max(obs.inf)
      max_obs_rem = max(obs.rem)
      rem.inf.ratio = max_obs_rem/max_obs_inf
      obs.rem = obs.rem/rem.inf.ratio
      sim.rem = sim.rem/rem.inf.ratio
      
      # #test to separately rescale removed simulated curve
      # max_sim_inf = max(sim.inf)
      # max_sim_rem = max(sim.rem)
      # rem.inf.ratio = max_sim_rem/max_sim_inf
      # sim.rem = sim.rem/rem.inf.ratio
      
      # Calculate differences after normalization
      diff.inf = (sim.inf - obs.inf)
      diff.rem = (sim.rem - obs.rem)
      
      #minimize using sum of absolute differences
      sum_squared_diff_I = sum(sum(abs(diff.inf))) #can multiply this by 2 or similar to weight it extra
      sum_squared_diff_R = sum(sum(abs(diff.rem)))
      # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
      sum_squared_diff = sum_squared_diff_R
      
      # #minimize using sum of squared differences
      # sum_squared_diff_I = sum(sum(diff.inf^2))
      # sum_squared_diff_R = sum(sum(diff.rem^2))
      # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
      
      return(sum_squared_diff) #minimized error
    }
    ############################## OPTIMIZE PARAMETERS ############################################################
    
    # uniform
    lower_bounds.tiss = c(0, 0, lambda.modifier)  #lower bounds for beta, gamma and lambda
    # lower_bounds.tiss = c(0, 0.6, 1)  #lower bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(1, 1, 0.3)  #upper bounds for beta, gamma and lambda
    upper_bounds.tiss = c(4, 4, lambda.modifier)  #upper bounds for beta, gamma and lambda
    # upper_bounds.tiss = c(5, 5, 1)  #upper bounds for beta, gamma and lambda
    
    
    control = list(itermax = 100)  # Maximum number of iterations. 200 is default
    
    # Run the optimization
    result.tiss = DEoptim(fn = objective_function, lower = lower_bounds.tiss, upper = upper_bounds.tiss,
                          data = coraldata.tiss, time = time, initial_state = initial_state.tiss,
                          control = control) # 411: turn this ON for tissue!
    
    # Extract the optimized parameters
    optimized_params.tiss = result.tiss$optim$bestmem  # 411: turn this ON for tissue!
    
    # Print the optimized parameters
    min.beta.tiss = as.numeric(optimized_params.tiss[1])
    min.gamma.tiss = as.numeric(optimized_params.tiss[2])
    min.lambda.tiss = as.numeric(optimized_params.tiss[3])
    min.beta.tiss.adj = min.beta.tiss * (1 / (1 + exp(-lambda.modifier * (cover.site))) + offset)
    # min.beta.tiss.adj = min.beta.tiss * (min.lambda.tiss * (1-exp(-130*(cover.site))))
    # min.beta.tiss.adj = (min.beta.tiss * (1 + min.lambda.tiss * sqrt(cover.site)))
    # min.beta.tiss.adj = (min.beta.tiss * (min.lambda.tiss * sqrt(cover.site)))
    R0 = min.beta.tiss.adj / min.gamma.tiss
    cat("Optimized Tissue Model Parameters for", site, " site:\n")
    cat("Beta:", min.beta.tiss, "\n")
    cat("Cover-adjusted beta:", min.beta.tiss.adj, "\n")
    cat("Gamma:", min.gamma.tiss, "\n")
    cat("Lambda:", min.lambda.tiss, "\n")
    cat("R0:", R0, '\n')
    
    params[[i]] = c(min.beta.tiss, min.beta.tiss.adj, min.gamma.tiss, min.lambda.tiss, R0, cover.site)
    
    #simulation using initial state variables and best-fit beta/gamma parameters
    SIR.out.tiss = data.frame(ode(c(S = initial_state.tiss[1], I = initial_state.tiss[2], R = initial_state.tiss[3]),
                                  time, SIR, c(b = min.beta.tiss, g = min.gamma.tiss,
                                               N = initial_state.tiss[4],
                                               l = min.lambda.tiss,
                                               C = initial_state.tiss[5])))
    my.SIRS[[i]] = SIR.out.tiss
    
  }
  
  # # Save/load workspace
  # save.image(file = "output_basic_COVER_1.0_workspace.RData")
  