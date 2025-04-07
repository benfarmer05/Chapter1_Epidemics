  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(DEoptim)
  library(deSolve)
  
  ################################## Set-up ##################################
  
  #import workspace from upstream script
  load(here("output/plots_basic_workspace.RData"))
  
  #all modeling output in this script is for a multi-host SIR model
  curr.host = 'Multi-host'
  
  susceptible_ref = susceptible_ref %>%
    mutate(Site = case_when(
      Site == "Offshore" ~ "off",
      Site == 'Midchannel' ~ 'mid',
      Site == 'Nearshore' ~ 'near'
    ))
  
  # # NOTE - 19 Feb 2025 - something I was testing to change the starting values of infected tissue
  # # Find the smallest nonzero tissue value across all categories within each Site
  # obs.model_updated <- obs.model %>%
  #   filter(Compartment == "Infected") %>%
  #   group_by(Site) %>%
  #   mutate(min_tissue = min(tissue[tissue > 0], na.rm = TRUE)) %>%
  #   ungroup()
  # 
  # # Apply changes only to "Total" and "High"
  # obs.model_updated <- obs.model_updated %>%
  #   arrange(Site, date) %>%
  #   group_by(Site, Category) %>%  # Ensure min(which(tissue > 0)) is calculated within each Category
  #   mutate(tissue = case_when(
  #     Category %in% c("Total", "High") & row_number() == min(which(tissue > 0)) ~ min_tissue / case_when(
  #       Site == "near" ~ polyp_SA.minimizer.nearshore,
  #       Site == "mid" ~ polyp_SA.minimizer.midchannel,
  #       Site == "off" ~ polyp_SA.minimizer.offshore
  #     ),
  #     TRUE ~ tissue
  #   )) %>%
  #   ungroup() %>%
  #   select(-min_tissue)
  # 
  # # Update obs.model in-place without changing row count
  # obs.model <- obs.model %>%
  #   rows_update(obs.model_updated, by = c("Site", "date", "Compartment", "Category"))
  
  
  
  # Update Error_eval with metrics and thresholds - pulled values from output_basic upstream
  error_eval <<- error_eval %>%
    mutate(
      SST_threshold = if_else(host == 'Multi-host', SST_threshold_value, SST_threshold),
      date_thresh = if_else(host == 'Multi-host', date_threshold, date_thresh)
    )
  
  # # Scenario 1 [maximum transmission modifier of 1.0, with 100% coral cover; uniform lambda].
  #lambda of 0.0: no effect of cover; overpredicted outbreaks. [Nearshore -> Offshore/Midchannel]
  #lambda of 0.3: essentially no effect of cover; overpredicted outbreaks. [Nearshore -> Offshore/Midchannel]
  #lambda of 0.6 / 0.7: essentially no effect of cover; overpredicted outbreaks. [Nearshore -> Offshore/Midchannel]
  #lambda of 3.0: essentially no effect of cover; overpredicted outbreaks. [Nearshore -> Offshore/Midchannel]
  #lambda of 7.0: super overpredicted LS/MS loss. the outbreak happens too early and too strong. [Offshore -> Nearshore/Midchannel]
  #lambda of 8.0: looks really good, at least going from Nearshore to Offshore. now trying reverse
  #lambda of 20.0: no outbreak!!
  
  # # Scenario 2 [maximum transmission modifier of 1.0, with 100% coral cover; varying lambda]
  #lambdas of 4.0, 5.0, 10.0: the amplifying effect of relative cover is different between groups. might be getting somewhere?
  #lambdas of 12.0, 12.0, 8.0: effect of cover at LS & MS is too great (going both directions)
  #lambdas of 3.0, 3.0, 8.0: effect of cover at LS & MS is too great (going both directions)
  
  # lambda = as.numeric(8.0) #0.7 worked well for basic SIR, 2.0 for multi-group [other script] - but likely needs to be higher here
  # offset = 1 - 1 / (1 + exp(-lambda * 1.0))
  
  lambda = 8.0 # NOTE - return to this; will be important to clearly define lambda parameters
  
  # # 4/5/10 experiment
  # #   - focusing on the ability to predict removal here - not necessarily infection [and that is what is described in inline comments]
  # #   - "overprediction" means relative to the 2/2/5 simulation. closeness to reality is relative to 2/2/5 simulation as well
  # lambda.LS = 4.0 # overpredicted off.to.near (and farther from reality; underpredicted near.to.off (and farther from reality)
  # lambda.MS = 5.0 # overpredicted off.to.near (and farther from reality; underpredicted near.to.off (and farther from reality)
  # lambda.HS = 10.0 # underpredicted off.to.near (but closer to reality); underpredicted near.to.off (and farther from reality)
  # #so, LS should go even lower, MS should go even lower, and HS should go ... I don't know? maybe stay the same?
  
  # # 2/2/5 experiment - the "original one"
  # lambda.LS = 2.0 # overpredicted off.to.near (by a lot); underpredicted near.to.off (by a lot)
  # lambda.MS = 2.5 # overpredicted off.to.near (by a little); underpredicted near.to.off (by a lot)
  # lambda.HS = 5.0 # underpredicted off.to.near (by a lot); underpredicted near.to.off (by a little)
  
  # # 0/0/0 experiment
  # #   - haven't tested these yet. NOTE / STOPPING POINT - 8 OCT 2024:
  # #       what exactly are the null conditions for effect of coral cover? should figure that out
  # lambda.LS = 0.0
  # lambda.MS = 0.0
  # lambda.HS = 0.0
  
  # # 05/125/5 experiment
  # # comparing to 2/2/5 experiment above, trying something new
  # lambda.LS = 0.5 # overpredicted off.to.near (still kind of bad but better); no change near.to.off (still quite bad)
  # lambda.MS = 1.25 # quite good off.to.near (closer to reality); no change near.to.off (still quite bad)
  # lambda.HS = 5.0 # underpredicted off.to.near (even worse); about right near.to.off (closer to reality)
  
  # # 0/05/20 experiment
  # # comparing to 01/125/15 experiment above
  # lambda.LS = 0 # overpredicted off.to.near (still bad and got worse); no change near.to.off (still quite bad)
  # lambda.MS = 0.5 # overpredicted off.to.near (not pretty good but got worse); underpredicted near.to.of (still bad but slightly better)
  # lambda.HS = 20.0 # underpredicted off.to.near (off but still getting better); shape changed near.to.off (still pretty good)
  
  # 01/125/15 experiment
  # comparing to 05/125/5 experiment above
  lambda.LS = 0.1 # overpredicted off.to.near (got worse); no change near.to.off (still quite bad)
  lambda.MS = 1.25 # overpredicted off.to.near (a bit worse but still good); underpredicted near.to.off (bad but better)
  lambda.HS = 15.0 # underpredicted off.to.near (still bad but a lot better); shape changed near.to.off (still pretty good)
  
  # # experiment based off of post-run optimization
  # # comparing to 05/125/5 experiment above
  # lambda.LS = 0 # overpredicted off.to.near (got worse); no change near.to.off (still quite bad)
  # lambda.MS = 0.2 # overpredicted off.to.near (a bit worse but still good); underpredicted near.to.off (bad but better)
  # lambda.HS = 13 # underpredicted off.to.near (still bad but a lot better); shape changed near.to.off (still pretty good)
  
  offset.LS = 1 - 1 / (1 + exp(-lambda.LS * 1.0))
  offset.MS = 1 - 1 / (1 + exp(-lambda.MS * 1.0))
  offset.HS = 1 - 1 / (1 + exp(-lambda.HS * 1.0))
  
  
  # ### testing
  # 
  # #null conditions
  # lambda.LS = 1
  # lambda.MS = 1
  # lambda.HS = 1
  # offset.LS = 0
  # offset.MS = 0
  # offset.HS = 0
  # 
  # 
  # ### testing
  
  
  
  ################################## Model: multi-host ##################################
  SIR.multi = function(t,y,p){
    {
      S.LS = y[1]
      I.LS = y[2]
      R.LS = y[3]
  
      S.MS = y[4]
      I.MS = y[5]
      R.MS = y[6]
  
      S.HS = y[7]
      I.HS = y[8]
      R.HS = y[9]
      
      # P = y[10]
    }
    with(as.list(p),{
        P = (I.LS + I.MS + I.HS)

      #null conditions
      transmission_modifier.LS = 1
      transmission_modifier.MS = 1
      transmission_modifier.HS = 1
      
      # #with effect of coral cover
      # transmission_modifier.LS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.LS)) / (1 - exp(-k_val)))
      # transmission_modifier.MS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.MS)) / (1 - exp(-k_val)))
      # transmission_modifier.HS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.HS)) / (1 - exp(-k_val)))
      
      #hybrid of frequency and density-dependent
      dS.LS.dt = -b.LS*S.LS*(P) / N.LS * transmission_modifier.LS
      dI.LS.dt = b.LS*S.LS*(P) / N.LS * transmission_modifier.LS - g.LS*I.LS
      dR.LS.dt = g.LS*I.LS

      dS.MS.dt = -b.MS*S.MS*(P) / N.MS * transmission_modifier.MS
      dI.MS.dt = b.MS*S.MS*(P) / N.MS * transmission_modifier.MS - g.MS*I.MS
      dR.MS.dt = g.MS*I.MS

      dS.HS.dt = -b.HS*S.HS*(P) / N.HS * transmission_modifier.HS
      dI.HS.dt = b.HS*S.HS*(P) / N.HS * transmission_modifier.HS - g.HS*I.HS
      dR.HS.dt = g.HS*I.HS
      
      # #more frequency-dependent
      # dS.LS.dt = -b.LS*S.LS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.LS
      # dI.LS.dt = b.LS*S.LS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.LS - g.LS*I.LS
      # dR.LS.dt = g.LS*I.LS
      # 
      # dS.MS.dt = -b.MS*S.MS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.MS
      # dI.MS.dt = b.MS*S.MS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.MS - g.MS*I.MS
      # dR.MS.dt = g.MS*I.MS
      # 
      # dS.HS.dt = -b.HS*S.HS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.HS
      # dI.HS.dt = b.HS*S.HS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.HS - g.HS*I.HS
      # dR.HS.dt = g.HS*I.HS
      
      # #more density-dependent
      # dS.LS.dt = -b.LS*S.LS*(P) * transmission_modifier.LS
      # dI.LS.dt = b.LS*S.LS*(P) * transmission_modifier.LS - g.LS*I.LS
      # dR.LS.dt = g.LS*I.LS
      # 
      # dS.MS.dt = -b.MS*S.MS*(P) * transmission_modifier.MS
      # dI.MS.dt = b.MS*S.MS*(P) * transmission_modifier.MS - g.MS*I.MS
      # dR.MS.dt = g.MS*I.MS
      # 
      # dS.HS.dt = -b.HS*S.HS*(P) * transmission_modifier.HS
      # dI.HS.dt = b.HS*S.HS*(P) * transmission_modifier.HS - g.HS*I.HS
      # dR.HS.dt = g.HS*I.HS
      
      return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, dS.HS.dt, dI.HS.dt, dR.HS.dt), P = P))
    })
  }
  
  ################################## Optimize multi-host model ##################################
  my.SIRS.multi = vector('list', length(sites))
  params.multi = vector('list', length(sites))
  curr.type = 'Fitted' #the below is for a fitting model for multi-host transmission (no DHW or projection)
  
  for(i in 1:length(sites)){
    
    site.loop = sites[i]
    
    # site.loop = "mid" #for testing purposes
    # i = 1
    # site.loop = "near" #for testing purposes
    # i = 2
    # site.loop = "off" #for testing purposes
    # i = 3
    
    days = summary %>%
      # drop_na(days.survey) %>% # NOTE - area to return to after fixing backtracking with patient zero corals. we do want to drop this
      filter(site == site.loop) %>%
      pull(days.inf.site)
    
    #adjust days to chop off observations after onset of Degree Heating Weeks
    days = days[1:(length(days) - DHW.modifier)]
    
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
    
    N.LS.site = susceptible_ref %>%
      filter(Site == site.loop, Category == 'Low') %>%
      pull(tissue_ref)
    N.MS.site = susceptible_ref %>%
      filter(Site == site.loop, Category == 'Moderate') %>%
      pull(tissue_ref)
    N.HS.site = susceptible_ref %>%
      filter(Site == site.loop, Category == 'High') %>%
      pull(tissue_ref)
    
    cover.LS.site = susceptible_ref %>%
      filter(Site == site.loop, Category == 'Low') %>%
      pull(cover_ref)
    cover.MS.site = susceptible_ref %>%
      filter(Site == site.loop, Category == 'Moderate') %>%
      pull(cover_ref)
    cover.HS.site = susceptible_ref %>%
      filter(Site == site.loop, Category == 'High') %>%
      pull(cover_ref)
    
    #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
    LS.inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Low", Compartment == "Infected") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    MS.inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Moderate", Compartment == "Infected") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    HS.inftiss = obs.model %>%
      filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    inftiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
      slice(head(row_number(), n()-DHW.modifier)) %>%
      pull(tissue)    
    
    LS.remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Low", Compartment == "Dead") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    MS.remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Moderate", Compartment == "Dead") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    HS.remtiss = obs.model %>%
      filter(Site == site.loop, Category == "High", Compartment == "Dead") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)    
    remtiss = obs.model %>%
      filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)   
    
    #trim
    LS.inftiss <- LS.inftiss[first_valid_idx:length(LS.inftiss)]
    MS.inftiss <- MS.inftiss[first_valid_idx:length(MS.inftiss)]
    HS.inftiss <- HS.inftiss[first_valid_idx:length(HS.inftiss)]
    LS.remtiss <- LS.remtiss[first_valid_idx:length(LS.remtiss)]
    MS.remtiss <- MS.remtiss[first_valid_idx:length(MS.remtiss)]
    HS.remtiss <- HS.remtiss[first_valid_idx:length(HS.remtiss)]
    inftiss <- inftiss[first_valid_idx:length(inftiss)]
    remtiss <- remtiss[first_valid_idx:length(remtiss)]
    
    #initial conditions
    I.LS.tiss = 0
    S.LS.tiss = N.LS.site - I.LS.tiss
    R.LS.tiss = 0
    
    I.MS.tiss = 0
    S.MS.tiss = N.MS.site - I.MS.tiss
    R.MS.tiss = 0
    
    I.HS.tiss = HS.inftiss[1] #first non-NA & non-zero infection entry
    # I.HS.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
    S.HS.tiss = N.HS.site - I.HS.tiss 
    R.HS.tiss = 0
    
    P.tiss = I.HS.tiss
    
    ################################## Optimize multi-host model ##################################
    # Set up the data and initial conditions
    coraldata.tiss = list(LS.inftiss, MS.inftiss, HS.inftiss, LS.remtiss, MS.remtiss, HS.remtiss)
    initial_state.tiss = c(S.LS.tiss, I.LS.tiss, R.LS.tiss, S.MS.tiss, I.MS.tiss, R.MS.tiss, S.HS.tiss, I.HS.tiss, R.HS.tiss,
                           # P.tiss,
                           N.LS.site, N.MS.site, N.HS.site,
                           cover.site,
                           cover.LS.site, cover.MS.site, cover.HS.site,
                           lambda = lambda)
    
    # Define the objective function for optimization
    objective_function = function(params, data, time, initial_state){
      
      # #testing - for mid
      # betas.LS = 0.614603
      # gammas.LS = 1.774520
      # betas.MS = 0.014265
      # gammas.MS = 2.748830
      # betas.HS = 1.111736
      # gammas.HS = 0.805012
      # initial_state = initial_state.tiss
      # time = days.model
      # data = coraldata.tiss

      betas.LS = params[1]
      gammas.LS = params[2]
      betas.MS = params[3]
      gammas.MS = params[4]
      betas.HS = params[5]
      gammas.HS = params[6]
      
      SIR.out = data.frame(ode(c(S.LS = initial_state[1], I.LS = initial_state[2], R.LS = initial_state[3],
                                 S.MS = initial_state[4], I.MS = initial_state[5], R.MS = initial_state[6],
                                 S.HS = initial_state[7], I.HS = initial_state[8], R.HS = initial_state[9]),
                               time, SIR.multi, c(b.LS = betas.LS, g.LS = gammas.LS,
                                            b.MS = betas.MS, g.MS = gammas.MS,
                                            b.HS = betas.HS, g.HS = gammas.HS,
                                            N.LS = initial_state[10], N.MS = initial_state[11], N.HS = initial_state[12],
                                            C = initial_state[13],
                                            C.LS = initial_state[14], C.MS = initial_state[15], C.HS = initial_state[16],
                                            l = as.numeric(initial_state[17]))))
      
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
      
      #extract simulated values
      sim.LS.inf = SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'I.LS')]
      sim.MS.inf = SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'I.MS')]
      sim.HS.inf = SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'I.HS')]
      sim.LS.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I.LS')]
      sim.MS.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I.MS')]
      sim.HS.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I.HS')]
      
      sim.LS.rem = SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'R.LS')]
      sim.MS.rem = SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'R.MS')]
      sim.HS.rem = SIR.out[which(SIR.out$time %in% days.obs_trimmed), which(colnames(SIR.out) %in% 'R.HS')]
      sim.LS.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R.LS')]
      sim.MS.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R.MS')]
      sim.HS.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R.HS')]
      
      #extract observed values
      obs.LS.inf = unlist(data[[1]])[days.obs < time_cutoff] # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.MS.inf = unlist(data[[2]])[days.obs < time_cutoff] # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.HS.inf = unlist(data[[3]])[days.obs < time_cutoff] # NOTE - this is a bit hard-coded; refer to this line if there are bugs
      obs.LS.inf.total = unlist(data[1])
      obs.MS.inf.total = unlist(data[2])
      obs.HS.inf.total = unlist(data[3])
      
      obs.LS.rem = unlist(data[[4]])[days.obs < time_cutoff]
      obs.MS.rem = unlist(data[[5]])[days.obs < time_cutoff]
      obs.HS.rem = unlist(data[[6]])[days.obs < time_cutoff]
      obs.LS.rem.total = unlist(data[4])
      obs.MS.rem.total = unlist(data[5])
      obs.HS.rem.total = unlist(data[6])
      
      #bind for model evaluation and comparison with single-host version
      sim.inf = cbind(sim.LS.inf, sim.MS.inf, sim.HS.inf)
      sim.rem = cbind(sim.LS.rem, sim.MS.rem, sim.HS.rem)
      obs.inf = cbind(obs.LS.inf, obs.MS.inf, obs.HS.inf)
      obs.rem = cbind(obs.LS.rem, obs.MS.rem, obs.HS.rem)

      sim.inf.total = cbind(sim.LS.inf.total, sim.MS.inf.total, sim.HS.inf.total)
      sim.rem.total = cbind(sim.LS.rem.total, sim.MS.rem.total, sim.HS.rem.total)
      obs.inf.total = cbind(obs.LS.inf.total, obs.MS.inf.total, obs.HS.inf.total)
      obs.rem.total = cbind(obs.LS.rem.total, obs.MS.rem.total, obs.HS.rem.total)
      
      # # Sum for model evaluation and comparison with single-host version
      # #   NOTE - this may result in very different output than treating residuals as within-group
      # sim.inf = rowSums(cbind(sim.LS.inf, sim.MS.inf, sim.HS.inf), na.rm = TRUE)
      # sim.rem = rowSums(cbind(sim.LS.rem, sim.MS.rem, sim.HS.rem), na.rm = TRUE)
      # obs.inf = rowSums(cbind(obs.LS.inf, obs.MS.inf, obs.HS.inf), na.rm = TRUE)
      # obs.rem = rowSums(cbind(obs.LS.rem, obs.MS.rem, obs.HS.rem), na.rm = TRUE)
      # 
      # sim.inf.total = rowSums(cbind(sim.LS.inf.total, sim.MS.inf.total, sim.HS.inf.total), na.rm = TRUE)
      # sim.rem.total = rowSums(cbind(sim.LS.rem.total, sim.MS.rem.total, sim.HS.rem.total), na.rm = TRUE)
      # obs.inf.total = rowSums(cbind(obs.LS.inf.total, obs.MS.inf.total, obs.HS.inf.total), na.rm = TRUE)
      # obs.rem.total = rowSums(cbind(obs.LS.rem.total, obs.MS.rem.total, obs.HS.rem.total), na.rm = TRUE)
      
      # # NOTE - this as a version where rescaling was required, because infections were part of the fitting process (not *just* removal)
      # #use ratio of maximum removed value to maximum infected value to rescale removed. fits the 'I' curve better this way
      # # NOTE - doing this within each susceptibility category. should it be global?
      # max.obs.LS.inf = max(obs.LS.inf)
      # max.obs.MS.inf = max(obs.MS.inf)
      # max.obs.HS.inf = max(obs.HS.inf)
      # 
      # max.obs.LS.rem = max(obs.LS.rem)
      # max.obs.MS.rem = max(obs.MS.rem)
      # max.obs.HS.rem = max(obs.HS.rem)
      # 
      # rem.inf.LS.ratio = max.obs.LS.rem/max.obs.LS.inf
      # rem.inf.MS.ratio = max.obs.MS.rem/max.obs.MS.inf
      # rem.inf.HS.ratio = max.obs.HS.rem/max.obs.HS.inf
      # 
      # obs.LS.rem = obs.LS.rem / rem.inf.LS.ratio
      # obs.MS.rem = obs.MS.rem / rem.inf.MS.ratio
      # obs.HS.rem = obs.HS.rem / rem.inf.HS.ratio
      # 
      # sim.LS.rem = sim.LS.rem / rem.inf.LS.ratio
      # sim.MS.rem = sim.MS.rem / rem.inf.MS.ratio
      # sim.HS.rem = sim.HS.rem / rem.inf.HS.ratio
      # 
      # # #global version of rescaling
      # # max.obs.inf = max(obs.LS.inf, obs.MS.inf, obs.HS.inf)
      # # max.obs.rem = max(obs.LS.rem, obs.MS.rem, obs.HS.rem)
      # # 
      # # rem.inf.ratio = max.obs.rem/max.obs.inf
      # # 
      # # obs.LS.rem = obs.LS.rem / rem.inf.ratio
      # # obs.MS.rem = obs.MS.rem / rem.inf.ratio
      # # obs.HS.rem = obs.HS.rem / rem.inf.ratio
      # # 
      # # sim.LS.rem = sim.LS.rem / rem.inf.ratio
      # # sim.MS.rem = sim.MS.rem / rem.inf.ratio
      # # sim.HS.rem = sim.HS.rem / rem.inf.ratio
      # 
      # # # Calculate differences after normalization
      # # diff.LS.inf = (sim.LS.inf - obs.LS.inf)
      # # diff.MS.inf = (sim.MS.inf - obs.MS.inf)
      # # diff.HS.inf = (sim.HS.inf - obs.HS.inf)
      # # diff.LS.rem = (sim.LS.rem - obs.LS.rem)
      # # diff.MS.rem = (sim.MS.rem - obs.MS.rem)
      # # diff.HS.rem = (sim.HS.rem - obs.HS.rem)
      # # 
      # # #minimize using sum of absolute differences
      # # sum_squared_diff_I = sum(sum(abs(diff.LS.inf)) + sum(abs(diff.MS.inf)) +
      # #                            sum(abs(diff.HS.inf))) #can multiply this by 2 or similar to weight it extra
      # # 
      # # sum_squared_diff_R = sum(sum(abs(diff.LS.rem)) + sum(abs(diff.MS.rem)) +
      # #                            sum(abs(diff.HS.rem)))
      # # 
      # # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
      
      # Calculate residuals
      diff.inf = sim.inf - obs.inf
      diff.inf.total = sim.inf.total - obs.inf.total
      diff.rem = sim.rem - obs.rem
      diff.rem.total = sim.rem.total - obs.rem.total
      
      #Version using absolute residuals - can reference if having issues with sum-of-squares
      # #version that was constrained to pre-thermal stress onset
      # # Minimize using sum of absolute residuals
      # # sum_absolute_diff_I.abs = sum(abs(diff_inf))
      # sum_absolute_diff_R.abs = sum(abs(diff_rem))
      # # sum_diff.abs = sum_absolute_diff_I + sum_absolute_diff_R
      # sum_diff.abs = sum_absolute_diff_R.abs
      sum_absolute_diff_I.abs.total = sum(abs(diff.inf.total))
      sum_absolute_diff_R.abs.total = sum(abs(diff.rem.total))
      # sum_diff.abs.total = sum_absolute_diff_R.abs.total
      sum_diff.abs.total = sum_absolute_diff_I.abs.total + sum_absolute_diff_R.abs.total #test to better constrain fit for nearshore
      
      # #Version where I was including the infected compartment in the fit, and also summing squares within groups before global sum. return to this if any issues
      # #minimize using sum of squared residuals
      # sum_squared_diff_I = sum(sum(diff.LS.inf^2) + sum(diff.MS.inf^2) +
      #                            sum(diff.HS.inf^2))
      # sum_squared_diff_R = sum(sum(diff.LS.rem^2) + sum(diff.MS.rem^2) +
      #                            sum(diff.HS.rem^2))
      # sum_diff = sum_squared_diff_I + sum_squared_diff_R
      
      #minimize using sum of squared residuals
      sum_diff = sum(diff.rem^2)
      sum_diff.total = sum(diff.rem.total^2)
      
      # # NOTE - see Kalizhanova et al. 2024 (TB SIR) for other error assessments - including mean absolute error (MAE)
      # # Total Sum of Squares (TSS) for removal only
      # mean_obs_rem <- obs.rem %>% mean(na.rm = TRUE) # Compute mean of observed removals
      # tss_rem = sum((obs.rem - mean_obs_rem)^2)   # Sum of squared differences from mean
      # mean_obs_rem.total <- obs.rem.total %>% mean(na.rm = TRUE)
      # tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
      # 
      # #R-squared
      # r_squared_rem = 1 - (sum_diff / tss_rem)
      # r_squared_rem.total = 1 - (sum_diff.total / tss_rem.total)
      # 
      # error_eval <<- error_eval %>%
      #   mutate(
      #     SSR = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', sum_diff.total, SSR),
      #     TSS = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', tss_rem.total, TSS),
      #     R_squared = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', r_squared_rem.total, R_squared),
      #     
      #     # Update list-columns with vectors
      #     sim_inf = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(sim.inf.total), sim_inf),
      #     sim_rem = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(sim.rem.total), sim_rem),
      #     obs_inf = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(obs.inf.total), obs_inf),
      #     obs_rem = if_else(site == site.loop & host == curr.host & type == curr.type & wave == 'Full', list(obs.rem.total), obs_rem)
      #   )
      
      # # Debugging checks
      # print(table(error_eval$wave))
      # print(any(error_eval$wave == 'Full'))  # Should return TRUE
      # print(str(sim.inf.total))  # Should show a valid list/vector
      # print(str(sim.rem.total))  # Should show a valid list/vector
      
      # # Updating error_eval with proper handling of list-columns
      # error_eval <<- error_eval %>%
      #   mutate(
      #     SSR = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ sum_diff.total,
      #       TRUE ~ SSR
      #     ),
      #     TSS = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ tss_rem.total,
      #       TRUE ~ TSS
      #     ),
      #     R_squared = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ r_squared_rem.total,
      #       TRUE ~ R_squared
      #     ),
      #     
      #     # Update list-columns with vectors
      #     sim_inf = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ list(sim.inf.total),
      #       TRUE ~ sim_inf
      #     ),
      #     sim_rem = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ list(sim.rem.total),
      #       TRUE ~ sim_rem
      #     ),
      #     obs_inf = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ list(obs.inf.total),
      #       TRUE ~ obs_inf
      #     ),
      #     obs_rem = case_when(
      #       site == site.loop & host == curr.host & type == curr.type & wave == 'Full' ~ list(obs.rem.total),
      #       TRUE ~ obs_rem
      #     )
      #   )
          
      return(sum_diff.abs.total) #return only the residual metric for the epidemic wave being fit to
      # return(sum_diff.total) #return only the residual metric for the epidemic wave being fit to
    }
    
    ############################## OPTIMIZE PARAMETERS ############################################################
    
    # NOTE - TEMP
    # uniform or no?
    lower_bounds.tiss = c(0, 0, 0, 0, 0, 0)  # Lower bounds for betas and gammas - maybe more relaxed?
    # upper_bounds.tiss = c(0.09/N.LS.site, 0.11, 0.15/N.MS.site, 0.16, 2.0/N.HS.site, 1.0)  # Upper bounds for betas and gammas
    # upper_bounds.tiss = c(0.15/N.LS.site, 0.10, 0.15/N.MS.site, 0.10, 1.5/N.HS.site, 1.5)  # Upper bounds for betas and gammas
    # upper_bounds.tiss = c(0.003/N.LS.site, 0.01, 0.01/N.MS.site, 0.02, 0.15/N.HS.site, 0.10)  # Upper bounds for betas and gammas

    # upper_bounds.tiss = c(1/N.LS.site, 1, 1/N.MS.site, 1, 1.5/N.HS.site, 1.5)  # Upper bounds for betas and gammas
    upper_bounds.tiss = c(4, 4, 4, 4, 4, 4)  # Upper bounds for betas and gammas

    # # NOTE - TEMP
    # STOPPING POINT - 12 feb 2025; not sure what on earth is going wrong but it seems like sum of squares (sum_diff) is not porting between functions correctly...need to make sure it is, so I can properly record R-squared values. and multi-host model should actually have a quite good R-squared so something is wonky
    #                                 - this likely applies to the single-host model as well unfortunately
    # lower_bounds.tiss = c(0.06, 0.08, 0.29, 2.57, 3.13, 3.49)  # Lower bounds for betas and gammas - maybe more relaxed?
    # upper_bounds.tiss = c(0.06, 0.08, 0.29, 2.57, 3.13, 3.49)  # Upper bounds for betas and gammas
    
    control = list(itermax = 300)  # Maximum number of iterations. 200 is default
    
    # Run the optimization
    result.tiss = DEoptim(fn = objective_function, lower = lower_bounds.tiss, upper = upper_bounds.tiss,
                          data = coraldata.tiss, time = days.model, initial_state = initial_state.tiss,
                          control = control)
    
    # #test based on 'SIR-testing.R' script which produced successful outbreak
    # output.tester = data.frame(ode(c(S.LS = S.LS.tiss, I.LS = I.LS.tiss, R.LS = R.LS.tiss,
    #                           S.MS = S.MS.tiss, I.MS = I.MS.tiss, R.MS = R.MS.tiss,
    #                           S.HS = S.HS.tiss, I.HS = I.HS.tiss, R.HS = R.HS.tiss),
    #                         time, SIR, c(b.LS = 0.0001149425, g.LS = 0.01,
    #                                      b.MS = 0.00009009009, g.MS = 0.02,
    #                                      b.HS = 0.004237288, g.HS = 0.1,
    #                                      N.LS = N.LS.site, N.MS = N.MS.site, N.HS = N.HS.site,
    #                                      C = 0.2,
    #                                      l = lambda)))
  
    # 
    # output.tester = pivot_longer(output.tester, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>%
    #   mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    #   mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
    #                               ifelse(Compartment == 'I.', 'Infected',
    #                                      ifelse(Compartment == 'R.', 'Dead', Compartment))))
    # colnames(output.tester)[1] = 'days'
    # colnames(output.tester)[4] = 'prop'
    # 
    # p.fit.tester = ggplot(data = output.tester, aes(days, prop, colour = Compartment, linetype = Category)) +
    #   xlab("Day of observation period") +
    #   ylab("Proportion of tissue") +
    #   geom_line() +
    #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    #   theme_classic()
    # p.fit.tester
    
    # Extract the optimized parameters
    optimized_params.tiss = result.tiss$optim$bestmem  # 411: turn this ON for tissue!
    
    # Print the optimized parameters
    min.beta.LS.tiss = as.numeric(optimized_params.tiss[1])
    min.gamma.LS.tiss = as.numeric(optimized_params.tiss[2])
    min.beta.LS.tiss.adj = min.beta.LS.tiss * (1 / (1 + exp(-lambda.LS * (cover.LS.site))) + offset.LS)
    R0.LS = min.beta.LS.tiss.adj / min.gamma.LS.tiss
    
    min.beta.MS.tiss = as.numeric(optimized_params.tiss[3])
    min.gamma.MS.tiss = as.numeric(optimized_params.tiss[4])
    min.beta.MS.tiss.adj = min.beta.MS.tiss * (1 / (1 + exp(-lambda.MS * (cover.MS.site))) + offset.MS)
    R0.MS = min.beta.MS.tiss.adj / min.gamma.MS.tiss
    
    min.beta.HS.tiss = as.numeric(optimized_params.tiss[5])
    min.gamma.HS.tiss = as.numeric(optimized_params.tiss[6])
    min.beta.HS.tiss.adj = min.beta.HS.tiss * (1 / (1 + exp(-lambda.HS * (cover.HS.site))) + offset.HS)
    R0.HS = min.beta.HS.tiss.adj / min.gamma.HS.tiss
    
    cat("Optimized Tissue Model Parameters for", site.loop, " site:\n\n")
    
    # LS (Low Susceptibility) parameters
    cat("Low Susceptibility (LS) Parameters:\n")
    cat("  Beta LS:", min.beta.LS.tiss, "\n")
    cat("  Cover-adjusted beta LS:", min.beta.LS.tiss.adj, "\n")
    cat("  Gamma LS:", min.gamma.LS.tiss, "\n\n")
    
    # MS (Medium Susceptibility) parameters
    cat("Medium Susceptibility (MS) Parameters:\n")
    cat("  Beta MS:", min.beta.MS.tiss, "\n")
    cat("  Cover-adjusted beta MS:", min.beta.MS.tiss.adj, "\n")
    cat("  Gamma MS:", min.gamma.MS.tiss, "\n\n")
    
    # HS (High Susceptibility) parameters
    cat("High Susceptibility (HS) Parameters:\n")
    cat("  Beta HS:", min.beta.HS.tiss, "\n")
    cat("  Cover-adjusted beta HS:", min.beta.HS.tiss.adj, "\n")
    cat("  Gamma HS:", min.gamma.HS.tiss, "\n\n")
    
    # R0
    cat("Basic Reproduction Number (R0) Parameters:\n")
    cat("  R0 LS:", R0.LS, "\n")
    cat("  R0 MS:", R0.MS, "\n")
    cat("  R0 HS:", R0.HS, "\n")
    
    params.multi[[i]] = c(min.beta.LS.tiss, min.beta.LS.tiss.adj, min.gamma.LS.tiss,
                          min.beta.MS.tiss, min.beta.MS.tiss.adj, min.gamma.MS.tiss,
                          min.beta.HS.tiss, min.beta.HS.tiss.adj, min.gamma.HS.tiss,
                          R0.LS, R0.MS, R0.HS,
                          cover.site,
                          cover.LS.site, cover.MS.site, cover.HS.site)
    
    #simulation using initial state variables and best-fit beta/gamma parameters
    SIR.out.tiss = data.frame(ode(c(S.LS = initial_state.tiss[1], I.LS = initial_state.tiss[2], R.LS = initial_state.tiss[3],
                                    S.MS = initial_state.tiss[4], I.MS = initial_state.tiss[5], R.MS = initial_state.tiss[6],
                                    S.HS = initial_state.tiss[7], I.HS = initial_state.tiss[8], R.HS = initial_state.tiss[9]),
                                    # S.HS = initial_state[7], I.HS = initial_state[8], R.HS = initial_state[9],
                                    # P = initial_state[10]),
                                    days.model, SIR.multi, c(b.LS = min.beta.LS.tiss, g.LS = min.gamma.LS.tiss,
                                               b.MS = min.beta.MS.tiss, g.MS = min.gamma.MS.tiss,
                                               b.HS = min.beta.HS.tiss, g.HS = min.gamma.HS.tiss,
                                               N.LS = initial_state.tiss[10], N.MS = initial_state.tiss[11], N.HS = initial_state.tiss[12],
                                               C = initial_state.tiss[13],
                                               C.LS = initial_state.tiss[14], C.MS = initial_state.tiss[15], C.HS = initial_state.tiss[16],
                                               l = as.numeric(initial_state.tiss[17]))))
    
    my.SIRS.multi[[i]] = SIR.out.tiss
  }
  
  ################################## Save output ##################################
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "multi_SIR_workspace.RData"))