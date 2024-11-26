  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(DEoptim)
  library(deSolve)
  
  #import workspace from upstream script
  load(here("output/plots_basic_workspace.RData"))
  
  susceptible_ref = susceptible_ref %>%
    mutate(Site = case_when(
      Site == "Offshore" ~ "off",
      Site == 'Midchannel' ~ 'mid',
      Site == 'Nearshore' ~ 'near'
    ))
  
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
  
  # 0/0/0 experiment
  #   - haven't tested these yet. NOTE / STOPPING POINT - 8 OCT 2024:
  #       what exactly are the null conditions for effect of coral cover? should figure that out
  lambda.LS = 0.0
  lambda.MS = 0.0
  lambda.HS = 0.0
  
  # # 05/125/5 experiment
  # # comparing to 2/2/5 experiment above, trying something new
  # lambda.LS = 0.5 # overpredicted off.to.near (still kind of bad but better); no change near.to.off (still quite bad)
  # lambda.MS = 1.25 # quite good off.to.near (closer to reality); no change near.to.off (still quite bad)
  # lambda.HS = 5.0 # underpredicted off.to.near (even worse); about right near.to.off (closer to reality)
  
  # # 01/125/15 experiment
  # # comparing to 05/125/5 experiment above
  # lambda.LS = 0.1 # overpredicted off.to.near (got worse); no change near.to.off (still quite bad)
  # lambda.MS = 1.25 # overpredicted off.to.near (a bit worse but still good); underpredicted near.to.off (bad but better)
  # lambda.HS = 15.0 # underpredicted off.to.near (still bad but a lot better); shape changed near.to.off (still pretty good)
  
  # # 0/05/20 experiment
  # # comparing to 01/125/15 experiment above
  # lambda.LS = 0 # overpredicted off.to.near (still bad and got worse); no change near.to.off (still quite bad)
  # lambda.MS = 0.5 # overpredicted off.to.near (not pretty good but got worse); underpredicted near.to.of (still bad but slightly better)
  # lambda.HS = 20.0 # underpredicted off.to.near (off but still getting better); shape changed near.to.off (still pretty good)
  
  offset.LS = 1 - 1 / (1 + exp(-lambda.LS * 1.0))
  offset.MS = 1 - 1 / (1 + exp(-lambda.MS * 1.0))
  offset.HS = 1 - 1 / (1 + exp(-lambda.HS * 1.0))
  
  SIR = function(t,y,p){
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
      
      P = y[10]
    }
    with(as.list(p),{
      # P = (I.LS + I.MS + I.HS)
      # dP = (I.LS*prop.LS + I.MS*prop.MS + I.HS*prop.HS)
      
      # transmission_modifier.LS = (1 / (1 + exp(-l * (C.LS))) + offset)
      # transmission_modifier.MS = (1 / (1 + exp(-l * (C.MS))) + offset)
      # transmission_modifier.HS = (1 / (1 + exp(-l * (C.HS))) + offset)
      transmission_modifier.LS = (1 / (1 + exp(-lambda.LS * (C.LS))) + offset.LS)
      transmission_modifier.MS = (1 / (1 + exp(-lambda.MS * (C.MS))) + offset.MS)
      transmission_modifier.HS = (1 / (1 + exp(-lambda.HS * (C.HS))) + offset.HS)
      
      dS.LS.dt = -b.LS*S.LS*(P) / N.LS * transmission_modifier.LS
      dI.LS.dt = b.LS*S.LS*(P) / N.LS * transmission_modifier.LS - g.LS*I.LS
      dR.LS.dt = g.LS*I.LS
  
      dS.MS.dt = -b.MS*S.MS*(P) / N.MS * transmission_modifier.MS
      dI.MS.dt = b.MS*S.MS*(P) / N.MS * transmission_modifier.MS - g.MS*I.MS
      dR.MS.dt = g.MS*I.MS
  
      dS.HS.dt = -b.HS*S.HS*(P) / N.HS * transmission_modifier.HS
      dI.HS.dt = b.HS*S.HS*(P) / N.HS * transmission_modifier.HS - g.HS*I.HS
      dR.HS.dt = g.HS*I.HS
      
      if (any((S.LS + I.LS + R.LS > N.LS) | 
              (S.MS + I.MS + R.MS > N.MS) | 
              (S.HS + I.HS + R.HS > N.HS))) {
        warning("Sum of S, I, and R exceeds N for at least one group during fitting process.")
      }
      
      # Update dP using the sum of infected individuals
      dP.dt = (I.LS + I.MS + I.HS) - P
  
  return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, dS.HS.dt, dI.HS.dt, dR.HS.dt, dP.dt)))
      # return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, dS.HS.dt, dI.HS.dt, dR.HS.dt)))
    })
  }
  
  my.SIRS.multi = vector('list', length(sites))
  params.multi = vector('list', length(sites))
  
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
    
    ############################## OPTIMIZE PARAMETERS ############################################################
    # Set up the data and initial conditions
    coraldata.tiss = list(LS.inftiss, MS.inftiss, HS.inftiss, LS.remtiss, MS.remtiss, HS.remtiss)
    initial_state.tiss = c(S.LS.tiss, I.LS.tiss, R.LS.tiss, S.MS.tiss, I.MS.tiss, R.MS.tiss, S.HS.tiss, I.HS.tiss, R.HS.tiss,
                           P.tiss,
                           N.LS.site, N.MS.site, N.HS.site,
                           cover.site,
                           cover.LS.site, cover.MS.site, cover.HS.site,
                           lambda = lambda)
    
    # Define the objective function for optimization
    objective_function = function(params, data, time, initial_state){
      
      # #testing
      # betas.LS = 2.0
      # gammas.LS = 1.0
      # betas.MS = 2.0
      # gammas.MS = 1.0
      # betas.HS = 2.0
      # gammas.HS = 1.0
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
                                 S.HS = initial_state[7], I.HS = initial_state[8], R.HS = initial_state[9],
                                 P = initial_state[10]),
                               time, SIR, c(b.LS = betas.LS, g.LS = gammas.LS,
                                            b.MS = betas.MS, g.MS = gammas.MS,
                                            b.HS = betas.HS, g.HS = gammas.HS,
                                            N.LS = initial_state[11], N.MS = initial_state[12], N.HS = initial_state[13],
                                            C = initial_state[14],
                                            C.LS = initial_state[15], C.MS = initial_state[16], C.HS = initial_state[17],
                                            l = as.numeric(initial_state[18]))))
      
      #extract simulated values at time points matching observations
      sim.LS.inf = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I.LS')]
      sim.MS.inf = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I.MS')]
      sim.HS.inf = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I.HS')]
  
      sim.LS.rem = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R.LS')]
      sim.MS.rem = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R.MS')]
      sim.HS.rem = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R.HS')]
      
      #extract observed values
      obs.LS.inf = unlist(data[1])
      obs.MS.inf = unlist(data[2])
      obs.HS.inf = unlist(data[3])
  
      obs.LS.rem = unlist(data[4])
      obs.MS.rem = unlist(data[5])
      obs.HS.rem = unlist(data[6])
  
      # NOTE - if there is NO rescaling done, the fit to removal is actually excellent (when I set
      #   the 'polyp_SA' to 1e-4 at least, as well as DHW modifier to 6 and initial conditions as below). but it's the fit to infection
      #   that will be really interesting to get right as well.
      #     initial conditions mentioned:
      #       upper_bounds.tiss = c(0.15/N.LS.site, 0.10, 0.15/N.MS.site, 0.10, 1/N.HS.site, 1)  # Upper bounds for betas and gammas
  
      #use ratio of maximum removed value to maximum infected value to rescale removed. fits the 'I' curve better this way
      # NOTE - doing this within each susceptibility category. should it be global?
      max.obs.LS.inf = max(obs.LS.inf)
      max.obs.MS.inf = max(obs.MS.inf)
      max.obs.HS.inf = max(obs.HS.inf)
  
      max.obs.LS.rem = max(obs.LS.rem)
      max.obs.MS.rem = max(obs.MS.rem)
      max.obs.HS.rem = max(obs.HS.rem)
  
      rem.inf.LS.ratio = max.obs.LS.rem/max.obs.LS.inf
      rem.inf.MS.ratio = max.obs.MS.rem/max.obs.MS.inf
      rem.inf.HS.ratio = max.obs.HS.rem/max.obs.HS.inf
  
      obs.LS.rem = obs.LS.rem / rem.inf.LS.ratio
      obs.MS.rem = obs.MS.rem / rem.inf.MS.ratio
      obs.HS.rem = obs.HS.rem / rem.inf.HS.ratio
  
      sim.LS.rem = sim.LS.rem / rem.inf.LS.ratio
      sim.MS.rem = sim.MS.rem / rem.inf.MS.ratio
      sim.HS.rem = sim.HS.rem / rem.inf.HS.ratio
      
      # #global version of rescaling
      # max.obs.inf = max(obs.LS.inf, obs.MS.inf, obs.HS.inf)
      # max.obs.rem = max(obs.LS.rem, obs.MS.rem, obs.HS.rem)
      # 
      # rem.inf.ratio = max.obs.rem/max.obs.inf
      # 
      # obs.LS.rem = obs.LS.rem / rem.inf.ratio
      # obs.MS.rem = obs.MS.rem / rem.inf.ratio
      # obs.HS.rem = obs.HS.rem / rem.inf.ratio
      # 
      # sim.LS.rem = sim.LS.rem / rem.inf.ratio
      # sim.MS.rem = sim.MS.rem / rem.inf.ratio
      # sim.HS.rem = sim.HS.rem / rem.inf.ratio
      
      # # Calculate differences after normalization
      # diff.LS.inf = (sim.LS.inf - obs.LS.inf)
      # diff.MS.inf = (sim.MS.inf - obs.MS.inf)
      # diff.HS.inf = (sim.HS.inf - obs.HS.inf)
      # diff.LS.rem = (sim.LS.rem - obs.LS.rem)
      # diff.MS.rem = (sim.MS.rem - obs.MS.rem)
      # diff.HS.rem = (sim.HS.rem - obs.HS.rem)
      # 
      # #minimize using sum of absolute differences
      # sum_squared_diff_I = sum(sum(abs(diff.LS.inf)) + sum(abs(diff.MS.inf)) +
      #                            sum(abs(diff.HS.inf))) #can multiply this by 2 or similar to weight it extra
      # 
      # sum_squared_diff_R = sum(sum(abs(diff.LS.rem)) + sum(abs(diff.MS.rem)) +
      #                            sum(abs(diff.HS.rem)))
      # 
      # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
      
      # Vectorized calculation of differences after normalization
      diff_inf = cbind(
        (sim.LS.inf - obs.LS.inf),
        (sim.MS.inf - obs.MS.inf),
        (sim.HS.inf - obs.HS.inf)
      )
  
      diff_rem = cbind(
        (sim.LS.rem - obs.LS.rem),
        (sim.MS.rem - obs.MS.rem),
        (sim.HS.rem - obs.HS.rem)
      )
  
      # Minimize using sum of absolute differences
      sum_squared_diff_I = sum(abs(diff_inf))
      sum_squared_diff_R = sum(abs(diff_rem))
      sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
  
      # #minimize using sum of squared differences
      # sum_squared_diff_I = sum(sum(diff.LS.inf^2) + sum(diff.MS.inf^2) +
      #                            sum(diff.HS.inf^2))
      # sum_squared_diff_R = sum(sum(diff.LS.rem^2) + sum(diff.MS.rem^2) +
      #                            sum(diff.HS.rem^2))
      # sum_squared_diff = sum_squared_diff_I + sum_squared_diff_R
  
      return(sum_squared_diff) #minimized error
    }
  
    ############################## OPTIMIZE PARAMETERS ############################################################
    
    # uniform or no?
    lower_bounds.tiss = c(0, 0, 0, 0, 0, 0)  # Lower bounds for betas and gammas - maybe more relaxed?
    # upper_bounds.tiss = c(0.09/N.LS.site, 0.11, 0.15/N.MS.site, 0.16, 2.0/N.HS.site, 1.0)  # Upper bounds for betas and gammas
    # upper_bounds.tiss = c(0.15/N.LS.site, 0.10, 0.15/N.MS.site, 0.10, 1.5/N.HS.site, 1.5)  # Upper bounds for betas and gammas
    # upper_bounds.tiss = c(0.003/N.LS.site, 0.01, 0.01/N.MS.site, 0.02, 0.15/N.HS.site, 0.10)  # Upper bounds for betas and gammas
    
    # upper_bounds.tiss = c(1/N.LS.site, 1, 1/N.MS.site, 1, 1.5/N.HS.site, 1.5)  # Upper bounds for betas and gammas
    upper_bounds.tiss = c(4, 4, 4, 4, 4, 4)  # Upper bounds for betas and gammas
    
    
    control = list(itermax = 200)  # Maximum number of iterations. 200 is default
    
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
                                    S.HS = initial_state.tiss[7], I.HS = initial_state.tiss[8], R.HS = initial_state.tiss[9],
                                    P = initial_state.tiss[10]),
                                  days.model, SIR, c(b.LS = min.beta.LS.tiss, g.LS = min.gamma.LS.tiss,
                                               b.MS = min.beta.MS.tiss, g.MS = min.gamma.MS.tiss,
                                               b.HS = min.beta.HS.tiss, g.HS = min.gamma.HS.tiss,
                                               N.LS = initial_state.tiss[11], N.MS = initial_state.tiss[12], N.HS = initial_state.tiss[13],
                                               C = initial_state.tiss[14],
                                               C.LS = initial_state.tiss[15], C.MS = initial_state.tiss[16], C.HS = initial_state.tiss[17],
                                               l = as.numeric(initial_state.tiss[18]))))  
    
    my.SIRS.multi[[i]] = SIR.out.tiss
    
  }
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "multi_SIR_workspace.RData"))
  