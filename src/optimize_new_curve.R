  
  
################################## alpha / k parameter plotting sandbox ##################################

  # CC = seq(0.001,1,.001)
  # lmbds = c(1,10,100,1000)
  # 
  # plot(CC,CC*CC, type = 'n')
  # for (l in 1:length(lmbds)) {
  #   mods.new <- 1-exp(-lmbds[l]*CC)
  #   lines(CC, mods.new, type = 'l')
  # }
  
  CC = seq(0.001,1,.001)
  a = seq(0,1,0.01) #alpha (weight of coral cover)
  # a = 0 #best_alpha
  k = 100 #shape of curve?
  
  plot(CC,CC*CC, type = 'n')
  for (i in 1:length(a)) {
    mods.1 <- (1-a[i]) + a[i]*((1-exp(-k*CC))/(1-exp(-k)))
    lines(CC, mods.1, type = 'l', col = 'red', lwd = 1)
  }
  
  # CC = seq(0.001,1,.001)
  # a = seq(0,1,0.1)
  # k = 5
  # 
  # plot(CC,CC*CC, type = 'n')
  # for (i in 1:length(a)) {
  #   mods.2 <- (1-a[i])*CC + a[i]*((1-exp(-k*CC))/(1-exp(-k)))
  #   lines(CC, mods.2, type = 'l', col = 'red', lwd = 1)
  # }

  
  
  ################################## single-host projection optimization ##################################
  library(dplyr)
  
  SIR_project = function(t,y,p){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p),{

      transmission_modifier = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C)) / (1 - exp(-k_val)))
      
      dS.dt = -b * S * I / N * transmission_modifier
      dI.dt = b * S * I / N * transmission_modifier - g * I
      dR.dt = g * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  # Define the search space
  alpha_values <- seq(0, 1, length.out = 1000)
  k_val = 3
  
  best_r_squared <- -Inf
  best_alpha <- NA
  
  for (i in alpha_values) {
    
    alpha_val = alpha_values[i] # NOTE - bad and hard-coded way to feed parameter into the SIR function - should look at this further
    
    # Run the model
    output.basic.offshore.transfer = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                                    days.model.offshore, SIR_project, c(b = beta.nearshore, g = gamma.nearshore,
                                                                                N = N.offshore,
                                                                                C = cover.offshore,
                                                                                l = lambda)))
    
    
    # Compute R-squared
    sim.rem.total = output.basic.offshore.transfer[which(output.basic.offshore.transfer$time %in% days.obs),
                                                   which(colnames(output.basic.offshore.transfer) %in% c('R'))]
    
    if (length(obs.rem.total) > length(sim.rem.total)) {
      obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
    }
    
    diff.rem.total = (sim.rem.total - obs.rem.total)
    sum_diff.total = sum(diff.rem.total^2, na.rm = TRUE) # NOTE - had to account for NAs when fit is really bad
    mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
    tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
    r_squared.near.to.off.basic = 1 - (sum_diff.total / tss_rem.total)
    
    # Check if this is the best so far
    if (r_squared.near.to.off.basic > best_r_squared) {
      best_r_squared <- r_squared.near.to.off.basic
      best_alpha <- i
    }
  }
  
  # Print the best parameters
  print(best_alpha)
  print(best_r_squared)
  
  
  
  
  
  
  
  ################################## multi-host projection optimization ##################################
  
  library(dplyr)
  
  SIR.multi_project = function(t,y,p){
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
    }
    with(as.list(p),{
      P = (I.LS + I.MS + I.HS)
      
      transmission_modifier.LS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.LS)) / (1 - exp(-k_val)))
      transmission_modifier.MS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.MS)) / (1 - exp(-k_val)))
      transmission_modifier.HS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.HS)) / (1 - exp(-k_val)))
      
      dS.LS.dt = -b.LS*S.LS*(P) / N.LS * transmission_modifier.LS
      dI.LS.dt = b.LS*S.LS*(P) / N.LS * transmission_modifier.LS - g.LS*I.LS
      dR.LS.dt = g.LS*I.LS
      
      dS.MS.dt = -b.MS*S.MS*(P) / N.MS * transmission_modifier.MS
      dI.MS.dt = b.MS*S.MS*(P) / N.MS * transmission_modifier.MS - g.MS*I.MS
      dR.MS.dt = g.MS*I.MS
      
      dS.HS.dt = -b.HS*S.HS*(P) / N.HS * transmission_modifier.HS
      dI.HS.dt = b.HS*S.HS*(P) / N.HS * transmission_modifier.HS - g.HS*I.HS
      dR.HS.dt = g.HS*I.HS
      
      return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, dS.HS.dt, dI.HS.dt, dR.HS.dt), P = P))
    })
  }
  
  
  # Define the search space
  alpha_values <- seq(0.1, 1, length.out = 1000)
  # alpha_values <- seq(0, 1, length.out = 1000)
  # alpha_values <- seq(0, 0, length.out = 10)
  k_val = 3
  
  best_r_squared <- -Inf
  best_lambda <- c(lambda.LS = NA, lambda.MS = NA, lambda.HS = NA)
  
  site.loop = 'Offshore'
  curr.site = 'off'
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  
  for (i in alpha_values) {
    
    alpha_val = i # NOTE - bad and hard-coded way to feed parameter into the SIR function - should look at this further
    
    # Run the model
    output.near.to.off.multi = data.frame(ode(c(S.LS = S.LS.offshore, I.LS = I.LS.offshore, R.LS = R.LS.offshore,
                                                S.MS = S.MS.offshore, I.MS = I.MS.offshore, R.MS = R.MS.offshore,
                                                S.HS = S.HS.offshore, I.HS = I.HS.offshore, R.HS = R.HS.offshore),
                                              days.model.offshore, SIR.multi_project, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                                                                b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                                                                b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                                                                N.LS = N.LS.offshore, N.MS = N.MS.offshore, N.HS = N.HS.offshore,
                                                                                C = cover.offshore,
                                                                                C.LS = cover.offshore.LS, C.MS = cover.offshore.MS, C.HS = cover.offshore.HS,
                                                                                l = lambda)))
    
    # Compute R-squared
    sim.rem.total = rowSums(output.near.to.off.multi[which(output.near.to.off.multi$time %in% days.obs),
                                                     which(colnames(output.near.to.off.multi) %in% c('R.LS', 'R.MS', 'R.HS'))],
                            na.rm = TRUE)
    
    if (length(obs.rem.total) > length(sim.rem.total)) {
      obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
    }
    
    diff.rem.total = (sim.rem.total - obs.rem.total)
    sum_diff.total = sum(diff.rem.total^2, na.rm = TRUE) # NOTE - had to account for NAs when fit is really bad
    mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
    tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
    r_squared.near.to.off.multi = 1 - (sum_diff.total / tss_rem.total)
    
    # Check if this is the best so far
    if (r_squared.near.to.off.multi > best_r_squared) {
      best_r_squared <- r_squared.near.to.off.multi
      best_alpha <- i
    }
  }
  
  # Print the best parameters
  print(best_alpha)
  print(best_r_squared)
  