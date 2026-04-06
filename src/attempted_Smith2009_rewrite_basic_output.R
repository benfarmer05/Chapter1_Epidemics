
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

################################## define cover params ##################################

#with effect of coral cover
alpha_val = 0.13
# 0.01 produces too much outbreak at offshore
# 0.05 too I think
# 0.1 produces just a little too much outbreak, but R-squared is very bad (-0.57)...interesting
# 0.13 looks really good!! It is slightly out of phase but about on the money
# 0.15 produces not quite enough outbreak (R-squared = 0.28)
# 0.2 produces no outbreak at all. R-squared is -0.1 (technically better than with a = 0.1, but only because values near 0 are technically closer than a slightly poorly predicted actual outbreak)
K_val = 172

q_val = 0.62
# q_val = 0

################################## Model: single-host ##################################
SIR = function(t,y,p){ # 'p' is parameters or params
  {
    S = y[1]
    I = y[2]
    R = y[3]
  }
  with(as.list(p),{
    
    # #null conditions
    # transmission_modifier = 1
    
    contact_rate = (K^q_val) * S * I / (N^q_val)
    
    #nonlinear scaling range as in Smith et al. (2009)
    dS.dt = -b * contact_rate
    dI.dt = b * contact_rate - g * I
    dR.dt = g * I
    
    # #frequency-dependent (and can be w/ NL dens contact rate)
    # dS.dt = -b * S * I / N * transmission_modifier
    # dI.dt = b * S * I / N * transmission_modifier - g * I
    # dR.dt = g * I
    
    # #density-dependent with NL dens contact rate
    # dS.dt = -b * S * I * transmission_modifier
    # dI.dt = b * S * I * transmission_modifier - g * I
    # dR.dt = g * I
    
    # #density-dependent
    # dS.dt = -b * S * I
    # dI.dt = b * S * I - g * I
    # dR.dt = g * I
    
    return(list(c(dS.dt, dI.dt, dR.dt)))
  })
}

################################## Optimize full outbreak ##################################
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
                                  q = q_val,
                                  K = K_val)))
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
    diff.inf.total = (sim.inf.total - obs.inf.total)
    diff.rem.total = (sim.rem.total - obs.rem.total)
    
    #aggregate sum of absolute residuals
    # NOTE - this is not sum of squares and should be clearly stated/defended in the manuscript if used
    sum_abs_diff_I.total = sum(sum(abs(diff.inf.total))) #can multiply this by 2 or similar to weight it extra
    sum_abs_diff_R.total = sum(sum(abs(diff.rem.total)))
    # sum_diff = sum_abs_diff_I + sum_abs_diff_R #this is the version where infections AND removal are fitted, not *just* removal
    sum_diff.abs.total = sum_abs_diff_R.total
    
    #minimize using sum of squared residuals
    sum_diff.total = sum(diff.rem.total^2)
    
    # # NOTE - see Kalizhanova et al. 2024 (TB SIR) for other error assessments - including mean absolute error (MAE)
    # # Total Sum of Squares (TSS) for removal only
    # mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
    # tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
    # 
    # #R-squared
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
    
    return(sum_diff.abs.total) #return only the residual metric for the epidemic wave being fit to
  }
  
  # uniform
  
  # # TEST
  # lower_bounds.tiss = c(0, 0, lambda.modifier)  #lower bounds for beta, gamma and lambda
  # upper_bounds.tiss = c(0.1, 3, lambda.modifier)  #upper bounds for beta, gamma and lambda
  # # TEST
  
  lower_bounds.tiss = c(0, 0, lambda.modifier)  #lower bounds for beta, gamma and lambda
  upper_bounds.tiss = c(4, 4, lambda.modifier)  #upper bounds for beta, gamma and lambda
  
  control = list(itermax = 200)  # Maximum number of iterations. 200 is default
  
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
                                                   q = q_val,
                                                   K = K_val)))
  my.SIRS.basic.full[[i]] = SIR.out.tiss
}

################################## Save output ##################################

# #pass workspace to downstream script
# save.image(file = here("output", "basic_SIR_workspace.RData"))