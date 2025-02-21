    
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(dplyr)
  library(purrr)
  
  ################################## Set-up ##################################
  
  #import workspace from upstream script
  load(here("output/plots_multi_workspace.RData"))
  # load(here("output/plots_multi_workspace_lower_start.RData"))
  # load(here("output/plots_multi_workspace_betterproj.RData"))
  
  # Function to calculate R-squared, RMSE, NRMSE, and sMAPE
  calculate_metrics <- function(output, obs_data, site, host, type, wave) {
    
    # **New Change: Ensure `days.obs` is filtered for the current site**
    days.obs <- days_sites %>%
      filter(site == !!site) %>%  # Dynamically filter by site
      pull(days.obs) %>%
      unlist()
    
    # Sum 'Dead' tissue values for observed and simulated data
    sim.rem.total <- output %>%
      filter(Compartment == "Dead", days.model %in% days.obs) %>%
      group_by(days.model) %>%
      summarize(total_tissue = sum(tissue, na.rm = TRUE), .groups = "drop") %>%
      pull(total_tissue)
    
    obs.rem.total <- obs_data %>%
      filter(Site == site, Category == "Total", Compartment == "Dead") %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)
    
    # Ensure obs and sim vectors are aligned in length
    min_length <- min(length(obs.rem.total), length(sim.rem.total))
    obs.rem.total <- tail(obs.rem.total, min_length)
    sim.rem.total <- tail(sim.rem.total, min_length)
    
    # Calculate differences and error metrics
    diff.rem.total <- sim.rem.total - obs.rem.total
    sum_diff.total <- sum(diff.rem.total^2, na.rm = TRUE)
    mean_obs_rem.total <- mean(obs.rem.total, na.rm = TRUE)
    tss_rem.total <- sum((obs.rem.total - mean_obs_rem.total)^2, na.rm = TRUE)
    
    # R-squared and RMSE
    r_squared <- 1 - (sum_diff.total / tss_rem.total)
    rmse <- sqrt(mean(diff.rem.total^2, na.rm = TRUE))
    
    # NRMSE calculations
    range_obs <- max(obs.rem.total, na.rm = TRUE) - min(obs.rem.total, na.rm = TRUE)
    mean_obs <- mean(obs.rem.total, na.rm = TRUE)
    nrmse_range <- rmse / range_obs
    nrmse_mean <- rmse / mean_obs
    
    # sMAPE calculation
    smape <- mean(2 * abs(sim.rem.total - obs.rem.total) / (abs(obs.rem.total) + abs(sim.rem.total)), na.rm = TRUE) * 100
    
    return(tibble(site = site, host = host, type = type, wave = wave, 
                  R_squared = r_squared, RMSE = rmse, NRMSE_range = nrmse_range, 
                  NRMSE_mean = nrmse_mean, sMAPE = smape))
  }
  
  
  # List of all models with corresponding metadata
  model_list <- list(
    # Nearshore models
    list(output = output.basic.nearshore, site = "near", host = "Single", type = "DHW", wave = "Pre-heat"),
    list(output = output.basic.nearshore.full, site = "near", host = "Single", type = "Fitted", wave = "Full"),
    list(output = output.basic.nearshore.DHW, site = "near", host = "Single", type = "DHW", wave = "Full"),
    list(output = output.nearshore, site = "near", host = "Multi", type = "Fitted", wave = "Full"),
    
    # Midchannel models
    list(output = output.basic.midchannel, site = "mid", host = "Single", type = "DHW", wave = "Pre-heat"),
    list(output = output.basic.midchannel.full, site = "mid", host = "Single", type = "Fitted", wave = "Full"),
    list(output = output.basic.midchannel.DHW, site = "mid", host = "Single", type = "DHW", wave = "Full"),
    list(output = output.midchannel, site = "mid", host = "Multi", type = "Fitted", wave = "Full"),
    
    # Offshore models
    list(output = output.basic.offshore, site = "off", host = "Single", type = "DHW", wave = "Pre-heat"),
    list(output = output.basic.offshore.full, site = "off", host = "Single", type = "Fitted", wave = "Full"),
    list(output = output.basic.offshore.DHW, site = "off", host = "Single", type = "DHW", wave = "Full"),
    list(output = output.offshore, site = "off", host = "Multi", type = "Fitted", wave = "Full"),
    
    #Projected models
    list(output = output.basic.offshore.transfer, site = "off", host = "Single", type = "Projected", wave = "Full"),
    list(output = output.near.to.off.multi, site = "off", host = "Multi", type = "Projected", wave = "Full"),
    list(output = output.basic.nearshore.transfer, site = "near", host = "Single", type = "Projected", wave = "Full"),
    list(output = output.off.to.near.multi, site = "near", host = "Multi", type = "Projected", wave = "Full")
    
  )
  
  # Apply function to all models and combine results
  error_metrics <- map_dfr(model_list, ~calculate_metrics(.x$output, obs.model, .x$site, .x$host, .x$type, .x$wave))
  
  # # Update error_eval table
  # error_eval2 <- error_eval %>%
  #   left_join(error_metrics, by = c("site", "host", "type", "wave")) %>%
  #   mutate(R_squared = coalesce(R_squared.y, R_squared.x),
  #          RMSE = coalesce(RMSE, 0)) %>% 
  #   select(-R_squared.x, -R_squared.y)
  
  # # Update error_eval table
  # error_eval2 <- error_eval %>%
  #   mutate(
  #     R_squared = error_metrics$R_squared,
  #     RMSE = error_metrics$RMSE,
  #     NRMSE_range = error_metrics$NRMSE_range,
  #     NRMSE_mean = error_metrics$NRMSE_mean,
  #     sMAPE = error_metrics$sMAPE
  #   )
  
  ################################## Save output ##################################
  
  save.image(file = here("output", "error_eval_workspace.RData"))
  # save.image(file = here("output", "error_eval_workspace_lower_start.RData"))
  