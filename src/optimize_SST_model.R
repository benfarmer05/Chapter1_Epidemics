  
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
  library(ggnewscale)
  library(extrafont)
  library(scales)

  ################################## Set-up ##################################
  
  #import workspace from upstream script
  # load(here("output/plots_basic_workspace.RData"))
  load(here("output/SST_effect_workspace.RData"))
  
  ################################## More mature sigmoid & logistic behaviors of temperature ##################################
  # AKA - zeta / eta parameter plotting sandbox
  
  #making choice to pull boundaries of SST from entire dataset (1985 - 2024). could constrain
  T_min = DHW.CRW.full %>%
    # pull(SST.90th_HS) %>%
    pull(SST_90th_HS_smoothed) %>%
    min
  T_max = DHW.CRW.full %>%
    # pull(SST.90th_HS) %>%
    pull(SST_90th_HS_smoothed) %>%
    max
  
  # Define temperature range from 23 to 33
  T <- seq(T_min, T_max, length.out = 1000)
  
  
  #APPROACH 1: Inflection point at 30.5C (original)
  temp_effect_inflection <- function(temp, T_inflection, zeta_bounded) {

    # Convert bounded zeta to effective zeta
    # zeta_effective <- 1 * zeta_bounded/(1-1*zeta_bounded)
    zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)

    # Logistic function with inflection point at T_inflection
    modifier <- 1 / (1 + exp(-zeta_effective * (temp - T_inflection)))

    return(modifier)
  }

  #APPROACH 2: Original model clarified
  temp_effect_original <- function(temp, tmin, zeta_bounded) {

    # Convert bounded zeta to effective zeta

    # zeta_effective <- 5 * zeta_bounded/(1-0.7*zeta_bounded)
    zeta_effective <- 1 * zeta_bounded/(1 - 1*zeta_bounded)

    # Original formula clarified
    # This is a non-centered sigmoid that:
    # - Starts near 0 at tmin
    # - Increases asymptotically toward 1 as temp increases
    # - Has no inflection point assumption
    modifier <- (1 - exp(-zeta_effective * (temp - tmin))) / (1 + exp(-zeta_effective * (temp - tmin)))

    return(modifier)
  }

  #APPROACH 3: Decreasing function with temperature
  temp_effect_decreasing <- function(temp, tmax, tmin, zeta_bounded) {

    # Convert bounded zeta to effective zeta
    # zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    # zeta_effective <- 1 * zeta_bounded/(1-1*zeta_bounded)
    zeta_effective <- zeta_bounded / (1 - zeta_bounded)
    
    # Calculate midpoint of temperature range for centering the sigmoid
    T_mid <- 30.5  # (tmin + tmax) / 2
    
    # Decreasing sigmoid that starts near 1 (at low temps) and approaches 0 as temp increases
    # Note the positive sign before zeta_effective which makes the function decrease with temperature
    modifier <- 1 / (1 + exp(zeta_effective * (temp - T_mid)))

    return(modifier)
  }

  # APPROACH 4: Non-centered sigmoid decreasing with temperature (FIXED)
  temp_effect_noncentered_decreasing <- function(temp, tmin, tmax, zeta_bounded) {

    # Convert bounded zeta to effective zeta
    zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded) # apply more tuning to behavior
    # zeta_effective <- zeta_bounded/(1-zeta_bounded)

    # Modified non-centered sigmoid that:
    # - Equals 1 at tmin (23°C)
    # - Decreases toward 0 as temperature increases
    # - Maintains the non-centered sigmoid shape but inverted

    # Calculate the relative position in the temperature range
    rel_temp <- (temp - tmin) / (tmax - tmin)

    # Apply the decreasing non-centered sigmoid
    modifier <- 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))
    # modifier <- 1 - ((1 - exp(-zeta_bounded * rel_temp)) / (1 + exp(-zeta_bounded * rel_temp)))

    return(modifier)
  }

  # temp_effect_noncentered_decreasing <- function(temp, T_min = 23, T_max = 33, zeta_bounded = 0.5) {
  #   # Convert bounded zeta to effective zeta
  #   zeta_effective <- zeta_bounded/(1-zeta_bounded)
  #
  #   # Calculate the relative position in the temperature range
  #   rel_temp <- (temp - T_min) / (T_max - T_min)
  #
  #   # Most simplified form - standard logistic function
  #   modifier <- 1 / (1 + exp(zeta_effective * rel_temp))
  #
  #   return(modifier)
  # }

  # Define bounded zeta values from 0 to 1 (avoiding 1 which would create infinity)
  zeta_values <- c(0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

  # Create data for all approaches
  plot_data_inflection <- data.frame()
  plot_data_original <- data.frame()
  plot_data_decreasing <- data.frame()
  plot_data_noncentered_decreasing <- data.frame()

  for (zeta in zeta_values) {
    # Data for approach 1
    modifiers_inflection <- sapply(T, function(t) temp_effect_inflection(t, T_inflection = 30.5, zeta_bounded = zeta))
    temp_data_inflection <- data.frame(
      Temperature = T,
      Modifier = modifiers_inflection,
      Zeta = as.factor(zeta)
    )
    plot_data_inflection <- rbind(plot_data_inflection, temp_data_inflection)

    # Data for approach 2
    modifiers_original <- sapply(T, function(t) temp_effect_original(t, tmin = T_min, zeta_bounded = zeta))
    temp_data_original <- data.frame(
      Temperature = T,
      Modifier = modifiers_original,
      Zeta = as.factor(zeta)
    )
    plot_data_original <- rbind(plot_data_original, temp_data_original)
  
    # Data for approach 3
    modifiers_decreasing <- sapply(T, function(t) temp_effect_decreasing(t, tmin = T_min, tmax = T_max, zeta_bounded = zeta))
    temp_data_decreasing <- data.frame(
      Temperature = T,
      Modifier = modifiers_decreasing,
      Zeta = as.factor(zeta)
    )
    plot_data_decreasing <- rbind(plot_data_decreasing, temp_data_decreasing)

    # Data for approach a4
    modifiers_noncentered_decreasing <- sapply(T, function(t) temp_effect_noncentered_decreasing(t, tmin = T_min, tmax = T_max, zeta_bounded = zeta))
    temp_data_noncentered_decreasing <- data.frame(
      Temperature = T,
      Modifier = modifiers_noncentered_decreasing,
      Zeta = as.factor(zeta)
    )
    plot_data_noncentered_decreasing <- rbind(plot_data_noncentered_decreasing, temp_data_noncentered_decreasing)
  }

  # Plot 1: Original inflection point model
  p1 <- ggplot(plot_data_inflection, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(size = 1) +
    labs(
      title = "Temperature Threshold with Classic Logistic Function",
      x = "Temperature (°C)",
      y = "Transmission Modifier"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    scale_color_viridis_d(option = "inferno") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(T_min, T_max, 2)) +
    geom_vline(xintercept = 30.5, linetype = "dashed", alpha = 0.5)

  # Plot 2: Original increasing model
  p2 <- ggplot(plot_data_original, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(size = 1) +
    labs(
      title = "Temperature Effect with non-centered Sigmoid Function",
      x = "Temperature (°C)",
      y = "Transmission Modifier"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    scale_color_viridis_d(option = "inferno") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(T_min, T_max, 2)) +
    geom_vline(xintercept = 23, linetype = "dashed", alpha = 0.5)
  
  # Plot 3: Classic sigmoid decreasing model
  p3 <- ggplot(plot_data_decreasing, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(size = 1, show.legend = FALSE) +
    labs(
      x = "Temperature (°C)",
      y = "beta modifier"
    ) +
    theme_classic(base_family = 'Georgia') +
    theme(
      legend.position = "right",
      # legend.position = "inside",
      # legend.position.inside = c(0, 0),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    scale_color_viridis_d(option = "inferno") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(T_min, T_max, 2)) +
    geom_vline(xintercept = 30.5, color = 'grey40', linetype = 'dashed', size = 1)
    
  # # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # quartz(h = 3, w = 3)
  # 
  # p3
  # 
  # # #ggplot-export to image
  # # ggsave(filename = here("output", "p3.jpg"), device = "jpg", width = 3, height = 3, dpi = 1200)
  # 
  # # Close the Quartz device
  # dev.off()
  
  zeta_values_overlay <- c(0.069, 0.044)
  
  # Generate new data for the overlay curves, adding a column "Curve"
  new_data_overlay <- data.frame()
  for (zeta in zeta_values_overlay) {
    modifiers_overlay <- sapply(T, function(t) temp_effect_decreasing(t, tmin = T_min, tmax = T_max, zeta_bounded = zeta))
    # Label 0.069 as "zeta" (red) and 0.044 as "eta" (blue)
    curve_label <- ifelse(zeta == 0.069, "zeta", "eta")
    temp_data_overlay <- data.frame(
      Temperature = T,
      Modifier = modifiers_overlay,
      Zeta = as.factor(zeta),
      Curve = curve_label
    )
    new_data_overlay <- rbind(new_data_overlay, temp_data_overlay)
  }
  
  # Start a new color scale so that overlay lines can have a different mapping
  zeta_eta_fig <- p3 + new_scale_color()
  
  # Overlay the new curves with custom colors and a legend
  zeta_eta_fig <- zeta_eta_fig +
    geom_line(
      data = new_data_overlay,
      aes(x = Temperature, y = Modifier, group = Curve, color = Curve),
      size = 1 #1.5
    ) +
    scale_color_manual(
      values = c("zeta" = "red", "eta" = "blue"),
      name = "",
      labels = c("zeta (0.069)", "eta (0.044)")
    )
  
  zeta_eta_fig = zeta_eta_fig + theme(
    # legend.position = "inside", legend.position.inside =  c(0.92, .85)
    legend.position = "bottom",
    legend.margin = margin(-10, 0, , 0),  # Moves the legend closer to the plot
    plot.margin = margin(10, 10, 10, 10)  # Adjusts overall plot margins
  )
  
  # # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # # quartz(h = 5, w = 3.35)
  # # quartz(h = 6, w = 7.087)
  # quartz(h = 3, w = 3)
  # 
  # zeta_eta_fig
  # 
  # # #ggplot-export to image
  # # ggsave(filename = here("output", "zeta_eta.jpg"), device = "jpg", width = 3, height = 3, dpi = 1200)
  # 
  # # Close the Quartz device
  # dev.off()
  
  # Plot 4: Fixed non-centered sigmoid decreasing model
  p4 <- ggplot(plot_data_noncentered_decreasing, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(size = 1) +
    labs(
      title = "Temperature Effect with Inverted Non-centered Sigmoid",
      x = "Temperature (°C)",
      y = "Transmission Modifier"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    scale_color_viridis_d(option = "inferno") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(T_min, T_max, 2)) +
    geom_vline(xintercept = 23, linetype = "dashed", alpha = 0.5)

  # Combine plots - using a 2x2 grid for better viewing
  (p1 + p2) / (p3 + p4)

  # Function explanations for mathematical clarity
  cat("Approach 1: Inflection Point at 30.5°C\n")
  cat("Formula: modifier = 1 / (1 + exp(-zeta_effective * (temp - T_inflection)))\n")
  cat("- Classic logistic function centered at T_inflection\n")
  cat("- At T_inflection (30.5°C), modifier = 0.5\n")
  cat("- Assumes symmetric behavior around inflection point\n\n")

  cat("Approach 2: Original Model Clarified\n")
  cat("Formula: modifier = (1 - exp(-zeta_effective * (temp - T_min))) / (1 + exp(-zeta_effective * (temp - T_min)))\n")
  cat("- At T_min (23°C), modifier ≈ 0\n")
  cat("- No assumption of inflection point\n")
  cat("- Asymmetric behavior that rises more steeply at first\n\n")

  cat("Approach 3: Decreasing Function with Temperature\n")
  cat("Formula: modifier = 1 / (1 + exp(zeta_effective * (temp - T_mid)))\n")
  cat("- At low temperatures (23°C), modifier ≈ 1\n")
  cat("- At high temperatures (33°C), modifier approaches 0\n")
  cat("- Higher zeta values create steeper drop in transmission\n")
  cat("- Centered at T_mid (28°C) where modifier = 0.5\n\n")

  cat("Approach 4: Inverted Non-centered Sigmoid with Temperature\n")
  cat("Formula: modifier = 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))\n")
  cat("- At T_min (23°C), modifier = 1 exactly for all zeta values\n")
  cat("- Decreases toward 0 as temperature increases\n")
  cat("- Higher zeta values create steeper drop in transmission\n")
  cat("- Preserves the asymmetric shape but inverted from approach 2\n")

  ################################## SST set-up ##################################

  # # Filter DHW.CRW.full for dates from the last date in SST_sites to the end of 2020
  # last_date_in_SST_sites <- max(SST_sites$date)
  # SST_DHW_extended <- DHW.CRW.full %>%
  #   filter(date > last_date_in_SST_sites & date <= as.Date("2020-12-31")) %>%
  #   mutate(site = "mid",  # Assuming the same site as in original data
  #          days.inf.site = NA,  # You may want to adjust this if needed
  #          time = as.numeric(difftime(date, min(date), units = "days")))
  # 
  # # Combine the original SST_sites with the extended data
  # SST_sites_extended <- bind_rows(SST_sites,
  #                                 select(SST_DHW_extended,
  #                                        site,
  #                                        date,
  #                                        days.inf.site,
  #                                        time,
  #                                        SST = SST.90th_HS))
  # 
  # # Sort the dataframe by site, then date
  # SST_sites_extended <- SST_sites_extended %>%
  #   arrange(site, date)
  
  # Find the last date in the original SST_sites
  last_date_in_SST_sites <- max(SST_sites$date)
  
  # Filter DHW.CRW.full for dates from the last date in SST_sites to the end of 2020
  SST_DHW_extended <- DHW.CRW.full %>%
    filter(date > last_date_in_SST_sites & date <= as.Date("2020-12-31")) %>%
    # mutate(site = "mid",
    #        days.inf.site = NA,
    #        time = NA)  # Assuming the same site as in original data
    slice(rep(1:n(), each = 3)) %>%
    mutate(site = rep(c("mid", "off", "near"), times = n()/3),
           days.inf.site = NA,
           time = NA)  
  
  # Combine the original SST_sites with the extended data
  SST_sites_extended <- bind_rows(SST_sites, 
                                  select(SST_DHW_extended, 
                                         site, 
                                         date, 
                                         days.inf.site, 
                                         SST = SST.90th_HS))
  
  DHW_sites_extended <- bind_rows(DHW_sites, 
                                  select(SST_DHW_extended, 
                                         site, 
                                         date, 
                                         days.inf.site, 
                                         DHW = DHW_from_90th_HS.1))
  
  # Reset time sequence for each site
  SST_sites_extended <- SST_sites_extended %>%
    group_by(site) %>%
    arrange(date) %>%
    mutate(time = row_number() - 1) %>%  # Start at 0
    ungroup() %>%
    arrange(site, date)
  
  DHW_sites_extended <- DHW_sites_extended %>%
    group_by(site) %>%
    arrange(date) %>%
    mutate(time = row_number() - 1) %>%  # Start at 0
    ungroup() %>%
    arrange(site, date)
  
  #create moving average of SST data to smooth it and reduce wiggliness of the data
  #  note - this is not needed for DHW since it represents accumulation of stress, rather than real-time
  SST_sites_extended <- SST_sites_extended %>%
    group_by(site) %>%
    mutate(SST_smoothed = smooth(SST, kind = "3R")) %>%
    ungroup()  
  
  #making choice to pull boundaries of SST from entire dataset (1985 - 2024). could constrain
  T_min = DHW.CRW.full %>%
    # pull(SST.90th_HS) %>%
    pull(SST_90th_HS_smoothed) %>%
    min
  T_max = DHW.CRW.full %>%
    # pull(SST.90th_HS) %>%
    pull(SST_90th_HS_smoothed) %>%
    max
  
  susceptible_ref_SST = susceptible_ref %>%
    mutate(Site = case_when(
      Site == "Offshore" ~ "off",
      Site == 'Midchannel' ~ 'mid',
      Site == 'Nearshore' ~ 'near'
    ))
  
  obs.total.figures_SST = obs.total.figures_SST %>%
    mutate(Site = case_when(
      Site == "Offshore" ~ "off",
      Site == 'Midchannel' ~ 'mid',
      Site == 'Nearshore' ~ 'near'
    ))
  
  sites_SST = c('off', 'mid', 'near')
  
  #update days_sites to include 2020 (manually added "observations")
  days_sites_SST <- days_sites %>%
    mutate(
      days = map2(site, days, function(site_name, existing_days) {
        new_days <- obs.total.figures_SST %>%
          filter(Site == site_name, Compartment == 'Infected') %>%
          pull(days.inf.site)
        c(existing_days, new_days[!new_days %in% existing_days])
      })
    )
  
  ################################## Model: single-host with SST ##################################
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
      current_SST = SST$SST_smoothed[closest_index]
      # current_DHW = DHW$DHW[closest_index]
      
      #host density-null conditions
      # transmission_modifier = 1
      # removal_rate = g
      # transmission_rate = b * transmission_modifier
      removal_rate = g
      transmission_rate = b
      
      #STOPPING POINT - 27 MAR 2025
      #  - I think a real issue here is that my equation behaviors are pretty weird
      
      # # APPROACH 4: no threshold
      # # NOTE - issue is that low temperatures actually have the decreasing effect here. high ones should instead
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
      
      #APPROACH 3: threshold
      # NOTE - issue here is that there is no meaningful inflection point - maybe?
      #         - or something else like boundary conditions
      zeta_effective <- z / (1 - z)
      eta_effective <- e / (1 - e)
      beta_scaling_factor <- 1 / (1 + exp(zeta_effective * (current_SST - 30.5)))
      gamma_scaling_factor <- 1 / (1 + exp(eta_effective * (current_SST - 30.5)))
      transmission_rate <- transmission_rate * beta_scaling_factor
      removal_rate <- removal_rate * gamma_scaling_factor
      
      # #APPROACH 0: halting transmission/removal
      # if(current_SST > 30.5){
      #   removal_rate = 0.01
      # }
      
      dS.dt = -transmission_rate * S * I / N
      dI.dt = transmission_rate * S * I / N - removal_rate * I
      dR.dt = removal_rate * I
  
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  # ################################## SST optimization ##################################
  # my.SIRS.SST = vector('list', length(sites_SST))
  # params.SST = vector('list', length(sites_SST))
  # curr.type = 'Fitted' #the below is for a fitting model for multi-host transmission (no DHW or projection)
  # 
  # # for(i in 1:length(sites_SST)){
  # #   
  # #   site.loop = sites_SST[i]
  #   
  #   site.loop = "mid" #for testing purposes
  #   i = 1
  #   # site.loop = "near" #for testing purposes
  #   # i = 2
  #   # site.loop = "off" #for testing purposes
  #   # i = 3
  #   
  #   days <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
  #     filter(site == site.loop) %>%
  #     pull(days) %>%
  #     unlist()
  #   
  #   days.obs = na.omit(days) %>%
  #     as.numeric()
  #   
  #   days.model = SST_sites_extended %>%
  #     filter(site == site.loop) %>%
  #     pull(time) %>%
  #     max() %>%
  #     seq(from = 0, to = .)
  #   
  #   SST_df = SST_sites_extended %>%
  #     filter(site == site.loop) %>%
  #     select(date, time, SST_smoothed)
  #   
  #   DHW_df <- DHW_sites_extended %>%
  #     filter(site == site.loop) %>%
  #     select(date, time, DHW)
  #   
  #   SST_values = SST_df %>%
  #     pull(SST_smoothed)
  #   
  #   DHW_values = DHW_df %>%
  #     pull(DHW)
  #   
  #   # Ensure that the lengths of SST_values and DHW_values match the length of days.model
  #   #   NOTE - this error checker really should be verifying dates, not the sequence - would likely require comparing to summary
  #   if (length(SST_values) != length(days.model)) {
  #     stop("Length of SST_values does not match length of days.model.")
  #   }
  #   if (length(DHW_values) != length(days.model)) {
  #     stop("Length of DHW_values does not match length of days.model.")
  #   }
  #   
  #   N.site = susceptible_ref_SST %>%
  #     filter(Site == site.loop) %>%
  #     slice(1) %>% #all values for N.site are the same between categories, so slice first row
  #     pull(N.site)
  #   
  #   cover.site = susceptible_ref_SST %>%
  #     filter(Site == site.loop) %>%
  #     slice(1) %>% #same as above
  #     pull(cover.site)
  #   
  #   #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
  #   inftiss = obs.total.figures_SST %>%
  #     filter(Site == site.loop, Compartment == "Infected") %>%
  #     slice(head(row_number(), n()-DHW.modifier)) %>%
  #     pull(tissue)    
  #   
  #   remtiss = obs.total.figures_SST %>%
  #     filter(Site == site.loop, Compartment == "Recovered") %>%
  #     slice(head(row_number(), n()-DHW.modifier)) %>%
  #     pull(tissue)
  #   
  #   # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
  #   first_valid_idx <- which(!is.na(days))[1] #find the first non-NA index
  #   inftiss <- inftiss[first_valid_idx:length(inftiss)]
  #   remtiss <- remtiss[first_valid_idx:length(remtiss)]
  #   
  #   #initial conditions
  #   I.tiss = inftiss[1] #first non-NA & non-zero infection entry
  #   # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  #   S.tiss = N.site - I.tiss
  #   R.tiss = 0
  #   
  #   # Set up the data and initial conditions
  #   coraldata.tiss = list(inftiss, remtiss)
  #   # coraldata.tiss = list(inftiss) # NOTE - attempting to fit using just infections
  #   initial_state.tiss = c(S.tiss, I.tiss, R.tiss, N.site, cover.site)
  #   
  #   # Define the objective function for optimization
  #   objective_function = function(params, data, time, initial_state){
  #     
  #     # # Print parameters to track what values are being tested
  #     # cat("Testing params:", params, "\n")
  #     
  #     # #testing
  #     # betas = 3.05 #for midchannel
  #     # gammas = 3.01 #for midchannel
  #     # lambdas = 1.0
  #     # zetas = 0.89
  #     # etas = 0.70
  #     # initial_state = initial_state.tiss
  #     # time = days.model
  #     # data = coraldata.tiss
  #     
  #     betas = params[1]
  #     gammas = params[2]
  #     lambdas = params[3]
  #     zetas = params[4]
  #     etas = params[5]
  #     
  #     # # sanity check
  #     # print(paste("SST_df time range:", min(SST_df$time), "to", max(SST_df$time)))
  #     # print(paste("days.model range:", min(days.model), "to", max(days.model)))
  #     # print(paste("Initial state:", initial_state.tiss))
  #     # print(paste("Parameters:", params[1], params[2], params[3], params[4], params[5]))
  #     
  #     # Check if SST_df and DHW_df are valid
  #     if(any(is.na(SST_df$time)) || any(is.na(SST_df$SST_smoothed)) || 
  #        any(is.na(DHW_df$time)) || any(is.na(DHW_df$DHW))) {
  #       cat("Warning: NAs in SST or DHW data\n")
  #     }
  #     
  #     # Expand tryCatch to provide more information
  #     SIR.out = tryCatch({
  #       result <- ode(c(S = initial_state[1], I = initial_state[2], R = initial_state[3]),
  #                     time, SIR_project, c(b = betas, g = gammas,
  #                                          N = initial_state[4],
  #                                          z = zetas,
  #                                          e = etas,
  #                                          l = lambdas,
  #                                          C = initial_state[5]),
  #                     SST = SST_df,
  #                     DHW = DHW_df)
  #       # Check if ODE output contains NAs
  #       if(any(is.na(result))) {
  #         cat("Warning: ODE solution contains NAs\n")
  #       }
  #       data.frame(result)
  #     }, error = function(e) {
  #       cat("Error in ODE:", conditionMessage(e), "\n")
  #       return(NA)
  #     })
  #     
  #     # Check if ODE solving failed
  #     if(all(is.na(SIR.out))) {
  #       cat("ODE solving failed completely\n")
  #       return(Inf)  # Return Inf instead of NaN to avoid optimizer error
  #     }
  #     
  #     # Extract simulated values
  #     sim.inf.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'I')]
  #     sim.rem.total = SIR.out[which(SIR.out$time %in% days.obs), which(colnames(SIR.out) %in% 'R')]
  #     
  #     # Extract observed values
  #     obs.inf.total = unlist(data[[1]])
  #     obs.rem.total = unlist(data[[2]])
  #     
  #     # # NOTE - this is only needed if you are NOT manually adding in data for 2020, because when doing that,
  #     # #         the last-timepoint issue is already corrected
  #     # # Remove last infected timepoint, since the last "observation" is logged as zero due to nature of 
  #     # #   disease backtracking upstream in the code, and not meant to signify an actual datapoint
  #     # sim.inf.total <- sim.inf.total %>% head(-1)
  #     # obs.inf.total <- obs.inf.total %>% head(-1)
  #     
  #     # #NOTE - this is only needed if you ARE manually adding in data for 2020, since only infection
  #     # #          was added in, not removal
  #     # sim.rem.total <- sim.rem.total %>%
  #     #   head(length(obs.rem.total))
  #     
  #     # #NOTE - test to try and get the optimizer working
  #     # sim.inf.total <- sim.inf.total %>%
  #     #   head(12)
  #     # obs.inf.total = obs.inf.total %>%
  #     #   head(12)
  #     
  #     # Calculate residuals
  #     diff.inf.total = (sim.inf.total - obs.inf.total)
  #     diff.rem.total = (sim.rem.total - obs.rem.total)
  #     
  #     #aggregate sum of absolute residuals
  #     # NOTE - this is not sum of squares and should be clearly stated/defended in the manuscript if used
  #     sum_abs_diff_I.total = sum(sum(abs(diff.inf.total)))
  #     sum_abs_diff_R.total = sum(sum(abs(diff.rem.total)))
  #     sum_diff.abs.total = sum_abs_diff_R.total #version where only removal is fitted to
  #     # sum_diff.abs.total = sum_abs_diff_I.total #version where only infections are fitted to
  #     
  #     # #minimize using sum of squared residuals
  #     # sum_diff.total = sum(diff.rem.total^2)
  #     
  #   return(sum_diff.abs.total) #return only the residual metric for the epidemic wave being fit to
  #   }
  #   
  #   #boundary conditions
  #   # lower_bounds.tiss = c(0.5, 0.4, lambda.modifier, 1, 1)  # Lower bounds for betas, gammas, lambdas, zetas, etas
  #   # upper_bounds.tiss = c(1.5, 2, lambda.modifier, 1, 1)  # Upper bounds for betas, gammas, lambdas, zetas, etas
  #   lower_bounds.tiss = c(0, 0, lambda.modifier, 0, 0)  # Lower bounds for betas, gammas, lambdas, zetas, etas
  #   upper_bounds.tiss = c(2, 2, lambda.modifier, 0.2, 0.2)  # Upper bounds for betas, gammas, lambdas, zetas, etas
  #   
  #   control = list(itermax = 200)  # number of optimizer iterations. 200 is default
  #   
  #   # Run the optimization
  #   result.tiss = DEoptim(fn = objective_function, lower = lower_bounds.tiss, upper = upper_bounds.tiss,
  #                         data = coraldata.tiss, time = days.model, initial_state = initial_state.tiss,
  #                         control = control)
  #   
  #   # Extract the optimized parameters
  #   optimized_params.tiss = result.tiss$optim$bestmem
  #   
  #   # Print the optimized parameters
  #   min.beta.tiss = as.numeric(optimized_params.tiss[1])
  #   min.gamma.tiss = as.numeric(optimized_params.tiss[2])
  #   min.zeta.tiss = as.numeric(optimized_params.tiss[4])
  #   min.eta.tiss = as.numeric(optimized_params.tiss[5])
  #   
  #   params.SST[[i]] = c(min.beta.tiss, min.gamma.tiss, min.zeta.tiss, min.eta.tiss)
  #   
  #   #simulation using initial state variables and best-fit beta/gamma/zeta/eta parameters
  #   SIR.out.tiss = data.frame(ode(c(S = initial_state.tiss[1], I = initial_state.tiss[2], R = initial_state.tiss[3]),
  #                                 days.model, SIR_project, c(b = min.beta.tiss, g = min.gamma.tiss,
  #                                                        N = initial_state.tiss[4],
  #                                                        z = min.zeta.tiss,
  #                                                        e = min.eta.tiss,
  #                                                        l = lambdas,
  #                                                        C = initial_state.tiss[5]),
  #                                 SST = SST_df,
  #                                 DHW = DHW_df))
  #   
  #   my.SIRS.SST[[i]] = SIR.out.tiss
  # # }
  # 
  # ################################## Plots of epidemics w/ SST ##################################
  # 
  # obs.total.figures.old = obs.total.figures %>%
  #   mutate(Site = case_when(
  #     Site == "Offshore" ~ "off",
  #     Site == 'Midchannel' ~ 'mid',
  #     Site == 'Nearshore' ~ 'near'
  #   )) %>%
  #   mutate(Compartment = case_when(
  #     Compartment == "Recovered" ~ "Dead",
  #     TRUE ~ Compartment #keep all other strings unchanged
  #   ))
  # 
  # obs.total.figures_SST = obs.total.figures_SST %>%
  #   mutate(Compartment = case_when(
  #     Compartment == "Recovered" ~ "Dead",
  #     TRUE ~ Compartment #keep all other strings unchanged
  #   ))
  # 
  # #Nearshore
  # site.loop = 'Nearshore'
  # curr.site = 'near'
  # days.obs <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
  #   filter(site == curr.site) %>%
  #   pull(days) %>%
  #   unlist()
  # days.obs = na.omit(days.obs) %>%
  #   as.numeric()
  # 
  # order = 3
  # 
  # output.basic.nearshore.SST = my.SIRS.SST[[order]]
  # params.basic.nearshore.SST = params.SST[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  # beta.nearshore.SST = params.basic.nearshore.SST[1]
  # gamma.nearshore.SST = params.basic.nearshore.SST[2]
  # zeta.nearshore.SST = params.basic.nearshore.SST[3]
  # eta.nearshore.SST = params.basic.nearshore.SST[4]
  # 
  # tab.nearshore = tibble(round(beta.nearshore.SST, 2), round(gamma.nearshore.SST, 2),
  #                        round(zeta.nearshore.SST, 2), round(eta.nearshore.SST, 2))
  # names(tab.nearshore) = c('beta', 'gamma', 'zeta', 'eta')
  # 
  # #calculate R-squared and update error table
  # # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  # curr.type = 'DHW'
  # curr.wave = 'Full'
  # 
  # # sim.rem.total = output.basic.nearshore.SST[which(output.basic.nearshore.SST$time %in% days.obs), which(colnames(output.basic.nearshore.SST) %in% 'R')]
  # # obs.rem.total = obs.model %>%
  # #   filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
  # #   slice(head(row_number(), n()-DHW.modifier)) %>%
  # #   pull(tissue)
  # #
  # # # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  # # if (length(obs.rem.total) > length(sim.rem.total)) {
  # #   obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  # # }
  # 
  # # diff.rem.total = (sim.rem.total - obs.rem.total)
  # # sum_diff.total = sum(diff.rem.total^2)
  # # mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  # # tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  # # r_squared.basic.nearshore.SST = 1 - (sum_diff.total / tss_rem.total)
  # 
  # # error_eval <- error_eval %>%
  # #   mutate(R_squared = ifelse(site == curr.site &
  # #                               host == curr.host &
  # #                               type == curr.type &
  # #                               wave == curr.wave,
  # #                             r_squared.basic.nearshore.SST,
  # #                             R_squared))
  # 
  # output.basic.nearshore.SST = pivot_longer(output.basic.nearshore.SST, cols = -1, names_to = c("Compartment")) %>%
  #   mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
  #   mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
  #                               ifelse(Compartment == 'I', 'Infected',
  #                                      ifelse(Compartment == 'R', 'Dead', Compartment))))
  # 
  # colnames(output.basic.nearshore.SST)[1] = 'days.model'
  # colnames(output.basic.nearshore.SST)[3] = 'tissue'
  # 
  # p.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST, aes(days.model, tissue, colour = Compartment)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of tissue (m2)") +
  #   ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures.old %>% filter(Site == curr.site), aes(days.inf.site, tissue, colour = Compartment)) +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site & Compartment == 'Dead'), aes(days.inf.site, tissue, colour = Compartment)) +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   annotate(geom = "table", x = min(output.basic.nearshore.SST$days.model), y = min(output.basic.nearshore.SST$tissue)*0.7, label = list(tab.nearshore),
  #            vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
  #   theme_classic(base_family = 'Georgia') +
  #   theme(panel.background = element_rect(fill = "gray90"))
  # 
  # p.S.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of susceptible tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.I.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of infected tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Infected"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.D.fit.nearshore.basic.SST = ggplot(data = output.basic.nearshore.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of dead tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Dead"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # #Midchannel
  # site.loop = 'Midchannel'
  # curr.site = 'mid'
  # days.obs <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
  #   filter(site == curr.site) %>%
  #   pull(days) %>%
  #   unlist()
  # days.obs = na.omit(days.obs) %>%
  #   as.numeric()
  # 
  # order = 2
  # # order = 1
  # 
  # output.basic.midchannel.SST = my.SIRS.SST[[order]]
  # params.basic.midchannel.SST = params.SST[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  # beta.midchannel.SST = params.basic.midchannel.SST[1]
  # gamma.midchannel.SST = params.basic.midchannel.SST[2]
  # zeta.midchannel.SST = params.basic.midchannel.SST[3]
  # eta.midchannel.SST = params.basic.midchannel.SST[4]
  # 
  # tab.midchannel = tibble(round(beta.midchannel.SST, 2), round(gamma.midchannel.SST, 2),
  #                        round(zeta.midchannel.SST, 2), round(eta.midchannel.SST, 2))
  # names(tab.midchannel) = c('beta', 'gamma', 'zeta', 'eta')
  # 
  # #calculate R-squared and update error table
  # # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  # curr.type = 'DHW'
  # curr.wave = 'Full'
  # 
  # # sim.rem.total = output.basic.midchannel.SST[which(output.basic.midchannel.SST$time %in% days.obs), which(colnames(output.basic.midchannel.SST) %in% 'R')]
  # # obs.rem.total = obs.model %>%
  # #   filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
  # #   slice(head(row_number(), n()-DHW.modifier)) %>%
  # #   pull(tissue)
  # # 
  # # # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  # # if (length(obs.rem.total) > length(sim.rem.total)) {
  # #   obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  # # }
  # 
  # # diff.rem.total = (sim.rem.total - obs.rem.total)
  # # sum_diff.total = sum(diff.rem.total^2)
  # # mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  # # tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  # # r_squared.basic.midchannel.SST = 1 - (sum_diff.total / tss_rem.total)
  # 
  # # error_eval <- error_eval %>%
  # #   mutate(R_squared = ifelse(site == curr.site & 
  # #                               host == curr.host & 
  # #                               type == curr.type & 
  # #                               wave == curr.wave,
  # #                             r_squared.basic.midchannel.SST,
  # #                             R_squared))
  # 
  # output.basic.midchannel.SST = pivot_longer(output.basic.midchannel.SST, cols = -1, names_to = c("Compartment")) %>%
  #   mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
  #   mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
  #                               ifelse(Compartment == 'I', 'Infected',
  #                                      ifelse(Compartment == 'R', 'Dead', Compartment))))
  # 
  # colnames(output.basic.midchannel.SST)[1] = 'days.model'
  # colnames(output.basic.midchannel.SST)[3] = 'tissue'
  # 
  # p.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST, aes(days.model, tissue, colour = Compartment)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of tissue (m2)") +
  #   ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site), aes(days.inf.site, tissue, colour = Compartment)) +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   annotate(geom = "table", x = min(output.basic.midchannel.SST$days.model), y = min(output.basic.midchannel.SST$tissue)*0.7, label = list(tab.midchannel),
  #            vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
  #   theme_classic(base_family = 'Georgia') +
  #   theme(panel.background = element_rect(fill = "gray90"))
  # 
  # p.S.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of susceptible tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.I.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of infected tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Infected"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.D.fit.midchannel.basic.SST = ggplot(data = output.basic.midchannel.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of dead tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Dead"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # #Offshore
  # site.loop = 'Offshore'
  # curr.site = 'off'
  # days.obs <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
  #   filter(site == curr.site) %>%
  #   pull(days) %>%
  #   unlist()
  # days.obs = na.omit(days.obs) %>%
  #   as.numeric()
  # 
  # order = 1
  # 
  # output.basic.offshore.SST = my.SIRS.SST[[order]]
  # params.basic.offshore.SST = params.SST[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  # beta.offshore.SST = params.basic.offshore.SST[1]
  # gamma.offshore.SST = params.basic.offshore.SST[2]
  # zeta.offshore.SST = params.basic.offshore.SST[3]
  # eta.offshore.SST = params.basic.offshore.SST[4]
  # 
  # tab.offshore = tibble(round(beta.offshore.SST, 2), round(gamma.offshore.SST, 2),
  #                         round(zeta.offshore.SST, 2), round(eta.offshore.SST, 2))
  # names(tab.offshore) = c('beta', 'gamma', 'zeta', 'eta')
  # 
  # #calculate R-squared and update error table
  # # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  # curr.type = 'DHW'
  # curr.wave = 'Full'
  # 
  # # sim.rem.total = output.basic.offshore.SST[which(output.basic.offshore.SST$time %in% days.obs), which(colnames(output.basic.offshore.SST) %in% 'R')]
  # # obs.rem.total = obs.model %>%
  # #   filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
  # #   slice(head(row_number(), n()-DHW.modifier)) %>%
  # #   pull(tissue)
  # #
  # # # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  # # if (length(obs.rem.total) > length(sim.rem.total)) {
  # #   obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  # # }
  # 
  # # diff.rem.total = (sim.rem.total - obs.rem.total)
  # # sum_diff.total = sum(diff.rem.total^2)
  # # mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  # # tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  # # r_squared.basic.offshore.SST = 1 - (sum_diff.total / tss_rem.total)
  # 
  # # error_eval <- error_eval %>%
  # #   mutate(R_squared = ifelse(site == curr.site &
  # #                               host == curr.host &
  # #                               type == curr.type &
  # #                               wave == curr.wave,
  # #                             r_squared.basic.offshore.SST,
  # #                             R_squared))
  # 
  # output.basic.offshore.SST = pivot_longer(output.basic.offshore.SST, cols = -1, names_to = c("Compartment")) %>%
  #   mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
  #   mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
  #                               ifelse(Compartment == 'I', 'Infected',
  #                                      ifelse(Compartment == 'R', 'Dead', Compartment))))
  # 
  # colnames(output.basic.offshore.SST)[1] = 'days.model'
  # colnames(output.basic.offshore.SST)[3] = 'tissue'
  # 
  # p.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST, aes(days.model, tissue, colour = Compartment)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of tissue (m2)") +
  #   ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site), aes(days.inf.site, tissue, colour = Compartment)) +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   annotate(geom = "table", x = min(output.basic.offshore.SST$days.model), y = min(output.basic.offshore.SST$tissue)*0.7, label = list(tab.offshore),
  #            vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
  #   theme_classic(base_family = 'Georgia') +
  #   theme(panel.background = element_rect(fill = "gray90"))
  # 
  # p.S.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of susceptible tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.I.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of infected tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Infected"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.D.fit.offshore.basic.SST = ggplot(data = output.basic.offshore.SST %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of dead tissue") +
  #   ggtitle(paste(c("", curr.site), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total.figures_SST %>% filter(Site == curr.site, Compartment == "Dead"), aes(days.inf.site, tissue)) +
  #   theme_classic(base_family = 'Georgia')
  # 
  # p.I.fit.nearshore.basic.SST
  # p.D.fit.nearshore.basic.SST
  # 
  # p.I.fit.midchannel.basic.SST
  # p.D.fit.midchannel.basic.SST
  # 
  # p.I.fit.offshore.basic.SST
  # p.D.fit.offshore.basic.SST
  # 
  # p.fit.nearshore.basic.SST
  # p.fit.midchannel.basic.SST
  # p.fit.offshore.basic.SST
  # 
  ############################## sandbox for manually tweaking zeta/eta ##################################
  
  # 27 march 2025
  #   - these inputs were interesting for Midchannel (approach 3), with 0 to 1 b/g and 0 to 0.1 z/e bounds:
  #   - > tab.midchannel
  # # A tibble: 1 × 4
  # beta gamma  zeta   eta
  # <dbl> <dbl> <dbl> <dbl>
  #   1  0.99  0.99   0.1  0.08
  
  # try just playing around with these values ... but also could I just try fitting a straight up scalar for z/e ?? how would that incorporate temperature? maybe try literally rescaling temperature b/t 23 and 33 to 0 to 1 and fit that ?
  
  
  # params.SST[[1]] #off
  # params.SST[[2]] #mid
  # params.SST[[3]] #near
  
  
  #sandbox conditions
  # mid
  # beta.sand.SST.midchannel = 0.647 #these work well for midchannel
  # gamma.sand.SST.midchannel = 0.686
  # zeta.sand.SST.midchannel = 0.077
  # eta.sand.SST.midchannel = 0.00383
  # beta.sand.SST.midchannel = 0.99 #these work well for midchannel
  # gamma.sand.SST.midchannel = 0.99
  # zeta.sand.SST.midchannel = 0.1
  # eta.sand.SST.midchannel = 0.08
  beta.sand.SST.midchannel = 0.98629746 #these work well for midchannel. VERY sensitive to rounding
  gamma.sand.SST.midchannel = 0.98478648
  zeta.sand.SST.midchannel = 0.06909903
  eta.sand.SST.midchannel = 0.04377232
  # off
  # beta.sand.SST.offshore = 0.647 #these do not translate well from midchannel
  # gamma.sand.SST.offshore = 0.686
  # zeta.sand.SST.offshore = 0.077
  # eta.sand.SST.offshore = 0.00383
  beta.sand.SST.offshore = 1.6767164 #these work well for offshore. sensitive to rounding
  gamma.sand.SST.offshore = 1.5338912
  zeta.sand.SST.offshore = 0.1312599
  eta.sand.SST.offshore = 0.1614401
  # near
  # beta.sand.SST.nearshore = 0.647 #these do not translate well from midchannel
  # gamma.sand.SST.nearshore = 0.686
  # zeta.sand.SST.nearshore = 0.077
  # eta.sand.SST.nearshore = 0.00383
  beta.sand.SST.nearshore = 1.1704596 #these are good for nearshore - sort of. not really waves. sensitive to rounding
  gamma.sand.SST.nearshore = 0.9287771
  zeta.sand.SST.nearshore = 0.1232244
  eta.sand.SST.nearshore = 0.1933079
  
  SIR.sand.no_cover.SST = function(t,y,p,SST,DHW){ # 'p' is parameters or params
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
      current_SST = SST$SST_smoothed[closest_index]
      # current_DHW = DHW$DHW[closest_index]
      
      #host density-null conditions
      # transmission_modifier = 1
      # removal_rate = g
      # transmission_rate = b * transmission_modifier
      removal_rate = g
      transmission_rate = b
      
      #STOPPING POINT - 27 MAR 2025
      #  - I think a real issue here is that my equation behaviors are pretty weird
      
      # # APPROACH 4: no threshold
      # # NOTE - issue is that low temperatures actually have the decreasing effect here. high ones should instead
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
      
      #APPROACH 3: threshold
      # NOTE - issue here is that there is no meaningful inflection point - maybe?
      #         - or something else like boundary conditions
      zeta_effective <- z / (1 - z)
      eta_effective <- e / (1 - e)
      beta_scaling_factor <- 1 / (1 + exp(zeta_effective * (current_SST - 30.5)))
      gamma_scaling_factor <- 1 / (1 + exp(eta_effective * (current_SST - 30.5)))
      transmission_rate <- transmission_rate * beta_scaling_factor
      removal_rate <- removal_rate * gamma_scaling_factor
      
      # #APPROACH 0: halting transmission/removal
      # if(current_SST > 30.5){
      #   removal_rate = 0.01
      # }
      
      dS.dt = -transmission_rate * S * I / N
      dI.dt = transmission_rate * S * I / N - removal_rate * I
      dR.dt = removal_rate * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  #midchannel
  site.loop = 'mid'
  days.sand.SST.midchannel <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
    filter(site == site.loop) %>%
    pull(days) %>%
    unlist()
  
  days.obs.sand.SST.midchannel = na.omit(days.sand.SST.midchannel) %>%
    as.numeric()
  
  days.model.sand.SST.midchannel = SST_sites_extended %>%
    filter(site == site.loop) %>%
    pull(time) %>%
    max() %>%
    seq(from = 0, to = .)
  
  SST_df.sand.SST.midchannel = SST_sites_extended %>%
    filter(site == site.loop) %>%
    select(date, time, SST_smoothed)
  
  DHW_df.sand.SST.midchannel <- DHW_sites_extended %>%
    filter(site == site.loop) %>%
    select(date, time, DHW)
  
  SST_values.sand.SST.midchannel = SST_df.sand.SST.midchannel %>%
    pull(SST_smoothed)
  
  DHW_values.sand.SST.midchannel = DHW_df.sand.SST.midchannel %>%
    pull(DHW)
  
  N.sand.SST.midchannel = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  cover.sand.SST.midchannel = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #same as above
    pull(cover.site)
  
  #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
  inftiss.sand.SST.midchannel = obs.total.figures_SST %>%
    filter(Site == site.loop, Compartment == "Infected") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)    
  
  remtiss.sand.SST.midchannel = obs.total.figures_SST %>%
    filter(Site == site.loop, Compartment == "Recovered") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
  first_valid_idx <- which(!is.na(days.sand.SST.midchannel))[1] #find the first non-NA index
  inftiss.sand.SST.midchannel <- inftiss.sand.SST.midchannel[first_valid_idx:length(inftiss.sand.SST.midchannel)]
  remtiss.sand.SST.midchannel <- remtiss.sand.SST.midchannel[first_valid_idx:length(remtiss.sand.SST.midchannel)]
  
  #initial conditions
  I.tiss = inftiss.sand.SST.midchannel[1] #first non-NA & non-zero infection entry
  # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.tiss = N.sand.SST.midchannel - I.tiss
  R.tiss = 0
  cover.sand.SST.midchannel = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #same as above
    pull(cover.site)
  
  initial_state.tiss.sand.SST.midchannel = c(S.tiss, I.tiss, R.tiss, N.sand.SST.midchannel, cover.sand.SST.midchannel)
  
  #run midchannel model
  SIR.out.tiss.sand.SST.midchannel = data.frame(ode(c(S = initial_state.tiss.sand.SST.midchannel[1], I = initial_state.tiss.sand.SST.midchannel[2], R = initial_state.tiss.sand.SST.midchannel[3]),
                                days.model.sand.SST.midchannel, SIR_project, c(b = beta.sand.SST.midchannel, g = gamma.sand.SST.midchannel,
                                                           N = initial_state.tiss.sand.SST.midchannel[4],
                                                           z = zeta.sand.SST.midchannel,
                                                           e = eta.sand.SST.midchannel,
                                                           l = lambda,
                                                           C = initial_state.tiss.sand.SST.midchannel[5]),
                                SST = SST_df.sand.SST.midchannel,
                                DHW = DHW_df.sand.SST.midchannel))

  output.basic.sand.SST.midchannel = SIR.out.tiss.sand.SST.midchannel
  
  output.basic.sand.SST.midchannel = pivot_longer(output.basic.sand.SST.midchannel, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.sand.SST.midchannel)[1] = 'days.model'
  colnames(output.basic.sand.SST.midchannel)[3] = 'tissue'
  
  tab.sand.SST.midchannel = tibble(round(beta.sand.SST.midchannel, 2), round(gamma.sand.SST.midchannel, 2),
                          round(zeta.sand.SST.midchannel, 2), round(eta.sand.SST.midchannel, 2))
  names(tab.sand.SST.midchannel) = c('beta', 'gamma', 'zeta', 'eta')
  
  p.fit.sand.SST.midchannel = ggplot(data = output.basic.sand.SST.midchannel, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total.figures_SST %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.sand.SST.midchannel$days.model), y = min(output.basic.sand.SST.midchannel$tissue)*0.7, label = list(tab.sand.SST.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.sand.SST.midchannel = ggplot(data = output.basic.sand.SST.midchannel %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.sand.SST.midchannel = ggplot(data = output.basic.sand.SST.midchannel %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.sand.SST.midchannel = ggplot(data = output.basic.sand.SST.midchannel %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  #offshore
  site.loop = 'off'
  days.sand.SST.offshore <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
    filter(site == site.loop) %>%
    pull(days) %>%
    unlist()
  
  days.obs.sand.SST.offshore = na.omit(days.sand.SST.offshore) %>%
    as.numeric()
  
  days.model.sand.SST.offshore = SST_sites_extended %>%
    filter(site == site.loop) %>%
    pull(time) %>%
    max() %>%
    seq(from = 0, to = .)
  
  SST_df.sand.SST.offshore = SST_sites_extended %>%
    filter(site == site.loop) %>%
    select(date, time, SST_smoothed)
  
  DHW_df.sand.SST.offshore <- DHW_sites_extended %>%
    filter(site == site.loop) %>%
    select(date, time, DHW)
  
  SST_values.sand.SST.offshore = SST_df.sand.SST.offshore %>%
    pull(SST_smoothed)
  
  DHW_values.sand.SST.offshore = DHW_df.sand.SST.offshore %>%
    pull(DHW)
  
  N.sand.SST.offshore = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  cover.sand.SST.offshore = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #same as above
    pull(cover.site)
  
  #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
  inftiss.sand.SST.offshore = obs.total.figures_SST %>%
    filter(Site == site.loop, Compartment == "Infected") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)    
  
  remtiss.sand.SST.offshore = obs.total.figures_SST %>%
    filter(Site == site.loop, Compartment == "Recovered") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
  first_valid_idx <- which(!is.na(days.sand.SST.offshore))[1] #find the first non-NA index
  inftiss.sand.SST.offshore <- inftiss.sand.SST.offshore[first_valid_idx:length(inftiss.sand.SST.offshore)]
  remtiss.sand.SST.offshore <- remtiss.sand.SST.offshore[first_valid_idx:length(remtiss.sand.SST.offshore)]
  
  #initial conditions
  I.tiss = inftiss.sand.SST.offshore[1] #first non-NA & non-zero infection entry
  # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.tiss = N.sand.SST.offshore - I.tiss
  R.tiss = 0
  cover.sand.SST.offshore = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #same as above
    pull(cover.site)
  
  initial_state.tiss.sand.SST.offshore = c(S.tiss, I.tiss, R.tiss, N.sand.SST.offshore, cover.sand.SST.offshore)
  
  #run offshore model
  SIR.out.tiss.sand.SST.offshore = data.frame(ode(c(S = initial_state.tiss.sand.SST.offshore[1], I = initial_state.tiss.sand.SST.offshore[2], R = initial_state.tiss.sand.SST.offshore[3]),
                                                  days.model.sand.SST.offshore, SIR_project, c(b = beta.sand.SST.offshore, g = gamma.sand.SST.offshore,
                                                                                               N = initial_state.tiss.sand.SST.offshore[4],
                                                                                               z = zeta.sand.SST.offshore,
                                                                                               e = eta.sand.SST.offshore,
                                                                                               l = lambda,
                                                                                               C = initial_state.tiss.sand.SST.offshore[5]),
                                                  SST = SST_df.sand.SST.offshore,
                                                  DHW = DHW_df.sand.SST.offshore))
  
  output.basic.sand.SST.offshore = SIR.out.tiss.sand.SST.offshore
  
  output.basic.sand.SST.offshore = pivot_longer(output.basic.sand.SST.offshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.sand.SST.offshore)[1] = 'days.model'
  colnames(output.basic.sand.SST.offshore)[3] = 'tissue'
  
  tab.sand.SST.offshore = tibble(round(beta.sand.SST.offshore, 2), round(gamma.sand.SST.offshore, 2),
                                 round(zeta.sand.SST.offshore, 2), round(eta.sand.SST.offshore, 2))
  names(tab.sand.SST.offshore) = c('beta', 'gamma', 'zeta', 'eta')
  
  p.fit.sand.SST.offshore = ggplot(data = output.basic.sand.SST.offshore, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total.figures_SST %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.sand.SST.offshore$days.model), y = min(output.basic.sand.SST.offshore$tissue)*0.7, label = list(tab.sand.SST.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.sand.SST.offshore = ggplot(data = output.basic.sand.SST.offshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.sand.SST.offshore = ggplot(data = output.basic.sand.SST.offshore %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.sand.SST.offshore = ggplot(data = output.basic.sand.SST.offshore %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Nearshore code follows the same pattern
  site.loop = 'near'
  days.sand.SST.nearshore <- days_sites_SST %>% # NOTE - make sure this is working right with backtracked patient zero corals
    filter(site == site.loop) %>%
    pull(days) %>%
    unlist()
  
  days.obs.sand.SST.nearshore = na.omit(days.sand.SST.nearshore) %>%
    as.numeric()
  
  days.model.sand.SST.nearshore = SST_sites_extended %>%
    filter(site == site.loop) %>%
    pull(time) %>%
    max() %>%
    seq(from = 0, to = .)
  
  SST_df.sand.SST.nearshore = SST_sites_extended %>%
    filter(site == site.loop) %>%
    select(date, time, SST_smoothed)
  
  DHW_df.sand.SST.nearshore <- DHW_sites_extended %>%
    filter(site == site.loop) %>%
    select(date, time, DHW)
  
  SST_values.sand.SST.nearshore = SST_df.sand.SST.nearshore %>%
    pull(SST_smoothed)
  
  DHW_values.sand.SST.nearshore = DHW_df.sand.SST.nearshore %>%
    pull(DHW)
  
  N.sand.SST.nearshore = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  cover.sand.SST.nearshore = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #same as above
    pull(cover.site)
  
  #sequence of infected & removed SA's for each SusCat within site. remove timepoints before the first infection (zero omit)
  inftiss.sand.SST.nearshore = obs.total.figures_SST %>%
    filter(Site == site.loop, Compartment == "Infected") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)    
  
  remtiss.sand.SST.nearshore = obs.total.figures_SST %>%
    filter(Site == site.loop, Compartment == "Recovered") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim 'inftiss' and 'remtiss' using the same indices as 'days'
  first_valid_idx <- which(!is.na(days.sand.SST.nearshore))[1] #find the first non-NA index
  inftiss.sand.SST.nearshore <- inftiss.sand.SST.nearshore[first_valid_idx:length(inftiss.sand.SST.nearshore)]
  remtiss.sand.SST.nearshore <- remtiss.sand.SST.nearshore[first_valid_idx:length(remtiss.sand.SST.nearshore)]
  
  #initial conditions
  I.tiss = inftiss.sand.SST.nearshore[1] #first non-NA & non-zero infection entry
  # I.tiss = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.tiss = N.sand.SST.nearshore - I.tiss
  R.tiss = 0
  cover.sand.SST.nearshore = susceptible_ref_SST %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #same as above
    pull(cover.site)
  
  initial_state.tiss.sand.SST.nearshore = c(S.tiss, I.tiss, R.tiss, N.sand.SST.nearshore, cover.sand.SST.nearshore)
  
  #run nearshore model
  SIR.out.tiss.sand.SST.nearshore = data.frame(ode(c(S = initial_state.tiss.sand.SST.nearshore[1], I = initial_state.tiss.sand.SST.nearshore[2], R = initial_state.tiss.sand.SST.nearshore[3]),
                                                   days.model.sand.SST.nearshore, SIR_project, c(b = beta.sand.SST.nearshore, g = gamma.sand.SST.nearshore,
                                                                                                 N = initial_state.tiss.sand.SST.nearshore[4],
                                                                                                 z = zeta.sand.SST.nearshore,
                                                                                                 e = eta.sand.SST.nearshore,
                                                                                                 l = lambda,
                                                                                                 C = initial_state.tiss.sand.SST.nearshore[5]),
                                                   SST = SST_df.sand.SST.nearshore,
                                                   DHW = DHW_df.sand.SST.nearshore))
  
  output.basic.sand.SST.nearshore = SIR.out.tiss.sand.SST.nearshore
  
  output.basic.sand.SST.nearshore = pivot_longer(output.basic.sand.SST.nearshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.sand.SST.nearshore)[1] = 'days.model'
  colnames(output.basic.sand.SST.nearshore)[3] = 'tissue'
  
  tab.sand.SST.nearshore = tibble(round(beta.sand.SST.nearshore, 2), round(gamma.sand.SST.nearshore, 2),
                                  round(zeta.sand.SST.nearshore, 2), round(eta.sand.SST.nearshore, 2))
  names(tab.sand.SST.nearshore) = c('beta', 'gamma', 'zeta', 'eta')
  
  p.fit.sand.SST.nearshore = ggplot(data = output.basic.sand.SST.nearshore, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total.figures_SST %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.sand.SST.nearshore$days.model), y = min(output.basic.sand.SST.nearshore$tissue)*0.7, label = list(tab.sand.SST.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.sand.SST.nearshore = ggplot(data = output.basic.sand.SST.nearshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.sand.SST.nearshore = ggplot(data = output.basic.sand.SST.nearshore %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.sand.SST.nearshore = ggplot(data = output.basic.sand.SST.nearshore %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    # geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  ############################## Plots of sandbox ##################################
  
  plot_data_decreasing <- data.frame()
  
  # zeta_values <- seq(0, 1, length.out = 50)
  zeta_values <- (seq(0, 1, length.out = 30))^2
  
  T <- seq(T_min, T_max, length.out = 1000)
  
  for (zeta in zeta_values) {
    modifiers_decreasing <- sapply(T, function(t) temp_effect_decreasing(t, tmin = T_min, tmax = T_max, zeta_bounded = zeta))
    temp_data_decreasing <- data.frame(
      Temperature = T,
      Modifier = modifiers_decreasing,
      Zeta = zeta
    )
    plot_data_decreasing <- rbind(plot_data_decreasing, temp_data_decreasing)
  }
  
  # Manually add fitted zeta and eta
  zeta_select <- sapply(T, function(t) temp_effect_decreasing(t, tmin = T_min, tmax = T_max, zeta_bounded = zeta.sand.SST.midchannel))
  temp <- data.frame(
    Temperature = T,
    Modifier = zeta_select,
    Zeta = zeta.sand.SST.midchannel
  )
  plot_data_decreasing_select <- rbind(plot_data_decreasing, temp)
  
  eta_select <- sapply(T, function(t) temp_effect_decreasing(t, tmin = T_min, tmax = T_max, zeta_bounded = eta.sand.SST.midchannel))
  temp <- data.frame(
    Temperature = T,
    Modifier = eta_select,
    Zeta = eta.sand.SST.midchannel
  )
  plot_data_decreasing_select <- rbind(plot_data_decreasing_select, temp)
  
  # Plot 3: Classic sigmoid decreasing model
  fig2b = ggplot(plot_data_decreasing_select, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(data = plot_data_decreasing, size = 0.6, show.legend = TRUE) +
    geom_line(
      data = subset(plot_data_decreasing_select, Zeta == zeta.sand.SST.midchannel),
      aes(x = Temperature, y = Modifier),
      color = "darkorange2", linetype = 'dashed', size = 0.6, inherit.aes = FALSE
    ) +
    geom_line(
      data = subset(plot_data_decreasing_select, Zeta == eta.sand.SST.midchannel),
      aes(x = Temperature, y = Modifier),
      color = "orange", linetype = 'dashed', size = 0.6, inherit.aes = FALSE
    ) +
    labs(
      x = "Temperature (°C)",
      y = "Temperature scalar"
    ) +
    theme_classic(base_family = 'Georgia') +
    theme(
      # legend.position = "right",
      # legend.position = "bottom",
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.8),
      legend.key.size = unit(0.5, "lines"),  # smaller boxes
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7, color = 'black'),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
      legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    # scale_color_viridis_c(name = expression(zeta~"or"~eta), option = "inferno") + #option - 'D'
    scale_color_viridis_c(name = expression("ζ or η"), option = "inferno") + #option - 'D'
    
    scale_y_continuous(, limits = c(0, 1),
                       expand = c(0, 0),
                       breaks = seq(0, 1, 0.2)) +
    
    scale_x_continuous(labels = label_number(accuracy = 0.1),
                       expand = c(0, 0),
                       breaks = seq(T_min, T_max, 2)) +
    
    # geom_vline(xintercept = 30.5, color = 'grey40', linetype = 'dashed', size = 0.6) +
    # annotate("segment", x = 24.2, xend = 26, y = 0.2, yend = 0.2, 
    #          color = "darkorange2", linetype = "dashed", size = 1) +
    # annotate("segment", x = 24.2, xend = 26, y = 0.1, yend = 0.1, 
    #          color = "orange", linetype = "dashed", size = 1) +
    # annotate("text", x = 22, y = 0.2, label = "ζ = 0.07", hjust = 0, family = 'Georgia', size = 3) +
    # annotate("text", x = 22, y = 0.1, label = "η = 0.04", hjust = 0, family = 'Georgia', size = 3)
  
  geom_vline(xintercept = 30.5, color = 'grey40', linetype = 'dashed', size = 0.6) +
    annotate("segment", x = 25.2, xend = 27, y = 0.25, yend = 0.25, 
             color = "darkorange2", linetype = "dashed", size = 1) +
    annotate("segment", x = 25.2, xend = 27, y = 0.15, yend = 0.15, 
             color = "orange", linetype = "dashed", size = 1) +
    annotate("text", x = 23, y = 0.25, label = "ζ = 0.07", hjust = 0, family = 'Georgia', size = 3) +
    annotate("text", x = 23, y = 0.15, label = "η = 0.04", hjust = 0, family = 'Georgia', size = 3)
  
  
  ################################## export Figure 2b ##################################
  
  # #max dimensions are 7.087 in. wide by 9.45 in. tall (3.35 inches preferred)
  # quartz(h = 3, w = 3.35)
  # 
  # fig2b
  # 
  # # # Save the Quartz output directly as a PDF
  # # quartz.save(file = here("output", "fig2b.pdf"), type = "pdf")
  # #
  # # #ggplot-export to image
  # # ggsave(filename = here("output", "fig2b.png"), device = "png", width = 3.35, height = 3, dpi = 1200)
  # 
  # # Close the Quartz device
  # dev.off()
  # 
  # saveRDS(fig2b, file = here("output", "fig2b.rds"))
  
  ############################## Draft figures of SST for paper ##################################

  # Define styling parameters
  linewidths = 0.75 #0.4 #0.75 is roughly 1 pt. ggplot measures these in mm, not points
  symbsizes = 1.3
  titlesize = 10 #9  #text sizes in ggplot are actually in units of points, when specified using element_text
  textsize = 9
  
  # Data preparation for DEAD tissue plot
  output.basic.sand.SST.midchannel.SST <- output.basic.sand.SST.midchannel %>%
    mutate(date = as.Date("2018-11-09") + days.model)
  
  compartment_dead = "Dead"
  
  # Filter and prepare the data for dead tissue
  SST_data <- DHW.CRW.full
  output.basic.sand.SST.midchannel.SST.scaled.dead <- output.basic.sand.SST.midchannel.SST %>%
    filter(Compartment == compartment_dead) %>%
    mutate(date = as.Date(date), tissue_scaled = tissue)
  
  dead_data <- obs.total[obs.total$Site == 'Midchannel' & obs.total$Compartment == "Dead", ]
  dead_dates <- as.Date(dead_data$date)
  
  # Calculate scaling factor for secondary y-axis (dead)
  tissue_max_dead <- max(output.basic.sand.SST.midchannel.SST.scaled.dead$tissue_scaled)
  sst_min <- 23
  sst_max <- 32.5
  scale_factor_dead <- (sst_max - sst_min) / tissue_max_dead
  
  # Data preparation for INFECTED tissue plot
  compartment_infected = "Infected"
  
  # Filter and prepare the data for infected tissue
  output.basic.sand.SST.midchannel.SST.scaled.infected <- output.basic.sand.SST.midchannel.SST %>%
    filter(Compartment == compartment_infected) %>%
    mutate(date = as.Date(date), tissue_scaled = tissue)
  
  infected_data <- obs.total[obs.total$Site == 'Midchannel' & obs.total$Compartment == "Infected", ]
  infected_dates <- as.Date(infected_data$date)
  
  # Calculate scaling factor for secondary y-axis (infected)
  tissue_max_infected <- max(output.basic.sand.SST.midchannel.SST.scaled.infected$tissue_scaled)
  scale_factor_infected <- (sst_max - sst_min) / tissue_max_infected
  
  # Create the INFECTED tissue plot (Panel A)
  plot_A <- ggplot() +
    # SST line
    geom_line(data = SST_data, 
              aes(x = date, y = SST.90th_HS), 
              color = "#E69F00", 
              linewidth = linewidths) +
    
    # Tissue line (scaled to primary y-axis)
    geom_line(data = output.basic.sand.SST.midchannel.SST.scaled.infected,
              aes(x = date, y = sst_min + tissue_scaled * scale_factor_infected),
              color = "black",
              linewidth = linewidths) +
    
    # Infected tissue points (scaled to primary y-axis)
    geom_point(data = data.frame(date = infected_dates, tissue = infected_data$tissue),
               aes(x = date, y = sst_min + tissue * scale_factor_infected),
               color = "black",
               shape = 17,  # Triangle
               size = symbsizes) +
    
    # Vertical dashed line
    geom_vline(xintercept = as.Date("2019-12-06"), 
               color = "grey40", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Horizontal dashed line for SST threshold
    geom_hline(yintercept = SST_threshold, 
               color = "red", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Set scales
    scale_x_date(limits = c(as.Date("2018-01-01"), as.Date("2020-12-30")),
                 date_labels = "%Y",
                 date_breaks = "1 year") +
    
    scale_y_continuous(
      name = "SST (°C)",
      limits = c(23, 32.5),
      breaks = seq(24, 32, by = 2),
      # Secondary y-axis for infected tissue
      sec.axis = sec_axis(trans = ~ (. - sst_min) / scale_factor_infected,
                          name = "Infected tissue (m²)",
                          labels = function(x) formatC(x, format = "f", digits = 3))
    ) +
    
    # Labels and theme
    labs(x = "Year") +
    
    theme_classic(base_size = 14, base_family = "Georgia") +
    
    theme(
      axis.title = element_text(size = titlesize),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.title.y.right = element_text(angle = 90),
      legend.text = element_text(size = textsize),
      legend.title = element_text(size = titlesize),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(t = 8, r = 15, b = 8, l = 8),
      # Position legends inside the plot
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.65),
      legend.box = "vertical",
      legend.margin = margin(5, 5, 5, 5),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.3)
    ) +
    
    # Create proper legend in top left (move left and make much wider)
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(fill = alpha("white", 0.8), col = "black")),
      xmin = as.Date("2017-12-01"), xmax = as.Date("2019-04-15"),
      ymin = 30.8, ymax = 32.5
    ) +
    
    # Legend text and lines (fixed segment positioning)
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 32.1, yend = 32.1, color = "#E69F00", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 32.1, 
             label = "SST (°C)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0) +
    
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 31.3, yend = 31.3, color = "black", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 31.3, 
             label = "Tissue (m²)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0)
  
  # Create the DEAD tissue plot (Panel B)
  plot_B <- ggplot() +
    # SST line
    geom_line(data = SST_data, 
              aes(x = date, y = SST.90th_HS), 
              color = "#E69F00", 
              linewidth = linewidths) +
    
    # Tissue line (scaled to primary y-axis)
    geom_line(data = output.basic.sand.SST.midchannel.SST.scaled.dead,
              aes(x = date, y = sst_min + tissue_scaled * scale_factor_dead),
              color = "black",
              linewidth = linewidths) +
    
    # Tissue points (scaled to primary y-axis)
    geom_point(data = data.frame(date = dead_dates, tissue = dead_data$tissue),
               aes(x = date, y = sst_min + tissue * scale_factor_dead),
               color = "black",
               shape = 15,  # Square
               size = symbsizes) +
    
    # Vertical dashed line
    geom_vline(xintercept = as.Date("2019-12-06"), 
               color = "grey40", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Horizontal dashed line for SST threshold
    geom_hline(yintercept = SST_threshold, 
               color = "red", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Set scales
    scale_x_date(limits = c(as.Date("2018-01-01"), as.Date("2020-12-30")),
                 date_labels = "%Y",
                 date_breaks = "1 year") +
    
    scale_y_continuous(
      name = "SST (°C)",
      limits = c(23, 32.5),
      breaks = seq(24, 32, by = 2),
      # Secondary y-axis for dead tissue
      sec.axis = sec_axis(trans = ~ (. - sst_min) / scale_factor_dead,
                          name = "Removed tissue (m²)",
                          labels = function(x) formatC(x, format = "f", digits = 2))
    ) +
    
    # Labels and theme
    labs(x = "Year") +
    
    theme_classic(base_size = 14, base_family = "Georgia") +
    
    theme(
      axis.title = element_text(size = titlesize),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.title.y.right = element_text(angle = 90),
      legend.text = element_text(size = textsize),
      legend.title = element_text(size = titlesize),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(t = 8, r = 15, b = 8, l = 8),
      # Position legends inside the plot
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.65),
      legend.box = "vertical",
      legend.margin = margin(5, 5, 5, 5),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.3)
    ) +
    
    # Create proper legend in top left (move left and make much wider)
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(fill = alpha("white", 0.8), col = "black")),
      xmin = as.Date("2017-12-01"), xmax = as.Date("2019-04-15"),
      ymin = 30.8, ymax = 32.5
    ) +
    
    # Legend text and lines (fixed segment positioning)
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 32.1, yend = 32.1, color = "#E69F00", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 32.1, 
             label = "SST (°C)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0) +
    
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 31.3, yend = 31.3, color = "black", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 31.3, 
             label = "Tissue (m²)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0)
  
    # # Version where the panel B legend box is bottom right
    # # Create proper legend in bottom right (move right and center text better)
    # annotation_custom(
    #   grob = grid::rectGrob(gp = grid::gpar(fill = alpha("white", 0.8), col = "black")),
    #   xmin = as.Date("2019-09-01"), xmax = as.Date("2021-01-15"),
    #   ymin = 22.9, ymax = 24.6 #22.9
    # ) +
    # 
    # # Legend text and lines (better centered in box)
    # annotate("segment", x = as.Date("2019-10-01"), xend = as.Date("2019-12-01"),
    #          y = 24.15, yend = 24.15, color = "#E69F00", linewidth = linewidths) +
    # annotate("text", x = as.Date("2019-12-15"), y = 24.15,
    #          label = "SST (°C)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0) +
    # 
    # annotate("segment", x = as.Date("2019-10-01"), xend = as.Date("2019-12-01"),
    #          y = 23.35, yend = 23.35, color = "black", linewidth = linewidths) +
    # annotate("text", x = as.Date("2019-12-15"), y = 23.35,
    #          label = "Tissue (m²)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0)
  
  # Combine plots using patchwork (Infected first, then Dead)
  combined_plot_fig4 <- plot_A + plot_B + plot_annotation(tag_levels = 'A')
  
  
  
  
  load(here("output/tables_figures_workspace.RData"))
  
  # Prepare the fitted data with proper date conversion
  data_fig3_midchannel <- data_fig3 %>%
    filter(Site == "Midchannel", Host == "Single-host", Type == "Fitted") %>%
    mutate(date = as.Date("2018-11-09") + days.model)
  
  # Create Panel A with proper layer order (fitted lines behind points)
  plot_A_with_fitted <- ggplot() +
    # SST line
    geom_line(data = SST_data, 
              aes(x = date, y = SST.90th_HS), 
              color = "#E69F00", 
              linewidth = linewidths) +
    
    # Tissue line (scaled to primary y-axis)
    geom_line(data = output.basic.sand.SST.midchannel.SST.scaled.infected,
              aes(x = date, y = sst_min + tissue_scaled * scale_factor_infected),
              color = "black",
              linewidth = linewidths) +
    
    # Add gray fitted line for Infected compartment (BEHIND points)
    geom_line(data = data_fig3_midchannel %>% filter(Compartment == "Infected"),
              aes(x = date, y = sst_min + tissue * scale_factor_infected),
              color = "gray40", linewidth = linewidths) +
    
    # Infected tissue points (ABOVE fitted lines)
    geom_point(data = data.frame(date = infected_dates, tissue = infected_data$tissue),
               aes(x = date, y = sst_min + tissue * scale_factor_infected),
               color = "black",
               shape = 17,  # Triangle
               size = symbsizes) +
    
    # Vertical dashed line
    geom_vline(xintercept = as.Date("2019-12-06"), 
               color = "grey40", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Horizontal dashed line for SST threshold
    geom_hline(yintercept = SST_threshold, 
               color = "red", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Set scales
    scale_x_date(limits = c(as.Date("2018-01-01"), as.Date("2020-12-30")),
                 date_labels = "%Y",
                 date_breaks = "1 year") +
    
    scale_y_continuous(
      name = "SST (°C)",
      limits = c(23, 32.5),
      breaks = seq(24, 32, by = 2),
      # Secondary y-axis for infected tissue
      sec.axis = sec_axis(trans = ~ (. - sst_min) / scale_factor_infected,
                          name = "Infected tissue (m²)",
                          labels = function(x) formatC(x, format = "f", digits = 3))
    ) +
    
    # Labels and theme
    labs(x = "Year") +
    
    theme_classic(base_size = 14, base_family = "Georgia") +
    
    theme(
      axis.title = element_text(size = titlesize),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.title.y.right = element_text(angle = 90),
      legend.text = element_text(size = textsize),
      legend.title = element_text(size = titlesize),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(t = 8, r = 15, b = 8, l = 8),
      # Position legends inside the plot
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.65),
      legend.box = "vertical",
      legend.margin = margin(5, 5, 5, 5),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.3)
    ) +
    
    # Create proper legend in top left (your updated positioning)
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(fill = alpha("white", 0.8), col = "black")),
      xmin = as.Date("2017-12-01"), xmax = as.Date("2019-04-15"),
      ymin = 30.8, ymax = 32.5
    ) +
    
    # Legend text and lines (your updated positioning)
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 32.1, yend = 32.1, color = "#E69F00", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 32.1, 
             label = "SST (°C)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0) +
    
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 31.3, yend = 31.3, color = "black", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 31.3, 
             label = "Tissue (m²)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0)
  
  # Create Panel B with proper layer order (fitted lines behind points)
  plot_B_with_fitted <- ggplot() +
    # SST line
    geom_line(data = SST_data, 
              aes(x = date, y = SST.90th_HS), 
              color = "#E69F00", 
              linewidth = linewidths) +
    
    # Tissue line (scaled to primary y-axis)
    geom_line(data = output.basic.sand.SST.midchannel.SST.scaled.dead,
              aes(x = date, y = sst_min + tissue_scaled * scale_factor_dead),
              color = "black",
              linewidth = linewidths) +
    
    # Add gray fitted line for Removed compartment (BEHIND points)
    geom_line(data = data_fig3_midchannel %>% filter(Compartment == "Removed"),
              aes(x = date, y = sst_min + tissue * scale_factor_dead),
              color = "gray40", linewidth = linewidths) +
    
    # Tissue points (ABOVE fitted lines)
    geom_point(data = data.frame(date = dead_dates, tissue = dead_data$tissue),
               aes(x = date, y = sst_min + tissue * scale_factor_dead),
               color = "black",
               shape = 15,  # Square
               size = symbsizes) +
    
    # Vertical dashed line
    geom_vline(xintercept = as.Date("2019-12-06"), 
               color = "grey40", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Horizontal dashed line for SST threshold
    geom_hline(yintercept = SST_threshold, 
               color = "red", 
               linetype = "dashed", 
               linewidth = 1) +
    
    # Set scales
    scale_x_date(limits = c(as.Date("2018-01-01"), as.Date("2020-12-30")),
                 date_labels = "%Y",
                 date_breaks = "1 year") +
    
    scale_y_continuous(
      name = "SST (°C)",
      limits = c(23, 32.5),
      breaks = seq(24, 32, by = 2),
      # Secondary y-axis for dead tissue
      sec.axis = sec_axis(trans = ~ (. - sst_min) / scale_factor_dead,
                          name = "Removed tissue (m²)",
                          labels = function(x) formatC(x, format = "f", digits = 2))
    ) +
    
    # Labels and theme
    labs(x = "Year") +
    
    theme_classic(base_size = 14, base_family = "Georgia") +
    
    theme(
      axis.title = element_text(size = titlesize),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.title.y.right = element_text(angle = 90),
      legend.text = element_text(size = textsize),
      legend.title = element_text(size = titlesize),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(t = 8, r = 15, b = 8, l = 8),
      # Position legends inside the plot
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.65),
      legend.box = "vertical",
      legend.margin = margin(5, 5, 5, 5),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.3)
    ) +
    
    # Create proper legend in top left (matching Panel A)
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(fill = alpha("white", 0.8), col = "black")),
      xmin = as.Date("2017-12-01"), xmax = as.Date("2019-04-15"),
      ymin = 30.8, ymax = 32.5
    ) +
    
    # Legend text and lines (matching Panel A)
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 32.1, yend = 32.1, color = "#E69F00", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 32.1, 
             label = "SST (°C)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0) +
    
    annotate("segment", x = as.Date("2018-01-01"), xend = as.Date("2018-03-01"), 
             y = 31.3, yend = 31.3, color = "black", linewidth = linewidths) +
    annotate("text", x = as.Date("2018-03-15"), y = 31.3, 
             label = "Tissue (m²)", color = "black", size = textsize/ggplot2::.pt, family = "Georgia", hjust = 0)
  
  # Update the combined plot with fitted lines
  combined_plot_fig4_final <- plot_A_with_fitted + plot_B_with_fitted + plot_annotation(tag_levels = 'A')
  
  # Display the updated plot
  print(combined_plot_fig4_final)  
  
  ############################## Export Figure 4 ##################################
  
  #max dimensions are 7.087 in. wide by 9.45 in. tall (3.35 inches preferred)
  quartz(h = 3, w = 7.087)
  
  # combined_plot
  combined_plot_fig4_final

  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "fig4.pdf"), type = "pdf")

  #ggplot-export to image
  ggsave(filename = here("output", "fig4.png"), device = "png", width = 7.087, height = 3, dpi = 1200)

  # Close the Quartz device
  dev.off()
  
  # saveRDS(combined_plot_fig4, file = here("output", "fig4.rds"))
  
  ################################## Save output ##################################
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "optimize_SST_model_workspace.RData"))