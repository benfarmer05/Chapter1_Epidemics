  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(ggplot2)
  library(dplyr)

  ################################## Set-up ##################################
  
  #import workspace from upstream script
  # load(here("output", "plots_multi_workspace.RData"))
  load(here("output/tables_figures_workspace.RData"))
  
  ################################## Full range of SST ##################################
  
  # # Define parameters
  # SST_threshold <- 30.5  # Temperature threshold
  # z_values <- c(0.5, 1, 2.5, 3)  # Different z values for transmission effect
  # e_values <- c(0.0001, 1, 2, 3)  # Different e values for mortality effect
  # SST_range <- seq(26, 32, by = 0.1)  # SST from 29 to 32 degrees
  # 
  # # Set up plot
  # plot(NULL, xlim = range(SST_range), ylim = c(0, 1), 
  #      xlab = "Sea Surface Temperature (°C)", 
  #      ylab = "Scaling Factor", 
  #      main = "Effect of SST on Transmission and Mortality Rates")
  # 
  # # Define colors and line types
  # colors <- c("blue", "red", "green", "purple")
  # linetypes <- c(1, 2)  # Solid for transmission, dashed for mortality
  # 
  # # Plot transmission rate curves
  # for (i in seq_along(z_values)) {
  #   z <- z_values[i]
  #   scaling_factor_z <- 1 / (1 + exp(-z * (SST_range - SST_threshold)))
  #   lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2, lty = linetypes[1])
  # }
  # 
  # # Plot mortality rate curves
  # for (i in seq_along(e_values)) {
  #   e <- e_values[i]
  #   scaling_factor_e <- 1 / (1 + exp(-e * (SST_range - SST_threshold)))
  #   lines(SST_range, scaling_factor_e, col = colors[i], lwd = 2, lty = linetypes[2])
  # }
  # 
  # # Add legend
  # legend("topleft", legend = c(paste("Transmission (z) =", z_values), 
  #                              paste("Mortality (e) =", e_values)), 
  #        col = rep(colors, 2), lwd = 2, lty = rep(linetypes, each = length(z_values)))
  
  ################################## SST above threshold ##################################
  
  # # Define parameters
  # SST_threshold <- 30.5  # Temperature threshold
  # z_values <- c(0.63, 2.98, 3)  # Different z values for transmission effect
  # e_values <- c(0.0001, 1.89, 2.57)  # Different e values for mortality effect
  # SST_range <- seq(30.5, 32, by = 0.1)  # SST from 30.5°C onward
  # 
  # # Set up plot
  # plot(NULL, xlim = range(SST_range), ylim = c(0, 1), 
  #      xlab = "Sea Surface Temperature (°C)", 
  #      ylab = "Scaling Factor", 
  #      main = "Effect of SST on Transmission and Mortality Rates")
  # 
  # # Define colors and line types
  # colors <- c("blue", "red", "green")
  # linetypes <- c(1, 2)  # Solid for transmission, dashed for mortality
  # 
  # # Plot transmission rate curves (only for SST > 30.5)
  # for (i in seq_along(z_values)) {
  #   z <- z_values[i]
  #   scaling_factor_z <- 1 / (1 + exp(-z * (SST_range - SST_threshold)))
  #   lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2, lty = linetypes[1])
  # }
  # 
  # # Plot mortality rate curves (only for SST > 30.5)
  # for (i in seq_along(e_values)) {
  #   e <- e_values[i]
  #   scaling_factor_e <- 1 / (1 + exp(-e * (SST_range - SST_threshold)))
  #   lines(SST_range, scaling_factor_e, col = colors[i], lwd = 2, lty = linetypes[2])
  # }
  # 
  # # Add legend
  # legend("topleft", legend = c(paste("Transmission (z) =", z_values),
  #                              paste("Mortality (e) =", e_values)),
  #        col = rep(colors, 2), lwd = 2, lty = rep(linetypes, each = length(z_values)))
  
  ################################## Wide range ##################################
  
  # # Define parameters
  # SST_threshold <- 30.5  # Temperature threshold
  # z_values <- seq(0.001, 10, length.out = 8)  # 8 different z values in range 1 to 3
  # e_values <- seq(0.001, 10, length.out = 8)  # 8 different e values in range 1 to 3
  # SST_range <- seq(30.5, 34, by = 0.01)  # SST from 30.5°C onward
  # 
  # # Set up plot
  # plot(NULL, xlim = range(SST_range), ylim = c(0, 1), 
  #      xlab = "Sea Surface Temperature (°C)", 
  #      ylab = "Scaling Factor", 
  #      main = "Effect of SST on Transmission and Mortality Rates")
  # 
  # # Define colors and line types
  # colors <- colorRampPalette(c("blue", "red"))(length(z_values))  # Gradient color scheme
  # linetypes <- c(1, 2)  # Solid for transmission, dashed for mortality
  # 
  # # Plot transmission rate curves (only for SST > 30.5)
  # for (i in seq_along(z_values)) {
  #   z <- z_values[i]
  #   scaling_factor_z <- 1 / (1 + exp(-z * (SST_range - SST_threshold)))
  #   lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2, lty = linetypes[1])
  # }
  # 
  # # Plot mortality rate curves (only for SST > 30.5)
  # for (i in seq_along(e_values)) {
  #   e <- e_values[i]
  #   scaling_factor_e <- 1 / (1 + exp(-e * (SST_range - SST_threshold)))
  #   lines(SST_range, scaling_factor_e, col = colors[i], lwd = 2, lty = linetypes[2])
  # }
  # 
  # # Add legend
  # legend("bottomright", legend = paste("Rates (z, e) =", round(z_values, 2)), 
  #        col = colors, lwd = 2, lty = 1)

  ################################## Simplified function, 0 to 1 ##################################
  
  # # Define parameters
  # SST_threshold <- 30.5  # Threshold temperature
  # zeta_values <- seq(0, 1, by = 0.25)  # Varying zeta from 0 to 1
  # SST_range <- seq(30.5, 32, by = 0.1)  # SST only above threshold
  # 
  # # Set up plot
  # plot(NULL, xlim = range(SST_range), ylim = c(0, 1), 
  #      xlab = "Sea Surface Temperature (°C)", 
  #      ylab = "Scaling Factor", 
  #      main = "Effect of SST on Transmission & Mortality Rates")
  # 
  # # Colors for different zeta values
  # colors <- rainbow(length(zeta_values))
  # 
  # # Loop through different zeta values
  # for (i in seq_along(zeta_values)) {
  #   zeta <- zeta_values[i]
  #   scaling_factor <- 1 - zeta * (1 / (1 + exp(-10 * (SST_range - SST_threshold))))
  #   lines(SST_range, scaling_factor, col = colors[i], lwd = 2)
  # }
  # 
  # # Add legend
  # legend("topright", legend = paste("zeta =", zeta_values), 
  #        col = colors, lwd = 2)

  ################################## Maybe better simplified function, 0 to 1 ##################################
  
  # NOTE - 13 march 2025
  # - i still don't think this is quite right. Even at high z values, as temperature rises, the scaling factor decreases LESS quickly - and it never actually hits 0 (full depression of epidemic rates)
  
  # Define parameters
  SST_threshold <- 30.5  # Temperature threshold
  z_values <- seq(0, 1, length.out = 8)  # Zeta values including 0
  SST_range <- seq(30.5, 34, by = 0.01)  # SST range
  
  # Set up plot
  plot(NULL, xlim = range(SST_range), ylim = c(0, 1), 
       xlab = "Sea Surface Temperature (°C)", 
       ylab = "Scaling Factor", 
       main = "Effect of SST on Transmission and Mortality Rates")
  
  # Define colors
  colors <- colorRampPalette(c("blue", "red"))(length(z_values)) 
  
  # Adjusted scaling function
  for (i in seq_along(z_values)) {
    z <- z_values[i]
    scaling_factor_z <- exp(-z * (SST_range - SST_threshold))  # Exponential decay
    lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2)
  }
  
  # Add legend
  legend("bottomright", legend = paste("Rates (z, e) =", round(z_values, 2)),
         col = colors, lwd = 2)
  
  
  
  
  ################################## Behavior with no threshold ##################################
  
  # # Define parameters
  # SST_threshold <- 0  # Temperature threshold
  # z_values <- seq(0, 1, length.out = 1000)  # Zeta values including 0
  # SST_range <- seq(23, 33, by = 0.01)  # SST range
  # 
  # # Set up plot
  # plot(NULL, xlim = range(SST_range), ylim = c(0, 1), 
  #      xlab = "Sea Surface Temperature (°C)", 
  #      ylab = "Scaling Factor", 
  #      main = "Effect of SST on Transmission and Mortality Rates")
  # 
  # # Define colors
  # colors <- colorRampPalette(c("blue", "red"))(length(z_values)) 
  # 
  # # Adjusted scaling function
  # for (i in seq_along(z_values)) {
  #   z <- z_values[i]
  #   scaling_factor_z <- exp(-z * (SST_range - SST_threshold))  # Exponential decay
  #   lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2)
  # }
  
  # Define parameters
  tempmin = 23
  tempmax = 33
  SST_range = tempmax - tempmin
  SST_vals <- seq(23, 33, by = 0.01)  # SST range from 23°C to 33°C
  z_values <- seq(0, 1, length.out = 8)  # Zeta values including 0
  
  # Set up plot
  plot(NULL, xlim = range(SST_vals), ylim = c(0, 1),
  # plot(NULL, xlim = c(23, 35), ylim = c(0, 1), 
       xlab = "Sea Surface Temperature (°C)", 
       ylab = "Scaling Factor", 
       main = "Effect of SST on Transmission and Mortality Rates (No Threshold)")
  
  # Define colors
  colors <- colorRampPalette(c("blue", "red"))(length(z_values))
  
  # Adjusted scaling function with gradual temperature modulation
  for (i in seq_along(z_values)) {
    z <- z_values[i]
    scaling_factor_z <- exp(-z * (SST_vals - tempmin) / (SST_range))  # Gradual modulation
    # scaling_factor_z[scaling_factor_z > 1] <- 1  # Cap at 1 to avoid overmodulation
    lines(SST_vals, scaling_factor_z, col = colors[i], lwd = 2)
  }
  
  # # Add legend
  # legend("bottomright", legend = paste("Rates (z, e) =", round(z_values, 2)),
  #        col = colors, lwd = 2)
  
  
  
  
  ################################## Plot SST & DHW through time ##################################
  
  # Plot SST over time
  ggplot(SST_df, aes(x = date, y = SST)) +
    geom_line(color = "blue") +
    labs(title = "Sea Surface Temperature (SST) Over Time",
         x = "Date",
         y = "SST (°C)") +
    theme_minimal()
  
  # Merge the two datasets by date
  SST_DHW_df <- full_join(SST_df, DHW_df, by = c("date", "time"))
  
  # Plot SST and DHW over time
  ggplot(SST_DHW_df, aes(x = date)) +
    geom_line(aes(y = SST, color = "SST"), size = 1) +
    geom_line(aes(y = DHW, color = "DHW"), size = 1) +
    scale_color_manual(values = c("SST" = "#E69F00", "DHW" = "#0072B2"), 
                       name = "Variable", labels = c("SST (°C)", "DHW (°C-weeks)")) +
    labs(title = "Sea Surface Temperature (SST) and Degree Heating Weeks (DHW) Over Time",
         x = "Date",
         y = "Value") +
    theme_minimal()
  
  
  ################################## Attempt at unified plot ##################################
  
  # Merge the datasets
  SST_DHW_df <- full_join(SST_df, DHW_df, by = c("date", "time"))
  SST_DHW_df$date <- as.Date(SST_DHW_df$date)
  
  ggplot(SST_DHW_df, aes(x = date)) +
    geom_line(aes(y = SST, color = "SST"), size = 1) +
    geom_line(aes(y = DHW * 3, color = "DHW"), size = 1) +  # Scale DHW for visibility
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    scale_y_continuous(
      name = "SST (°C)", 
      sec.axis = sec_axis(~ . / 3, name = "DHW (°C-weeks)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("SST" = "#E69F00", "DHW" = "#0072B2"),  # Explicitly matching colors
      breaks = c("SST", "DHW"),  # Ensure correct order in legend
      labels = c("SST (°C)", "DHW (°C-weeks)")
    ) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
  
  ################################## Plot multiple years ##################################
  
  #removed tissue and DHW / SST
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = SST.90th_HS, color = "SST"), size = 1) +
    geom_line(aes(y = DHW_from_90th_HS.1 * 2, color = "DHW"), size = 1) +  # Scale DHW for visibility
    geom_line(data = obs.total.figures %>%
                 filter(Compartment == "Recovered") %>%
                 mutate(date = as.Date(date)) %>%
                 group_by(date) %>%
                 summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
               aes(x = date, y = tissue_sum), color = "red") +
    geom_point(data = obs.total.figures %>%
                 filter(Compartment == "Recovered") %>%
                 mutate(date = as.Date(date)) %>%
                 group_by(date) %>%
                 summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
               aes(x = date, y = tissue_sum), color = "black") +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    scale_y_continuous(
      name = "SST (°C)", 
      sec.axis = sec_axis(~ . / 2, name = "DHW (°C-weeks)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("SST" = "#E69F00", "DHW" = "#0072B2"),  # Explicitly matching colors
      breaks = c("SST", "DHW"),  # Ensure correct order in legend
      labels = c("SST (°C)", "DHW (°C-weeks)")
    ) +
    xlim(as.Date(c("2018-01-01", "2024-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
    
  
  #SST and infections
  
  # sst_threshold = 30.5
  sst_threshold <- SST_threshold
  sst_periods <- DHW.CRW.full %>%
    filter(SST.90th_HS >= SST_threshold) %>%
    mutate(next_date = lead(date)) %>%
    filter(!is.na(next_date))
  
  rescaler = 75
  compartment = "Infected"
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = SST.90th_HS, color = "SST"), size = 1) +
    geom_line(data = obs.total.figures %>%
                filter(Compartment == compartment) %>%
                mutate(date = as.Date(date)) %>%
                group_by(date) %>%
                summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
              aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_point(data = obs.total.figures %>%
                 filter(Compartment == compartment) %>%
                 mutate(date = as.Date(date)) %>%
                 group_by(date) %>%
                 summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
               aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    geom_hline(yintercept = SST_threshold, linetype = "dashed", color = "red", size = 1) +
    scale_y_continuous(
      name = "SST (°C)", 
      sec.axis = sec_axis(~ . / rescaler, name = "Infected tissue (m2)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("SST" = "#E69F00", "Tissue" = "#0072B2"),  # Explicitly matching colors
      breaks = c("SST", "Tissue"),  # Ensure correct order in legend
      labels = c("SST (°C)", "Infected tissue (m2)")
    ) +
    xlim(as.Date(c("2018-01-01", "2020-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
  
  #SST and removal
  rescaler = 1.5
  compartment = "Recovered"
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = SST.90th_HS, color = "SST"), size = 1) +
    geom_line(data = obs.total.figures %>%
                filter(Compartment == compartment) %>%
                mutate(date = as.Date(date)) %>%
                group_by(date) %>%
                summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
              aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_point(data = obs.total.figures %>%
                 filter(Compartment == compartment) %>%
                 mutate(date = as.Date(date)) %>%
                 group_by(date) %>%
                 summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
               aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    scale_y_continuous(
      name = "SST (°C)", 
      sec.axis = sec_axis(~ . / rescaler, name = "Removed tissue (m2)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("SST" = "#E69F00", "Tissue" = "#0072B2"),  # Explicitly matching colors
      breaks = c("SST", "Tissue"),  # Ensure correct order in legend
      labels = c("SST (°C)", "Removed tissue (m2)")
    ) +
    xlim(as.Date(c("2018-01-01", "2020-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
  
  
  #DHW and infections
  rescaler = 45
  compartment = "Infected"
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = DHW_from_90th_HS.1, color = "DHW"), size = 1) +
    geom_line(data = obs.total.figures %>%
                filter(Compartment == compartment) %>%
                mutate(date = as.Date(date)) %>%
                group_by(date) %>%
                summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
              aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_point(data = obs.total.figures %>%
                 filter(Compartment == compartment) %>%
                 mutate(date = as.Date(date)) %>%
                 group_by(date) %>%
                 summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
               aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    scale_y_continuous(
      name = "DHW (°C-weeks)", 
      sec.axis = sec_axis(~ . / rescaler, name = "Infected tissue (m2)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("DHW" = "#E69F00", "Tissue" = "#0072B2"),  # Explicitly matching colors
      breaks = c("DHW", "Tissue"),  # Ensure correct order in legend
      labels = c("DHW (°C-weeks)", "Infected tissue (m2)")
    ) +
    xlim(as.Date(c("2018-01-01", "2020-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
  
  #DHW and removal
  rescaler = 0.5
  compartment = "Recovered"
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = DHW_from_90th_HS.1, color = "DHW"), size = 1) +
    geom_line(data = obs.total.figures %>%
                filter(Compartment == compartment) %>%
                mutate(date = as.Date(date)) %>%
                group_by(date) %>%
                summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
              aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_point(data = obs.total.figures %>%
                 filter(Compartment == compartment) %>%
                 mutate(date = as.Date(date)) %>%
                 group_by(date) %>%
                 summarize(tissue_sum = sum(tissue, na.rm = TRUE), .groups = "drop"),  # Summing tissue values across sites
               aes(x = date, y = tissue_sum * rescaler, color = "Tissue")) +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    scale_y_continuous(
      name = "DHW (°C-weeks)", 
      sec.axis = sec_axis(~ . / rescaler, name = "Removed tissue (m2)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("DHW" = "#E69F00", "Tissue" = "#0072B2"),  # Explicitly matching colors
      breaks = c("DHW", "Tissue"),  # Ensure correct order in legend
      labels = c("DHW (°C-week)", "Removed tissue (m2)")
    ) +
    xlim(as.Date(c("2018-01-01", "2020-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
  
  
  ################################## Focus on Nearshore ##################################
  
  #SST and infections
  rescaler = 75
  compartment = "Infected"
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = SST.90th_HS, color = "SST"), size = 1) +
    geom_line(data = obs.total.figures %>%
                filter(Compartment == compartment, Site == "Nearshore") %>%
                mutate(date = as.Date(date)),  # Summing tissue values across sites
              aes(x = date, y = tissue * rescaler, color = "Tissue")) +
    geom_point(data = obs.total.figures %>%
                filter(Compartment == compartment, Site == "Nearshore") %>%
                mutate(date = as.Date(date)),  # Summing tissue values across sites
              aes(x = date, y = tissue * rescaler, color = "Tissue")) +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    geom_hline(yintercept = SST_threshold, linetype = "dashed", color = "red", size = 1) +
    scale_y_continuous(
      name = "SST (°C)", 
      sec.axis = sec_axis(~ . / rescaler, name = "Infected tissue (m2)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("SST" = "#E69F00", "Tissue" = "#0072B2"),  # Explicitly matching colors
      breaks = c("SST", "Tissue"),  # Ensure correct order in legend
      labels = c("SST (°C)", "Infected tissue (m2)")
    ) +
    xlim(as.Date(c("2018-01-01", "2020-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()
  
  #SST and removal
  rescaler = 1.5
  compartment = "Recovered"
  ggplot(DHW.CRW.full, aes(x = date)) +
    geom_line(aes(y = SST.90th_HS, color = "SST"), size = 1) +
    geom_line(data = obs.total.figures %>%
                filter(Compartment == compartment, Site == "Nearshore") %>%
                mutate(date = as.Date(date)),  # Summing tissue values across sites
              aes(x = date, y = tissue * rescaler, color = "Tissue")) +
    geom_point(data = obs.total.figures %>%
                 filter(Compartment == compartment, Site == "Nearshore") %>%
                 mutate(date = as.Date(date)),  # Summing tissue values across sites
               aes(x = date, y = tissue * rescaler, color = "Tissue")) +
    geom_vline(xintercept = as.Date("2019-11-16"), color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.Date("2019-12-06"), color = "black", linetype = "dashed") +
    geom_hline(yintercept = SST_threshold, linetype = "dashed", color = "red", size = 1) +
    scale_y_continuous(
      name = "SST (°C)", 
      sec.axis = sec_axis(~ . / rescaler, name = "Removed tissue (m2)")  # Reverse the scaling
    ) +
    scale_color_manual(
      name = "Variable",
      values = c("SST" = "#E69F00", "Tissue" = "#0072B2"),  # Explicitly matching colors
      breaks = c("SST", "Tissue"),  # Ensure correct order in legend
      labels = c("SST (°C)", "Removed tissue (m2)")
    ) +
    xlim(as.Date(c("2018-01-01", "2020-12-30"))) +
    labs(x = "Date", color = "Variable") +
    theme_minimal()

  ################################## Plot in base R to avoid scaling issues ##################################
  
  #infected tissue
  # Sample data preparation
  rescaler = 1.5
  compartment = "Infected"
  
  # Filter and prepare the data
  SST_data <- DHW.CRW.full
  tissue_data <- obs.total.figures %>%
    filter(Compartment == compartment, Site == "Nearshore") %>%
    mutate(date = as.Date(date), tissue_scaled = tissue * rescaler)
  
  # Plot the SST data
  plot(SST_data$date, SST_data$SST.90th_HS, type = "l", col = "#E69F00", 
       xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
       ylim = c(23, 32), xlab = "Date", ylab = "SST (°C)", 
       main = "SST and Infected Tissue Over Time", lwd = 2)
  
  # Add vertical dashed lines
  abline(v = as.Date(c("2019-11-16", "2019-12-06")), col = "black", lty = 2)
  
  # Add horizontal dashed line for SST threshold
  abline(h = SST_threshold, col = "red", lty = 2, lwd = 1)
  
  # Overlay tissue data
  par(new = TRUE)
  
  # Plot the tissue data with a secondary axis
  plot(tissue_data$date, tissue_data$tissue_scaled, type = "b", col = "#0072B2", 
       xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
       ylim = c(0, max(tissue_data$tissue_scaled)), 
       xlab = "", ylab = "", axes = FALSE)  # No y-axis, to avoid double labeling
  
  # Add the secondary y-axis (right axis)
  axis(side = 4, at = pretty(tissue_data$tissue_scaled), 
       labels = pretty(tissue_data$tissue_scaled) / rescaler)  # Adjust to match tissue scaling
  
  # Add labels to the secondary axis
  mtext("Removed tissue (m2)", side = 4, line = 3)
  
  # Add a legend
  legend("topright", legend = c("SST (°C)", "Infected tissue (m2)"), 
         col = c("#E69F00", "#0072B2"), lwd = 2, pch = c(NA, 16), bty = "n")
  
  #removed tissue
  # Sample data preparation
  rescaler = 1.5
  compartment = "Recovered"
  
  # Filter and prepare the data
  SST_data <- DHW.CRW.full
  tissue_data <- obs.total.figures %>%
    filter(Compartment == compartment, Site == "Nearshore") %>%
    mutate(date = as.Date(date), tissue_scaled = tissue * rescaler)
  
  # Plot the SST data
  plot(SST_data$date, SST_data$SST.90th_HS, type = "l", col = "#E69F00", 
       xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
       ylim = c(23, 32), xlab = "Date", ylab = "SST (°C)", 
       main = "SST and Removed Tissue Over Time", lwd = 2)
  
  # Add vertical dashed lines
  abline(v = as.Date(c("2019-11-16", "2019-12-06")), col = "black", lty = 2)
  
  # Add horizontal dashed line for SST threshold
  abline(h = SST_threshold, col = "red", lty = 2, lwd = 1)
  
  # Overlay tissue data
  par(new = TRUE)
  
  # Plot the tissue data with a secondary axis
  plot(tissue_data$date, tissue_data$tissue_scaled, type = "b", col = "#0072B2", 
       xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
       ylim = c(0, max(tissue_data$tissue_scaled)), 
       xlab = "", ylab = "", axes = FALSE)  # No y-axis, to avoid double labeling
  
  # Add the secondary y-axis (right axis)
  axis(side = 4, at = pretty(tissue_data$tissue_scaled), 
       labels = pretty(tissue_data$tissue_scaled) / rescaler)  # Adjust to match tissue scaling
  
  # Add labels to the secondary axis
  mtext("Removed tissue (m2)", side = 4, line = 3)
  
  # Add a legend
  legend("topright", legend = c("SST (°C)", "Removed tissue (m2)"), 
         col = c("#E69F00", "#0072B2"), lwd = 2, pch = c(NA, 16), bty = "n")
  
  
  ################################## More mature sigmoid & logistic behaviors of temperature ##################################
  
  # Load necessary libraries
  library(ggplot2)
  library(viridis)
  library(patchwork)  # For arranging multiple plots
  
  # Define temperature range from 23 to 33
  T <- seq(23, 33, length.out = 100)
  
  # ---- APPROACH 1: Inflection point at 30.5C (original) ----
  temp_effect_inflection <- function(temp, T_inflection = 30.5, zeta_bounded = 0.5) {
    # Convert bounded zeta to effective zeta
    zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    
    # Logistic function with inflection point at T_inflection
    modifier <- 1 / (1 + exp(-zeta_effective * (temp - T_inflection)))
    
    return(modifier)
  }
  
  # ---- APPROACH 2: Original model clarified ----
  temp_effect_original <- function(temp, T_min = 23, zeta_bounded = 0.5) {
    
    # Convert bounded zeta to effective zeta
    
    zeta_effective <- 5 * zeta_bounded/(1-0.7*zeta_bounded)
    # zeta_effective <- zeta_bounded/(1 - zeta_bounded)
    
    # Original formula clarified
    # This is a non-centered sigmoid that:
    # - Starts near 0 at T_min
    # - Increases asymptotically toward 1 as temp increases
    # - Has no inflection point assumption
    modifier <- (1 - exp(-zeta_effective * (temp - T_min))) / (1 + exp(-zeta_effective * (temp - T_min)))
    
    return(modifier)
  }
  
  # ---- APPROACH 3: Decreasing function with temperature ----
  temp_effect_decreasing <- function(temp, T_max = 33, T_min = 23, zeta_bounded = 0.5) {
    # Convert bounded zeta to effective zeta
    zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    
    # Calculate midpoint of temperature range for centering the sigmoid
    T_mid <- (T_min + T_max) / 2
    
    # Decreasing sigmoid that starts near 1 (at low temps) and approaches 0 as temp increases
    # Note the positive sign before zeta_effective which makes the function decrease with temperature
    modifier <- 1 / (1 + exp(zeta_effective * (temp - T_mid)))
    
    return(modifier)
  }
  
  # ---- APPROACH 4: Non-centered sigmoid decreasing with temperature (FIXED) ----
  temp_effect_noncentered_decreasing <- function(temp, T_min = 23, T_max = 33, zeta_bounded = 0.5) {
    # Convert bounded zeta to effective zeta
    zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    
    # Modified non-centered sigmoid that:
    # - Equals 1 at T_min (23°C)
    # - Decreases toward 0 as temperature increases
    # - Maintains the non-centered sigmoid shape but inverted
    
    # Calculate the relative position in the temperature range
    rel_temp <- (temp - T_min) / (T_max - T_min)
    
    # Apply the decreasing non-centered sigmoid
    modifier <- 1 - ((1 - exp(-zeta_effective * rel_temp)) / (1 + exp(-zeta_effective * rel_temp)))
    
    return(modifier)
  }
  
  # Define bounded zeta values from 0 to 1 (avoiding 1 which would create infinity)
  zeta_values <- c(0, 0.01, 0.05, 0.1, 0.3, 0.5, 1)
  
  # Create data for all approaches
  plot_data_inflection <- data.frame()
  plot_data_original <- data.frame()
  plot_data_decreasing <- data.frame()
  plot_data_noncentered_decreasing <- data.frame()
  
  for (zeta in zeta_values) {
    # Data for approach 1
    modifiers_inflection <- sapply(T, function(t) temp_effect_inflection(t, zeta_bounded = zeta))
    temp_data_inflection <- data.frame(
      Temperature = T,
      Modifier = modifiers_inflection,
      Zeta = as.factor(zeta)
    )
    plot_data_inflection <- rbind(plot_data_inflection, temp_data_inflection)
    
    # Data for approach 2
    modifiers_original <- sapply(T, function(t) temp_effect_original(t, zeta_bounded = zeta))
    temp_data_original <- data.frame(
      Temperature = T,
      Modifier = modifiers_original,
      Zeta = as.factor(zeta)
    )
    plot_data_original <- rbind(plot_data_original, temp_data_original)
    
    # Data for approach 3
    modifiers_decreasing <- sapply(T, function(t) temp_effect_decreasing(t, zeta_bounded = zeta))
    temp_data_decreasing <- data.frame(
      Temperature = T,
      Modifier = modifiers_decreasing,
      Zeta = as.factor(zeta)
    )
    plot_data_decreasing <- rbind(plot_data_decreasing, temp_data_decreasing)
    
    # Data for approach a4
    modifiers_noncentered_decreasing <- sapply(T, function(t) temp_effect_noncentered_decreasing(t, zeta_bounded = zeta))
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
    scale_x_continuous(breaks = seq(23, 33, 2)) +
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
    scale_x_continuous(breaks = seq(23, 33, 2)) +
    geom_vline(xintercept = 23, linetype = "dashed", alpha = 0.5)
  
  # Plot 3: Classic sigmoid decreasing model
  p3 <- ggplot(plot_data_decreasing, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(size = 1) +
    labs(
      title = "Temperature Effect with Decreasing Sigmoid Function",
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
    scale_x_continuous(breaks = seq(23, 33, 2)) +
    geom_vline(xintercept = 28, linetype = "dashed", alpha = 0.5)
  
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
    scale_x_continuous(breaks = seq(23, 33, 2)) +
    geom_vline(xintercept = 23, linetype = "dashed", alpha = 0.5)
  
  # Combine plots - using a 2x2 grid for better viewing
  (p1 + p2) / (p3 + p4)
  
  # # Let's focus only on the fourth plot
  # p4
  
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
  