  
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

  # Define parameters
  SST_threshold <- 30.5  # Temperature threshold
  z_values <- c(0.5, 1, 2.5, 3)  # Different z values for transmission effect
  e_values <- c(0.0001, 1, 2, 3)  # Different e values for mortality effect
  SST_range <- seq(26, 32, by = 0.1)  # SST from 29 to 32 degrees

  # Set up plot
  plot(NULL, xlim = range(SST_range), ylim = c(0, 1),
       xlab = "Sea Surface Temperature (°C)",
       ylab = "Scaling Factor",
       main = "Effect of SST on Transmission and Mortality Rates")

  # Define colors and line types
  colors <- c("blue", "red", "green", "purple")
  linetypes <- c(1, 2)  # Solid for transmission, dashed for mortality

  # Plot transmission rate curves
  for (i in seq_along(z_values)) {
    z <- z_values[i]
    scaling_factor_z <- 1 / (1 + exp(-z * (SST_range - SST_threshold)))
    lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2, lty = linetypes[1])
  }

  # Plot mortality rate curves
  for (i in seq_along(e_values)) {
    e <- e_values[i]
    scaling_factor_e <- 1 / (1 + exp(-e * (SST_range - SST_threshold)))
    lines(SST_range, scaling_factor_e, col = colors[i], lwd = 2, lty = linetypes[2])
  }

  # Add legend
  legend("topleft", legend = c(paste("Transmission (z) =", z_values),
                               paste("Mortality (e) =", e_values)),
         col = rep(colors, 2), lwd = 2, lty = rep(linetypes, each = length(z_values)))

  ################################## SST above threshold ##################################
  
  # Define parameters
  SST_threshold <- 30.5  # Temperature threshold
  z_values <- c(0.63, 2.98, 3)  # Different z values for transmission effect
  e_values <- c(0.0001, 1.89, 2.57)  # Different e values for mortality effect
  SST_range <- seq(30.5, 32, by = 0.1)  # SST from 30.5°C onward

  # Set up plot
  plot(NULL, xlim = range(SST_range), ylim = c(0, 1),
       xlab = "Sea Surface Temperature (°C)",
       ylab = "Scaling Factor",
       main = "Effect of SST on Transmission and Mortality Rates")

  # Define colors and line types
  colors <- c("blue", "red", "green")
  linetypes <- c(1, 2)  # Solid for transmission, dashed for mortality

  # Plot transmission rate curves (only for SST > 30.5)
  for (i in seq_along(z_values)) {
    z <- z_values[i]
    scaling_factor_z <- 1 / (1 + exp(-z * (SST_range - SST_threshold)))
    lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2, lty = linetypes[1])
  }

  # Plot mortality rate curves (only for SST > 30.5)
  for (i in seq_along(e_values)) {
    e <- e_values[i]
    scaling_factor_e <- 1 / (1 + exp(-e * (SST_range - SST_threshold)))
    lines(SST_range, scaling_factor_e, col = colors[i], lwd = 2, lty = linetypes[2])
  }

  # Add legend
  legend("topleft", legend = c(paste("Transmission (z) =", z_values),
                               paste("Mortality (e) =", e_values)),
         col = rep(colors, 2), lwd = 2, lty = rep(linetypes, each = length(z_values)))
  
  ################################ Wide range ##################################
  
  # Define parameters
  SST_threshold <- 30.5  # Temperature threshold
  z_values <- seq(0.001, 10, length.out = 8)  # 8 different z values in range 1 to 3
  e_values <- seq(0.001, 10, length.out = 8)  # 8 different e values in range 1 to 3
  SST_range <- seq(30.5, 34, by = 0.01)  # SST from 30.5°C onward

  # Set up plot
  plot(NULL, xlim = range(SST_range), ylim = c(0, 1),
       xlab = "Sea Surface Temperature (°C)",
       ylab = "Scaling Factor",
       main = "Effect of SST on Transmission and Mortality Rates")

  # Define colors and line types
  colors <- colorRampPalette(c("blue", "red"))(length(z_values))  # Gradient color scheme
  linetypes <- c(1, 2)  # Solid for transmission, dashed for mortality

  # Plot transmission rate curves (only for SST > 30.5)
  for (i in seq_along(z_values)) {
    z <- z_values[i]
    scaling_factor_z <- 1 / (1 + exp(-z * (SST_range - SST_threshold)))
    lines(SST_range, scaling_factor_z, col = colors[i], lwd = 2, lty = linetypes[1])
  }

  # Plot mortality rate curves (only for SST > 30.5)
  for (i in seq_along(e_values)) {
    e <- e_values[i]
    scaling_factor_e <- 1 / (1 + exp(-e * (SST_range - SST_threshold)))
    lines(SST_range, scaling_factor_e, col = colors[i], lwd = 2, lty = linetypes[2])
  }

  # Add legend
  legend("bottomright", legend = paste("Rates (z, e) =", round(z_values, 2)),
         col = colors, lwd = 2, lty = 1)

  ################################## Simplified function, 0 to 1 ##################################

  # Define parameters
  SST_threshold <- 30.5  # Threshold temperature
  zeta_values <- seq(0, 1, by = 0.25)  # Varying zeta from 0 to 1
  SST_range <- seq(30.5, 32, by = 0.1)  # SST only above threshold

  # Set up plot
  plot(NULL, xlim = range(SST_range), ylim = c(0, 1),
       xlab = "Sea Surface Temperature (°C)",
       ylab = "Scaling Factor",
       main = "Effect of SST on Transmission & Mortality Rates")

  # Colors for different zeta values
  colors <- rainbow(length(zeta_values))

  # Loop through different zeta values
  for (i in seq_along(zeta_values)) {
    zeta <- zeta_values[i]
    scaling_factor <- 1 - zeta * (1 / (1 + exp(-10 * (SST_range - SST_threshold))))
    lines(SST_range, scaling_factor, col = colors[i], lwd = 2)
  }

  # Add legend
  legend("topright", legend = paste("zeta =", zeta_values),
         col = colors, lwd = 2)

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

  ################################## Plot Nearshore in base R to avoid scaling issues ##################################
  
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
  
  
  ################################## Plot Midchannel in base R to avoid scaling issues ##################################
  
  #infected tissue
  # Sample data preparation
  rescaler = 1.5
  compartment = "Infected"
  
  # Filter and prepare the data
  SST_data <- DHW.CRW.full
  tissue_data <- obs.total.figures %>%
    filter(Compartment == compartment, Site == "Midchannel") %>%
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
    filter(Compartment == compartment, Site == "Midchannel") %>%
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
  
  
  
  
  ################################## Plot Offshore in base R to avoid scaling issues ##################################
  
  #infected tissue
  # Sample data preparation
  rescaler = 1.5
  compartment = "Infected"
  
  # Filter and prepare the data
  SST_data <- DHW.CRW.full
  tissue_data <- obs.total.figures %>%
    filter(Compartment == compartment, Site == "Offshore") %>%
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
    filter(Compartment == compartment, Site == "Offshore") %>%
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
  
  
  ################################## Infection extension 2020 ##################################
  
  # Prepare the base data for synthetic extension - Nearshore
  infected_base_data_nearshore <- obs.total.figures %>%
    filter(Compartment == "Infected", Site == "Nearshore") %>%
    mutate(date = as.Date(date))
  
  # Prepare the base data for synthetic extension - Offshore
  infected_base_data_offshore <- obs.total.figures %>%
    filter(Compartment == "Infected", Site == "Offshore") %>%
    mutate(date = as.Date(date))
  
  # Prepare the base data for synthetic extension - Midchannel
  infected_base_data_midchannel <- obs.total.figures %>%
    filter(Compartment == "Infected", Site == "Midchannel") %>%
    mutate(date = as.Date(date))
  
  # Synthetic data for Nearshore
  infected_synthetic_data_nearshore <- tibble(
    Compartment = "Infected",
    Site = "Nearshore",
    date = as.Date(c(
      "2019-12-06",
      "2020-01-15", 
      "2020-02-15", 
      "2020-03-15", 
      "2020-04-15", 
      "2020-05-15", 
      "2020-06-15", 
      "2020-07-15", 
      "2020-08-15", 
      "2020-09-15", 
      "2020-10-15", 
      "2020-11-15", 
      "2020-12-15"
    )),
    tissue = c(
      # Start rising immediately from the 2019-11-12 trajectory
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 1.2,   # January - continuing rise
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 1.4,   # January - continuing rise
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 1.8,   # February - continued increase
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 1.3,   # March - further rise
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 1.1,   # April - continued growth
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0.7,   # May - peak
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0.03,   # June - slight decline
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0.02,   # July - continued decline
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0.01,   # August - further reduction
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0,   # September - lower levels
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0,   # October - stabilizing
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0.2,   # November - slight decrease
      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 0.6    # December - further reduction
    )
  )
  
  # Synthetic data for Offshore
  infected_synthetic_data_offshore <- tibble(
    Compartment = "Infected",
    Site = "Offshore",
    date = as.Date(c(
      "2019-12-06",
      "2020-01-15", 
      "2020-02-15", 
      "2020-03-15", 
      "2020-04-15", 
      "2020-05-15", 
      "2020-06-15", 
      "2020-07-15", 
      "2020-08-15", 
      "2020-09-15", 
      "2020-10-15", 
      "2020-11-15", 
      "2020-12-15"
    )),
    tissue = c(
      # Start rising immediately from the 2019-11-12 trajectory
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 1.2,   # January - continuing rise
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 1.4,   # January - continuing rise
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 1.8,   # February - continued increase
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 1.1,   # March - further rise
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.5,   # April - continued growth
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.2,   # May - peak
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.03,   # June - slight decline
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.02,   # July - continued decline
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.01,   # August - further reduction
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0,   # September - lower levels
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0,   # October - stabilizing
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.2,   # November - slight decrease
      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 0.6    # December - further reduction
    )
  )
  
  # Synthetic data for Midchannel
  infected_synthetic_data_midchannel <- tibble(
    Compartment = "Infected",
    Site = "Midchannel",
    date = as.Date(c(
      "2019-12-06",
      "2020-01-15", 
      "2020-02-15", 
      "2020-03-15", 
      "2020-04-15", 
      "2020-05-15", 
      "2020-06-15", 
      "2020-07-15", 
      "2020-08-15", 
      "2020-09-15", 
      "2020-10-15", 
      "2020-11-15", 
      "2020-12-15"
    )),
    tissue = c(
      # If no base data exists, use a trend similar to other sites or a constant value
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 1,   # January - moderate rise
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.7,   # January - continued rise
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.9,   # February - continued increase
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.6,   # March - further rise
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.5,   # April - moderate growth
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.3,   # May - peak
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.01,   # June - slight decline
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.005,   # July - continued decline
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.001,   # August - further reduction
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0,   # September - lower levels
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0,   # October - stabilizing
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.1,   # November - slight decrease
      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 0.3    # December - further reduction
    )
  )
  
  # Combine original and synthetic data for Nearshore, Offshore, and Midchannel
  infected_combined_data_nearshore <- bind_rows(
    # Modify the 2019-12-06 row to be synthetic
    mutate(infected_base_data_nearshore, synthetic = ifelse(date == "2019-12-06", TRUE, FALSE)),
    mutate(infected_synthetic_data_nearshore, synthetic = TRUE)
  )
  
  infected_combined_data_offshore <- bind_rows(
    # Modify the 2019-12-06 row to be synthetic
    mutate(infected_base_data_offshore, synthetic = ifelse(date == "2019-12-06", TRUE, FALSE)),
    mutate(infected_synthetic_data_offshore, synthetic = TRUE)
  )
  
  infected_combined_data_midchannel <- bind_rows(
    # Modify the 2019-12-06 row to be synthetic
    mutate(infected_base_data_midchannel, synthetic = ifelse(date == "2019-12-06", TRUE, FALSE)),
    mutate(infected_synthetic_data_midchannel, synthetic = TRUE)
  )
  
  # Remove rows where date is 2019-12-06 and days.survey is NA
  infected_combined_data_nearshore <- infected_combined_data_nearshore %>%
    filter(!(date == as.Date("2019-12-06") & is.na(days.survey)))
  
  infected_combined_data_offshore <- infected_combined_data_offshore %>%
    filter(!(date == as.Date("2019-12-06") & is.na(days.survey)))
  
  infected_combined_data_midchannel <- infected_combined_data_midchannel %>%
    filter(!(date == as.Date("2019-12-06") & is.na(days.survey)))
  
  # Identify the last known values from 2019-12-06 for each site
  last_known_nearshore <- infected_combined_data_nearshore %>%
    filter(date == as.Date("2019-12-06")) %>%
    select(days.survey, days.inf.site, tissue)
  
  last_known_offshore <- infected_combined_data_offshore %>%
    filter(date == as.Date("2019-12-06")) %>%
    select(days.survey, days.inf.site, tissue)
  
  last_known_midchannel <- infected_combined_data_midchannel %>%
    filter(date == as.Date("2019-12-06")) %>%
    select(days.survey, days.inf.site, tissue)
  
  # Update the existing 2019-12-06 rows
  infected_combined_data_nearshore <- infected_combined_data_nearshore %>%
    mutate(
      tissue = ifelse(date == as.Date("2019-12-06"), 
                      infected_base_data_nearshore$tissue[which(infected_base_data_nearshore$date == "2019-11-12")] * 1.2, 
                      tissue)
    )
  
  infected_combined_data_offshore <- infected_combined_data_offshore %>%
    mutate(
      tissue = ifelse(date == as.Date("2019-12-06"), 
                      infected_base_data_offshore$tissue[which(infected_base_data_offshore$date == "2019-11-12")] * 1.2, 
                      tissue)
    )
  
  infected_combined_data_midchannel <- infected_combined_data_midchannel %>%
    mutate(
      tissue = ifelse(date == as.Date("2019-12-06"), 
                      infected_base_data_midchannel$tissue[which(infected_base_data_midchannel$date == "2019-11-12")] * 1.2, 
                      tissue)
    )
  
  # Define synthetic rows and compute new days.survey and days.inf.site for each site
  infected_combined_data_nearshore <- infected_combined_data_nearshore %>%
    arrange(date) %>%
    mutate(
      days.survey = ifelse(synthetic, 
                           last_known_nearshore$days.survey + as.numeric(date - as.Date("2019-12-06")), 
                           days.survey),
      days.inf.site = ifelse(synthetic, 
                             last_known_nearshore$days.inf.site + as.numeric(date - as.Date("2019-12-06")), 
                             days.inf.site)
    )
  
  infected_combined_data_offshore <- infected_combined_data_offshore %>%
    arrange(date) %>%
    mutate(
      days.survey = ifelse(synthetic, 
                           last_known_offshore$days.survey + as.numeric(date - as.Date("2019-12-06")), 
                           days.survey),
      days.inf.site = ifelse(synthetic, 
                             last_known_offshore$days.inf.site + as.numeric(date - as.Date("2019-12-06")), 
                             days.inf.site)
    )
  
  infected_combined_data_midchannel <- infected_combined_data_midchannel %>%
    arrange(date) %>%
    mutate(
      days.survey = ifelse(synthetic, 
                           last_known_midchannel$days.survey + as.numeric(date - as.Date("2019-12-06")), 
                           days.survey),
      days.inf.site = ifelse(synthetic, 
                             last_known_midchannel$days.inf.site + as.numeric(date - as.Date("2019-12-06")), 
                             days.inf.site)
    )
  
  # Create the new obs.total.figures_SST by removing existing site Infected data 
  # and adding the new combined data
  obs.total.figures_SST.infected <- obs.total.figures %>%
    filter(!(Compartment == "Infected" & Site %in% c("Nearshore", "Offshore", "Midchannel"))) %>%
    bind_rows(infected_combined_data_nearshore) %>%
    bind_rows(infected_combined_data_offshore) %>%
    bind_rows(infected_combined_data_midchannel)
  
  # Plotting section
  # Sample data preparation
  rescaler = 1.5
  compartment = "Infected"
  
  #create moving average of SST data to smooth it and reduce wiggliness of the data
  # Simple moving average smoothing
  DHW.CRW.full <- DHW.CRW.full %>%
    mutate(SST_90th_HS_smoothed = smooth(SST.90th_HS, kind = "3R"))
  
  # Prepare SST data
  SST_data <- DHW.CRW.full
  
  # Set up a 3x1 plot layout
  par(mfrow=c(3,1), mar=c(4, 4, 2, 4))
  
  # List of sites to plot
  sites <- c("Nearshore", "Offshore", "Midchannel")
  
  # Loop through each site to create plots
  for (site in sites) {
    # Filter tissue data for the current site
    tissue_data <- obs.total.figures_SST.infected %>%
      filter(Compartment == compartment, Site == site) %>%
      mutate(date = as.Date(date), tissue_scaled = tissue * rescaler)
    
    # Plot the SST data
    # plot(SST_data$date, SST_data$SST.90th_HS, type = "l", col = "#E69F00", 
    plot(SST_data$date, SST_data$SST_90th_HS_smoothed, type = "l", col = "#E69F00",
         xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
         ylim = c(23, 32), xlab = "Date", ylab = "SST (°C)", 
         main = paste("SST and Infected Tissue Over Time (", site, ")", sep = ""), lwd = 2)
    
    # Add vertical dashed lines
    abline(v = as.Date(c("2019-11-16", "2019-12-06")), col = "black", lty = 2)
    
    # Add horizontal dashed line for SST threshold
    abline(h = SST_threshold, col = "red", lty = 2, lwd = 1)
    
    # Overlay tissue data
    par(new = TRUE)
    
    # Plot the tissue data with a secondary axis
    plot(tissue_data$date, tissue_data$tissue_scaled, type = "b", 
         col = ifelse(tissue_data$synthetic, "red", "#0072B2"), 
         xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
         ylim = c(0, max(tissue_data$tissue_scaled)), 
         xlab = "", ylab = "", axes = FALSE)
    
    # Add the secondary y-axis (right axis)
    axis(side = 4, at = pretty(tissue_data$tissue_scaled), 
         labels = pretty(tissue_data$tissue_scaled) / rescaler)
    
    # Add labels to the secondary axis
    mtext("Removed tissue (m2)", side = 4, line = 3)
    
    # Add a legend
    legend("topright", legend = c("SST (°C)", "Infected tissue (m2)"), 
           col = c("#E69F00", "#0072B2"), lwd = 2, pch = c(NA, 16), bty = "n")
  }
  
  
  ################################## Removal extension 2020 ##################################
  
  # Prepare the base data for synthetic extension - Nearshore
  removed_base_data_nearshore <- obs.total.figures %>%
    filter(Compartment == "Recovered", Site == "Nearshore") %>%
    mutate(date = as.Date(date))
  
  # Prepare the base data for synthetic extension - Offshore
  removed_base_data_offshore <- obs.total.figures %>%
    filter(Compartment == "Recovered", Site == "Offshore") %>%
    mutate(date = as.Date(date))
  
  # Prepare the base data for synthetic extension - Midchannel
  removed_base_data_midchannel <- obs.total.figures %>%
    filter(Compartment == "Recovered", Site == "Midchannel") %>%
    mutate(date = as.Date(date))
  
  #scalars for simulating removed tissue in 2020
  increase_factors_nearshore <- c(
    1.04,   # January 2020
    1.10,  # February 2020
    1.13,   # March 2020
    1.14,  # April 2020
    1.15,   # May 2020 (plateau starts)
    1.15,   # June 2020 (plateau)
    1.15,   # July 2020 (plateau)
    1.15,   # August 2020 (plateau)
    1.155,  # September 2020 (resumes rising)
    1.16,   # October 2020
    1.17,  # November 2020
    1.20    # December 2020
  )
  increase_factors_offshore <- c(
    1.2,   # January 2020
    1.3,  # February 2020
    1.35,   # March 2020
    1.375,  # April 2020
    1.40,   # May 2020 (plateau starts)
    1.40,   # June 2020 (plateau)
    1.40,   # July 2020 (plateau)
    1.40,   # August 2020 (plateau)
    1.41,  # September 2020 (resumes rising)
    1.415,   # October 2020
    1.45,  # November 2020
    1.5    # December 2020
  )
  increase_factors_midchannel <- c(
    1.3,   # January 2020
    1.4,  # February 2020
    1.45,   # March 2020
    1.475,  # April 2020
    1.50,   # May 2020 (plateau starts)
    1.50,   # June 2020 (plateau)
    1.50,   # July 2020 (plateau)
    1.50,   # August 2020 (plateau)
    1.51,  # September 2020 (resumes rising)
    1.515,   # October 2020
    1.55,  # November 2020
    1.6    # December 2020
  )
  
  # Synthetic data for Nearshore
  removed_synthetic_data_nearshore <- tibble(
    Compartment = "Recovered",
    Site = "Nearshore",
    date = as.Date(c(
      "2020-01-15", 
      "2020-02-15", 
      "2020-03-15", 
      "2020-04-15", 
      "2020-05-15", 
      "2020-06-15", 
      "2020-07-15", 
      "2020-08-15", 
      "2020-09-15", 
      "2020-10-15", 
      "2020-11-15", 
      "2020-12-15"
    )),
    tissue = removed_base_data_nearshore$tissue[which(removed_base_data_nearshore$date == "2019-12-06")] * increase_factors_nearshore
  )
  
  # Synthetic data for Offshore
  removed_synthetic_data_offshore <- tibble(
    Compartment = "Recovered",
    Site = "Offshore",
    date = as.Date(c(
      "2020-01-15", 
      "2020-02-15", 
      "2020-03-15", 
      "2020-04-15", 
      "2020-05-15", 
      "2020-06-15", 
      "2020-07-15", 
      "2020-08-15", 
      "2020-09-15", 
      "2020-10-15", 
      "2020-11-15", 
      "2020-12-15"
    )),
    tissue = removed_base_data_offshore$tissue[which(removed_base_data_offshore$date == "2019-12-06")] * increase_factors_offshore
  )  
  
  # Synthetic data for Midchannel
  removed_synthetic_data_midchannel <- tibble(
    Compartment = "Recovered",
    Site = "Midchannel",
    date = as.Date(c(
      "2020-01-15", 
      "2020-02-15", 
      "2020-03-15", 
      "2020-04-15", 
      "2020-05-15", 
      "2020-06-15", 
      "2020-07-15", 
      "2020-08-15", 
      "2020-09-15", 
      "2020-10-15", 
      "2020-11-15", 
      "2020-12-15"
    )),
    tissue = removed_base_data_midchannel$tissue[which(removed_base_data_midchannel$date == "2019-12-06")] * increase_factors_midchannel
  )  
  
  # Combine original and synthetic data for Nearshore, Offshore, and Midchannel
  removed_combined_data_nearshore <- bind_rows(
    mutate(removed_base_data_nearshore, synthetic = FALSE),
    mutate(removed_synthetic_data_nearshore, synthetic = TRUE)
  )
  
  removed_combined_data_offshore <- bind_rows(
    mutate(removed_base_data_offshore, synthetic = FALSE),
    mutate(removed_synthetic_data_offshore, synthetic = TRUE)
  )
  
  removed_combined_data_midchannel <- bind_rows(
    mutate(removed_base_data_midchannel, synthetic = FALSE),
    mutate(removed_synthetic_data_midchannel, synthetic = TRUE)
  )
  
  # Identify the last known values from 2019-12-06 for each site
  last_known_nearshore <- removed_combined_data_nearshore %>%
    filter(date == as.Date("2019-12-06")) %>%
    select(days.survey, days.inf.site, tissue)
  
  last_known_offshore <- removed_combined_data_offshore %>%
    filter(date == as.Date("2019-12-06")) %>%
    select(days.survey, days.inf.site, tissue)
  
  last_known_midchannel <- removed_combined_data_midchannel %>%
    filter(date == as.Date("2019-12-06")) %>%
    select(days.survey, days.inf.site, tissue)
  
  # Define synthetic rows and compute new days.survey and days.inf.site for each site
  removed_combined_data_nearshore <- removed_combined_data_nearshore %>%
    arrange(date) %>%
    mutate(
      days.survey = ifelse(synthetic, 
                           last_known_nearshore$days.survey + as.numeric(date - as.Date("2019-12-06")), 
                           days.survey),
      days.inf.site = ifelse(synthetic, 
                             last_known_nearshore$days.inf.site + as.numeric(date - as.Date("2019-12-06")), 
                             days.inf.site)
    )
  
  removed_combined_data_offshore <- removed_combined_data_offshore %>%
    arrange(date) %>%
    mutate(
      days.survey = ifelse(synthetic, 
                           last_known_offshore$days.survey + as.numeric(date - as.Date("2019-12-06")), 
                           days.survey),
      days.inf.site = ifelse(synthetic, 
                             last_known_offshore$days.inf.site + as.numeric(date - as.Date("2019-12-06")), 
                             days.inf.site)
    )
  
  removed_combined_data_midchannel <- removed_combined_data_midchannel %>%
    arrange(date) %>%
    mutate(
      days.survey = ifelse(synthetic, 
                           last_known_midchannel$days.survey + as.numeric(date - as.Date("2019-12-06")), 
                           days.survey),
      days.inf.site = ifelse(synthetic, 
                             last_known_midchannel$days.inf.site + as.numeric(date - as.Date("2019-12-06")), 
                             days.inf.site)
    )
  
  # Create the new obs.total.figures_SST by removing existing site Recovered data 
  # and adding the new combined data
  obs.total.figures_SST.removal <- obs.total.figures %>%
    filter(!(Compartment == "Recovered" & Site %in% c("Nearshore", "Offshore", "Midchannel"))) %>%
    bind_rows(removed_combined_data_nearshore) %>%
    bind_rows(removed_combined_data_offshore) %>%
    bind_rows(removed_combined_data_midchannel)
  
  obs.total.figures_SST <- obs.total.figures %>%
    filter(!(Compartment == "Recovered" | Compartment == "Infected" & Site %in% c("Nearshore", "Offshore", "Midchannel"))) %>%
    bind_rows(removed_combined_data_nearshore) %>%
    bind_rows(removed_combined_data_offshore) %>%
    bind_rows(removed_combined_data_midchannel) %>%
    bind_rows(infected_combined_data_nearshore) %>%
    bind_rows(infected_combined_data_offshore) %>%
    bind_rows(infected_combined_data_midchannel) 
  
  # Plotting section
  # Sample data preparation
  rescaler = 1.5
  compartment = "Recovered"
  
  #create moving average of SST data to smooth it and reduce wiggliness of the data
  # Simple moving average smoothing
  DHW.CRW.full <- DHW.CRW.full %>%
    mutate(SST_90th_HS_smoothed = smooth(SST.90th_HS, kind = "3R"))
  
  # Prepare SST data
  SST_data <- DHW.CRW.full
  
  # Set up a 3x1 plot layout
  par(mfrow=c(3,1), mar=c(4, 4, 2, 4))
  
  # List of sites to plot
  sites <- c("Nearshore", "Offshore", "Midchannel")
  
  # Loop through each site to create plots
  for (site in sites) {
    # Filter tissue data for the current site
    tissue_data <- obs.total.figures_SST.removal %>%
      filter(Compartment == compartment, Site == site) %>%
      mutate(date = as.Date(date), tissue_scaled = tissue * rescaler)
    
    # Plot the SST data
    # plot(SST_data$date, SST_data$SST.90th_HS, type = "l", col = "#E69F00", 
    plot(SST_data$date, SST_data$SST_90th_HS_smoothed, type = "l", col = "#E69F00",
         xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
         ylim = c(23, 32), xlab = "Date", ylab = "SST (°C)", 
         main = paste("SST and Recovered Tissue Over Time (", site, ")", sep = ""), lwd = 2)
    
    # Add vertical dashed lines
    abline(v = as.Date(c("2019-11-16", "2019-12-06")), col = "black", lty = 2)
    
    # Add horizontal dashed line for SST threshold
    abline(h = SST_threshold, col = "red", lty = 2, lwd = 1)
    
    # Overlay tissue data
    par(new = TRUE)
    
    # Plot the tissue data with a secondary axis
    plot(tissue_data$date, tissue_data$tissue_scaled, type = "b", 
         col = ifelse(tissue_data$synthetic, "red", "#0072B2"), 
         xlim = c(as.Date("2018-01-01"), as.Date("2020-12-30")), 
         ylim = c(0, max(tissue_data$tissue_scaled)), 
         xlab = "", ylab = "", axes = FALSE)
    
    # Add the secondary y-axis (right axis)
    axis(side = 4, at = pretty(tissue_data$tissue_scaled), 
         labels = pretty(tissue_data$tissue_scaled) / rescaler)
    
    # Add labels to the secondary axis
    mtext("Recovered tissue (m2)", side = 4, line = 3)
    
    # Add a legend
    legend("topright", legend = c("SST (°C)", "Recovered tissue (m2)"), 
           col = c("#E69F00", "#0072B2"), lwd = 2, pch = c(NA, 16), bty = "n")
  }
  
  ################################## More mature sigmoid & logistic behaviors of temperature ##################################
  
  # Load necessary libraries
  library(ggplot2)
  library(viridis)
  library(patchwork)  # For arranging multiple plots
  
  # Define temperature range from 23 to 33
  T <- seq(23, 33, length.out = 100)
  
  #APPROACH 1: Inflection point at 30.5C (original)
  temp_effect_inflection <- function(temp, T_inflection = 30.5, zeta_bounded = 0.5) {
    # Convert bounded zeta to effective zeta
    zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    
    # Logistic function with inflection point at T_inflection
    modifier <- 1 / (1 + exp(-zeta_effective * (temp - T_inflection)))
    
    return(modifier)
  }
  
  #APPROACH 2: Original model clarified
  temp_effect_original <- function(temp, T_min = 23, zeta_bounded = 0.5) {
    
    # Convert bounded zeta to effective zeta
    
    # zeta_effective <- 5 * zeta_bounded/(1-0.7*zeta_bounded)
    zeta_effective <- zeta_bounded/(1 - zeta_bounded)
    
    # Original formula clarified
    # This is a non-centered sigmoid that:
    # - Starts near 0 at T_min
    # - Increases asymptotically toward 1 as temp increases
    # - Has no inflection point assumption
    modifier <- (1 - exp(-zeta_effective * (temp - T_min))) / (1 + exp(-zeta_effective * (temp - T_min)))
    
    return(modifier)
  }
  
  #APPROACH 3: Decreasing function with temperature
  temp_effect_decreasing <- function(temp, T_max = 33, T_min = 23, zeta_bounded = 0.5) {
    
    # Convert bounded zeta to effective zeta
    # zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    zeta_effective <- zeta_bounded / (1 - zeta_bounded)
    
    # Calculate midpoint of temperature range for centering the sigmoid
    # T_mid <- (T_min + T_max) / 2
    T_mid = 30.5
    
    # Decreasing sigmoid that starts near 1 (at low temps) and approaches 0 as temp increases
    # Note the positive sign before zeta_effective which makes the function decrease with temperature
    modifier <- 1 / (1 + exp(zeta_effective * (temp - T_mid)))
    
    return(modifier)
  }
  
  #APPROACH 4: Non-centered sigmoid decreasing with temperature (FIXED)
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
  zeta_values <- c(0, 0.01, 0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  
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
      title = "APPROACH 1: Temperature Threshold with Classic Logistic Function",
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
      title = "APPROACH 2: Temperature Effect with non-centered Sigmoid Function",
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
      title = "APPROACH 3: Temperature Effect with Decreasing Sigmoid Function",
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
      title = "APPROACH 4: Temperature Effect with Inverted Non-centered Sigmoid",
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
  
  
  ################################## Possibly convoluted power function ##################################
  #APPROACH 4: Modified non-centered sigmoid decreasing with temperature (NORTHEAST BOWING)
  temp_effect_noncentered_decreasing_northeast <- function(temp, T_min = 23, T_max = 33, zeta_bounded = 0.5) {
    
    # Convert bounded zeta to effective zeta
    zeta_effective <- 1 * zeta_bounded/(1-1*zeta_bounded)
    # zeta_effective <- 5 * zeta_bounded/(1-0.8*zeta_bounded)
    
    
    # Calculate the relative position in the temperature range
    rel_temp <- (temp - T_min) / (T_max - T_min)
    
    # Apply a modified decreasing function that bows northeast
    # Using a power function to create the convex shape
    # Higher powers of rel_temp create more convex curves
    
    # The power is determined by zeta - higher zeta means higher power (more convex)
    power <- 1 + 2 * zeta_effective
    
    # Decreasing function that equals 1 at T_min and approaches 0 at T_max
    # with convex (northeast bowing) shape
    modifier <- 1 - rel_temp^power
    
    return(modifier)
  }
  
  # Add data generation for the new approach
  plot_data_northeast_bowing <- data.frame()
  
  for (zeta in zeta_values) {
    # Data for new approach
    modifiers_northeast <- sapply(T, function(t) temp_effect_noncentered_decreasing_northeast(t, zeta_bounded = zeta))
    temp_data_northeast <- data.frame(
      Temperature = T,
      Modifier = modifiers_northeast,
      Zeta = as.factor(zeta)
    )
    plot_data_northeast_bowing <- rbind(plot_data_northeast_bowing, temp_data_northeast)
  }
  
  # Plot for the new approach
  p5 <- ggplot(plot_data_northeast_bowing, aes(x = Temperature, y = Modifier, color = Zeta, group = Zeta)) +
    geom_line(size = 1) +
    labs(
      title = "Temperature Effect with Northeast Bowing Curves",
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
  
  # To view just this plot
  p5
  
  # Explanation of the new approach
  cat("Approach 5: Northeast Bowing Decreasing Function\n")
  cat("Formula: modifier = 1 - rel_temp^power where power = 1 + 2 * zeta_effective\n")
  cat("- At T_min (23°C), modifier = 1 exactly for all zeta values\n")
  cat("- Decreases toward 0 as temperature increases\n")
  cat("- Higher zeta values create more convex curves (bowling northeast)\n")
  cat("- As zeta increases, the drop in modifier becomes more delayed then steeper\n")
  ################################## Save output ##################################

  # #pass workspace to downstream script
  # save.image(file = here("output", "SST_effect_workspace.RData"))
  