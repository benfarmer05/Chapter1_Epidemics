    
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  ################################## SST logic ##################################

  
  # if (t >= thermal_onset_time) {
  #   
  #   # SST-based logic for increasing or decreasing rates
  #   if (current_SST < SST_threshold) {
  #     # SST below threshold - increase rates
  #     # transmission_rate <- transmission_rate * (1 + z * (SST_threshold - current_SST))
  #     # removal_rate <- removal_rate * (1 + e * (SST_threshold - current_SST))
  #     transmission_rate <- transmission_rate
  #     removal_rate <- removal_rate
  #   } else {
  #     # SST above threshold - decrease rates
  #     transmission_rate <- transmission_rate * (1 / (1 + exp(-z * (current_SST - SST_threshold))))
  #     removal_rate <- removal_rate * (1 / (1 + exp(-e * (current_SST - SST_threshold))))
  #   }
  # }
  
  
  ################################## Full range of SST ##################################

  # Define parameters
  SST_threshold <- 30.5  # Temperature threshold
  z_values <- c(0.5, 1, 2.5, 3)  # Different z values for transmission effect
  e_values <- c(0.0001, 1, 2, 3)  # Different e values for mortality effect
  SST_range <- seq(29, 32, by = 0.1)  # SST from 29 to 32 degrees
  
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
  
  
  
  
  
  ################################## Wide range ##################################
  
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
  