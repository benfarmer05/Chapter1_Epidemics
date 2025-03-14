  # Load necessary libraries
  library(ggplot2)
  library(viridis)  # For colorblind-friendly palettes
  
  # Define temperature range from 23 to 33
  T <- seq(23, 33, length.out = 100)
  
  # Function to compute y based on zeta and temperature T
  compute_y <- function(T, zeta) {
    # When zeta=0, this will return 0
    # As zeta increases, the curve will asymptote to 1 as temperature increases
    y <- (1 - exp(-zeta * (T - 23))) / (1 + exp(-zeta * (T - 23)))
    return(y)
  }
  
  # Define a broader range of zeta values
  # zeta_values <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0)
  zeta_values = seq(0, 2, length.out = 1000)
  
  # Create a data frame to store the results
  data <- data.frame()
  for (zeta in zeta_values) {
    temp_data <- data.frame(
      Temperature = T,
      Zeta = rep(zeta, length(T)),
      Y = compute_y(T, zeta)
    )
    data <- rbind(data, temp_data)
  }
  
  # Convert Zeta to factor for better legend display
  data$Zeta_factor <- factor(data$Zeta, levels = zeta_values)
  
  # # Plot the results with viridis color palette (colorblind-friendly)
  # ggplot(data, aes(x = Temperature, y = Y, color = Zeta_factor, group = Zeta_factor)) +
  #   geom_line(size = 1) +
  #   labs(
  #     title = "Effect of Zeta on Y vs Temperature",
  #     subtitle = "Higher zeta values result in steeper transitions",
  #     x = "Temperature (°C)", 
  #     y = "Y",
  #     color = "Zeta Value"
  #   ) +
  #   theme_minimal(base_size = 12) +
  #   theme(
  #     legend.position = "right",
  #     panel.grid.minor = element_blank(),
  #     plot.title = element_text(face = "bold"),
  #     legend.title = element_text(face = "bold")
  #   ) +
  #   scale_color_viridis_d(option = "plasma") +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  #   scale_x_continuous(breaks = seq(23, 33, 1))
  
  # Plot the results with viridis color palette (colorblind-friendly)
  ggplot(data, aes(x = Temperature, y = Y)) +
    geom_line(size = 1) +
    scale_color_viridis_d(option = "plasma") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(breaks = seq(23, 33, 1))
  