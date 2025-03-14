

# NOTE - run this while working within plots_basic.R



################################## More mature sigmoid & logistic behaviors of temperature ##################################
# AKA - zeta / eta parameter plotting sandbox

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
