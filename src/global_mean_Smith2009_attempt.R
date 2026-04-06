# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(tidyverse)
library(DEoptim)
library(deSolve)

################################## Set-up ##################################

load(here("output", "plots_obs_workspace.RData"))

curr.host = 'Single-host'

obs.model <- obs %>%
  mutate(
    Site = case_when(
      Site == "Offshore"   ~ "off",
      Site == "Midchannel" ~ "mid",
      Site == "Nearshore"  ~ "near",
      TRUE ~ Site
    )
  )

DHW.onset.date = as.Date("2019-07-01")
DHW.modifier   = 0

SST_threshold_value = 30.5
DHW_threshold_value = NA

min_sst_date <- DHW.CRW %>%
  filter(SST.90th_HS == min(SST.90th_HS, na.rm = TRUE)) %>%
  arrange(date) %>%
  slice(1) %>%
  pull(date)

date_threshold <- DHW.CRW %>%
  filter(date > min_sst_date, SST.90th_HS >= SST_threshold_value) %>%
  arrange(date) %>%
  slice(1) %>%
  pull(date) - 1

error_eval <- error_eval %>%
  mutate(
    SST_threshold = if_else(host == 'Single-host', SST_threshold_value, SST_threshold),
    date_thresh   = if_else(host == 'Single-host', date_threshold, date_thresh)
  )

lambda.modifier = 1.0
offset          = 1 - 1 / (1 + exp(-lambda.modifier))

sites = unique(summary$site)

################################## q toggle ##################################

# Set fit_q = TRUE to optimize q alongside beta and gamma.
# Set fit_q = FALSE to fix q at q_val_fixed.
fit_q       = FALSE   # <-- TOGGLE HERE
q_val_fixed = 1.0     # used only when fit_q = FALSE

################################## define cover params ##################################

alpha_val = 0.13
K_val     = 172

################################## Model: single-host ##################################

SIR = function(t, y, p) {
  S = y[1]; I = y[2]; R = y[3]
  with(as.list(p), {
    contact_rate = (K^q) * S * I / (N^q)
    dS.dt = -b * contact_rate
    dI.dt =  b * contact_rate - g * I
    dR.dt =  g * I
    return(list(c(dS.dt, dI.dt, dR.dt)))
  })
}

################################## Data prep ##################################

days_sites <- summary %>%
  group_by(site) %>%
  summarize(
    days = list(days.inf.site[1:length(days.inf.site)]),
    days.obs = list({
      days_cleaned <- na.omit(days.inf.site)
      days_cleaned[which(!is.na(days_cleaned))[1]:length(days_cleaned)]
    }),
    days.model = list(seq(from = min(na.omit(days.inf.site)),
                          to   = max(na.omit(days.inf.site)), by = 1))
  )

site_data_list = lapply(seq_along(sites), function(i) {
  
  site.loop = sites[i]
  
  days.obs   <- days_sites %>% filter(site == site.loop) %>% pull(days.obs)   %>% unlist()
  days.model <- days_sites %>% filter(site == site.loop) %>% pull(days.model) %>% unlist()
  
  N.site = susceptible_ref %>% filter(Site == site.loop) %>% slice(1) %>% pull(N.site)
  
  inftiss = obs.model %>%
    filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
    pull(tissue)
  
  remtiss = obs.model %>%
    filter(Site == site.loop, Category == "Total", Compartment == "Dead") %>%
    pull(tissue)
  
  first_valid_idx <- which(!is.na(
    days_sites %>% filter(site == site.loop) %>% pull(days) %>% unlist()
  ))[1]
  
  inftiss <- inftiss[first_valid_idx:length(inftiss)]
  remtiss <- remtiss[first_valid_idx:length(remtiss)]
  
  list(
    site       = site.loop,
    days.obs   = days.obs,
    days.model = days.model,
    N.site     = N.site,
    inftiss    = inftiss,
    remtiss    = remtiss,
    S0         = N.site - inftiss[1],
    I0         = inftiss[1],
    R0_init    = 0
  )
})

################################## Objective function ##################################

objective_function_global = function(params, site_data_list, K, q_fixed = NULL) {
  
  beta  = params[1]
  gamma = params[2]
  # if q_fixed is NULL, q is being fitted as params[3]
  q     = if (is.null(q_fixed)) params[3] else q_fixed
  
  total_residual = 0
  
  for (i in seq_along(site_data_list)) {
    
    sd = site_data_list[[i]]
    
    SIR.out = tryCatch({
      data.frame(ode(
        y     = c(S = sd$S0, I = sd$I0, R = sd$R0_init),
        times = sd$days.model,
        func  = SIR,
        parms = list(b = beta, g = gamma, N = sd$N.site, q = q, K = K)
      ))
    }, error = function(e) return(NULL))
    
    if (is.null(SIR.out)) {
      total_residual = total_residual + 1e9
      next
    }
    
    sim_R = SIR.out[which(SIR.out$time %in% sd$days.obs), "R"]
    obs_R = sd$remtiss
    n     = min(length(sim_R), length(obs_R))
    
    # site_residual = sum(abs(sim_R[1:n] - obs_R[1:n]))
    site_residual = sum(abs(sim_R[1:n] - obs_R[1:n])) / sd$N.site
    # site_residual = sum(abs(sim_R[1:n] - obs_R[1:n])) / max(obs_R, na.rm = TRUE)
    
    total_residual = total_residual + site_residual
  }
  
  return(total_residual)
}

################################## Optimization ##################################

control = list(itermax = 200)

if (fit_q) {
  cat("Fitting beta, gamma, AND q...\n")
  result_global = DEoptim(
    fn             = objective_function_global,
    lower          = c(0, 0, 0),
    upper          = c(4, 4, 2),
    site_data_list = site_data_list,
    K              = K_val,
    q_fixed        = NULL,      # q is free
    control        = control
  )
  beta_global  = result_global$optim$bestmem[1]
  gamma_global = result_global$optim$bestmem[2]
  q_global     = result_global$optim$bestmem[3]
} else {
  cat(sprintf("Fitting beta and gamma with q fixed at %.2f...\n", q_val_fixed))
  result_global = DEoptim(
    fn             = objective_function_global,
    lower          = c(0, 0),
    upper          = c(4, 4),
    site_data_list = site_data_list,
    K              = K_val,
    q_fixed        = q_val_fixed,  # q is fixed
    control        = control
  )
  beta_global  = result_global$optim$bestmem[1]
  gamma_global = result_global$optim$bestmem[2]
  q_global     = q_val_fixed
}

R0_global = beta_global * K_val / gamma_global  # valid at q = 1; see note below

cat("\nGlobal fitted parameters:\n")
cat(sprintf("  Beta  (global): %.6f\n", beta_global))
cat(sprintf("  Gamma (global): %.6f\n", gamma_global))
cat(sprintf("  q     (%s): %.4f\n",
            if (fit_q) "fitted" else "fixed ", q_global))
cat(sprintf("  R0 (beta*K/gamma, exact at q=1): %.4f\n", R0_global))

################################## Per-site diagnostics ##################################

global_SIR_output = lapply(seq_along(sites), function(i) {
  
  sd = site_data_list[[i]]
  
  SIR.out = data.frame(ode(
    y     = c(S = sd$S0, I = sd$I0, R = sd$R0_init),
    times = sd$days.model,
    func  = SIR,
    parms = list(b = beta_global, g = gamma_global,
                 N = sd$N.site, q = q_global, K = K_val)
  ))
  
  sim_R = SIR.out[which(SIR.out$time %in% sd$days.obs), "R"]
  obs_R = sd$remtiss
  n     = min(length(sim_R), length(obs_R))
  
  site_residual = sum(abs(sim_R[1:n] - obs_R[1:n]))
  
  cat("Site:", sd$site, "| Residual:", round(site_residual, 3),
      "| N:", sd$N.site, "\n")
  
  list(site = sd$site, SIR.out = SIR.out, sim_R = sim_R, obs_R = obs_R)
})

################################## Visualize global fit ##################################

site_labels <- c(near = "Nearshore", mid = "Midchannel", off = "Offshore")
site_colors <- c(Nearshore = "orange", Midchannel = "purple", Offshore = "magenta")

global_sims <- lapply(seq_along(sites), function(i) {
  sd <- site_data_list[[i]]
  data.frame(ode(
    y     = c(S = sd$S0, I = sd$I0, R = sd$R0_init),
    times = sd$days.model,
    func  = SIR,
    parms = list(b = beta_global, g = gamma_global,
                 N = sd$N.site, q = q_global, K = K_val)
  )) %>% mutate(site = sd$site)
}) %>%
  bind_rows() %>%
  pivot_longer(c(S, I, R), names_to = "Compartment", values_to = "tissue") %>%
  mutate(
    Compartment = case_match(Compartment,
                             "S" ~ "Susceptible", "I" ~ "Infected", "R" ~ "Dead"),
    site_label  = site_labels[site]
  )

obs_global <- obs.model %>%
  filter(Category == "Total") %>%
  mutate(site_label = case_match(Site,
                                 "near" ~ "Nearshore",
                                 "mid"  ~ "Midchannel",
                                 "off"  ~ "Offshore"))

q_label = if (fit_q) sprintf("q = %.3f (fitted)", q_global) else sprintf("q = %.2f (fixed)", q_global)

p_global_fit <- ggplot(global_sims,
                       aes(time, tissue,
                           color = site_label,
                           group = interaction(site_label, Compartment))) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_global,
             aes(days.inf.site, tissue, shape = Compartment, color = site_label),
             size = 1.8, inherit.aes = FALSE) +
  scale_color_manual(values = site_colors, name = "Site") +
  scale_shape_manual(values = c(Susceptible = 16, Infected = 17, Dead = 15), name = NULL) +
  facet_wrap(~ site_label, ncol = 1, scales = "free_y") +
  labs(x     = "Day of outbreak",
       y     = "Surface area of tissue (m²)",
       title = sprintf("Global SIR fit  |  β = %.4f, γ = %.4f, R₀ = %.3f  |  %s",
                       beta_global, gamma_global, R0_global, q_label)) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p_global_fit)





# ── Smith 2009 Fig. 2A recreation ─────────────────────────────────────────────

N_seq  <- seq(0, K_val, length.out = 500)
q_vals <- c(1, 0, -1, -5)  # reference curves as in Smith Fig. 2A

# Build reference curves
smith_curves <- expand.grid(N = N_seq, q = q_vals) %>%
  mutate(
    # contact = N^(1 - q),
    contact = N^(1 - q) * K_val^(q - 1),
    q_label = factor(q)
  )

# Build fitted q curve (or fixed q, whichever was used)
# fitted_curve <- data.frame(N = N_seq) %>%
#   mutate(contact = N^(1 - q_global))

fitted_curve <- data.frame(N = N_seq) %>%
  mutate(contact = N^(1 - q_global) / K_val^(1 - q_global))

# Label positions at N = K_val for each reference curve
label_pts <- smith_curves %>%
  group_by(q) %>%
  slice_max(N, n = 1) %>%
  mutate(q_label_txt = paste0("q=", q))

p_smith_A <- ggplot() +
  geom_line(data = smith_curves,
            aes(N, contact, group = q_label),
            color = "gray50", linewidth = 0.6) +
  geom_line(data = fitted_curve,
            aes(N, contact),
            color = "black", linewidth = 1.2) +
  geom_text(data = label_pts,
            aes(N, contact, label = q_label_txt),
            hjust = -0.1, size = 3, color = "gray30") +
  annotate("text", x = 10, y = max(fitted_curve$contact) * 0.85,
           label = sprintf("fitted q = %.3f", q_global),
           hjust = 0, size = 3.2, fontface = "bold") +
  scale_x_continuous("Population size, N",
                     limits = c(0, K_val * 1.15), expand = c(0, 0)) +
  scale_y_continuous("Scaled contact rate",
                     expand = c(0, 0)) +
  theme_classic() +
  ggtitle("A")

print(p_smith_A)
