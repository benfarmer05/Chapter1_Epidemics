# .rs.restartR(clean = TRUE)
rm(list = ls())

library(here)
library(tidyverse)
library(DEoptim)
library(deSolve)

################################## Set-up ##################################

load(here("output/multi_SIR_workspace.RData"))

curr.host <- "Multi-host"

obs.model <- obs %>%
  mutate(Site = case_when(
    Site == "Offshore"   ~ "off",
    Site == "Midchannel" ~ "mid",
    Site == "Nearshore"  ~ "near",
    TRUE ~ Site
  ))

DHW.modifier <- 0

sites <- unique(summary$site)

################################## Data prep ##################################

days_sites <- summary %>%
  group_by(site) %>%
  summarize(
    days.obs = list({
      d <- na.omit(days.inf.site)
      d[which(!is.na(d))[1]:length(d)]
    }),
    days.model = list(seq(
      from = min(na.omit(days.inf.site)),
      to   = max(na.omit(days.inf.site)),
      by   = 1
    )),
    .groups = "drop"
  )

site_data_list <- lapply(sites, function(s) {
  
  days.obs   <- days_sites %>% filter(site == s) %>% pull(days.obs)   %>% unlist()
  days.model <- days_sites %>% filter(site == s) %>% pull(days.model) %>% unlist()
  
  N.site <- susceptible_ref %>% filter(Site == s) %>% slice(1) %>% pull(N.site)
  N.LS   <- susceptible_ref %>% filter(Site == s, Category == "Low")      %>% pull(tissue_ref)
  N.MS   <- susceptible_ref %>% filter(Site == s, Category == "Moderate") %>% pull(tissue_ref)
  N.HS   <- susceptible_ref %>% filter(Site == s, Category == "High")     %>% pull(tissue_ref)
  
  pull_tiss <- function(cat, comp) {
    obs.model %>%
      filter(Site == s, Category == cat, Compartment == comp) %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)
  }
  
  HS.inf.raw <- pull_tiss("High", "Infected")
  first_idx  <- which(!is.na(HS.inf.raw) & HS.inf.raw != 0)[1]
  trim       <- function(x) x[first_idx:length(x)]
  
  list(
    site       = s,
    days.obs   = days.obs,
    days.model = days.model,
    N.site     = N.site,
    N.LS       = N.LS,
    N.MS       = N.MS,
    N.HS       = N.HS,
    LS.remtiss = trim(pull_tiss("Low",      "Dead")),
    MS.remtiss = trim(pull_tiss("Moderate", "Dead")),
    HS.remtiss = trim(pull_tiss("High",     "Dead")),
    LS.inftiss = trim(pull_tiss("Low",      "Infected")),
    MS.inftiss = trim(pull_tiss("Moderate", "Infected")),
    HS.inftiss = trim(pull_tiss("High",     "Infected")),
    init = c(
      S.LS = N.LS, I.LS = 0, R.LS = 0,
      S.MS = N.MS, I.MS = 0, R.MS = 0,
      S.HS = N.HS - HS.inf.raw[first_idx], I.HS = HS.inf.raw[first_idx], R.HS = 0
    )
  )
})

# K fixed at largest observed total community size across all sites,
# consistent with single-host approach: one community-level reference
K <- max(sapply(site_data_list, function(sd) sd$N.site))

cat(sprintf("K: %.4f\n", K))

################################## q toggle ##################################

# Set fit_q = TRUE to optimize q alongside beta and gamma.
# Set fit_q = FALSE to fix q at q_val_fixed.
#
# q controls the scaling of the contact rate with population size,
# following Smith et al. (2009):
#   q = 1  -> pure frequency dependence  (contact rate independent of N0)
#   q = 0  -> pure density dependence    (contact rate proportional to N0)
#   q > 1  -> super-frequency dependence (contact rate declines with N0)
#   q < 0  -> super-density dependence   (contact rate accelerates with N0)

fit_q       <- FALSE  # <-- TOGGLE HERE
q_val_fixed <- 1.0    # used only when fit_q = FALSE

################################## Model ##################################

SIR.multi <- function(t, y, p) {
  S.LS <- y[1]; I.LS <- y[2]; R.LS <- y[3]
  S.MS <- y[4]; I.MS <- y[5]; R.MS <- y[6]
  S.HS <- y[7]; I.HS <- y[8]; R.HS <- y[9]
  
  with(as.list(p), {
    
    # Total initial community size at this site — denominator for
    # community-level frequency-dependent scaling, following Smith et al. (2009)
    N0 <- N.LS + N.MS + N.HS
    
    # Pooled infectious pressure (raw count, no normalization)
    P <- I.LS + I.MS + I.HS
    
    # Contact rates: (K_i / N0)^q * S_i * P
    # mirrors single-host: (K / N)^q * S * I
    # q = 1 -> (K_i/N0)^1: frequency-dependent, portable across sites
    # q = 0 -> (K_i/N0)^0 = 1: density-dependent, no population scaling
    # contact.LS <- (K^q) * S.LS * P / (N0^q) # NOTE - community freq.
    # contact.MS <- (K^q) * S.MS * P / (N0^q) #q2...
    # contact.HS <- (K^q) * S.HS * P / (N0^q)
    # contact.LS <- (K^q.LS) * S.LS * P / (N0^q.LS) # NOTE - community freq.
    # contact.MS <- (K^q.MS) * S.MS * P / (N0^q.MS) #q2...
    # contact.HS <- (K^q.HS) * S.HS * P / (N0^q.HS)
    contact.LS <- (K^q.LS) * S.LS * P / (N.LS^q.LS) # NOTE - within-group freq.
    contact.MS <- (K^q.MS) * S.MS * P / (N.MS^q.MS)
    contact.HS <- (K^q.HS) * S.HS * P / (N.HS^q.HS)
    
    # Within-group scaling (disabled; left for reference)
    # contact.LS <- (K.LS / N.LS)^q * S.LS * P
    # contact.MS <- (K.MS / N.MS)^q * S.MS * P
    # contact.HS <- (K.HS / N.HS)^q * S.HS * P
    
    dS.LS <- -b.LS * contact.LS
    dI.LS <-  b.LS * contact.LS - g.LS * I.LS
    dR.LS <-  g.LS * I.LS
    
    dS.MS <- -b.MS * contact.MS
    dI.MS <-  b.MS * contact.MS - g.MS * I.MS
    dR.MS <-  g.MS * I.MS
    
    dS.HS <- -b.HS * contact.HS
    dI.HS <-  b.HS * contact.HS - g.HS * I.HS
    dR.HS <-  g.HS * I.HS
    
    list(c(dS.LS, dI.LS, dR.LS,
           dS.MS, dI.MS, dR.MS,
           dS.HS, dI.HS, dR.HS),
         P = P)
  })
}

################################## Objective function ##################################

objective_function_global <- function(params, site_data_list, q_fixed = NULL) {
  
  b.LS <- params[1]; g.LS <- params[2]
  b.MS <- params[3]; g.MS <- params[4]
  b.HS <- params[5]; g.HS <- params[6]
  q    <- if (is.null(q_fixed)) params[7] else q_fixed
  
  total_residual <- 0
  
  for (sd in site_data_list) {
    
    SIR.out <- tryCatch({
      data.frame(ode(
        y     = sd$init,
        times = sd$days.model,
        func  = SIR.multi,
        parms = list(b.LS = b.LS, g.LS = g.LS,
                     b.MS = b.MS, g.MS = g.MS,
                     b.HS = b.HS, g.HS = g.HS,
                     N.LS = sd$N.LS, N.MS = sd$N.MS, N.HS = sd$N.HS,
                     K    = K,
                     q    = q)
      ))
    }, error = function(e) NULL)
    
    if (is.null(SIR.out)) {
      total_residual <- total_residual + 1e9
      next
    }
    
    rows <- which(SIR.out$time %in% sd$days.obs)
    
    sim.LS.rem <- SIR.out[rows, "R.LS"]
    sim.MS.rem <- SIR.out[rows, "R.MS"]
    sim.HS.rem <- SIR.out[rows, "R.HS"]
    
    sim.LS.inf <- SIR.out[rows, "I.LS"]
    sim.MS.inf <- SIR.out[rows, "I.MS"]
    sim.HS.inf <- SIR.out[rows, "I.HS"]
    
    n.LS <- min(length(sim.LS.rem), length(sd$LS.remtiss))
    n.MS <- min(length(sim.MS.rem), length(sd$MS.remtiss))
    n.HS <- min(length(sim.HS.rem), length(sd$HS.remtiss))
    
    n.LS.i <- min(length(sim.LS.inf), length(sd$LS.inftiss))
    n.MS.i <- min(length(sim.MS.inf), length(sd$MS.inftiss))
    n.HS.i <- min(length(sim.HS.inf), length(sd$HS.inftiss))
    
    site_residual <- (
      sum(abs(sim.LS.rem[1:n.LS] - sd$LS.remtiss[1:n.LS])) +
        sum(abs(sim.MS.rem[1:n.MS] - sd$MS.remtiss[1:n.MS])) +
        sum(abs(sim.HS.rem[1:n.HS] - sd$HS.remtiss[1:n.HS])) +
        sum(abs(sim.LS.inf[1:n.LS.i] - sd$LS.inftiss[1:n.LS.i])) +
        sum(abs(sim.MS.inf[1:n.MS.i] - sd$MS.inftiss[1:n.MS.i])) +
        sum(abs(sim.HS.inf[1:n.HS.i] - sd$HS.inftiss[1:n.HS.i]))
    ) / sd$N.site
    
    total_residual <- total_residual + site_residual
  }
  
  total_residual
}

# ################################## Convert site-specific betas to K-scaled betas ##################################
# 
# # NOTE - only run this if setting q = 1 and not fitting q
# 
# q.for.conversion <- q_val_fixed  # or a specific numeric value e.g. 1.0
# 
# # Prior betas from site-specific script (no K scaling, pure frequency-dependent)
# # Replace these with your actual fitted values
# prior_b.LS <- c(near = 0.03, mid = 0.12, off = 0.22)
# prior_b.MS <- c(near = 0.14, mid = 0.10, off = 0.28)
# prior_b.HS <- c(near = 2.08, mid = 1.94, off = 0.48)
# 
# # Convert: beta_new = beta_old / (K / N)^q
# # i.e. undo the (K/N)^q multiplier that the new script applies
# # This gives the beta the global script needs to reproduce the same effective transmission rate
# converted_b.LS <- sapply(site_data_list, function(sd)
#   prior_b.LS[sd$site] / (K / sd$N.LS)^q.for.conversion)
# converted_b.MS <- sapply(site_data_list, function(sd)
#   prior_b.MS[sd$site] / (K / sd$N.MS)^q.for.conversion)
# converted_b.HS <- sapply(site_data_list, function(sd)
#   prior_b.HS[sd$site] / (K / sd$N.HS)^q.for.conversion)
# 
# cat("Converted betas (K-scaled, per site):\n")
# cat(sprintf("  LS: near=%.6f  mid=%.6f  off=%.6f\n",
#             converted_b.LS[1], converted_b.LS[2], converted_b.LS[3]))
# cat(sprintf("  MS: near=%.6f  mid=%.6f  off=%.6f\n",
#             converted_b.MS[1], converted_b.MS[2], converted_b.MS[3]))
# cat(sprintf("  HS: near=%.6f  mid=%.6f  off=%.6f\n",
#             converted_b.HS[1], converted_b.HS[2], converted_b.HS[3]))
# 
# # Summary values to center bounds around — using mean across sites
# # since the global script fits one beta per group across all sites
# center_b.LS <- mean(converted_b.LS)
# center_b.MS <- mean(converted_b.MS)
# center_b.HS <- mean(converted_b.HS)
# 
# tol <- 0.5  # additive tolerance around converted priors — adjust as needed
# 
# lower_constrained <- c(max(0, center_b.LS - tol), 0,
#            max(0, center_b.MS - tol), 0,
#            max(0, center_b.HS - tol), 0)
# upper_constrained <- c(center_b.LS + tol, 4,
#            center_b.MS + tol, 4,
#            center_b.HS + tol, 4)
# 
# 
################################## Optimization ##################################

control <- list(itermax = 300)

if (fit_q) {
  cat("Fitting beta, gamma, AND q...\n")
  result_global <- DEoptim(
    fn             = objective_function_global,
    lower          = c(0, 0, 0, 0, 0, 0,  0.5),
    upper          = c(4, 4, 4, 4, 4, 4,  1.0),    site_data_list = site_data_list,
    q_fixed        = NULL,
    control        = control
  )
  q.global <- result_global$optim$bestmem[7]
} else {
  cat(sprintf("Fitting beta and gamma with q fixed at %.2f...\n", q_val_fixed))
  result_global <- DEoptim(
    fn             = objective_function_global,
    lower          = c(0, 0, 0, 0, 0, 0),
    upper          = c(4, 4, 4, 4, 4, 4),
    # lower = lower_constrained, # NOTE - only run if setting q = 1 and not fitting q
    # upper = upper_constrained,
    site_data_list = site_data_list,
    q_fixed        = q_val_fixed,
    control        = control
  )
  q.global <- q_val_fixed
}

b.LS.global <- result_global$optim$bestmem[1]
g.LS.global <- result_global$optim$bestmem[2]
b.MS.global <- result_global$optim$bestmem[3]
g.MS.global <- result_global$optim$bestmem[4]
b.HS.global <- result_global$optim$bestmem[5]
g.HS.global <- result_global$optim$bestmem[6]

cat("\nGlobal fitted parameters:\n")
cat(sprintf("  beta  LS: %.6f   gamma LS: %.6f   R0 LS: %.4f\n",
            b.LS.global, g.LS.global, b.LS.global / g.LS.global))
cat(sprintf("  beta  MS: %.6f   gamma MS: %.6f   R0 MS: %.4f\n",
            b.MS.global, g.MS.global, b.MS.global / g.MS.global))
cat(sprintf("  beta  HS: %.6f   gamma HS: %.6f   R0 HS: %.4f\n",
            b.HS.global, g.HS.global, b.HS.global / g.HS.global))
cat(sprintf("  q: %.4f  (%s)\n", q.global,
            if (fit_q) "fitted" else "fixed"))

################################## Diagnostics ##################################

SUS_COLORS  <- c(Low = "#1E90FF", Moderate = "#FFD700", High = "#FF1493")
site_labels <- c(near = "Nearshore", mid = "Midchannel", off = "Offshore")

global_sims <- lapply(site_data_list, function(sd) {
  data.frame(ode(
    y     = sd$init,
    times = sd$days.model,
    func  = SIR.multi,
    parms = list(b.LS = b.LS.global, g.LS = g.LS.global,
                 b.MS = b.MS.global, g.MS = g.MS.global,
                 b.HS = b.HS.global, g.HS = g.HS.global,
                 N.LS = sd$N.LS, N.MS = sd$N.MS, N.HS = sd$N.HS,
                 K    = K,
                 q    = q.global)
  )) %>%
    select(-any_of("P")) %>%
    mutate(site = sd$site)
}) %>%
  bind_rows() %>%
  pivot_longer(cols = -c(time, site),
               names_pattern = "(.*)\\.(..)",
               names_to = c("Compartment", "Susceptibility")) %>%
  mutate(
    Compartment    = case_match(Compartment,
                                "S" ~ "Susceptible", "I" ~ "Infected", "R" ~ "Dead"),
    Susceptibility = case_match(Susceptibility,
                                "LS" ~ "Low", "MS" ~ "Moderate", "HS" ~ "High"),
    site_label     = site_labels[site]
  )

obs_plot <- obs.model %>%
  filter(Category != "Total") %>%
  rename(Susceptibility = Category) %>%
  mutate(site_label = site_labels[Site])

p_global_fit <- ggplot(global_sims %>% filter(Compartment == "Dead"),
                       aes(time, value, color = Susceptibility)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_plot %>% filter(Compartment == "Dead"),
             aes(days.inf.site, tissue, color = Susceptibility),
             shape = 15, size = 1.5) +
  scale_color_manual(values = SUS_COLORS) +
  facet_wrap(~ site_label, ncol = 1, scales = "free_y") +
  labs(x     = "Day of outbreak",
       y     = "Surface area of tissue (m²)",
       title = sprintf(
         "Global multi-host SIR fit (q = %.3f, %s)\nβ_LS=%.4f γ_LS=%.4f | β_MS=%.4f γ_MS=%.4f | β_HS=%.4f γ_HS=%.4f",
         q.global, if (fit_q) "fitted" else "fixed",
         b.LS.global, g.LS.global,
         b.MS.global, g.MS.global,
         b.HS.global, g.HS.global
       )) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p_global_fit)