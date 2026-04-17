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
    init = c(
      S.LS = N.LS, I.LS = 0, R.LS = 0,
      S.MS = N.MS, I.MS = 0, R.MS = 0,
      S.HS = N.HS - HS.inf.raw[first_idx], I.HS = HS.inf.raw[first_idx], R.HS = 0
    )
  )
})

# K_sg fixed at largest observed N_sg across all sites
K.LS <- max(sapply(site_data_list, function(sd) sd$N.LS))
K.MS <- max(sapply(site_data_list, function(sd) sd$N.MS))
K.HS <- max(sapply(site_data_list, function(sd) sd$N.HS))

cat(sprintf("K.LS: %.4f  K.MS: %.4f  K.HS: %.4f\n", K.LS, K.MS, K.HS))

################################## Model ##################################

SIR.multi <- function(t, y, p) {
  S.LS <- y[1]; I.LS <- y[2]; R.LS <- y[3]
  S.MS <- y[4]; I.MS <- y[5]; R.MS <- y[6]
  S.HS <- y[7]; I.HS <- y[8]; R.HS <- y[9]
  
  with(as.list(p), {
    P <- I.LS + I.MS + I.HS
    
    dS.LS <- -b.LS * S.LS * P * (K.LS / N.LS)^q
    dI.LS <-  b.LS * S.LS * P * (K.LS / N.LS)^q - g.LS * I.LS
    dR.LS <-  g.LS * I.LS
    
    dS.MS <- -b.MS * S.MS * P * (K.MS / N.MS)^q
    dI.MS <-  b.MS * S.MS * P * (K.MS / N.MS)^q - g.MS * I.MS
    dR.MS <-  g.MS * I.MS
    
    dS.HS <- -b.HS * S.HS * P * (K.HS / N.HS)^q
    dI.HS <-  b.HS * S.HS * P * (K.HS / N.HS)^q - g.HS * I.HS
    dR.HS <-  g.HS * I.HS
    
    list(c(dS.LS, dI.LS, dR.LS,
           dS.MS, dI.MS, dR.MS,
           dS.HS, dI.HS, dR.HS),
         P = P)
  })
}

################################## Objective function ##################################

objective_function_global <- function(params, site_data_list) {
  
  b.LS <- params[1]; g.LS <- params[2]
  b.MS <- params[3]; g.MS <- params[4]
  b.HS <- params[5]; g.HS <- params[6]
  q    <- params[7]
  
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
                     K.LS = K.LS,   K.MS = K.MS,   K.HS = K.HS,
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
    
    n.LS <- min(length(sim.LS.rem), length(sd$LS.remtiss))
    n.MS <- min(length(sim.MS.rem), length(sd$MS.remtiss))
    n.HS <- min(length(sim.HS.rem), length(sd$HS.remtiss))
    
    site_residual <- (
      sum(abs(sim.LS.rem[1:n.LS] - sd$LS.remtiss[1:n.LS])) +
        sum(abs(sim.MS.rem[1:n.MS] - sd$MS.remtiss[1:n.MS])) +
        sum(abs(sim.HS.rem[1:n.HS] - sd$HS.remtiss[1:n.HS]))
    ) / sd$N.site
    
    total_residual <- total_residual + site_residual
  }
  
  total_residual
}

################################## Optimization ##################################

control <- list(itermax = 300)

result_global <- DEoptim(
  fn             = objective_function_global,
  lower          = c(0,   0,   0,   0,   0,   0,   0.5),
  upper          = c(4,   4,   4,   4,   4,   4,   1.0),
  site_data_list = site_data_list,
  control        = control
)

b.LS.global <- result_global$optim$bestmem[1]
g.LS.global <- result_global$optim$bestmem[2]
b.MS.global <- result_global$optim$bestmem[3]
g.MS.global <- result_global$optim$bestmem[4]
b.HS.global <- result_global$optim$bestmem[5]
g.HS.global <- result_global$optim$bestmem[6]
q.global    <- result_global$optim$bestmem[7]

cat("\nGlobal fitted parameters:\n")
cat(sprintf("  beta  LS: %.6f   gamma LS: %.6f   R0 LS: %.4f\n",
            b.LS.global, g.LS.global, b.LS.global / g.LS.global))
cat(sprintf("  beta  MS: %.6f   gamma MS: %.6f   R0 MS: %.4f\n",
            b.MS.global, g.MS.global, b.MS.global / g.MS.global))
cat(sprintf("  beta  HS: %.6f   gamma HS: %.6f   R0 HS: %.4f\n",
            b.HS.global, g.HS.global, b.HS.global / g.HS.global))
cat(sprintf("  q: %.4f\n", q.global))

################################## Diagnostics ##################################

SUS_COLORS <- c(Low = "#1E90FF", Moderate = "#FFD700", High = "#FF1493")
site_labels <- c(near = "Nearshore", mid = "Midchannel", off = "Offshore")

global_sims <- lapply(site_data_list, function(sd) {
  data.frame(ode(
    y     = sd$init,
    times = sd$days.model,
    func  = SIR.multi,
    parms = list(b.LS = b.LS.global, g.LS = g.LS.global,
                 b.MS = b.MS.global, g.MS = g.MS.global,
                 b.HS = b.HS.global, g.HS = g.HS.global,
                 N.LS = sd$N.LS, N.MS = sd$N.MS, N.HS = sd$N.HS, # Note - would change this to try within-group freq.-dep. vs community freq.-dep.
                 K.LS = K.LS,   K.MS = K.MS,   K.HS = K.HS,
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
         "Global multi-host SIR fit (q = %.3f)\nβ_LS=%.4f γ_LS=%.4f | β_MS=%.4f γ_MS=%.4f | β_HS=%.4f γ_HS=%.4f",
         q.global,
         b.LS.global, g.LS.global,
         b.MS.global, g.MS.global,
         b.HS.global, g.HS.global
       )) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p_global_fit)