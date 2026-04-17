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
  N.LMS  <- susceptible_ref %>% filter(Site == s, Category %in% c("Low", "Moderate")) %>% pull(tissue_ref) %>% sum()
  N.HS   <- susceptible_ref %>% filter(Site == s, Category == "High") %>% pull(tissue_ref)
  
  pull_tiss <- function(cat, comp) {
    obs.model %>%
      filter(Site == s, Category == cat, Compartment == comp) %>%
      slice(head(row_number(), n() - DHW.modifier)) %>%
      pull(tissue)
  }
  
  HS.inf.raw <- pull_tiss("High", "Infected")
  first_idx  <- which(!is.na(HS.inf.raw) & HS.inf.raw != 0)[1]
  trim       <- function(x) x[first_idx:length(x)]
  
  LS.remtiss <- trim(pull_tiss("Low",      "Dead"))
  MS.remtiss <- trim(pull_tiss("Moderate", "Dead"))
  LS.inftiss <- trim(pull_tiss("Low",      "Infected"))
  MS.inftiss <- trim(pull_tiss("Moderate", "Infected"))
  
  list(
    site        = s,
    days.obs    = days.obs,
    days.model  = days.model,
    N.site      = N.site,
    N.LMS       = N.LMS,
    N.HS        = N.HS,
    LMS.remtiss = LS.remtiss + MS.remtiss,
    HS.remtiss  = trim(pull_tiss("High", "Dead")),
    LMS.inftiss = LS.inftiss + MS.inftiss,
    HS.inftiss  = trim(pull_tiss("High", "Infected")),
    init = c(
      S.LMS = N.LMS, I.LMS = 0, R.LMS = 0,
      S.HS  = N.HS - HS.inf.raw[first_idx], I.HS = HS.inf.raw[first_idx], R.HS = 0
    )
  )
})

# K fixed at largest observed total community size across all sites (always uses full list)
K <- max(sapply(site_data_list, function(sd) sd$N.site))

cat(sprintf("K: %.4f\n", K))

################################## Toggles ##################################

# --- Site selection ---
# Choose any subset of sites to include in optimization.
# Diagnostics will always show all three for visual comparison.
# Options: any combination of "near", "mid", "off"
# fit_sites <- c("near", "mid", "off")  # <-- TOGGLE HERE
fit_sites <- c("near", "off")  # <-- TOGGLE HERE

site_data_list_fit <- Filter(function(sd) sd$site %in% fit_sites, site_data_list)
cat(sprintf("Fitting to site(s): %s\n", paste(fit_sites, collapse = ", ")))

# --- q toggle ---
# q controls the scaling of the contact rate with population size,
# following Smith et al. (2009):
#   q = 1  -> pure frequency dependence  (contact rate independent of N0)
#   q = 0  -> pure density dependence    (contact rate proportional to N0)
#   q > 1  -> super-frequency dependence (contact rate declines with N0)
#   q < 0  -> super-density dependence   (contact rate accelerates with N0)
fit_q       <- FALSE  # <-- TOGGLE HERE
joint_q     <- FALSE   # TRUE = one shared q; FALSE = separate q.LMS and q.HS; ignored when fit_q = FALSE
q_val_fixed <- 1.0    # used only when fit_q = FALSE

# --- Gamma toggles ---
fix_g.LMS   <- FALSE; g.LMS.fixed <- 0.14
fix_g.HS    <- FALSE; g.HS.fixed  <- 2.49

# --- Beta upper bound toggle ---
cap_beta    <- TRUE; beta_upper <- 0.1

# --- Residual normalization toggle ---
# TRUE  = normalize by max observed R (dead tissue) at each site
# FALSE = normalize by N.site
norm_by_max_R <- FALSE

  # NOTES - 16 April 2026 [when weighting each site by N.site...could consider weighting by R.site. also, capping betas at 0.1; minimizing absolute total I + R]
  # - fixing q = 1 with community freq-dep. is interesting. it converges on nearshore, and midchannel/offshore take off much too slowly (but realistically)
  # - freeing q (but LMS == HS) with community freq-dep. does basically the same as above.
  # - freeing q (with LMS != HS) with community freq-dep. (BUT only considering Nearshore & Offshore) gets really cool. similar as above, and converges on nearshore and q's -> 1.0...but a bit closer to working. setting q's = 1.0 directly does not work.
  # - freeing q (with LMS != HS) with within-group freq-dep. (BUT only considering nearshore & offshore) is the best result so far, I think. Timing is right, and magnitude isn't bad at all...could consider normalizing by max R instead of max N next ? or making q in LMS == HS? or fixing q = 1
  # - fixing q = 1 with within-group freq-dep.  (only considering nearshore * offshore) is the result closest to what we already showed in paper...which makes some sense, but very weird that gamma in LMS is twice that of HS. kind of a quirk
  # - freeing q (with LMS == HS) with within-group freq-dep. (only considering nearshore * offshore) is very similar to just above. same thing with strangely much higher LMS gamma than HS
# - so do we need to consider within-group freq-dep. with 2 groups ?
#     - and/or let q.LMS != q.HS ?
#     - optimize globally across 2 instead of 3 sites ?
# - may consider using 'initialpop' term for DEoptim to help with convergence

# NOTES - 16 April 2026 [when weighting by N.group or max R.group. also, capping betas at 0.1]
# - freeing q (with LMS != HS) with within-group freq-dep. (BUT only considering nearshore & offshore); weighting by N.group. result is VERY good for nearshore and offshore HS, but predicted removal in LMS is literally almost zero. so, now attempting same thing but weighting by max R.group
# - freeing q (with LMS != HS) with within-group freq-dep. (BUT only considering nearshore & offshore); weighting by max R.group. result is okay. no "curve-switching", but generally timing and curves look real good; fitted parameters make sense. 
# - freeing q (with LMS != HS) with community freq-dep. (BUT only considering Nearshore & Offshore); weighting by max R.group. result is craperama
# - freeing q (with LMS != HS) with community freq-dep. (BUT only considering Nearshore & Offshore); weighting by N.group. result is SUPER COOL and new best result so far
# - fixing q with community freq-dep. and weighting by N.group (BUT only considering Nearshore & Offshore) is GREAT too; if anything it converges on offshore here, but generally shows why it's important to have nonlinear contact rate scaling
# - fixing q with within-group freq-dep. and weighting by N.group (BUT only considering Nearshore & Offshore) is ... actually about what we'd expect! "invariant", like has been shown before.


################################## Model ##################################

SIR.multi <- function(t, y, p) {
  S.LMS <- y[1]; I.LMS <- y[2]; R.LMS <- y[3]
  S.HS  <- y[4]; I.HS  <- y[5]; R.HS  <- y[6]
  
  with(as.list(p), {
    
    N0 <- N.LMS + N.HS
    P  <- I.LMS + I.HS
    
    # # NOTE - community freq-dep.
    # contact.LMS <- (K^q.LMS) * S.LMS * P / (N0^q.LMS)
    # contact.HS  <- (K^q.HS)  * S.HS  * P / (N0^q.HS)
    
    # NOTE - within-group freq-dep.
    contact.LMS <- (K^q.LMS) * S.LMS * P / (N.LMS^q.LMS)
    contact.HS <- (K^q.HS) * S.HS * P / (N.HS^q.HS)
    
    dS.LMS <- -b.LMS * contact.LMS
    dI.LMS <-  b.LMS * contact.LMS - g.LMS * I.LMS
    dR.LMS <-  g.LMS * I.LMS
    
    dS.HS <- -b.HS * contact.HS
    dI.HS <-  b.HS * contact.HS - g.HS * I.HS
    dR.HS <-  g.HS * I.HS
    
    list(c(dS.LMS, dI.LMS, dR.LMS,
           dS.HS,  dI.HS,  dR.HS),
         P = P)
  })
}

################################## Objective function ##################################

objective_function_global <- function(params, site_data_list,
                                      q_fixed = NULL, joint_q = TRUE,
                                      g.LMS.fixed = NULL, g.HS.fixed = NULL,
                                      norm_by_max_R = FALSE) {
  
  idx <- 1
  b.LMS <- params[idx]; idx <- idx + 1
  g.LMS <- if (!is.null(g.LMS.fixed)) g.LMS.fixed else { v <- params[idx]; idx <- idx + 1; v }
  b.HS  <- params[idx]; idx <- idx + 1
  g.HS  <- if (!is.null(g.HS.fixed))  g.HS.fixed  else { v <- params[idx]; idx <- idx + 1; v }
  
  if (!is.null(q_fixed)) {
    q.LMS <- q_fixed; q.HS <- q_fixed
  } else if (joint_q) {
    q.LMS <- params[idx]; q.HS <- params[idx]; idx <- idx + 1
  } else {
    q.LMS <- params[idx]; idx <- idx + 1; q.HS <- params[idx]; idx <- idx + 1
  }
  
  total_residual <- 0
  
  for (sd in site_data_list) {
    
    SIR.out <- tryCatch({
      data.frame(ode(
        y     = sd$init,
        times = sd$days.model,
        func  = SIR.multi,
        parms = list(b.LMS = b.LMS, g.LMS = g.LMS,
                     b.HS  = b.HS,  g.HS  = g.HS,
                     N.LMS = sd$N.LMS, N.HS = sd$N.HS,
                     K     = K,
                     q.LMS = q.LMS, q.HS = q.HS)
      ))
    }, error = function(e) NULL)
    
    if (is.null(SIR.out)) {
      total_residual <- total_residual + 1e9
      next
    }
    
    rows <- which(SIR.out$time %in% sd$days.obs)
    
    sim.LMS.rem <- SIR.out[rows, "R.LMS"]
    sim.HS.rem  <- SIR.out[rows, "R.HS"]
    sim.LMS.inf <- SIR.out[rows, "I.LMS"]
    sim.HS.inf  <- SIR.out[rows, "I.HS"]
    
    n.LMS   <- min(length(sim.LMS.rem), length(sd$LMS.remtiss))
    n.HS    <- min(length(sim.HS.rem),  length(sd$HS.remtiss))
    n.LMS.i <- min(length(sim.LMS.inf), length(sd$LMS.inftiss))
    n.HS.i  <- min(length(sim.HS.inf),  length(sd$HS.inftiss))
    
    norm <- if (norm_by_max_R) {
      max(c(sd$LMS.remtiss, sd$HS.remtiss), na.rm = TRUE)
    } else {
      sd$N.site
    }
    
    # # NOTE - normalize by site-level N or max R
    # site_residual <- (
    #   sum(abs(sim.LMS.rem[1:n.LMS]   - sd$LMS.remtiss[1:n.LMS]))   +
    #     sum(abs(sim.HS.rem[1:n.HS]     - sd$HS.remtiss[1:n.HS]))     +
    #     sum(abs(sim.LMS.inf[1:n.LMS.i] - sd$LMS.inftiss[1:n.LMS.i])) +
    #     sum(abs(sim.HS.inf[1:n.HS.i]   - sd$HS.inftiss[1:n.HS.i]))
    # ) / norm
    
    # NOTE - normalize by group-level N
    site_residual <-
      sum(abs(sim.LMS.rem - sd$LMS.remtiss)) / sd$N.LMS +
      sum(abs(sim.HS.rem - sd$HS.remtiss)) / sd$N.HS +
      sum(abs(sim.LMS.inf - sd$LMS.inftiss)) / sd$N.LMS +
      sum(abs(sim.HS.inf - sd$HS.inftiss)) / sd$N.HS
    
    # # NOTE - normalize by group-level max R
    # site_residual <-
    #   sum(abs(sim.LMS.rem - sd$LMS.remtiss)) / max(sd$LMS.remtiss) +
    #   sum(abs(sim.HS.rem  - sd$HS.remtiss))  / max(sd$HS.remtiss)  +
    #   sum(abs(sim.LMS.inf - sd$LMS.inftiss)) / max(sd$LMS.remtiss) +
    #   sum(abs(sim.HS.inf  - sd$HS.inftiss))  / max(sd$HS.remtiss)
    
    total_residual <- total_residual + site_residual
  }
  
  total_residual
}

################################## Optimization ##################################

control <- list(itermax = 300)

b_upper    <- if (cap_beta) beta_upper else 4
g.LMS.arg  <- if (fix_g.LMS) g.LMS.fixed else NULL
g.HS.arg   <- if (fix_g.HS)  g.HS.fixed  else NULL

lower_base <- c(0)
if (!fix_g.LMS) lower_base <- c(lower_base, 0)
lower_base <- c(lower_base, 0)
if (!fix_g.HS)  lower_base <- c(lower_base, 0)

upper_base <- c(b_upper)
if (!fix_g.LMS) upper_base <- c(upper_base, 4)
upper_base <- c(upper_base, b_upper)
if (!fix_g.HS)  upper_base <- c(upper_base, 4)

if (fit_q) {
  if (joint_q) {
    cat(sprintf("Fitting beta, %sAND shared q...\n",
                if (!fix_g.LMS || !fix_g.HS) "gamma, " else ""))
    result_global <- DEoptim(
      fn             = objective_function_global,
      lower          = c(lower_base, 0.5),
      upper          = c(upper_base, 1.0),
      site_data_list = site_data_list_fit,
      q_fixed        = NULL, joint_q = TRUE,
      g.LMS.fixed    = g.LMS.arg, g.HS.fixed = g.HS.arg,
      norm_by_max_R  = norm_by_max_R,
      control        = control
    )
    q.LMS.global <- tail(result_global$optim$bestmem, 1)
    q.HS.global  <- q.LMS.global
  } else {
    cat(sprintf("Fitting beta, %sAND separate q.LMS / q.HS...\n",
                if (!fix_g.LMS || !fix_g.HS) "gamma, " else ""))
    result_global <- DEoptim(
      fn             = objective_function_global,
      lower          = c(lower_base, 0.5, 0.5),
      upper          = c(upper_base, 1.0, 1.0),
      site_data_list = site_data_list_fit,
      q_fixed        = NULL, joint_q = FALSE,
      g.LMS.fixed    = g.LMS.arg, g.HS.fixed = g.HS.arg,
      norm_by_max_R  = norm_by_max_R,
      control        = control
    )
    q.LMS.global <- tail(result_global$optim$bestmem, 2)[1]
    q.HS.global  <- tail(result_global$optim$bestmem, 1)
  }
} else {
  cat(sprintf("Fitting beta, %swith q fixed at %.2f...\n",
              if (!fix_g.LMS || !fix_g.HS) "gamma, " else "", q_val_fixed))
  result_global <- DEoptim(
    fn             = objective_function_global,
    lower          = lower_base,
    upper          = upper_base,
    site_data_list = site_data_list_fit,
    q_fixed        = q_val_fixed, joint_q = TRUE,
    g.LMS.fixed    = g.LMS.arg, g.HS.fixed = g.HS.arg,
    norm_by_max_R  = norm_by_max_R,
    control        = control
  )
  q.LMS.global <- q_val_fixed
  q.HS.global  <- q_val_fixed
}

# Extract fitted parameters in order
res <- result_global$optim$bestmem
idx <- 1
b.LMS.global <- res[idx]; idx <- idx + 1
g.LMS.global <- if (fix_g.LMS) g.LMS.fixed else { v <- res[idx]; idx <- idx + 1; v }
b.HS.global  <- res[idx]; idx <- idx + 1
g.HS.global  <- if (fix_g.HS)  g.HS.fixed  else { v <- res[idx]; idx <- idx + 1; v }
if (fit_q) {
  q.LMS.global <- res[idx]; idx <- idx + 1
  q.HS.global  <- if (joint_q) q.LMS.global else { v <- res[idx]; idx <- idx + 1; v }
}

cat("\nGlobal fitted parameters:\n")
cat(sprintf("  beta  LMS: %.6f   gamma LMS: %.6f   R0 LMS: %.4f\n",
            b.LMS.global, g.LMS.global, b.LMS.global / g.LMS.global))
cat(sprintf("  beta  HS:  %.6f   gamma HS:  %.6f   R0 HS:  %.4f\n",
            b.HS.global, g.HS.global, b.HS.global / g.HS.global))
cat(sprintf("  q.LMS: %.4f   q.HS: %.4f  (%s)\n", q.LMS.global, q.HS.global,
            if (!fit_q) "fixed" else if (joint_q) "fitted jointly" else "fitted separately"))
cat(sprintf("  Fitted to: %s\n", paste(fit_sites, collapse = ", ")))

################################## Diagnostics ##################################

# Diagnostics always simulate all three sites for visual comparison,
# regardless of which sites were used in fitting.

SUS_COLORS  <- c("Low/Moderate" = "#1E90FF", High = "#FF1493")
site_labels <- c(near = "Nearshore", mid = "Midchannel", off = "Offshore")

global_sims <- lapply(site_data_list, function(sd) {
  data.frame(ode(
    y     = sd$init,
    times = sd$days.model,
    func  = SIR.multi,
    parms = list(b.LMS = b.LMS.global, g.LMS = g.LMS.global,
                 b.HS  = b.HS.global,  g.HS  = g.HS.global,
                 N.LMS = sd$N.LMS, N.HS = sd$N.HS,
                 K     = K,
                 q.LMS = q.LMS.global, q.HS = q.HS.global)
  )) %>%
    select(-any_of("P")) %>%
    mutate(
      site      = sd$site,
      fitted    = sd$site %in% fit_sites  # flag whether this site was in the fit
    )
}) %>%
  bind_rows() %>%
  pivot_longer(cols = -c(time, site, fitted),
               names_pattern = "(.*)\\.(.*)",
               names_to = c("Compartment", "Susceptibility")) %>%
  mutate(
    Compartment    = case_match(Compartment,
                                "S" ~ "Susceptible", "I" ~ "Infected", "R" ~ "Dead"),
    Susceptibility = case_match(Susceptibility,
                                "LMS" ~ "Low/Moderate", "HS" ~ "High"),
    site_label     = site_labels[site]
  )

obs_plot <- obs.model %>%
  filter(Category %in% c("Low", "Moderate", "High")) %>%
  mutate(Category = if_else(Category %in% c("Low", "Moderate"), "Low/Moderate", Category)) %>%
  group_by(Site, Category, Compartment, days.inf.site) %>%
  summarize(tissue = sum(tissue, na.rm = TRUE), .groups = "drop") %>%
  rename(Susceptibility = Category) %>%
  mutate(
    site_label = site_labels[Site],
    fitted     = Site %in% fit_sites
  )

# Build subtitle flagging which sites were held out
holdout_sites <- setdiff(c("near", "mid", "off"), fit_sites)
fit_label <- if (length(holdout_sites) == 0) {
  "all sites fitted"
} else {
  sprintf("fitted: %s | held out: %s",
          paste(site_labels[fit_sites],    collapse = ", "),
          paste(site_labels[holdout_sites], collapse = ", "))
}

p_global_fit <- ggplot(global_sims %>% filter(Compartment == "Dead"),
                       aes(time, value, color = Susceptibility,
                           linetype = fitted)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_plot %>% filter(Compartment == "Dead"),
             aes(days.inf.site, tissue, color = Susceptibility,
                 shape = fitted),
             size = 1.5) +
  scale_color_manual(values = SUS_COLORS) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                        labels = c("TRUE" = "Fitted", "FALSE" = "Held out"),
                        name   = NULL) +
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 1),
                     labels = c("TRUE" = "Fitted", "FALSE" = "Held out"),
                     name   = NULL) +
  facet_wrap(~ site_label, ncol = 1, scales = "free_y") +
  labs(x     = "Day of outbreak",
       y     = "Surface area of tissue (m²)",
       title = sprintf(
         "Global multi-host SIR fit (q.LMS=%.3f q.HS=%.3f, %s)\nβ_LMS=%.4f γ_LMS=%.4f | β_HS=%.4f γ_HS=%.4f",
         q.LMS.global, q.HS.global,
         if (!fit_q) "fixed" else if (joint_q) "fitted jointly" else "fitted separately",
         b.LMS.global, g.LMS.global,
         b.HS.global,  g.HS.global
       ),
       subtitle = fit_label) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p_global_fit)