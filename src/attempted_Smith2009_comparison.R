# .rs.restartR(clean = TRUE)
rm(list = ls())

library(here)
library(tidyverse)
library(DEoptim)
library(deSolve)
library(patchwork)

load(here("output", "plots_obs_workspace.RData"))

# ── Constants ────────────────────────────────────────────────────────────────
lambda.modifier <- 1.0
offset          <- 1 - 1 / (1 + exp(-lambda.modifier))
alpha_val       <- 0.13
k_val           <- 3
DHW.modifier    <- 0
site.fit        <- "near"

# Alpha-exponent model controls
q_val        <- 1      # transmission scaling exponent: 0 = density-dep, 1 = frequency-dep
use_freq_dep <- FALSE  # legacy flag, no longer used in SIR_alphaexp

obs.model <- obs %>%
  mutate(Site = case_when(
    Site == "Offshore"   ~ "off",
    Site == "Midchannel" ~ "mid",
    Site == "Nearshore"  ~ "near",
    TRUE ~ Site
  ))

# ── Shared data prep ─────────────────────────────────────────────────────────
days_sites <- summary %>%
  group_by(site) %>%
  summarize(
    days     = list(days.inf.site[1:length(days.inf.site)]),
    days.obs = list({
      d <- na.omit(days.inf.site)
      d[which(!is.na(d))[1]:length(d)]
    }),
    days.model = list(seq(min(na.omit(days.inf.site)),
                          max(na.omit(days.inf.site)), by = 1))
  )

get_site_data <- function(s) {
  days       <- days_sites %>% filter(site == s) %>% pull(days)      %>% unlist()
  days.obs   <- days_sites %>% filter(site == s) %>% pull(days.obs) %>% unlist()
  days.model <- days_sites %>% filter(site == s) %>% pull(days.model) %>% unlist()
  
  N.site     <- susceptible_ref %>% filter(Site == s) %>% slice(1) %>% pull(N.site)
  cover.site <- susceptible_ref %>% filter(Site == s) %>% slice(1) %>% pull(cover.site)
  
  first_idx <- which(!is.na(days))[1]
  
  inftiss <- obs.model %>%
    filter(Site == s, Category == "Total", Compartment == "Infected") %>%
    pull(tissue) %>% .[first_idx:length(.)]
  
  remtiss <- obs.model %>%
    filter(Site == s, Category == "Total", Compartment == "Dead") %>%
    pull(tissue) %>% .[first_idx:length(.)]
  
  list(days = days, days.obs = days.obs, days.model = days.model,
       N = N.site, cover = cover.site,
       inftiss = inftiss, remtiss = remtiss,
       I0 = inftiss[1], S0 = N.site - inftiss[1], R0_init = 0)
}

sd_near <- get_site_data("near")
sd_mid  <- get_site_data("mid")
sd_off  <- get_site_data("off")

# ── K: rescaling constant (Smith et al. 2009) ────────────────────────────────
# Set to maximum N across all sites so that beta units are invariant to q.
# (N/K) replaces raw N in the denominator, keeping beta in (tissue * time)^-1
# regardless of the exponent chosen.
K_scale <- max(sd_near$N, sd_mid$N, sd_off$N)
cat(sprintf("K_scale set to %.4f (max N across sites)\n", K_scale))

# ── ODE: cover modifier (unchanged) ─────────────────────────────────────────
SIR_cover <- function(t, y, p) {
  S <- y[1]; I <- y[2]; R <- y[3]
  with(as.list(p), {
    tm <- (1 - alpha_val) + alpha_val * ((1 - exp(-k_val * C)) / (1 - exp(-k_val)))
    dS <- -b * S * I / N * tm
    dI <-  b * S * I / N * tm - g * I
    dR <-  g * I
    list(c(dS, dI, dR))
  })
}

# # ── ODE: alpha-exponent (Smith et al. framework) ─────────────────────────────
# # Infection rate = b * S * I / N^q  (Smith et al. Eq. 4)
# #   q = 0  ->  b * S * I            (density-dependent, mass action)
# #   q = 1  ->  b * S * I / N        (frequency-dependent)
# #   0 < q < 1  -> saturating contact rate with density
# #
# # K_scale does NOT enter the ODE. It is used only after fitting to back-transform
# # beta to a common scale: beta_rescaled = beta * K^q  (see Smith et al. Eq. 4).
# # This makes beta comparable across runs that use different q values.
# SIR_alphaexp <- function(t, y, p) {
#   S <- y[1]; I <- y[2]; R <- y[3]
#   with(as.list(p), {
#     I_safe <- max(I, 1e-10)
#     N_safe <- max(N, 1e-10)
#     force  <- b * S * I_safe / N_safe^q
#     dS <- -force
#     dI <-  force - g * I
#     dR <-  g * I
#     list(c(dS, dI, dR))
#   })
# }

SIR_alphaexp <- function(t, y, p) {
  S <- y[1]; I <- y[2]; R <- y[3]
  with(as.list(p), {
    contact_rate <- (K^q_val) * S * I / (N^q_val)
    dS <- -b * contact_rate
    dI <-  b * contact_rate - g * I
    dR <-  g * I
    list(c(dS, dI, dR))
  })
}

# ── Generic fitter ────────────────────────────────────────────────────────────
fit_sir <- function(sir_fn, sir_type, sd) {
  
  obj <- function(params) {
    # cat("N =", sd$N, "days.obs length =", length(sd$days.obs), "remtiss length =", length(sd$remtiss), "\n")
    b <- params[1]
    g <- params[2]

    parms <- if (sir_type == "cover") {
      c(b = b, g = g, N = sd$N, C = sd$cover, l = lambda.modifier)
    } else {
      c(b = b, g = g, N = sd$N, K = K_scale, q = q_val)
    }

    out <- tryCatch(
      data.frame(ode(c(S = sd$S0, I = sd$I0, R = sd$R0_init),
                     sd$days.model, sir_fn, parms)),
      error = function(e) return(NULL)
    )
    if (is.null(out)) return(1e9)

    # sim.R <- out[which(out$time %in% sd$days.obs), "R"]
    # sim.I <- out[which(out$time %in% sd$days.obs), "I"]
    # obs.R <- sd$remtiss
    # obs.I <- sd$inftiss
    # w_R <- 2   # weight on R (dead tissue) — adjust as needed
    # w_I <- 1   # weight on I (infected tissue)
    # # w_R * sum(abs(sim.R - obs.R)) + w_I * sum(abs(sim.I - obs.I))

    obs_idx <- which(sd$days.obs %in% out$time)
    sim.R <- out[which(out$time %in% sd$days.obs), "R"]
    obs.R <- sd$remtiss[obs_idx]
    sim.I <- out[which(out$time %in% sd$days.obs), "I"]
    obs.I <- sd$inftiss[obs_idx]


    sum(abs(sim.R - obs.R))
  }
  
  # # version with weighted I fitting
  # obj <- function(params) {
  #   b <- params[1]
  #   g <- params[2]
  #   
  #   parms <- if (sir_type == "cover") {
  #     c(b = b, g = g, N = sd$N, C = sd$cover, l = lambda.modifier)
  #   } else {
  #     c(b = b, g = g, N = sd$N, K = K_scale, q = q_val)
  #   }
  #   
  #   out <- tryCatch(
  #     data.frame(ode(c(S = sd$S0, I = sd$I0, R = sd$R0_init),
  #                    sd$days.model, sir_fn, parms)),
  #     error = function(e) return(NULL)
  #   )
  #   if (is.null(out)) return(1e9)
  #   
  #   sim.R <- out[which(out$time %in% sd$days.obs), "R"]
  #   sim.I <- out[which(out$time %in% sd$days.obs), "I"]
  #   obs.R <- sd$remtiss
  #   obs.I <- sd$inftiss
  #   
  #   err.R <- sum(abs(sim.R - obs.R)) / max(obs.R)
  #   err.I <- sum(abs(sim.I - obs.I)) / max(obs.I)
  #   
  #   2 * err.R + 3 * err.I
  # }
  
  res <- DEoptim(fn = obj, lower = c(0, 0), upper = c(4, 4),
                 control = list(itermax = 200, trace = FALSE))
  
  p <- res$optim$bestmem
  names(p) <- c("beta", "gamma")
  p
}

# ── Fit both models for Nearshore ─────────────────────────────────────────────
cat("Fitting cover-modifier model...\n")
p_cover <- fit_sir(SIR_cover, sir_type = "cover", sd = sd_near)
p_cover["C"] <- sd_near$cover

cat("Fitting alpha-exponent model...\n")
p_alpha <- fit_sir(SIR_alphaexp, sir_type = "alpha", sd = sd_near)

cat("\nCover-modifier params:\n"); print(round(p_cover, 4))
cat("\nAlpha-exponent params (q =", q_val, "):\n")
print(round(p_alpha, 4))
# The estimated beta is already Smith's rescaled beta (beta = beta_qD / K^q),
# because K^q is embedded inside the ODE via (N/K)^q in the denominator.
# To recover the raw contact rate parameter beta_qD (units depend on q):
beta_qD <- unname(p_alpha["beta"]) * K_scale^q_val
cat(sprintf("  Estimated beta (rescaled, comparable across q): %.6f\n",
            unname(p_alpha["beta"])))
cat(sprintf("  beta_qD = beta * K^q (raw, units depend on q):  %.6f\n", beta_qD))
cat(sprintf("  K_scale = %.4f, q = %.2f\n", K_scale, q_val))

# ── Simulate both — Nearshore fit ─────────────────────────────────────────────
sim_cover <- data.frame(ode(
  c(S = sd_near$S0, I = sd_near$I0, R = sd_near$R0_init),
  sd_near$days.model, SIR_cover,
  c(b = unname(p_cover["beta"]), g = unname(p_cover["gamma"]),
    N = sd_near$N, C = sd_near$cover, l = lambda.modifier)
)) %>% mutate(method = "Cover modifier")

sim_alpha <- data.frame(ode(
  c(S = sd_near$S0, I = sd_near$I0, R = sd_near$R0_init),
  sd_near$days.model, SIR_alphaexp,
  c(b = unname(p_alpha["beta"]), g = unname(p_alpha["gamma"]),
    N = sd_near$N, K = K_scale, q = q_val)
)) %>% mutate(method = "Alpha exponent")

sims <- bind_rows(sim_cover, sim_alpha) %>%
  pivot_longer(c(S, I, R), names_to = "Compartment", values_to = "tissue") %>%
  mutate(Compartment = case_match(Compartment,
                                  "S" ~ "Susceptible", "I" ~ "Infected", "R" ~ "Dead"))

obs_near <- obs.model %>%
  filter(Site == site.fit, Category == "Total")

METHOD_COLORS    <- c("Cover modifier" = "steelblue", "Alpha exponent" = "tomato")
METHOD_LINETYPES <- c("Cover modifier" = "solid",     "Alpha exponent" = "dashed")
site_colors      <- c(Nearshore = "orange", Midchannel = "purple", Offshore = "magenta")

p_fit <- ggplot(sims, aes(time, tissue,
                          color    = method,
                          linetype = method,
                          group    = interaction(method, Compartment))) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_near,
             aes(days.inf.site, tissue, shape = Compartment),
             color = "black", size = 1.5, inherit.aes = FALSE) +
  scale_color_manual(values = METHOD_COLORS, name = NULL) +
  scale_linetype_manual(values = METHOD_LINETYPES, name = NULL) +
  scale_shape_manual(values = c(Susceptible = 16, Infected = 17, Dead = 15), name = NULL) +
  labs(x = "Day of outbreak", y = "Surface area of tissue (m²)",
       title = "Nearshore — fitted comparison") +
  theme_classic() +
  theme(legend.position = "bottom")

# ── Behavior panel A: cover modifier ─────────────────────────────────────────
CC     <- seq(0.001, 1, 0.001)
a_vals <- sort(c(seq(0, 1, 0.1), 0.13))

beh_cover <- expand.grid(CC = CC, alpha = a_vals) %>%
  mutate(scalar       = (1 - alpha) + alpha * ((1 - exp(-k_val * CC)) / (1 - exp(-k_val))),
         is_highlight = alpha == 0.13)

site_pts_cover <- data.frame(
  CC_pct = c(0.247, 0.0215) * 100,
  scalar = c((1-0.13) + 0.13*((1-exp(-k_val*0.247))/(1-exp(-k_val))),
             (1-0.13) + 0.13*((1-exp(-k_val*0.0215))/(1-exp(-k_val)))),
  site = c("Nearshore", "Offshore")
)

pA <- ggplot() +
  geom_line(data = filter(beh_cover, !is_highlight),
            aes(CC * 100, scalar, color = alpha, group = alpha), linewidth = 0.5) +
  geom_line(data = filter(beh_cover, is_highlight),
            aes(CC * 100, scalar), color = "blue", linewidth = 0.9) +
  geom_point(data = site_pts_cover,
             aes(CC_pct, scalar, fill = site),
             shape = 23, size = 2.5, color = "black") +
  scale_color_viridis_c(name = expression(alpha[val]), option = "inferno") +
  scale_fill_manual(values = c(Nearshore = "orange", Offshore = "magenta"), name = NULL) +
  scale_x_continuous("Coral cover (%)", limits = c(0, 100), expand = c(0,0)) +
  scale_y_continuous("Contact rate scalar", limits = c(0, 1), expand = c(0,0)) +
  theme_classic() +
  ggtitle("A: Cover modifier behavior")

# ── Behavior panel B: alpha-exponent (mirrors Smith et al. Fig. 2A) ──────────
# X-axis: population size N (0 to K_scale)
# Y-axis: scaled contact rate = (N/K)^(1-q)  (Smith Eq. 2, kappa = 1)
# At N = K, every curve evaluates to (K/K)^(1-q) = 1 regardless of q,
# producing exact convergence at the top right — this is a mathematical
# property of the equation, not a post-hoc normalization.
#   q = 1  -> flat line at 1     (frequency-dependent)
#   q = 0  -> linear rise        (density-dependent)
#   0 < q < 1 -> saturating      (intermediate)
#   q < 0  -> accelerating       (super-linear)

N_seq  <- seq(0.01, K_scale, length.out = 500)
q_vals <- sort(c(0, 1, -1, -5, q_val))

beh_alpha <- expand.grid(N = N_seq, q = q_vals) %>%
  mutate(
    contact      = (N / K_scale)^(1 - q),
    is_highlight = abs(q - q_val) < 1e-6,
    q_label      = factor(q)
  )

site_pts_B <- data.frame(
  N    = c(sd_near$N, sd_mid$N, sd_off$N),
  site = c("Nearshore", "Midchannel", "Offshore")
)

pB <- ggplot() +
  geom_vline(data = site_pts_B,
             aes(xintercept = N),
             color    = site_colors[site_pts_B$site],
             linetype = "dotted", linewidth = 0.6) +
  geom_line(data = filter(beh_alpha, !is_highlight),
            aes(N, contact, group = q_label), color = "gray40", linewidth = 0.5) +
  geom_line(data = filter(beh_alpha, is_highlight),
            aes(N, contact), color = "black", linewidth = 1.1) +
  geom_text(data = filter(beh_alpha, !is_highlight) %>%
              group_by(q) %>% slice_max(N, n = 1),
            aes(N, contact, label = paste0("q=", q)),
            hjust = -0.1, size = 2.8, color = "gray30") +
  geom_text(data = filter(beh_alpha, is_highlight) %>% slice_max(N, n = 1),
            aes(N, contact, label = paste0("q=", q_val, " (active)")),
            hjust = -0.1, size = 2.8, color = "black") +
  scale_x_continuous("Population size (N)",
                     limits = c(0, K_scale * 1.25), expand = c(0, 0)) +
  scale_y_continuous("Scaled contact rate  [(N/K)^(1-q)]",
                     limits = c(0, NA), expand = c(0, 0)) +
  theme_classic() +
  ggtitle("B: Alpha-exponent behavior")

p_behavior <- pA + pB + plot_layout(ncol = 2, guides = "keep")

print(p_fit)
print(p_behavior)

# ── Projections ───────────────────────────────────────────────────────────────
make_proj_cover <- function(sd, label) {
  data.frame(ode(
    c(S = sd$S0, I = sd$I0, R = sd$R0_init),
    sd$days.model, SIR_cover,
    c(b = unname(p_cover["beta"]), g = unname(p_cover["gamma"]),
      N = sd$N, C = sd$cover, l = lambda.modifier)
  )) %>% mutate(method = "Cover modifier", projection = label)
}

make_proj_alpha <- function(sd, label) {
  data.frame(ode(
    c(S = sd$S0, I = sd$I0, R = sd$R0_init),
    sd$days.model, SIR_alphaexp,
    c(b = unname(p_alpha["beta"]), g = unname(p_alpha["gamma"]),
      N = sd$N, K = K_scale, q = q_val)
  )) %>% mutate(method = "Alpha exponent", projection = label)
}

projs <- bind_rows(
  make_proj_cover(sd_mid, "near_to_mid"),
  make_proj_cover(sd_off, "near_to_off"),
  make_proj_alpha(sd_mid, "near_to_mid"),
  make_proj_alpha(sd_off, "near_to_off")
) %>%
  pivot_longer(c(S, I, R), names_to = "Compartment", values_to = "tissue") %>%
  mutate(Compartment = case_match(Compartment,
                                  "S" ~ "Susceptible", "I" ~ "Infected", "R" ~ "Dead"))

obs_projs <- bind_rows(
  obs.model %>% filter(Site == "mid", Category == "Total") %>% mutate(projection = "near_to_mid"),
  obs.model %>% filter(Site == "off", Category == "Total") %>% mutate(projection = "near_to_off")
) %>%
  mutate(Compartment = case_match(Compartment,
                                  "Susceptible" ~ "Susceptible",
                                  "Infected"    ~ "Infected",
                                  "Dead"        ~ "Dead"))

p_projs <- ggplot(projs,
                  aes(time, tissue,
                      color    = method,
                      linetype = method,
                      group    = interaction(method, Compartment))) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_projs,
             aes(days.inf.site, tissue, shape = Compartment),
             color = "black", size = 1.5, inherit.aes = FALSE) +
  scale_color_manual(values = METHOD_COLORS, name = NULL) +
  scale_linetype_manual(values = METHOD_LINETYPES, name = NULL) +
  scale_shape_manual(values = c(Susceptible = 16, Infected = 17, Dead = 15), name = NULL) +
  facet_wrap(~projection, ncol = 1,
             labeller = labeller(projection = c(
               near_to_mid = "Near → Mid (Projected)",
               near_to_off = "Near → Off (Projected)"))) +
  labs(x = "Day of outbreak", y = "Surface area of tissue (m²)",
       title = "Projections from Nearshore") +
  theme_classic() +
  theme(legend.position = "bottom")

p_fit + p_projs + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")