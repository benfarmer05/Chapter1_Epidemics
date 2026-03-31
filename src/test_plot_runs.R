  
  # .rs.restartR(clean = TRUE)
  rm(list = ls())
  
  library(here)
  library(tidyverse)
  library(deSolve)
  library(patchwork)
  
  # ── Upstream context (obs data, initial conditions, helper objects) ──────────
  load(here("output/multi_SIR_workspace.RData"))
  
  # ── Load iteration results ───────────────────────────────────────────────────
  OUTPUT_DIR <- here("output/withingroup-freq-multi_SIR_runs")
  
  all_runs <- lapply(
    sort(list.files(OUTPUT_DIR, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)),
    readRDS
  )
  N_ITER <- length(all_runs)
  message(N_ITER, " iterations loaded.")
  
  ################################## Toggles ##################################
  
  INSPECT_ITER <- "random"   # set to an integer to pin, e.g. 7
  
  ################################## Constants ################################
  
  # Index positions within each params.multi[[site]] vector (see plots_multi.R)
  PARAM_NAMES <- c(
    "beta.LS", "beta.LS.adj", "gamma.LS",
    "beta.MS", "beta.MS.adj", "gamma.MS",
    "beta.HS", "beta.HS.adj", "gamma.HS",
    "R0.LS",   "R0.MS",       "R0.HS",
    "cover",   "cover.LS",    "cover.MS", "cover.HS"
  )
  
  SITE_ORDER  <- c(mid = 1, near = 2, off = 3)   # order within my.SIRS.multi
  SITE_DISPLAY_ORDER = c("off", "mid", "near")
  # SITE_LABELS <- c(near = "Nearshore", mid = "Midchannel", off = "Offshore")
  SITE_LABELS <- c(off = "Offshore", mid = "Midchannel", near = "Nearshore")
  PROJ_LABELS <- c(
    fitted      = "Fitted",
    near_to_off = "Near → Off (Projected)",
    near_to_mid = "Near → Mid (Projected)",
    off_to_near = "Off → Near (Projected)",
    off_to_mid = "Off → Mid (Projected)",
    mid_to_near = "Mid → Near (Projected)",
    mid_to_off = "Mid → Off (Projected)"
  )
  SUS_COLORS  <- c(Low = "#1E90FF", Moderate = "#FFD700", High = "#FF1493")
  
  ################################## Helpers ##################################
  
  # pivot a raw ODE data frame to long format (mirrors plots_multi.R)
  pivot_ode <- function(df, iter, site) {
    df %>%
      select(-last_col()) %>%
      pivot_longer(cols = -1,
                   names_pattern = "(.*)(..)$",
                   names_to = c("Compartment", "Susceptibility")) %>%
      mutate(
        Compartment = case_match(Compartment,
                                 "S." ~ "Susceptible", "I." ~ "Infected", "R." ~ "Dead", .default = "value"),
        Susceptibility = case_match(Susceptibility,
                                    "LS" ~ "Low", "MS" ~ "Moderate", "HS" ~ "High"),
        iteration = iter, site = site
      ) %>%
      rename(days.model = 1, tissue = value)
  }
  
  # run the ODE for a projection (source params → target initial conditions)
  run_projection <- function(src_params, tgt_init, days_vec) {
    data.frame(ode(
      y     = tgt_init,
      times = days_vec,
      func  = SIR.multi,
      parms = c(
        b.LS = unname(src_params["beta.LS"]),
        g.LS = unname(src_params["gamma.LS"]),
        b.MS = unname(src_params["beta.MS"]),
        g.MS = unname(src_params["gamma.MS"]),
        b.HS = unname(src_params["beta.HS"]),
        g.HS = unname(src_params["gamma.HS"]),
        N.LS = unname(tgt_init["S.LS"] + tgt_init["I.LS"]),
        N.MS = unname(tgt_init["S.MS"] + tgt_init["I.MS"]),
        N.HS = unname(tgt_init["S.HS"] + tgt_init["I.HS"]),
        C    = unname(src_params["cover"]),
        C.LS = unname(src_params["cover.LS"]),
        C.MS = unname(src_params["cover.MS"]),
        C.HS = unname(src_params["cover.HS"]),
        l    = lambda
      )
    ))
  }
  
  ################################## Build params_long ########################
  
  params_long <- map_dfr(all_runs, function(run) {
    map_dfr(names(SITE_ORDER), function(s) {
      p <- setNames(run$params.multi[[SITE_ORDER[s]]], PARAM_NAMES)
      as_tibble_row(p) %>% mutate(iteration = run$iteration, site = s, .before = 1)
    })
  })
  
  params_long = params_long %>%
    mutate(site = factor(site, levels = SITE_DISPLAY_ORDER))
  
  ################################## Build sims_long (fitted) #################
  
  sims_fitted <- map_dfr(all_runs, function(run) {
    map_dfr(names(SITE_ORDER), function(s) {
      pivot_ode(run$my.SIRS.multi[[SITE_ORDER[s]]], run$iteration, s)
    })
  })
  
  sims_fitted = sims_fitted %>%
    mutate(site = factor(site, levels = SITE_DISPLAY_ORDER))
  
  ################################## Build sims_long (projected) ##############
  
  # Pull site-level initial conditions from the upstream workspace,
  # replicating the HS.indices trimming logic from plots_multi.R so that
  # LS and MS can also have nonzero starting infections where observed.
  make_init <- function(s) {
    ref        <- susceptible_ref %>% filter(Site == s)
    site_label <- SITE_LABELS[s]
    
    # Identify the first timepoint where HS infection is nonzero,
    # then apply that same index to LS and MS (mirrors HS.indices in plots_multi.R)
    pull_tiss <- function(cat) {
      obs %>%
        filter(Site == site_label, Category == cat, Compartment == "Infected") %>%
        pull(tissue)
    }
    
    hs_tiss  <- pull_tiss("High")
    hs_idx   <- cumsum(!is.na(hs_tiss) & hs_tiss != 0) > 0   # same as HS.indices
    first_i  <- which(hs_idx)[1]
    
    I.LS <- pull_tiss("Low")[first_i]
    I.MS <- pull_tiss("Moderate")[first_i]
    I.HS <- hs_tiss[first_i]
    
    # Replace any NA (species not yet infected at that timepoint) with 0
    I.LS <- if (is.na(I.LS)) 0 else I.LS
    I.MS <- if (is.na(I.MS)) 0 else I.MS
    
    N.LS <- ref %>% filter(Category == "Low")      %>% pull(tissue_ref)
    N.MS <- ref %>% filter(Category == "Moderate") %>% pull(tissue_ref)
    N.HS <- ref %>% filter(Category == "High")     %>% pull(tissue_ref)
    
    c(S.LS = N.LS - I.LS, I.LS = I.LS, R.LS = 0,
      S.MS = N.MS - I.MS, I.MS = I.MS, R.MS = 0,
      S.HS = N.HS - I.HS, I.HS = I.HS, R.HS = 0)
  }
  
  # Days vectors per site (mirrors optimization script)
  make_days <- function(s) {
    d <- summary %>% filter(site == s) %>% pull(days.inf.site)
    d <- d[1:(length(d) - DHW.modifier)]
    d <- d[which(!is.na(d))[1]:length(d)]
    seq(min(d), max(d) + 120, by = 1)
  }
  
  inits <- setNames(lapply(names(SITE_ORDER), make_init), names(SITE_ORDER))
  days  <- setNames(lapply(names(SITE_ORDER), make_days),  names(SITE_ORDER))
  
  #establish N for scaling
  site_N = tibble(
    site = names(inits),
    N = sapply(inits, function(x) sum(x)) #sum of S + I + R
  ) %>%
    mutate(site = factor(site, levels = SITE_DISPLAY_ORDER))
  
  # Projections: list of (source site, target site) pairs
  PROJECTIONS <- list(
    near_to_off = c(src = "near", tgt = "off"),
    near_to_mid = c(src = "near", tgt = "mid"),
    off_to_near = c(src = "off",  tgt = "near"),
    off_to_mid = c(src = "off", tgt = "mid"),
    mid_to_near = c(src = "mid", tgt = "near"),
    mid_to_off = c(src = "mid", tgt = "off")
  )
  
  sims_projected <- map_dfr(names(PROJECTIONS), function(proj_name) {
    proj <- PROJECTIONS[[proj_name]]
    map_dfr(all_runs, function(run) {
      src_p <- setNames(run$params.multi[[SITE_ORDER[proj["src"]]]], PARAM_NAMES)
      raw   <- run_projection(src_p, inits[[proj["tgt"]]], days[[proj["tgt"]]])
      pivot_ode(raw, run$iteration, proj["tgt"]) %>%
        mutate(projection = proj_name)
    })
  }) %>% mutate(projection = coalesce(projection, "fitted"))
  
  sims_fitted <- sims_fitted %>% mutate(projection = "fitted")
  
  sims_long <- bind_rows(sims_fitted, sims_projected)
  
  sims_long = sims_long %>%
    mutate(site = factor(site, levels = SITE_DISPLAY_ORDER))
  
  sims_long_scaled = sims_long %>%
    left_join(site_N, by = "site") %>%
    mutate(tissue = tissue / N) %>%
    select(-N)
  
  obs_multi_scaled <- obs.multi %>%
    left_join(site_N %>% mutate(Site = SITE_LABELS[as.character(site)]),
              by = "Site") %>%
    mutate(tissue = tissue / N) %>%
    select(-N, -site)
  
  obs_model_scaled = obs.model %>%
    left_join(site_N, by = c("Site" = "site")) %>%
    mutate(tissue = tissue / N) %>%
    select(-N)
  
  ################################## Parameter summary ########################
  
  param_summary <- params_long %>%
    group_by(site) %>%
    summarise(across(
      c(beta.LS, gamma.LS, R0.LS, beta.MS, gamma.MS, R0.MS, beta.HS, gamma.HS, R0.HS),
      list(median = median,
           lo95   = ~quantile(., 0.025),
           hi95   = ~quantile(., 0.975)),
      .names = "{.col}__{.fn}"
    )) %>%
    pivot_longer(-site, names_to = c("parameter", "stat"), names_sep = "__") %>%
    pivot_wider(names_from = stat, values_from = value)
  
  print(param_summary)
  
  ################################## Plot helpers #############################
  
  theme_set(theme_classic(base_family = "Arial"))
  
  # median ribbon layer
  med_line <- function(data) {
    med <- data %>%
      group_by(days.model, Compartment, Susceptibility) %>%
      summarise(tissue = median(tissue), .groups = "drop")
    geom_line(data = med, aes(x = days.model, y = tissue,
                              color = Susceptibility, group = Susceptibility),
              linewidth = 1, alpha = 1)
  }
  
  # observed points layer
  obs_points <- function(site_label, compartment = NULL) {
    d <- obs.multi %>%
      rename(Susceptibility = Category) %>%
      filter(Site == site_label)
    if (!is.null(compartment)) d <- filter(d, Compartment == compartment)
    geom_point(data = d,
               aes(x = days.inf.site, y = tissue,
                   color = Susceptibility, shape = Compartment),
               size = 1.5)
  }
  
  # core layered plot function
  layer_plot <- function(sim_data, site, proj = "fitted",
                         compartments = c("Susceptible", "Infected", "Dead"),
                         title = NULL,
                         obs_src = obs.multi,
                         y_lab   = "Surface area of tissue (m²)") {
    
    d <- sim_data %>%
      filter(site == !!site, projection == proj,
             Compartment %in% compartments) %>%
      mutate(group_id = paste(iteration, Susceptibility))
    
    med <- d %>%
      group_by(days.model, Compartment, Susceptibility) %>%
      summarise(tissue = median(tissue), .groups = "drop")
    
    site_label <- SITE_LABELS[site]
    obs_d <- obs_src %>%
      rename(Susceptibility = Category) %>%
      filter(Site == site_label, Compartment %in% compartments)
    
    ggplot() +
      geom_line(data = d,
                aes(x = days.model, y = tissue,
                    color = Susceptibility, group = group_id),
                alpha = 0.12, linewidth = 0.4) +
      geom_line(data = med,
                aes(x = days.model, y = tissue,
                    color = Susceptibility, group = Susceptibility),
                linewidth = 1) +
      geom_point(data = obs_d,
                 aes(x = days.inf.site, y = tissue,
                     color = Susceptibility, shape = Compartment),
                 size = 1.5) +
      scale_color_manual(values = SUS_COLORS) +
      scale_shape_manual(values = c(Susceptible = 16, Infected = 17, Dead = 15)) +
      labs(x = "Day of outbreak", y = y_lab,
           title = title %||% paste(site_label, "-", PROJ_LABELS[proj])) +
      theme(legend.position = "bottom")
  }
  
  ################################## Single-iteration inspection ##############
  
  inspect_iter <- if (identical(INSPECT_ITER, "random")) sample(seq_len(N_ITER), 1) else INSPECT_ITER
  message("Inspecting iteration: ", inspect_iter)
  
  plot_single_iter <- function(iter, site, proj = "fitted",
                               compartments = c("Susceptible", "Infected", "Dead")) {
    d <- sims_long %>%
      filter(iteration == iter, site == !!site, projection == proj,
             Compartment %in% compartments)
    
    site_label <- SITE_LABELS[site]
    
    ggplot(d, aes(days.model, tissue, color = Susceptibility, linetype = Compartment)) +
      geom_line(linewidth = 0.8) +
      obs_points(site_label) +
      scale_color_manual(values = SUS_COLORS) +
      scale_shape_manual(values = c(Susceptible = 16, Infected = 17, Dead = 15)) +
      labs(x = "Day of outbreak", y = "Surface area of tissue (m²)",
           title = paste0("Iter ", iter, " | ", site_label, " | ", PROJ_LABELS[proj])) +
      theme(legend.position = "bottom")
  }
  
  # Example: view one site at a time
  p_inspect <- wrap_plots(
    lapply(SITE_DISPLAY_ORDER, \(s) plot_single_iter(inspect_iter, s)),
    nrow = 1
  ) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  p_inspect
  
  ################################## Layered all-iteration plots ##############
  
  # Equivalent of fig3: fitted nearshore + fitted offshore + near→off projection
  # (single-column: Dead compartment only, mirroring p.D.fit.* panels)
  
  make_fig3_panel <- function(site, proj, compartment = "Dead") {
    layer_plot(sims_long, site, proj, compartments = compartment,
               title = paste(SITE_LABELS[site], "-", PROJ_LABELS[proj]))
  }
  
  fig3_layered <-
    make_fig3_panel("near", "fitted")    +
    make_fig3_panel("off",  "fitted")    +
    make_fig3_panel("off",  "near_to_off") +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  fig3_layered
  
  # ── Projection-only layered plot (colour by susceptibility) ─────────────────
  
  layer_plot_proj <- function(proj, compartment = "Dead", title = NULL) {
    tgt <- PROJECTIONS[[proj]]["tgt"]
    layer_plot(sims_long, site = tgt, proj = proj,
               compartments = compartment,
               title = title %||% PROJ_LABELS[proj])
  }
  
  fig3_projections <-
    layer_plot_proj("near_to_off") +
    layer_plot_proj("near_to_mid") +
    layer_plot_proj("off_to_near") +
    layer_plot_proj("off_to_mid") +
    layer_plot_proj("mid_to_near") +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  fig3_projections
  
  # fig3_projections_SIR <-
  #   layer_plot_proj("near_to_off", compartment = c("Susceptible", "Infected", "Dead")) +
  #   layer_plot_proj("near_to_mid", compartment = c("Susceptible", "Infected", "Dead")) +
  #   layer_plot_proj("off_to_near", compartment = c("Susceptible", "Infected", "Dead")) +
  #   plot_layout(ncol = 1, guides = "collect") &
  #   theme(legend.position = "bottom")
  # 
  # fig3_projections_SIR
  
  # ── Same but collapsed to a single black total line ──────────────────────────
  
  layer_plot_collapsed <- function(proj, compartment = "Dead", title = NULL,
                                   sim_data = sims_long,
                                   obs_src  = obs.model,
                                   y_lab    = "Surface area of tissue (m²)") {
    
    tgt        <- PROJECTIONS[[proj]]["tgt"]
    site_label <- SITE_LABELS[tgt]
    
    d <- sim_data %>%
      filter(site == tgt, projection == proj, Compartment == compartment) %>%
      group_by(iteration, days.model) %>%
      summarise(tissue = sum(tissue), .groups = "drop")
    
    med <- d %>%
      group_by(days.model) %>%
      summarise(tissue = median(tissue), .groups = "drop")
    
    obs_d <- obs_src %>%
      filter(Site == tgt, Category == "Total", Compartment == compartment)
    
    ggplot() +
      geom_line(data = d,
                aes(x = days.model, y = tissue, group = iteration),
                color = "black", alpha = 0.12, linewidth = 0.4) +
      geom_line(data = med,
                aes(x = days.model, y = tissue),
                color = "black", linewidth = 1) +
      geom_point(data = obs_d,
                 aes(x = days.inf.site, y = tissue),
                 color = "black", shape = 15, size = 1.5) +
      labs(x = "Day of outbreak", y = y_lab,
           title = title %||% PROJ_LABELS[proj]) +
      theme(legend.position = "none")
  }
  
  fig3_projections_collapsed <-
    layer_plot_collapsed("near_to_off") +
    layer_plot_collapsed("near_to_mid") +
    layer_plot_collapsed("off_to_near") +
    plot_layout(ncol = 1)
  
  fig3_projections_collapsed
  
  
  
  layer_plot_collapsed_SIR <- function(proj, title = NULL,
                                       sim_data = sims_long,
                                       obs_src  = obs.model,
                                       y_lab    = "Surface area of tissue (m²)") {
    
    tgt        <- PROJECTIONS[[proj]]["tgt"]
    site_label <- SITE_LABELS[tgt]
    
    d <- sim_data %>%
      filter(site == tgt, projection == proj,
             Compartment %in% c("Susceptible", "Infected", "Dead")) %>%
      group_by(iteration, days.model, Compartment) %>%
      summarise(tissue = sum(tissue), .groups = "drop")
    
    med <- d %>%
      group_by(days.model, Compartment) %>%
      summarise(tissue = median(tissue), .groups = "drop")
    
    obs_d <- obs_src %>%
      filter(Site == tgt, Category == "Total",
             Compartment %in% c("Susceptible", "Infected", "Dead"))
    
    ggplot() +
      geom_line(data = d,
                aes(x = days.model, y = tissue,
                    group = interaction(iteration, Compartment),
                    color = Compartment),
                alpha = 0.12, linewidth = 0.4) +
      geom_line(data = med,
                aes(x = days.model, y = tissue, color = Compartment),
                linewidth = 1) +
      geom_point(data = obs_d,
                 aes(x = days.inf.site, y = tissue, color = Compartment),
                 shape = 15, size = 1.5) +
      scale_color_manual(values = c(Susceptible = "#2196F3",
                                    Infected    = "#FF9800",
                                    Dead        = "#9C27B0")) +
      labs(x = "Day of outbreak", y = y_lab,
           title = title %||% PROJ_LABELS[proj]) +
      theme(legend.position = "bottom")
  }
  
  fig3_projections_SIR_collapsed <-
    layer_plot_collapsed_SIR("near_to_off") +
    layer_plot_collapsed_SIR("near_to_mid") +
    layer_plot_collapsed_SIR("off_to_near") +
    plot_layout(ncol = 1)
  
  fig3_projections_SIR_collapsed
  
  
  fig3_projections_I <-
    layer_plot_proj("near_to_off", compartment = "Infected") +
    layer_plot_proj("near_to_mid", compartment = "Infected") +
    layer_plot_proj("off_to_near", compartment = "Infected") +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  fig3_projections_I
  
  
  
  fig3_projections_I_collapsed <-
    layer_plot_collapsed("near_to_off", compartment = "Infected") +
    layer_plot_collapsed("near_to_mid", compartment = "Infected") +
    layer_plot_collapsed("off_to_near", compartment = "Infected") +
    plot_layout(ncol = 1)
  
  fig3_projections_I_collapsed
  
  
  ################################## scaled plots #############
  
  # # ── fig3_layered (scaled) ────────────────────────────────────────────────────
  # 
  # fig3_layered_scaled <-
  #   layer_plot(sims_long_scaled, "near", "fitted",
  #              obs_src = obs_multi_scaled,
  #              y_lab   = "Proportion of tissue") +
  #   layer_plot(sims_long_scaled, "off", "fitted",
  #              obs_src = obs_multi_scaled,
  #              y_lab   = "Proportion of tissue") +
  #   layer_plot(sims_long_scaled, "off", "near_to_off",
  #              obs_src = obs_multi_scaled,
  #              y_lab   = "Proportion of tissue") +
  #   plot_layout(ncol = 1, guides = "collect") &
  #   theme(legend.position = "bottom")
  # 
  # fig3_layered
  # fig3_layered_scaled
  # 
  # # ── fig3_projections (scaled) ────────────────────────────────────────────────
  # 
  # fig3_projections_scaled <-
  #   layer_plot(sims_long_scaled, "off",  "near_to_off",
  #              obs_src = obs_multi_scaled, y_lab = "Proportion of tissue") +
  #   layer_plot(sims_long_scaled, "mid",  "near_to_mid",
  #              obs_src = obs_multi_scaled, y_lab = "Proportion of tissue") +
  #   layer_plot(sims_long_scaled, "near", "off_to_near",
  #              obs_src = obs_multi_scaled, y_lab = "Proportion of tissue") +
  #   plot_layout(ncol = 1, guides = "collect") &
  #   theme(legend.position = "bottom")
  # 
  # fig3_projections_scaled
  
  # ── fig3_projections_collapsed (scaled) ──────────────────────────────────────
  
  fig3_projections_collapsed_scaled <-
    layer_plot_collapsed("near_to_off",
                         sim_data = sims_long_scaled,
                         obs_src  = obs_model_scaled,
                         y_lab    = "Proportion of tissue") +
    layer_plot_collapsed("near_to_mid",
                         sim_data = sims_long_scaled,
                         obs_src  = obs_model_scaled,
                         y_lab    = "Proportion of tissue") +
    layer_plot_collapsed("off_to_near",
                         sim_data = sims_long_scaled,
                         obs_src  = obs_model_scaled,
                         y_lab    = "Proportion of tissue") +
    plot_layout(ncol = 1) &
    coord_cartesian(ylim = c(0, 1))
  
  fig3_projections_collapsed_scaled
  
  # ── fig3_projections_SIR_collapsed (scaled) ──────────────────────────────────
  
  fig3_projections_SIR_collapsed_scaled <-
    layer_plot_collapsed_SIR("near_to_off",
                             sim_data = sims_long_scaled,
                             obs_src  = obs_model_scaled,
                             y_lab    = "Proportion of tissue") +
    layer_plot_collapsed_SIR("near_to_mid",
                             sim_data = sims_long_scaled,
                             obs_src  = obs_model_scaled,
                             y_lab    = "Proportion of tissue") +
    layer_plot_collapsed_SIR("off_to_near",
                             sim_data = sims_long_scaled,
                             obs_src  = obs_model_scaled,
                             y_lab    = "Proportion of tissue") +
    plot_layout(ncol = 1)
  
  fig3_projections_SIR_collapsed_scaled
  
  # ── fig3_projections_I (scaled) ──────────────────────────────────────────────
  
  fig3_projections_I_scaled <-
    layer_plot(sims_long_scaled, "off",  "near_to_off",
               compartments = "Infected",
               obs_src = obs_multi_scaled, y_lab = "Proportion of tissue") +
    layer_plot(sims_long_scaled, "mid",  "near_to_mid",
               compartments = "Infected",
               obs_src = obs_multi_scaled, y_lab = "Proportion of tissue") +
    layer_plot(sims_long_scaled, "near", "off_to_near",
               compartments = "Infected",
               obs_src = obs_multi_scaled, y_lab = "Proportion of tissue") +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  fig3_projections_I_scaled
  
  # ── fig3_projections_I_collapsed (scaled) ────────────────────────────────────
  
  fig3_projections_I_collapsed_scaled <-
    layer_plot_collapsed("near_to_off", compartment = "Infected",
                         sim_data = sims_long_scaled,
                         obs_src  = obs_model_scaled,
                         y_lab    = "Proportion of tissue") +
    layer_plot_collapsed("near_to_mid", compartment = "Infected",
                         sim_data = sims_long_scaled,
                         obs_src  = obs_model_scaled,
                         y_lab    = "Proportion of tissue") +
    layer_plot_collapsed("off_to_near", compartment = "Infected",
                         sim_data = sims_long_scaled,
                         obs_src  = obs_model_scaled,
                         y_lab    = "Proportion of tissue") +
    plot_layout(ncol = 1)
  
  fig3_projections_I_collapsed_scaled
  
  ################################## Parameter distribution plots #############
  
  p_params <- params_long %>%
    select(iteration, site, beta.LS, beta.MS, beta.HS, gamma.LS, gamma.MS, gamma.HS) %>%
    pivot_longer(-c(iteration, site)) %>%
    mutate(name = factor(name, levels = c(
      "beta.LS", "beta.MS", "beta.HS",
      "gamma.LS", "gamma.MS", "gamma.HS"
    ))) %>%
    mutate(site = factor(SITE_LABELS[site], levels = SITE_LABELS[SITE_DISPLAY_ORDER])) %>%
    ggplot(aes(value, fill = site)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    # facet_wrap(~name, scales = "free") + 
    # facet_wrap(~name, scales = "free_y") +
    facet_wrap(~name, scales = "fixed") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Count", fill = "Site") +
    theme(legend.position = "bottom")
  
  p_params
  
  
  
  
  
  
  p_params <- params_long %>%
    select(iteration, site, beta.LS, beta.MS, beta.HS, gamma.LS, gamma.MS, gamma.HS) %>%
    pivot_longer(-c(iteration, site)) %>%
    mutate(
      param = case_when(
        startsWith(name, "beta")  ~ "beta",
        startsWith(name, "gamma") ~ "gamma"
      ),
      susceptibility = case_when(
        endsWith(name, "LS") ~ "Low",
        endsWith(name, "MS") ~ "Moderate",
        endsWith(name, "HS") ~ "High"
      ),
      susceptibility = factor(susceptibility, levels = c("Low", "Moderate", "High")),
      param = factor(param, levels = c("beta", "gamma"),
                     labels = c("β (transmission)", "γ (recovery)")),
      site  = factor(SITE_LABELS[site], levels = SITE_LABELS[SITE_DISPLAY_ORDER])
    ) %>%
    ggplot(aes(value, fill = susceptibility)) +
    # geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_density(alpha = 0.4) +
    # facet_grid(site ~ param, scales = "free_x") +
    facet_grid(site ~ param, scales = "free") +
    scale_fill_manual(values = SUS_COLORS) +
    # labs(x = NULL, y = "Count", fill = "Susceptibility") +
    labs(x = NULL, y = "Density", fill = "Susceptibility") +
    theme(legend.position = "bottom")
  
  p_params
  
  
  
  library(gt)
  
  param_summary %>%
    filter(!parameter %in% c("R0.LS", "R0.MS", "R0.HS")) %>%
    mutate(
      site      = factor(SITE_LABELS[site], levels = SITE_LABELS[SITE_DISPLAY_ORDER]),
      parameter = factor(parameter, levels = c(
        "beta.LS", "beta.MS", "beta.HS",
        "gamma.LS", "gamma.MS", "gamma.HS",
        "R0.LS", "R0.MS", "R0.HS"
      ))
    ) %>%
    arrange(site, parameter) %>%
    gt(groupname_col = "site") %>%
    cols_label(
      parameter = "Parameter",
      median    = "Median",
      lo95      = "2.5%",
      hi95      = "97.5%"
    ) %>%
    fmt_number(columns = c(median, lo95, hi95), decimals = 4) %>%
    tab_spanner(label = "95% Credible Interval", columns = c(lo95, hi95)) %>%
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_row_groups()
    )
  