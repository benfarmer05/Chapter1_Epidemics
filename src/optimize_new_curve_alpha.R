  
  # NOTE:
  #   - the below script is only intended to understanding the behavior of alpha (effect of host density
  #       on SCTLD transmission). the final value of alpha actually used in publication was found by starting here,
  #       and then going back to output_basic.R and working toward a value that best projected Nearshore -> Offshore
  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(ggplot2)
  library(dplyr)
  library(here)
  library(deSolve)
  library(extrafont)
  
  ################################## Set-up ##################################
  # NOTE - run this while working within plots_basic.R and/or plots_multi.R if doing more than just the
  #         plotting sandbox section
  # load(here("output/plots_basic_workspace.RData"))
  load(here("output/plots_multi_workspace.RData"))
  
  ################################## sandbox: cover ##################################
  
  CC = seq(0.001, 1, .001)
  a = seq(0, 1, 0.05) #alpha (weight of coral cover)
  a = sort(c(a, 0.13))  # add 0.13 and sort the sequence
  k = 3 #shape of curve?
  
  # Create data frame for all curves
  plot_data <- expand.grid(CC = CC, alpha = a) %>%
    mutate(
      mods = (1 - alpha) + alpha * ((1 - exp(-3 * CC)) / (1 - exp(-3))),
      coral_cover_pct = CC * 100,
      is_special = alpha == 0.13
    )
  
  # Create data for special points
  special_points <- data.frame(
    coral_cover_pct = c(0.247 * 100, 0.0215 * 100),
    mods = c(
      (1 - 0.13) + 0.13 * ((1 - exp(-3 * 0.247)) / (1 - exp(-3))),
      (1 - 0.13) + 0.13 * ((1 - exp(-3 * 0.0215)) / (1 - exp(-3)))
    ),
    point_type = c("Nearshore", "Offshore")
  )
  
  # Create the plot
  p <- ggplot() +
    # Add gradient lines (excluding the special alpha = 0.13)
    geom_line(data = filter(plot_data, !is_special), 
              aes(x = coral_cover_pct, y = mods, color = alpha, group = alpha),
              linewidth = 0.6) +
    
    # Add special alpha = 0.13 line in blue
    geom_line(data = filter(plot_data, is_special),
              aes(x = coral_cover_pct, y = mods),
              color = "blue", linewidth = 0.6) + #linetype = 'dashed'
    
    # Add special points (using fill aesthetic to avoid conflict)
    geom_point(data = special_points,
               aes(x = coral_cover_pct, y = mods, fill = point_type),
               shape = 23, size = 2, stroke = 1, color = "black") +
    
    # Scales and colors
    scale_color_viridis_c(name = expression(alpha), option = "inferno") +
    scale_fill_manual(name = "", 
                      values = c("Nearshore" = "orange", "Offshore" = "magenta"),
                      guide = "none") +
    
    # Axes
    scale_x_continuous(name = "Coral cover (%)", 
                       limits = c(0, 100),
                       expand = c(0, 0)) +
    scale_y_continuous(name = "Scalar",
                       limits = c(0, 1),
                       expand = c(0, 0),
                       breaks = seq(0, 1, 0.2)) +
    
    
    
    # Theme
    theme_classic(base_family = 'Georgia') +
    theme(
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7, color = 'black'),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(t = 8, r = 15, b = 8, l = 8),
      # Position legends inside the plot
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.75), # moved to upper right
      legend.box = "vertical",
      legend.margin = margin(5, 5, 5, 5),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.3)
    ) +
    
    # Guides for legends
    guides(
      color = guide_colorbar(
        title = expression(alpha),
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 0.8,
        barheight = 3,
        order = 1
      )
    )
  
  # Create a custom legend using annotations
  fig2a <- p + 
    # Add manual legend box (moved to middle right area)
    annotate("rect", xmin = 70, xmax = 95, ymin = 0.25, ymax = 0.45, 
             fill = "white", color = NA, linewidth = 0.3) +
    
    # Alpha = 0.13 line sample and label
    annotate("segment", x = 72, xend = 76, y = 0.41, yend = 0.41,
             color = "blue", linewidth = 1) + # ,linetype = 'dashed'
    annotate("text", x = 77, y = 0.41, 
             label = expression(alpha == 0.13), 
             hjust = 0, size = 3, family = "Georgia") +
    
    # Nearshore point and label
    annotate("point", x = 74, y = 0.36, 
             shape = 23, fill = "orange", color = "black", size = 2, stroke = 1) +
    annotate("text", x = 77, y = 0.36, 
             label = "Nearshore", 
             hjust = 0, size = 3, family = "Georgia") +
    
    # Offshore point and label  
    annotate("point", x = 74, y = 0.31, 
             shape = 23, fill = "magenta", color = "black", size = 2, stroke = 1) +
    annotate("text", x = 77, y = 0.31, 
             label = "Offshore", 
             hjust = 0, size = 3, family = "Georgia")
  
  #max dimensions are 7.087 in. wide by 9.45 in. tall (3.35 inches preferred)
  quartz(h = 3, w = 3.35)
  
  fig2a
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "fig2a.pdf"), type = "pdf")
  
  #ggplot-export to image
  ggsave(filename = here("output", "fig2a.png"), device = "png", width = 3.35, height = 3, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
  
  
  # # IN BASE R
  # 
  # # COVER-BASED VERSION
  # 
  # CC = seq(0.001,1,.001)
  # # a = seq(0,1,0.01) #alpha (weight of coral cover)
  # a = seq(0,1,0.05) #alpha (weight of coral cover)
  # a = sort(c(a, 0.13))  # add 0.13 and sort the sequence (if 0.13 or other value of interest is already there, can comment this out)
  # k = 3 #shape of curve?
  # 
  # dpi = 1200
  # png(here("output/transmission_plot.png"), width = 3.35 * dpi, height = 3 * dpi, res = dpi)
  # # pdf(here("output/transmission_plot.pdf"), width = 3.35, height = 3, family = 'Georgia')
  # 
  # par(family = "Georgia")
  # # par(mar = c(4.5, 4.5, 1, 1)) #adjust margins (bottom, left, top, right)
  # par(mar = c(3, 3, 1, 1)) #adjust margins (bottom, left, top, right)
  # par(mgp = c(1.5, 0.5, 0))  # Reduce space between axis labels and titles
  # 
  # plot(NA, NA, 
  #      xlim = c(0, 100), 
  #      ylim = c(0, 1),
  #      xlab = "Coral cover (%)", 
  #      ylab = "Scalar",
  #      cex.lab = 0.75,
  #      cex.axis = 0.6
  # )
  # 
  # for (i in 1:length(a)) {
  #   mods.1 <- (1 - a[i]) + a[i] * ((1 - exp(-3 * CC)) / (1 - exp(-3)))
  #   
  #   col_choice <- ifelse(a[i] == 0.13, 'blue', adjustcolor('red', alpha.f = 1))  
  #   lwd_choice <- ifelse(a[i] == 0.13, 2, 1)  # Thicker for α = 0.13
  #   
  #   #plot in terms of 0 to 100 coral cover
  #   lines(CC * 100, mods.1, col = col_choice, lwd = lwd_choice)
  #   
  #   if (a[i] == 0.13) {
  #     points(0.247 * 100, (1 - a[i]) + a[i] * ((1 - exp(-3 * 0.247)) / (1 - exp(-3))), 
  #            pch = 8, col = "orange", cex = 0.75, lwd = 2)  # Nearshore (orange)
  #     
  #     points(0.0215 * 100, (1 - a[i]) + a[i] * ((1 - exp(-3 * 0.0215)) / (1 - exp(-3))), 
  #            pch = 8, col = "magenta", cex = 0.75, lwd = 2)  # Offshore (magenta)
  #   }
  #   
  # }
  # 
  # # Add legend slightly lower than the top right
  # legend("topright", 
  #        legend = c(expression(alpha == 0.13), "Nearshore", "Offshore"), 
  #        col = c("blue", "orange", "magenta"), 
  #        pch = c(NA, 8, 8), 
  #        lwd = c(2, NA, NA), 
  #        bty = "n",
  #        cex = 0.6,
  #        inset = c(0, 0.1)) # Move the legend slightly downward
  # 
  # dev.off()
  
  ################################## sandbox: N ##################################
  
  # N-BASED VERSION
  reef_area = 1000 #or could do 600, the 2-D area of reef sites, in m2
  N.offshore.scaled = N.offshore / reef_area
  N.nearshore.scaled = N.nearshore / reef_area
  N_range = seq(0,reef_area,1)
  # N_range = N_range / reef_area
  a = seq(0,1,0.01) #alpha (weight of N)
  k = 3 #shape of curve
  
  #self-scaled
  
  # Set font to Georgia
  par(family = "Georgia")  
  
  plot(NA, NA, 
       xlim = range(N_range), 
       ylim = c(0, 1),
       # xlab = "Population size (density of tissue SA per 1000 m²)", 
       xlab = "Population size (m²)", 
       ylab = "Transmission Modifier")
  
  for (i in 1:length(a)) {
    
    mods.1 <- (1 - a[i]) + a[i] * ((1 - exp(-3 * N_range * (1/max(N_range)))) / (1 - exp(-3)))
    
    col_choice <- ifelse(a[i] == 0.13, 'blue', adjustcolor('red', alpha.f = 0.5))  
    lwd_choice <- ifelse(a[i] == 0.13, 2, 0.5)  # Thicker for α = 0.13
    
    lines(N_range, mods.1, col = col_choice, lwd = lwd_choice)
    
    # Add site-specific stars only for α = 0.13
    if (a[i] == 0.13) {
      points(N.nearshore, (1 - a[i]) + a[i] * ((1 - exp(-3 * N.nearshore * (1/max(N_range)))) / (1 - exp(-3))), 
             pch = 8, col = "orange", cex = 1.5, lwd = 2)  # Nearshore (orange)
      
      points(N.offshore, (1 - a[i]) + a[i] * ((1 - exp(-3 * N.offshore * (1/max(N_range)))) / (1 - exp(-3))), 
             pch = 8, col = "magenta", cex = 1.5, lwd = 2)  # Offshore (magenta)
    }
  }
  
  # Add legend slightly lower than the top right
  legend("topright", 
         legend = c(expression(alpha == 0.13), "Nearshore", "Offshore"), 
         col = c("blue", "orange", "magenta"), 
         pch = c(NA, 8, 8), 
         lwd = c(2, NA, NA), 
         bty = "n",
         # y.intersp = 1.5,   # Increase spacing between legend items
         inset = c(0, 0.1)) # Move the legend slightly downward
  
  
  ################################## single-host projection optimization ##################################
  
  SIR_project = function(t,y,p){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p),{

      #with cover
      transmission_modifier = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C)) / (1 - exp(-k_val)))
      
      # #with N
      # transmission_modifier = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*N*(1/max(N_range)))) / (1 - exp(-k_val)))
      
      #frequency-dependent
      dS.dt = -b * S * I / N * transmission_modifier
      dI.dt = b * S * I / N * transmission_modifier - g * I
      dR.dt = g * I
      
      # #density-dependent
      # dS.dt = -b * S * I * transmission_modifier
      # dI.dt = b * S * I * transmission_modifier - g * I
      # dR.dt = g * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  # Define the search space
  alpha_values <- seq(0, 1, length.out = 1000)
  k_val = 3
  
  best_r_squared <- -Inf
  best_alpha <- NA
  
  for (i in alpha_values) {
    
    alpha_val = i # NOTE - bad and hard-coded way to feed parameter into the SIR function - should look at this further
    
    # Run the model
    output.basic.offshore.transfer = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                                    days.model.offshore, SIR_project, c(b = params.basic.nearshore.full[1], g = params.basic.nearshore.full[3],
                                                                                N = N.offshore,
                                                                                C = cover.offshore.full,
                                                                                l = lambda.modifier)))
    
    
    # Compute R-squared
    sim.rem.total = output.basic.offshore.transfer[which(output.basic.offshore.transfer$time %in% days.obs),
                                                   which(colnames(output.basic.offshore.transfer) %in% c('R'))]
    
    # NOTE - this assumes that obs.rem.total is pulling data for offshore correctly from the upstream environment. note if any issues
    if (length(obs.rem.total) > length(sim.rem.total)) {
      obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
    }
    
    diff.rem.total = (sim.rem.total - obs.rem.total)
    sum_diff.total = sum(diff.rem.total^2, na.rm = TRUE) # NOTE - had to account for NAs when fit is really bad
    mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
    tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
    r_squared.near.to.off.basic = 1 - (sum_diff.total / tss_rem.total)
    
    # Check if this is the best so far
    if (r_squared.near.to.off.basic > best_r_squared) {
      best_r_squared <- r_squared.near.to.off.basic
      best_alpha <- i
    }
  }
  
  # Print the best parameters
  print(best_alpha)
  print(best_r_squared)
  
  #best alpha value was: 0.07107107
  
  ################################## multi-host projection optimization ##################################
  
  SIR.multi_project = function(t,y,p){
    {
      S.LS = y[1]
      I.LS = y[2]
      R.LS = y[3]
      
      S.MS = y[4]
      I.MS = y[5]
      R.MS = y[6]
      
      S.HS = y[7]
      I.HS = y[8]
      R.HS = y[9]
    }
    with(as.list(p),{
      P = (I.LS + I.MS + I.HS)
      
      # ### TEST ###
      # transmission_modifier.LS = (1 - 0) + 0*((1 - exp(-3*0.01460644)) / (1 - exp(-3)))
      # transmission_modifier.MS = (1 - alpha_val) + alpha_val*((1 - exp(-3*0.003109595)) / (1 - exp(-3)))
      # transmission_modifier.HS = (1 - alpha_val) + alpha_val*((1 - exp(-3*0.0009200924)) / (1 - exp(-3)))
      # 
      # transmission_modifier.LS
      # transmission_modifier.MS
      # transmission_modifier.HS
      # ### TEST ###
      
      transmission_modifier.LS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.LS)) / (1 - exp(-k_val)))
      transmission_modifier.MS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.MS)) / (1 - exp(-k_val)))
      transmission_modifier.HS = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C.HS)) / (1 - exp(-k_val)))
      
      # #hybrid of frequency and density-dependent
      # dS.LS.dt = -b.LS*S.LS*(P) / N.LS * transmission_modifier.LS
      # dI.LS.dt = b.LS*S.LS*(P) / N.LS * transmission_modifier.LS - g.LS*I.LS
      # dR.LS.dt = g.LS*I.LS
      # 
      # dS.MS.dt = -b.MS*S.MS*(P) / N.MS * transmission_modifier.MS
      # dI.MS.dt = b.MS*S.MS*(P) / N.MS * transmission_modifier.MS - g.MS*I.MS
      # dR.MS.dt = g.MS*I.MS
      # 
      # dS.HS.dt = -b.HS*S.HS*(P) / N.HS * transmission_modifier.HS
      # dI.HS.dt = b.HS*S.HS*(P) / N.HS * transmission_modifier.HS - g.HS*I.HS
      # dR.HS.dt = g.HS*I.HS
      
      #more frequency-dependent
      dS.LS.dt = -b.LS*S.LS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.LS
      dI.LS.dt = b.LS*S.LS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.LS - g.LS*I.LS
      dR.LS.dt = g.LS*I.LS
      
      dS.MS.dt = -b.MS*S.MS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.MS
      dI.MS.dt = b.MS*S.MS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.MS - g.MS*I.MS
      dR.MS.dt = g.MS*I.MS
      
      dS.HS.dt = -b.HS*S.HS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.HS
      dI.HS.dt = b.HS*S.HS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.HS - g.HS*I.HS
      dR.HS.dt = g.HS*I.HS
      
      return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, dS.HS.dt, dI.HS.dt, dR.HS.dt), P = P))
    })
  }
  
  
  # Define the search space
  # alpha_values <- seq(0.1, 1, length.out = 1000)
  alpha_values <- seq(0, 1, length.out = 1000)
  # alpha_values <- seq(0, 0, length.out = 10)
  k_val = 3
  
  best_r_squared <- -Inf
  best_alpha <- c(alpha.LS = NA, alpha.MS = NA, alpha.HS = NA)
  
  site.loop = 'Offshore'
  curr.site = 'off'
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  for (i in alpha_values) {
    
    alpha_val = i # NOTE - bad and hard-coded way to feed parameter into the SIR function - should look at this further
    # alpha_val = 0 #test
    
    # Run the model
    output.near.to.off.multi = data.frame(ode(c(S.LS = S.LS.offshore, I.LS = I.LS.offshore, R.LS = R.LS.offshore,
                                                S.MS = S.MS.offshore, I.MS = I.MS.offshore, R.MS = R.MS.offshore,
                                                S.HS = S.HS.offshore, I.HS = I.HS.offshore, R.HS = R.HS.offshore),
                                              days.model.offshore, SIR.multi_project, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                                                                b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                                                                b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                                                                N.LS = N.LS.offshore, N.MS = N.MS.offshore, N.HS = N.HS.offshore,
                                                                                C = cover.offshore,
                                                                                C.LS = cover.offshore.LS, C.MS = cover.offshore.MS, C.HS = cover.offshore.HS,
                                                                                l = lambda)))
    
    # Compute R-squared
    sim.rem.total = rowSums(output.near.to.off.multi[which(output.near.to.off.multi$time %in% days.obs),
                                                     which(colnames(output.near.to.off.multi) %in% c('R.LS', 'R.MS', 'R.HS'))],
                            na.rm = TRUE)
    
    if (length(obs.rem.total) > length(sim.rem.total)) {
      obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
    }
    
    diff.rem.total = (sim.rem.total - obs.rem.total)
    sum_diff.total = sum(diff.rem.total^2, na.rm = TRUE) # NOTE - had to account for NAs when fit is really bad
    mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
    tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
    r_squared.near.to.off.multi = 1 - (sum_diff.total / tss_rem.total)
    
    # Check if this is the best so far
    if (r_squared.near.to.off.multi > best_r_squared) {
      best_r_squared <- r_squared.near.to.off.multi
      best_alpha <- i
    }
  }
  
  # Print the best parameters
  print(best_alpha)
  print(best_r_squared)
  
  
  
  #best alpha value was: 0.004004004
  
  
  
  
  
  ################################## sandbox for manually tweaking alpha ##################################
  
  #sandbox conditions
  beta.nearshore.sand = params.basic.nearshore.full[1] #0.65137 # 0.76 #0.64
  gamma.nearshore.sand = params.basic.nearshore.full[3] #0.5622153 # 0.56
  alpha_val = 0.07107107 # 0.03 # 0.07107107 # 0
  k_val = 3
  lambda_val = 1.0
  offset_val = 1 - 1 / (1 + exp(-lambda_val * 1.0))
  
  
  SIR.sand.no_cover = function(t,y,p){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p),{
      
      #null conditions
      transmission_modifier = 1
      
      # #with effect of coral cover
      # transmission_modifier = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C)) / (1 - exp(-k_val)))
      
      dS.dt = -b * S * I / N * transmission_modifier
      dI.dt = b * S * I / N * transmission_modifier - g * I
      dR.dt = g * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  SIR.sand.cover = function(t,y,p){ # 'p' is parameters or params
    {
      S = y[1]
      I = y[2]
      R = y[3]
    }
    with(as.list(p),{
      
      # #null conditions
      # transmission_modifier = 1
      
      #with effect of coral cover
      transmission_modifier = (1 - alpha_val) + alpha_val*((1 - exp(-k_val*C)) / (1 - exp(-k_val)))
      # transmission_modifier = (1 / (1 + exp(-lambda_val * (C))) + offset_val)
      
      dS.dt = -b * S * I / N * transmission_modifier
      dI.dt = b * S * I / N * transmission_modifier - g * I
      dR.dt = g * I
      
      return(list(c(dS.dt, dI.dt, dR.dt)))
    })
  }
  
  #run nearshore model with rates we have already fitted with null condition for coral cover effect 
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  output.basic.nearshore.sand = data.frame(ode(c(S = S.nearshore, I = I.nearshore, R = R.nearshore),
                                                  days.model.nearshore, SIR.sand.no_cover, c(b = beta.nearshore.sand, g = gamma.nearshore.sand,
                                                                              N = N.nearshore,
                                                                              l = lambda.nearshore,
                                                                              C = cover.nearshore)))
  
  sim.rem.total = output.basic.nearshore.sand[which(output.basic.nearshore.sand$time %in% days.obs), which(colnames(output.basic.nearshore.sand) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore.sand.no_cover = 1 - (sum_diff.total / tss_rem.total)
  
  #run nearshore model with rates we already have but also with new alpha values (effect of coral cover)
  output.basic.nearshore.sand = data.frame(ode(c(S = S.nearshore, I = I.nearshore, R = R.nearshore),
                                               days.model.nearshore, SIR.sand.cover, c(b = beta.nearshore.sand, g = gamma.nearshore.sand,
                                                                                          N = N.nearshore,
                                                                                          l = lambda.nearshore,
                                                                                          C = cover.nearshore)))
  
  sim.rem.total = output.basic.nearshore.sand[which(output.basic.nearshore.sand$time %in% days.obs), which(colnames(output.basic.nearshore.sand) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore.sand.cover = 1 - (sum_diff.total / tss_rem.total)
  
  #project model to offshore
  curr.site = 'off'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  output.basic.offshore.transfer.sand = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                                  days.model.offshore, SIR.sand.cover, c(b = beta.nearshore.sand, g = gamma.nearshore.sand,
                                                                              N = N.offshore,
                                                                              l = lambda.nearshore,
                                                                              C = cover.offshore)))
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = output.basic.offshore.transfer.sand[which(output.basic.offshore.transfer.sand$time %in% days.obs), which(colnames(output.basic.offshore.transfer.sand) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.near.to.off.basic.sand = 1 - (sum_diff.total / tss_rem.total)
  
  #project model to midchannel
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  output.basic.midchannel.transfer.sand = data.frame(ode(c(S = S.midchannel, I = I.midchannel, R = R.midchannel),
                                                       days.model.midchannel, SIR.sand.cover, c(b = beta.nearshore.sand, g = gamma.nearshore.sand,
                                                                                              N = N.midchannel,
                                                                                              l = lambda.nearshore,
                                                                                              C = cover.midchannel)))
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = output.basic.midchannel.transfer.sand[which(output.basic.midchannel.transfer.sand$time %in% days.obs), which(colnames(output.basic.midchannel.transfer.sand) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.near.to.mid.basic.sand = 1 - (sum_diff.total / tss_rem.total)
  
  r_squared.basic.nearshore.sand.no_cover
  r_squared.basic.nearshore.sand.cover
  r_squared.near.to.off.basic.sand
  r_squared.near.to.mid.basic.sand
  