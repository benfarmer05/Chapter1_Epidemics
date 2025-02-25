
  # .rs.restartR(clean = TRUE)
  rm(list=ls())

  library(here)
  library(tidyverse)
  library(ggplot2)
  library(ggpmisc)
  library(patchwork)
  library(deSolve)

  #import workspace from upstream script
  load(here("output/basic_SIR_workspace.RData"))
  # load(here("output/basic_SIR_workspace_abs.RData"))
  # load(here("output/basic_SIR_workspace_abs_correct.RData"))

  ################################## Set-up ##################################
  
  curr.host = 'Single-host'
  
  obs.total = obs %>%
    filter(Category == 'Total') %>%
    group_by(Site, Compartment, Category) %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    ungroup()
  obs = obs %>%
    filter(Category != "Total") %>%
    group_by(Site, Compartment, Category) %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    ungroup()
  
  #ensure all 'geom_text' portions of plots are in the same font as theme_classic ('Arial' in this case)
  theme_set(theme_classic(base_family = "Arial"))
  update_geom_defaults("text", list(colour = "black", family = theme_get()$text$family))
  # display.brewer.all(colorblindFriendly = TRUE)
  
  #for plotting parameter value output in a tibble table format
  val.hjust = 0
  val.vjust = -4
  
  susceptible_ref = susceptible_ref %>%
    mutate(Site = case_when(
      Site == "off" ~ "Offshore",
      Site == 'mid' ~ 'Midchannel',
      Site == 'near' ~ 'Nearshore'
    ))
  
  
  ################################## Fitted pre-thermal stress ##################################
  
  # Nearshore
  site.loop = 'Nearshore'
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 2
  
  N.nearshore = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.basic.nearshore = my.SIRS.basic[[order]]
  params.basic.nearshore = params.basic[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.nearshore = params.basic.nearshore[1]
  beta.nearshore.adj = params.basic.nearshore[2]
  gamma.nearshore = params.basic.nearshore[3]
  lambda.nearshore = params.basic.nearshore[4]
  R0.nearshore = params.basic.nearshore[5]
  cover.nearshore = params.basic.nearshore[6]
  
  tab.nearshore = tibble(round(beta.nearshore, 2), round(beta.nearshore.adj, 2), round(gamma.nearshore, 2),
                          round(R0.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  
  
  
  
  
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Pre-heat'
  
  sim.rem.total = output.basic.nearshore[which(output.basic.nearshore$time %in% days.obs), which(colnames(output.basic.nearshore) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.nearshore,
                              R_squared))
  
  output.basic.nearshore = pivot_longer(output.basic.nearshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.nearshore)[1] = 'days.model'
  colnames(output.basic.nearshore)[3] = 'tissue'
  
  p.fit.nearshore.basic = ggplot(data = output.basic.nearshore, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.nearshore$days.model), y = min(output.basic.nearshore$tissue)*0.7, label = list(tab.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia')

  p.S.fit.nearshore.basic = ggplot(data = output.basic.nearshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic = ggplot(data = output.basic.nearshore %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.basic = ggplot(data = output.basic.nearshore %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Midchannel
  site.loop = 'Midchannel'
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 1
  
  N.midchannel = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.basic.midchannel = my.SIRS.basic[[order]]
  params.basic.midchannel = params.basic[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.midchannel = params.basic.midchannel[1]
  beta.midchannel.adj = params.basic.midchannel[2]
  gamma.midchannel = params.basic.midchannel[3]
  lambda.midchannel = params.basic.midchannel[4]
  R0.midchannel = params.basic.midchannel[5]
  cover.midchannel = params.basic.midchannel[6]
  
  tab.midchannel = tibble(round(beta.midchannel, 2), round(beta.midchannel.adj, 2), round(gamma.midchannel, 2),
                          round(R0.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Pre-heat'
  
  sim.rem.total = output.basic.midchannel[which(output.basic.midchannel$time %in% days.obs), which(colnames(output.basic.midchannel) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.nearshore,
                              R_squared))
  
  
  
  
  
  output.basic.midchannel = pivot_longer(output.basic.midchannel, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.midchannel)[1] = 'days.model'
  colnames(output.basic.midchannel)[3] = 'tissue'
  
  p.fit.midchannel.basic = ggplot(data = output.basic.midchannel, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.midchannel$days.model), y = min(output.basic.midchannel$tissue)*0.7, label = list(tab.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel.basic = ggplot(data = output.basic.midchannel %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.basic = ggplot(data = output.basic.midchannel %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.basic = ggplot(data = output.basic.midchannel %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Offshore
  site.loop = 'Offshore'
  curr.site = 'off'
  days.obs = days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist()
  
  order = 3
  
  N.offshore = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.basic.offshore = my.SIRS.basic[[order]]
  params.basic.offshore = params.basic[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.offshore = params.basic.offshore[1]
  beta.offshore.adj = params.basic.offshore[2]
  gamma.offshore = params.basic.offshore[3]
  lambda.offshore = params.basic.offshore[4]
  R0.offshore = params.basic.offshore[5]
  cover.offshore = params.basic.offshore[6]
  
  tab.offshore = tibble(round(beta.offshore, 2), round(beta.offshore.adj, 2), round(gamma.offshore, 2),
                          round(R0.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  
  
  
  
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Pre-heat'
  
  sim.rem.total = output.basic.offshore[which(output.basic.offshore$time %in% days.obs), which(colnames(output.basic.offshore) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.offshore = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.offshore,
                              R_squared))
  
  output.basic.offshore = pivot_longer(output.basic.offshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.offshore)[1] = 'days.model'
  colnames(output.basic.offshore)[3] = 'tissue'
  
  p.fit.offshore.basic = ggplot(data = output.basic.offshore, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.offshore$days.model), y = min(output.basic.offshore$tissue)*0.7, label = list(tab.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore.basic = ggplot(data = output.basic.offshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.basic = ggplot(data = output.basic.offshore %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.basic = ggplot(data = output.basic.offshore %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  ################################## Fitted whole-outbreak ##################################
  # Nearshore
  site.loop = 'Nearshore'
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 2
  N.nearshore = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.basic.nearshore.full = my.SIRS.basic.full[[order]]
  params.basic.nearshore.full = params.basic.full[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.nearshore.full = params.basic.nearshore.full[1]
  beta.nearshore.adj.full = params.basic.nearshore.full[2]
  gamma.nearshore.full = params.basic.nearshore.full[3]
  lambda.nearshore.full = params.basic.nearshore.full[4]
  R0.nearshore.full = params.basic.nearshore.full[5]
  cover.nearshore.full = params.basic.nearshore.full[6]
  
  tab.nearshore.full = tibble(round(beta.nearshore.full, 2), round(beta.nearshore.adj.full, 2), round(gamma.nearshore.full, 2),
                         round(R0.nearshore.full, 2), round(cover.nearshore.full*100, 2))
  names(tab.nearshore.full) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'Fitted'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.nearshore.full[which(output.basic.nearshore.full$time %in% days.obs), which(colnames(output.basic.nearshore.full) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore.full = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.nearshore.full,
                              R_squared))
  
  output.basic.nearshore.full = pivot_longer(output.basic.nearshore.full, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.nearshore.full)[1] = 'days.model'
  colnames(output.basic.nearshore.full)[3] = 'tissue'
  
  p.fit.nearshore.basic.full = ggplot(data = output.basic.nearshore.full, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.nearshore.full$days.model), y = min(output.basic.nearshore.full$tissue)*0.7, label = list(tab.nearshore.full),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia')
  
  p.S.fit.nearshore.basic.full = ggplot(data = output.basic.nearshore.full %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic.full = ggplot(data = output.basic.nearshore.full %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.basic.full = ggplot(data = output.basic.nearshore.full %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Midchannel
  site.loop = 'Midchannel'
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 1
  
  N.midchannel = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.basic.midchannel.full = my.SIRS.basic.full[[order]]
  params.basic.midchannel.full = params.basic.full[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.midchannel.full = params.basic.midchannel.full[1]
  beta.midchannel.adj.full = params.basic.midchannel.full[2]
  gamma.midchannel.full = params.basic.midchannel.full[3]
  lambda.midchannel.full = params.basic.midchannel.full[4]
  R0.midchannel.full = params.basic.midchannel.full[5]
  cover.midchannel.full = params.basic.midchannel.full[6]
  
  tab.midchannel.full = tibble(round(beta.midchannel.full, 2), round(beta.midchannel.adj.full, 2), round(gamma.midchannel.full, 2),
                          round(R0.midchannel.full, 2), round(cover.midchannel.full*100, 2))
  names(tab.midchannel.full) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'Fitted'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.midchannel.full[which(output.basic.midchannel.full$time %in% days.obs), which(colnames(output.basic.midchannel.full) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.midchannel.full = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.midchannel.full,
                              R_squared))
  
  output.basic.midchannel.full = pivot_longer(output.basic.midchannel.full, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.midchannel.full)[1] = 'days.model'
  colnames(output.basic.midchannel.full)[3] = 'tissue'
  
  p.fit.midchannel.basic.full = ggplot(data = output.basic.midchannel.full, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.midchannel.full$days.model), y = min(output.basic.midchannel.full$tissue)*0.7, label = list(tab.midchannel.full),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel.basic.full = ggplot(data = output.basic.midchannel.full %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.basic.full = ggplot(data = output.basic.midchannel.full %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.basic.full = ggplot(data = output.basic.midchannel.full %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Offshore
  site.loop = 'Offshore'
  curr.site = 'off'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 3
  
  N.offshore = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.basic.offshore.full = my.SIRS.basic.full[[order]]
  params.basic.offshore.full = params.basic.full[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.offshore.full = params.basic.offshore.full[1]
  beta.offshore.adj.full = params.basic.offshore.full[2]
  gamma.offshore.full = params.basic.offshore.full[3]
  lambda.offshore.full = params.basic.offshore.full[4]
  R0.offshore.full = params.basic.offshore.full[5]
  cover.offshore.full = params.basic.offshore.full[6]
  
  tab.offshore.full = tibble(round(beta.offshore.full, 2), round(beta.offshore.adj.full, 2), round(gamma.offshore.full, 2),
                        round(R0.offshore.full, 2), round(cover.offshore.full*100, 2))
  names(tab.offshore.full) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'Fitted'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.offshore.full[which(output.basic.offshore.full$time %in% days.obs), which(colnames(output.basic.offshore.full) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.offshore.full = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.offshore.full,
                              R_squared))
  
  output.basic.offshore.full = pivot_longer(output.basic.offshore.full, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.offshore.full)[1] = 'days.model'
  colnames(output.basic.offshore.full)[3] = 'tissue'
  
  p.fit.offshore.basic.full = ggplot(data = output.basic.offshore.full, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.offshore.full$days.model), y = min(output.basic.offshore.full$tissue)*0.7, label = list(tab.offshore.full),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore.basic.full = ggplot(data = output.basic.offshore.full %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.basic.full = ggplot(data = output.basic.offshore.full %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.basic.full = ggplot(data = output.basic.offshore.full %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  ################################## Thermal stress  ##################################
  #Nearshore
  site.loop = 'Nearshore'
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 2
  
  output.basic.nearshore.DHW = my.SIRS.basic.DHW[[order]]
  params.basic.nearshore.DHW = params.basic.DHW[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.nearshore = params.basic.nearshore.DHW[1]
  eta.nearshore = params.basic.nearshore.DHW[2]
  
  tab.nearshore = tibble(round(beta.nearshore, 2), round(beta.nearshore.adj, 2), round(gamma.nearshore, 2),
                          round(R0.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.nearshore.DHW[which(output.basic.nearshore.DHW$time %in% days.obs), which(colnames(output.basic.nearshore.DHW) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.nearshore.DHW = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.nearshore.DHW,
                              R_squared))
  
  output.basic.nearshore.DHW = pivot_longer(output.basic.nearshore.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.nearshore.DHW)[1] = 'days.model'
  colnames(output.basic.nearshore.DHW)[3] = 'tissue'
  
  p.fit.nearshore.basic.DHW = ggplot(data = output.basic.nearshore.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.nearshore.DHW$days.model), y = min(output.basic.nearshore.DHW$tissue)*0.7, label = list(tab.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.nearshore.basic.DHW = ggplot(data = output.basic.nearshore.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic.DHW = ggplot(data = output.basic.nearshore.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.basic.DHW = ggplot(data = output.basic.nearshore.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #Midchannel
  site.loop = 'Midchannel'
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 1
  
  output.basic.midchannel.DHW = my.SIRS.basic.DHW[[order]]
  params.basic.midchannel.DHW = params.basic.DHW[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.midchannel = params.basic.midchannel.DHW[1]
  eta.midchannel = params.basic.midchannel.DHW[2]
  
  tab.midchannel = tibble(round(beta.midchannel, 2), round(beta.midchannel.adj, 2), round(gamma.midchannel, 2),
                          round(R0.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.midchannel.DHW[which(output.basic.midchannel.DHW$time %in% days.obs), which(colnames(output.basic.midchannel.DHW) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.midchannel.DHW = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.midchannel.DHW,
                              R_squared))
  
  output.basic.midchannel.DHW = pivot_longer(output.basic.midchannel.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.midchannel.DHW)[1] = 'days.model'
  colnames(output.basic.midchannel.DHW)[3] = 'tissue'
  
  p.fit.midchannel.basic.DHW = ggplot(data = output.basic.midchannel.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.midchannel.DHW$days.model), y = min(output.basic.midchannel.DHW$tissue)*0.7, label = list(tab.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel.basic.DHW = ggplot(data = output.basic.midchannel.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.basic.DHW = ggplot(data = output.basic.midchannel.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_line(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.basic.DHW = ggplot(data = output.basic.midchannel.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #Offshore
  site.loop = 'Offshore'
  curr.site = 'off'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 3
  
  output.basic.offshore.DHW = my.SIRS.basic.DHW[[order]]
  params.basic.offshore.DHW = params.basic.DHW[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.offshore = params.basic.offshore.DHW[1]
  eta.offshore = params.basic.offshore.DHW[2]
  
  tab.offshore = tibble(round(beta.offshore, 2), round(beta.offshore.adj, 2), round(gamma.offshore, 2),
                         round(R0.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  curr.type = 'DHW'
  curr.wave = 'Full'
  
  sim.rem.total = output.basic.offshore.DHW[which(output.basic.offshore.DHW$time %in% days.obs), which(colnames(output.basic.offshore.DHW) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops off zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.basic.offshore.DHW = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.basic.offshore.DHW,
                              R_squared))
  
  output.basic.offshore.DHW = pivot_longer(output.basic.offshore.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.basic.offshore.DHW)[1] = 'days.model'
  colnames(output.basic.offshore.DHW)[3] = 'tissue'
  
  p.fit.offshore.basic.DHW = ggplot(data = output.basic.offshore.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.basic.offshore.DHW$days.model), y = min(output.basic.offshore.DHW$tissue)*0.7, label = list(tab.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore.basic.DHW = ggplot(data = output.basic.offshore.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.basic.DHW = ggplot(data = output.basic.offshore.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.basic.DHW = ggplot(data = output.basic.offshore.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic
  p.I.fit.nearshore.basic.DHW
  p.D.fit.nearshore.basic
  p.D.fit.nearshore.basic.DHW
  
  p.I.fit.midchannel.basic
  p.I.fit.midchannel.basic.DHW
  p.D.fit.midchannel.basic
  p.D.fit.midchannel.basic.DHW
  
  p.I.fit.offshore.basic
  p.I.fit.offshore.basic.DHW
  p.D.fit.offshore.basic
  p.D.fit.offshore.basic.DHW
  
  p.fit.nearshore.basic
  p.fit.nearshore.basic.DHW
  p.fit.midchannel.basic
  p.fit.midchannel.basic.DHW
  p.fit.offshore.basic
  p.fit.offshore.basic.DHW
  
  ################################## Project outbreaks  ##################################
  curr.host <- "Single-host"
  curr.type <- "Projected"
  curr.wave = "Full"
  
  days.model.offshore = unique(output.basic.offshore$days.model)
  days.model.midchannel = unique(output.basic.midchannel$days.model)
  days.model.nearshore = unique(output.basic.nearshore$days.model)
  
  site.loop = 'Offshore'
  inftiss.offshore = obs.total %>%
    filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
    pull(tissue) %>%
    # {trimmed <- .[cumsum(!is.na(.) & . != 0) > 0]  # Trim leading NA or zero values
    # head(trimmed, length(trimmed) - DHW.modifier)}  # Apply DHW modifier
    { .[cumsum(!is.na(.) & . != 0) > 0] }
  
  site.loop = 'Midchannel'
  inftiss.midchannel = obs.total %>%
    filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
    pull(tissue) %>%
    { .[cumsum(!is.na(.) & . != 0) > 0] }
  
  site.loop = 'Nearshore'
  inftiss.nearshore = obs.total %>%
    filter(Site == site.loop, Category == "Total", Compartment == "Infected") %>%
    pull(tissue) %>%
    { .[cumsum(!is.na(.) & . != 0) > 0] }
  
  #run SIR for Offshore based on fit from Nearshore
  site.loop = 'Offshore'
  curr.site = 'off'
  curr.type = 'Projected'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  I.offshore = inftiss.offshore[1]
  # I.offshore = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.offshore = N.offshore - I.offshore
  R.offshore = 0
  
  
  
  
  # # TESTING

  lambda = 1.0

  offset = 1 - 1 / (1 + exp(-lambda * 1.0))
  # # TESTING
  
  
  
  
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.basic.offshore.transfer = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                            days.model.offshore, SIR, c(b = beta.nearshore, g = gamma.nearshore,# NOTE - g = gamma.nearshore or 0.60 / 0.66 to run a test, b = beta.nearshore or .82 / 0.89 to run a test 6 dec 2024
                                                N = N.offshore,
                                                l = lambda.nearshore,
                                                C = cover.offshore)))

  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = output.basic.offshore.transfer[which(output.basic.offshore.transfer$time %in% days.obs), which(colnames(output.basic.offshore.transfer) %in% 'R')]
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  first_valid_idx <- which(!is.na(days))[1] #find the first non-NA index
  obs.rem.total <- obs.rem.total[first_valid_idx:length(obs.rem.total)]
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.near.to.off.basic = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
  mutate(R_squared = ifelse(site == curr.site & 
                              host == curr.host & 
                              type == curr.type & 
                              wave == curr.wave,
                            r_squared.near.to.off.basic,
                            R_squared))
  
  output.basic.offshore.transfer = pivot_longer(output.basic.offshore.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.basic.offshore.transfer)[1] = 'days.model'
  colnames(output.basic.offshore.transfer)[3] = 'tissue'
  
  p.fit.near.to.off.basic = ggplot(data = output.basic.offshore.transfer, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia')
  
  beta.offshore.transfer = beta.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.offshore))) + offset)
  
  # #NOTE - test! 6 dec 2024
  # gamma.offshore.transfer = gamma.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.offshore))) + offset)
  # 
  # beta.offshore.transfer = .82 * (1 / (1 + exp(-lambda.modifier * (cover.offshore))) + offset)
  
  
  # beta.offshore.transfer = (beta.nearshore * (lambda.nearshore * (1-exp(-130*(cover.offshore)))))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.offshore)))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.nearshore)))
  R0.offshore.transfer = beta.offshore.transfer / gamma.nearshore
  tab.offshore.transfer = tibble(round(beta.nearshore, 2), round(beta.offshore.transfer, 2), round(gamma.nearshore, 2), round(R0.offshore.transfer, 2), round(cover.offshore*100, 2))
  names(tab.offshore.transfer) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')

  p.fit.near.to.off.basic = p.fit.near.to.off.basic +
    annotate(geom = "table", x = min(p.fit.near.to.off.basic$data$days.model), y = min(p.fit.near.to.off.basic$data$tissue)*0.7, label = list(tab.offshore.transfer),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia')
    # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.offshore.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)

  p.S.fit.near.to.off.basic = ggplot(data = output.basic.offshore.transfer %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  p.I.fit.near.to.off.basic = ggplot(data = output.basic.offshore.transfer %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  p.D.fit.near.to.off.basic = ggplot(data = output.basic.offshore.transfer %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  #run SIR for Midchannel based on fit from Nearshore
  site.loop = 'Midchannel'
  I.midchannel = inftiss.midchannel[1]
  # I.midchannel = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.midchannel = N.midchannel - I.midchannel
  R.midchannel = 0
  
  #simulation using initial state variables from naive site.loop and parameters from fitted site.loop
  output.basic.midchannel.transfer = data.frame(ode(c(S = S.midchannel, I = I.midchannel, R = R.midchannel),
                                            days.model.midchannel, SIR, c(b = beta.nearshore, g = gamma.nearshore,
                                                         N = N.midchannel,
                                                         l = lambda.nearshore,
                                                         C = cover.midchannel)))
  
  output.basic.midchannel.transfer = pivot_longer(output.basic.midchannel.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.basic.midchannel.transfer)[1] = 'days.model'
  colnames(output.basic.midchannel.transfer)[3] = 'tissue'
  
  p.fit.near.to.mid.basic = ggplot(data = output.basic.midchannel.transfer, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia') #+
    # xlim(0, 325)
  
  beta.midchannel.transfer = beta.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.midchannel))) + offset)
  # beta.midchannel.transfer = (beta.nearshore * (lambda.nearshore * (1-exp(-130*(cover.midchannel)))))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.midchannel)))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.nearshore)))
  R0.midchannel.transfer = beta.midchannel.transfer / gamma.nearshore
  tab.midchannel.transfer = tibble(round(beta.nearshore, 2), round(beta.midchannel.transfer, 2), round(gamma.nearshore, 2), round(R0.midchannel.transfer, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel.transfer) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')
  
  p.fit.near.to.mid.basic = p.fit.near.to.mid.basic +
    annotate(geom = "table", x = min(p.fit.near.to.mid.basic$data$days.model), y = min(p.fit.near.to.mid.basic$data$tissue)*0.7, label = list(tab.midchannel.transfer),
             family = 'Georgia', hjust = val.hjust, vjust = val.vjust)
  # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.midchannel.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.near.to.mid.basic = ggplot(data = output.basic.midchannel.transfer %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.near.to.mid.basic = ggplot(data = output.basic.midchannel.transfer %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.near.to.mid.basic = ggplot(data = output.basic.midchannel.transfer %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #run SIR for Nearshore based on fit from Offshore
  site.loop = 'Nearshore'
  I.nearshore = inftiss.nearshore[1]
  # I.nearshore = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.nearshore = N.nearshore - I.nearshore
  R.nearshore = 0
  
  #simulation using initial state variables from naive site.loop and parameters from fitted site.loop
  output.basic.nearshore.transfer = data.frame(ode(c(S = S.nearshore, I = I.nearshore, R = R.nearshore),
                                              days.model.nearshore, SIR, c(b = beta.offshore, g = gamma.offshore,
                                                                            N = N.nearshore,
                                                                            l = lambda.offshore,
                                                                            C = cover.nearshore)))
  
  output.basic.nearshore.transfer = pivot_longer(output.basic.nearshore.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.basic.nearshore.transfer)[1] = 'days.model'
  colnames(output.basic.nearshore.transfer)[3] = 'tissue'
  
  p.fit.off.to.near.basic = ggplot(data = output.basic.nearshore.transfer, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia') #+
  # xlim(0, 325)
  
  beta.nearshore.transfer = beta.offshore * (1 / (1 + exp(-lambda.modifier * (cover.nearshore))) + offset)
  # beta.nearshore.transfer = (beta.offshore * (lambda.offshore * (1-exp(-130*(cover.nearshore)))))
  # (beta.offshore * (1 + lambda.offshore * sqrt(cover.nearshore)))
  # (beta.offshore * (1 + lambda.offshore * sqrt(cover.offshore)))
  R0.nearshore.transfer = beta.nearshore.transfer / gamma.offshore
  tab.nearshore.transfer = tibble(round(beta.offshore, 2), round(beta.nearshore.transfer, 2), round(gamma.offshore, 2), round(R0.nearshore.transfer, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore.transfer) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')
  
  p.fit.off.to.near.basic = p.fit.off.to.near.basic +
    annotate(geom = "table", x = min(p.fit.off.to.near.basic$data$days.model), y = min(p.fit.off.to.near.basic$data$tissue)*0.7, label = list(tab.nearshore.transfer),
             family = 'Georgia', hjust = val.hjust, vjust = val.vjust)
  # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.nearshore.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.off.to.near.basic = ggplot(data = output.basic.nearshore.transfer %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.off.to.near.basic = ggplot(data = output.basic.nearshore.transfer %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.off.to.near.basic = ggplot(data = output.basic.nearshore.transfer %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Combine the plots
  # nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.fit.offshore.basic | p.fit.offshore.transfer.basic) + plot_layout(guides = "collect",
  #   axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  # 
  # nearshore.to.midchannel.basic = (p.fit.nearshore.basic | p.fit.midchannel.basic | p.fit.midchannel.transfer.basic) + plot_layout(guides = "collect",
  #   axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  # 
  # offshore.to.nearshore.basic = (p.fit.offshore.basic | p.fit.nearshore.basic | p.fit.nearshore.transfer.basic) + plot_layout(guides = "collect",
  #   axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  
  nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.offshore.basic | p.I.fit.offshore.basic | p.D.fit.offshore.basic) / (p.fit.near.to.off.basic | p.I.fit.near.to.off.basic | p.D.fit.near.to.off.basic) + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                                                         axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  # nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.offshore.basic | p.I.fit.offshore.basic | p.D.fit.offshore.basic) / (p.fit.near.to.off.basic | p.I.fit.near.to.off.basic | p.D.fit.near.to.off.basic) + plot_layout(guides = "collect") &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  
  nearshore.to.midchannel.basic = (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.midchannel.basic | p.I.fit.midchannel.basic | p.D.fit.midchannel.basic) / (p.fit.near.to.mid.basic | p.I.fit.near.to.mid.basic | p.D.fit.near.to.mid.basic)  + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                                                                  axis_titles = 'collect') &
    theme(legend.position = 'bottom') # & xlim(0, 325)
  
  offshore.to.nearshore.basic = (p.fit.offshore.basic | p.I.fit.offshore.basic | p.D.fit.offshore.basic) / (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.off.to.near.basic | p.I.fit.off.to.near.basic | p.D.fit.off.to.near.basic)  + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                                                          axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  ################################## Thermal projections  ##################################
  site.loop = 'Offshore'

  output.basic.offshore.transfer.DHW = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                                      days.model.offshore, SIR.DHW, c(b = beta.nearshore, g = gamma.nearshore,# NOTE - g = gamma.nearshore or 0.60 / 0.66 to run a test, b = beta.nearshore or .82 / 0.89 to run a test 6 dec 2024
                                                                                      z = zeta.nearshore,
                                                                                      e = eta.nearshore,
                                                                                      N = N.offshore,
                                                                                      l = lambda.nearshore,
                                                                                      C = cover.offshore,
                                                                                      SST_threshold = SST_threshold_value,
                                                                                      DHW_threshold = DHW_threshold_value),
                                                      SST = SST_df,
                                                      DHW = DHW_df))

  output.basic.offshore.transfer.DHW = pivot_longer(output.basic.offshore.transfer.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.basic.offshore.transfer.DHW)[1] = 'days.model'
  colnames(output.basic.offshore.transfer.DHW)[3] = 'tissue'

  p.fit.near.to.off.basic.DHW = ggplot(data = output.basic.offshore.transfer.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia')
  
  beta.offshore.transfer.DHW = beta.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.offshore))) + offset)
  
  R0.offshore.transfer.DHW = beta.offshore.transfer.DHW / gamma.nearshore
  tab.offshore.transfer.DHW = tibble(round(beta.nearshore, 2), round(beta.offshore.transfer.DHW, 2), round(gamma.nearshore, 2), round(R0.offshore.transfer.DHW, 2), round(cover.offshore*100, 2))
  names(tab.offshore.transfer.DHW) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')
  
  p.fit.near.to.off.basic.DHW = p.fit.near.to.off.basic.DHW +
    annotate(geom = "table", x = min(p.fit.near.to.off.basic.DHW$data$days.model), y = min(p.fit.near.to.off.basic.DHW$data$tissue)*0.7, label = list(tab.offshore.transfer.DHW ),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia')
  
  p.S.fit.near.to.off.basic.DHW = ggplot(data = output.basic.offshore.transfer.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.near.to.off.basic.DHW = ggplot(data = output.basic.offshore.transfer.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.near.to.off.basic.DHW = ggplot(data = output.basic.offshore.transfer.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #run SIR for Midchannel based on fit from Nearshore
  site.loop = 'Midchannel'
  
  output.basic.midchannel.transfer.DHW = data.frame(ode(c(S = S.midchannel, I = I.midchannel, R = R.midchannel),
                                                      days.model.midchannel, SIR.DHW, c(b = beta.nearshore, g = gamma.nearshore,# NOTE - g = gamma.nearshore or 0.60 / 0.66 to run a test, b = beta.nearshore or .82 / 0.89 to run a test 6 dec 2024
                                                                                      z = zeta.midchannel,
                                                                                      e = eta.midchannel,
                                                                                      N = N.midchannel,
                                                                                      l = lambda.nearshore,
                                                                                      C = cover.midchannel,
                                                                                      SST_threshold = SST_threshold_value,
                                                                                      DHW_threshold = DHW_threshold_value),
                                                      SST = SST_df,
                                                      DHW = DHW_df))
  
  output.basic.midchannel.transfer.DHW = pivot_longer(output.basic.midchannel.transfer.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.basic.midchannel.transfer.DHW)[1] = 'days.model'
  colnames(output.basic.midchannel.transfer.DHW)[3] = 'tissue'
  
  p.fit.near.to.mid.basic.DHW = ggplot(data = output.basic.midchannel.transfer.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia') #+
  # xlim(0, 325)
  
  beta.midchannel.transfer.DHW = beta.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.midchannel))) + offset)
  # beta.midchannel.transfer = (beta.nearshore * (lambda.nearshore * (1-exp(-130*(cover.midchannel)))))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.midchannel)))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.nearshore)))
  R0.midchannel.transfer.DHW = beta.midchannel.transfer / gamma.nearshore
  tab.midchannel.transfer.DHW = tibble(round(beta.nearshore, 2), round(beta.midchannel.transfer.DHW, 2), round(gamma.nearshore, 2), round(R0.midchannel.transfer.DHW, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel.transfer.DHW) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')
  
  p.fit.near.to.mid.basic.DHW = p.fit.near.to.mid.basic.DHW +
    annotate(geom = "table", x = min(p.fit.near.to.mid.basic.DHW$data$days.model), y = min(p.fit.near.to.mid.basic.DHW$data$tissue)*0.7, label = list(tab.midchannel.transfer.DHW),
             family = 'Georgia', hjust = val.hjust, vjust = val.vjust)
  # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.midchannel.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.near.to.mid.basic.DHW = ggplot(data = output.basic.midchannel.transfer.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.near.to.mid.basic.DHW = ggplot(data = output.basic.midchannel.transfer.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.near.to.mid.basic.DHW = ggplot(data = output.basic.midchannel.transfer.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  #run SIR for Nearshore based on fit from Offshore
  site.loop = 'Nearshore'
  
  output.basic.nearshore.transfer.DHW = data.frame(ode(c(S = S.nearshore, I = I.nearshore, R = R.nearshore),
                                                      days.model.nearshore, SIR.DHW, c(b = beta.offshore, g = gamma.offshore,# NOTE - g = gamma.offshore or 0.60 / 0.66 to run a test, b = beta.offshore or .82 / 0.89 to run a test 6 dec 2024
                                                                                      z = zeta.offshore,
                                                                                      e = eta.offshore,
                                                                                      N = N.nearshore,
                                                                                      l = lambda.offshore,
                                                                                      C = cover.nearshore,
                                                                                      SST_threshold = SST_threshold_value,
                                                                                      DHW_threshold = DHW_threshold_value),
                                                      SST = SST_df,
                                                      DHW = DHW_df))
  
  output.basic.nearshore.transfer.DHW = pivot_longer(output.basic.nearshore.transfer.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.basic.nearshore.transfer.DHW)[1] = 'days.model'
  colnames(output.basic.nearshore.transfer.DHW)[3] = 'tissue'
  
  
  # STOPPING POINT - 7 feb 2025
  #   - for some reason this transfer just doesn't seem to be working correctly.
  #   - but also, not sure it matters because these transfers don't appear to be useful when bringing in thermal stress
  #   - so, we probably leave it out of the paper or plot it and explain why it might not have worked well (math ? )
  #
  #   - the next problem is the multi-host models ... their fit really is not good when limiting it to pre-thermal stress
  #       - so, especially given thermal stress doesn't even work with the single-host projections, I say we simply go back and fit
  #         the multi-host model to the *entire* outbreak and remember to simply frame that as a separate exercise
  #         in the manuscript. it would be lovely if host type, projections, and DHWs all worked together, but that might just
  #         be asking too much and we can justify presenting things separately (just with a note that synthesizing the approaches
  #         had mixed results)
  #
  #   - so ... last steps really are:
  #       1.) making sure if I want to stick with R-squared or absolute residual optimization
  #       2.) whichever I choose for #1, make sure it looks right in the output
  #       3.) triple check that calculation of R-squared and optimization is being done correctly
  #       4.) choose how to present effect of thermal stress (within existing figures, also maybe separate fine detail of it into fig. 4)
  #       5.) as appropriate, embed R-squared values into figures
  #       6.) work with Gaby to wrap up fig. 1
  #       7.) work with Dan to finalize typesetting on figs. 2-3 (or 2-4)
  #       8.) update manuscript accordingly to all above. ensure any necessary references to SIR methods are included
  
  
  p.fit.off.to.near.basic.DHW = ggplot(data = output.basic.nearshore.transfer.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia')
  
  beta.nearshore.transfer.DHW = beta.offshore * (1 / (1 + exp(-lambda.modifier * (cover.nearshore))) + offset)
  
  R0.nearshore.transfer.DHW = beta.nearshore.transfer.DHW / gamma.offshore
  tab.nearshore.transfer.DHW = tibble(round(beta.offshore, 2), round(beta.nearshore.transfer.DHW, 2), round(gamma.offshore, 2), round(R0.nearshore.transfer.DHW, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore.transfer.DHW) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')
  
  p.fit.near.to.off.basic.DHW = p.fit.near.to.off.basic.DHW +
    annotate(geom = "table", x = min(p.fit.near.to.off.basic.DHW$data$days.model), y = min(p.fit.near.to.off.basic.DHW$data$tissue)*0.7, label = list(tab.nearshore.transfer.DHW ),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia')
  
  p.S.fit.near.to.off.basic.DHW = ggplot(data = output.basic.nearshore.transfer.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.near.to.off.basic.DHW = ggplot(data = output.basic.nearshore.transfer.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.near.to.off.basic.DHW = ggplot(data = output.basic.nearshore.transfer.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop, '- Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  # Combine the plots
  # nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.fit.offshore.basic | p.fit.offshore.transfer.basic) + plot_layout(guides = "collect",
  #   axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  #
  # nearshore.to.midchannel.basic = (p.fit.nearshore.basic | p.fit.midchannel.basic | p.fit.midchannel.transfer.basic) + plot_layout(guides = "collect",
  #   axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  #
  # offshore.to.nearshore.basic = (p.fit.offshore.basic | p.fit.nearshore.basic | p.fit.nearshore.transfer.basic) + plot_layout(guides = "collect",
  #   axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)

  nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.offshore.basic | p.I.fit.offshore.basic | p.D.fit.offshore.basic) / (p.fit.near.to.off.basic | p.I.fit.near.to.off.basic | p.D.fit.near.to.off.basic) + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                                                         axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  # nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.offshore.basic | p.I.fit.offshore.basic | p.D.fit.offshore.basic) / (p.fit.near.to.off.basic | p.I.fit.near.to.off.basic | p.D.fit.near.to.off.basic) + plot_layout(guides = "collect") &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)

  nearshore.to.midchannel.basic = (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.midchannel.basic | p.I.fit.midchannel.basic | p.D.fit.midchannel.basic) / (p.fit.near.to.mid.basic | p.I.fit.near.to.mid.basic | p.D.fit.near.to.mid.basic)  + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                                                                  axis_titles = 'collect') &
    theme(legend.position = 'bottom') # & xlim(0, 325)

  offshore.to.nearshore.basic = (p.fit.offshore.basic | p.I.fit.offshore.basic | p.D.fit.offshore.basic) / (p.fit.nearshore.basic | p.I.fit.nearshore.basic | p.D.fit.nearshore.basic) / (p.fit.off.to.near.basic | p.I.fit.off.to.near.basic | p.D.fit.off.to.near.basic)  + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                                                          axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  
  ################################## Plots ##################################
  # # Observations only [lines are observations]
  # #overlaid
  # ggplot(data = obs, aes(days.inf.site, tissue, colour = Compartment, linetype = Category, shape = Site)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of tissue") +
  #   ggtitle(paste(c("", 'All sites'), collapse="")) +
  #   geom_line() +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   theme_classic(base_family = 'Georgia')
  
  #facet-wrapped observations
  (p.SIR.offshore.basic | p.SIR.midchannel.basic | p.SIR.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom') #&
  # scale_color_brewer(name = 'Disease compartment', labels = c("High", "Low", "Medium"), palette = 'Dark2')
  
  # Infection observations only [lines are observations]
  (p.I.offshore.basic | p.I.midchannel.basic | p.I.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  #whole outbreak simulations
  (p.fit.offshore.basic | p.fit.midchannel.basic | p.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.fit.offshore.basic.full | p.fit.midchannel.basic.full | p.fit.nearshore.basic.full) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.fit.offshore.basic.DHW | p.fit.midchannel.basic.DHW | p.fit.nearshore.basic.DHW) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  #projection vs fitted
  (p.fit.offshore.basic.full | p.fit.near.to.off.basic)
  (p.fit.nearshore.basic.full | p.fit.off.to.near.basic)
  (p.fit.midchannel.basic.full | p.fit.near.to.mid.basic)
  
  #susceptible compartment
  (p.S.fit.offshore.basic | p.S.fit.midchannel.basic | p.S.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.S.fit.offshore.basic.full | p.S.fit.midchannel.basic.full | p.S.fit.nearshore.basic.full) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.S.fit.offshore.basic.DHW | p.S.fit.midchannel.basic.DHW | p.S.fit.nearshore.basic.DHW) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  #infected compartment
  (p.I.fit.offshore.basic | p.I.fit.midchannel.basic | p.I.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.I.fit.offshore.basic.full | p.I.fit.midchannel.basic.full | p.I.fit.nearshore.basic.full) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.I.fit.offshore.basic.DHW | p.I.fit.midchannel.basic.DHW | p.I.fit.nearshore.basic.DHW) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.I.fit.near.to.off.basic | p.I.fit.off.to.near.basic)
  
  #dead compartment
  (p.D.fit.offshore.basic | p.D.fit.midchannel.basic | p.D.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.D.fit.offshore.basic.full | p.D.fit.midchannel.basic.full | p.D.fit.nearshore.basic.full) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.D.fit.offshore.basic.DHW | p.D.fit.midchannel.basic.DHW | p.D.fit.nearshore.basic.DHW) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.D.fit.near.to.off.basic | p.D.fit.off.to.near.basic)
  
  # STOPPING POINT - 11 FEB 2025
  #   - need to add R-squared to projected fit. also, decide if want to project from pre-thermal stress, thermal stress, or whole-outbreak fit
  #   - then need to do the same for multi-host model
  
  ################################## Save workspace  ##################################
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "plots_basic_workspace.RData"))
  