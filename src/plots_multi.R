  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(ggplot2)
  library(tidyverse)
  library(ggpmisc)
  library(patchwork)
  library(deSolve)
  
  #import workspace from upstream script
  load(here("output/multi_SIR_workspace.RData"))
  # load(here("output/multi_SIR_workspace_lower_start.RData"))
  
  ################################## Set-up ##################################
  
  # load(here("output/pre_DHW_integration_and_cover_ratio_9Oct2024/multi_SIR_workspace_01_125_15.RData"))
  
  # obs.multi = obs.multi %>%
  #   mutate(Category = case_when(
  #     Category == "Low" ~ "LS",
  #     Category == 'Moderate' ~ 'MS',
  #     Category == 'High' ~ 'HS'
  #   ))
  
  # NOTE - should ideally go back upstream in the scripts to rename Category to Susceptibility. does not affect
  #         function of the code, though
  obs.multi = obs.multi %>%
    rename(Susceptibility = Category)
  
  suscat.names = c('Low', 'Moderate', 'High')
  
  susceptible_ref = susceptible_ref %>%
    mutate(Site = case_when(
      Site == "off" ~ "Offshore",
      Site == 'mid' ~ 'Midchannel',
      Site == 'near' ~ 'Nearshore'
    ))
  
  #ensure all 'geom_text' portions of plots are in the same font as theme_classic ('Arial' in this case)
  theme_set(theme_classic(base_family = "Arial"))
  update_geom_defaults("text", list(colour = "black", family = theme_get()$text$family))
  # display.brewer.all(colorblindFriendly = TRUE)
  
  ################################## Fitted whole-outbreak ##################################
  
  # Nearshore
  site.loop = 'Nearshore'
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 2
  
  N.LS.nearshore = susceptible_ref %>%
    filter(Site == site.loop, Category == 'Low') %>%
    pull(tissue_ref)
  N.MS.nearshore = susceptible_ref %>%
    filter(Site == site.loop, Category == 'Moderate') %>%
    pull(tissue_ref)
  N.HS.nearshore = susceptible_ref %>%
    filter(Site == site.loop, Category == 'High') %>%
    pull(tissue_ref)
  
  output.raw.nearshore = my.SIRS.multi[[order]]
  output.nearshore = output.raw.nearshore
  params.nearshore = params.multi[[order]]
  
  cover.nearshore.LS = params.nearshore[14]
  cover.nearshore.MS = params.nearshore[15]
  cover.nearshore.HS = params.nearshore[16]
  cover.nearshore = c(cover.nearshore.LS, cover.nearshore.MS, cover.nearshore.HS)
  
  beta.nearshore.LS = params.nearshore[1]
  beta.nearshore.LS.adj = params.nearshore[2]
  gamma.nearshore.LS = params.nearshore[3]
  R0.nearshore.LS = params.nearshore[10]
  
  beta.nearshore.MS = params.nearshore[4]
  beta.nearshore.MS.adj = params.nearshore[5]
  gamma.nearshore.MS = params.nearshore[6]
  R0.nearshore.MS = params.nearshore[11]
  
  beta.nearshore.HS = params.nearshore[7]
  beta.nearshore.HS.adj = params.nearshore[8]
  gamma.nearshore.HS = params.nearshore[9]
  R0.nearshore.HS = params.nearshore[12]
  
  betas.nearshore = c(beta.nearshore.LS, beta.nearshore.MS, beta.nearshore.HS)
  betas.nearshore.adj = c(beta.nearshore.LS.adj, beta.nearshore.MS.adj, beta.nearshore.HS.adj)
  gammas.nearshore = c(gamma.nearshore.LS, gamma.nearshore.MS, gamma.nearshore.HS)
  R0s.nearshore = c(R0.nearshore.LS, R0.nearshore.MS, R0.nearshore.HS)
  
  tab.nearshore = tibble(suscat.names, round(betas.nearshore, 2), round(betas.nearshore.adj, 2), round(gammas.nearshore, 2),
                         round(R0s.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = rowSums(output.nearshore[which(output.nearshore$time %in% days.obs), 
                                           which(colnames(output.nearshore) %in% c('R.LS', 'R.MS', 'R.HS'))], 
                          na.rm = TRUE)
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
  r_squared.multi.nearshore = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.multi.nearshore,
                              R_squared))
  
  output.nearshore = output.nearshore %>% select(-last_col())
  
  output.nearshore = pivot_longer(output.nearshore, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.nearshore)[1] = 'days.model'
  colnames(output.nearshore)[4] = 'tissue'
  
  p.fit.nearshore.multi = ggplot(data = output.nearshore, aes(days.model, tissue, colour = Susceptibility, linetype = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Susceptibility, shape = Compartment)) +
    scale_color_brewer(name = 'Susceptibility', palette = 'Set2') +
    annotate(geom = "table", x = max(output.nearshore$days.model), y = max(output.nearshore$tissue)*0.7, label = list(tab.nearshore),
             vjust = 1, hjust = 1) +
    theme_classic()
  
  p.S.fit.nearshore.multi = ggplot(data = output.nearshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of live tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.I.fit.nearshore.multi = ggplot(data = output.nearshore %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic() + 
    guides(color = 'none')
  
  p.D.fit.nearshore.multi = ggplot(data = output.nearshore %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic() +
    guides(color = 'none')
  
  # Midchannel
  site.loop = 'Midchannel'
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 1
  
  N.LS.midchannel = susceptible_ref %>%
    filter(Site == site.loop, Category == 'Low') %>%
    pull(tissue_ref)
  N.MS.midchannel = susceptible_ref %>%
    filter(Site == site.loop, Category == 'Moderate') %>%
    pull(tissue_ref)
  N.HS.midchannel = susceptible_ref %>%
    filter(Site == site.loop, Category == 'High') %>%
    pull(tissue_ref)
  
  output.raw.midchannel = my.SIRS.multi[[order]]
  output.midchannel = output.raw.midchannel
  params.midchannel = params.multi[[order]]
  
  cover.midchannel.LS = params.midchannel[14]
  cover.midchannel.MS = params.midchannel[15]
  cover.midchannel.HS = params.midchannel[16]
  cover.midchannel = c(cover.midchannel.LS, cover.midchannel.MS, cover.midchannel.HS)
  
  beta.midchannel.LS = params.midchannel[1]
  beta.midchannel.LS.adj = params.midchannel[2]
  gamma.midchannel.LS = params.midchannel[3]
  R0.midchannel.LS = params.midchannel[10]
  
  beta.midchannel.MS = params.midchannel[4]
  beta.midchannel.MS.adj = params.midchannel[5]
  gamma.midchannel.MS = params.midchannel[6]
  R0.midchannel.MS = params.midchannel[11]
  
  beta.midchannel.HS = params.midchannel[7]
  beta.midchannel.HS.adj = params.midchannel[8]
  gamma.midchannel.HS = params.midchannel[9]
  R0.midchannel.HS = params.midchannel[12]
  
  betas.midchannel = c(beta.midchannel.LS, beta.midchannel.MS, beta.midchannel.HS)
  betas.midchannel.adj = c(beta.midchannel.LS.adj, beta.midchannel.MS.adj, beta.midchannel.HS.adj)
  gammas.midchannel = c(gamma.midchannel.LS, gamma.midchannel.MS, gamma.midchannel.HS)
  R0s.midchannel = c(R0.midchannel.LS, R0.midchannel.MS, R0.midchannel.HS)
  
  tab.midchannel = tibble(suscat.names, round(betas.midchannel, 2), round(betas.midchannel.adj, 2), round(gammas.midchannel, 2),
                         round(R0s.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = rowSums(output.midchannel[which(output.midchannel$time %in% days.obs), 
                                           which(colnames(output.midchannel) %in% c('R.LS', 'R.MS', 'R.HS'))], 
                          na.rm = TRUE)
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
  r_squared.multi.midchannel = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.multi.midchannel,
                              R_squared))
  
  output.midchannel = output.midchannel %>% select(-last_col())
  
  output.midchannel = pivot_longer(output.midchannel, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.midchannel)[1] = 'days.model'
  colnames(output.midchannel)[4] = 'tissue'
  
  p.fit.midchannel.multi = ggplot(data = output.midchannel, aes(days.model, tissue, colour = Susceptibility, linetype = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Susceptibility, shape = Compartment)) +
    scale_color_brewer(name = 'Susceptibility', palette = 'Set2') +
    annotate(geom = "table", x = max(output.midchannel$days.model), y = max(output.midchannel$tissue)*0.7, label = list(tab.midchannel),
             vjust = 1, hjust = 1) +
    theme_classic()
  
  p.S.fit.midchannel.multi = ggplot(data = output.midchannel %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of live tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.I.fit.midchannel.multi = ggplot(data = output.midchannel %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic() + 
    guides(color = 'none')
  
  p.D.fit.midchannel.multi = ggplot(data = output.midchannel %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic() +
    guides(color = 'none')  
  
  # Offshore
  site.loop = 'Offshore'
  curr.site = 'off'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  order = 3
  
  N.LS.offshore = susceptible_ref %>%
    filter(Site == site.loop, Category == 'Low') %>%
    pull(tissue_ref)
  N.MS.offshore = susceptible_ref %>%
    filter(Site == site.loop, Category == 'Moderate') %>%
    pull(tissue_ref)
  N.HS.offshore = susceptible_ref %>%
    filter(Site == site.loop, Category == 'High') %>%
    pull(tissue_ref)
  
  output.raw.offshore = my.SIRS.multi[[order]]
  output.offshore = output.raw.offshore
  params.offshore = params.multi[[order]]
  
  cover.offshore.LS = params.offshore[14]
  cover.offshore.MS = params.offshore[15]
  cover.offshore.HS = params.offshore[16]
  cover.offshore = c(cover.offshore.LS, cover.offshore.MS, cover.offshore.HS)
  
  beta.offshore.LS = params.offshore[1]
  beta.offshore.LS.adj = params.offshore[2]
  gamma.offshore.LS = params.offshore[3]
  R0.offshore.LS = params.offshore[10]
  
  beta.offshore.MS = params.offshore[4]
  beta.offshore.MS.adj = params.offshore[5]
  gamma.offshore.MS = params.offshore[6]
  R0.offshore.MS = params.offshore[11]
  
  beta.offshore.HS = params.offshore[7]
  beta.offshore.HS.adj = params.offshore[8]
  gamma.offshore.HS = params.offshore[9]
  R0.offshore.HS = params.offshore[12]
  
  betas.offshore = c(beta.offshore.LS, beta.offshore.MS, beta.offshore.HS)
  betas.offshore.adj = c(beta.offshore.LS.adj, beta.offshore.MS.adj, beta.offshore.HS.adj)
  gammas.offshore = c(gamma.offshore.LS, gamma.offshore.MS, gamma.offshore.HS)
  R0s.offshore = c(R0.offshore.LS, R0.offshore.MS, R0.offshore.HS)
  
  tab.offshore = tibble(suscat.names, round(betas.offshore, 2), round(betas.offshore.adj, 2), round(gammas.offshore, 2),
                          round(R0s.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = rowSums(output.offshore[which(output.offshore$time %in% days.obs), 
                                            which(colnames(output.offshore) %in% c('R.LS', 'R.MS', 'R.HS'))], 
                          na.rm = TRUE)
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
  r_squared.multi.offshore = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.multi.offshore,
                              R_squared))
  
  
  output.offshore = output.offshore %>% select(-last_col())
  
  output.offshore = pivot_longer(output.offshore, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.offshore)[1] = 'days.model'
  colnames(output.offshore)[4] = 'tissue'
  
  p.fit.offshore.multi = ggplot(data = output.offshore, aes(days.model, tissue, colour = Susceptibility, linetype = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Susceptibility, shape = Compartment)) +
    scale_color_brewer(name = 'Susceptibility', palette = 'Set2') +
    annotate(geom = "table", x = max(output.offshore$days.model), y = max(output.offshore$tissue)*0.7, label = list(tab.offshore),
             vjust = 1, hjust = 1) +
    theme_classic()
  
  p.S.fit.offshore.multi = ggplot(data = output.offshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of live tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  # note - experimented with color mapping here
  p.I.fit.offshore.multi = ggplot(data = output.offshore %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    # geom_line(color = '#fc8e62') +
    # geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility), color = '#fc8e62') +
    theme_classic() + 
    guides(color = 'none')
  
  p.D.fit.offshore.multi = ggplot(data = output.offshore %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    # geom_line(color = '#66c2a4') +
    # geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility), color = '#66c2a4') +
    theme_classic() +
    guides(color = 'none')  
  
  # p.fit.offshore;p.I.fit.offshore;p.D.fit.offshore
  
  ################################## Project outbreaks  ##################################
  
  site.loop = 'Offshore'
  
  # Calculate the indices to keep based on HS.inftiss.offshore (the high-susceptibility corals kick off the epidemic at every site)
  HS.indices = obs %>%
    filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
    pull(tissue) %>%
    { cumsum(!is.na(.) & . != 0) > 0 }  # Logical vector of indices to keep
  
  # Apply the same indices to LS and MS
  HS.inftiss.offshore = obs %>%
    filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]
  LS.inftiss.offshore = obs %>%
    filter(Site == site.loop, Category == "Low", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]
  MS.inftiss.offshore = obs %>%
    filter(Site == site.loop, Category == "Moderate", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]  
  
  site.loop = 'Midchannel'
  HS.indices = obs %>%
    filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
    pull(tissue) %>%
    { cumsum(!is.na(.) & . != 0) > 0 }  # Logical vector of indices to keep
  HS.inftiss.midchannel = obs %>%
    filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]
  LS.inftiss.midchannel = obs %>%
    filter(Site == site.loop, Category == "Low", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]
  MS.inftiss.midchannel = obs %>%
    filter(Site == site.loop, Category == "Moderate", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]  
  
  site.loop = 'Nearshore'
  HS.indices = obs %>%
    filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
    pull(tissue) %>%
    { cumsum(!is.na(.) & . != 0) > 0 }  # Logical vector of indices to keep
  HS.inftiss.nearshore = obs %>%
    filter(Site == site.loop, Category == "High", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]
  LS.inftiss.nearshore = obs %>%
    filter(Site == site.loop, Category == "Low", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]
  MS.inftiss.nearshore = obs %>%
    filter(Site == site.loop, Category == "Moderate", Compartment == "Infected") %>%
    pull(tissue) %>%
    .[HS.indices]  
  
  #run SIR for Offshore based on fit from Nearshore
  site.loop = 'Offshore'
  curr.site = 'off'
  curr.type = 'Projected'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 
  
  # TESTING
  # HS.inftiss.offshore[1] = HS.inftiss.offshore[1] / 10
  # TESTING
  
  I.LS.offshore = LS.inftiss.offshore[1]
  S.LS.offshore = N.LS.offshore - I.LS.offshore
  R.LS.offshore = 0
  I.MS.offshore = MS.inftiss.offshore[1]
  S.MS.offshore = N.MS.offshore - I.MS.offshore
  R.MS.offshore = 0
  I.HS.offshore = HS.inftiss.offshore[1]
  S.HS.offshore = N.HS.offshore - I.HS.offshore 
  R.HS.offshore = 0
  
  P.offshore = I.LS.offshore + I.MS.offshore + I.HS.offshore
  
  # # TESTING
  # days.model.offshore <- c(days.model.offshore, max(days.model.offshore) + 1:max(days.model.offshore) + 100)
  # # TESTING
  
  # # TESTING
  # lambda.LS = 150 #0.1
  # lambda.MS = 50 #1.25
  # lambda.HS = 50 #15
  # lambda.LS = 0 #0.1
  # lambda.MS = 0 #1.25
  # lambda.HS = 0 #15
  # lambda.LS = 0.1
  # lambda.MS = 1.25
  # lambda.HS = 15
  # lambda.LS = 0
  # lambda.MS = 0
  # lambda.HS = 11
  lambda.LS = 0
  lambda.MS = 0.2
  lambda.HS = 13
  
  lambda.LS = 0
  lambda.MS = 0
  lambda.HS = 0
  
  offset.LS = 1 - 1 / (1 + exp(-lambda.LS * 1.0))
  offset.MS = 1 - 1 / (1 + exp(-lambda.MS * 1.0))
  offset.HS = 1 - 1 / (1 + exp(-lambda.HS * 1.0))
  
  # #null conditions
  # lambda.LS = 1
  # lambda.MS = 1
  # lambda.HS = 1
  # offset.LS = 0
  # offset.MS = 0
  # offset.HS = 0
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.near.to.off.multi = data.frame(ode(c(S.LS = S.LS.offshore, I.LS = I.LS.offshore, R.LS = R.LS.offshore,
                                                S.MS = S.MS.offshore, I.MS = I.MS.offshore, R.MS = R.MS.offshore,
                                                S.HS = S.HS.offshore, I.HS = I.HS.offshore, R.HS = R.HS.offshore),
                                              days.model.offshore, SIR.multi, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                                           b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                                           b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                                           N.LS = N.LS.offshore, N.MS = N.MS.offshore, N.HS = N.HS.offshore,
                                                           C = cover.offshore,
                                                           C.LS = cover.offshore.LS, C.MS = cover.offshore.MS, C.HS = cover.offshore.HS,
                                                           l = lambda)))
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = rowSums(output.near.to.off.multi[which(output.near.to.off.multi$time %in% days.obs), 
                                           which(colnames(output.near.to.off.multi) %in% c('R.LS', 'R.MS', 'R.HS'))], 
                          na.rm = TRUE)
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
  r_squared.near.to.off.multi = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.near.to.off.multi,
                              R_squared))
  
  
  # NOTE / STOPPING POINT - 8 OCT 2024
  #   - it may be a real problem that lambda is not passed as a parameter from the multi-group SIR. if we end up wanting to "fit" lambda, 
  #       refactoring will need to be done in the SIR loop and downstream. I guess for now, can focus on getting this running, and if I
  #       ever require lambda to be fitted, can cross the bridge later
  
  output.near.to.off.multi = output.near.to.off.multi %>% select(-last_col())
  
  output.near.to.off.multi = pivot_longer(output.near.to.off.multi, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.near.to.off.multi)[1] = 'days.model'
  colnames(output.near.to.off.multi)[4] = 'tissue'
  
  # # TESTING
  # days.model <- c(days.model, max(days.model) + 1:max(days.model) + 100)
  # # TESTING
  
  p.fit.near.to.off.multi = ggplot(data = output.near.to.off.multi, aes(days.model, tissue, colour = Compartment, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment, shape = Susceptibility)) +
    scale_color_brewer(name = 'Compartment', palette = 'Set2') +
    theme_classic()
  
  beta.offshore.transfer.LS = beta.nearshore.LS * (1 / (1 + exp(-lambda.LS * (cover.offshore.LS))) + offset.LS)
  R0.offshore.transfer.LS = beta.offshore.transfer.LS / gamma.nearshore.LS
  beta.offshore.transfer.MS = beta.nearshore.MS * (1 / (1 + exp(-lambda.MS * (cover.offshore.MS))) + offset.MS)
  R0.offshore.transfer.MS = beta.offshore.transfer.MS / gamma.nearshore.MS
  beta.offshore.transfer.HS = beta.nearshore.HS * (1 / (1 + exp(-lambda.HS * (cover.offshore.HS))) + offset.HS)
  R0.offshore.transfer.HS = beta.offshore.transfer.HS / gamma.nearshore.HS
  
  betas.offshore.transfer = c(beta.offshore.transfer.LS, beta.offshore.transfer.MS, beta.offshore.transfer.HS)
  R0s.offshore.transfer = c(R0.offshore.transfer.LS, R0.offshore.transfer.MS, R0.offshore.transfer.HS)
  tab.offshore.transfer = tibble(suscat.names, round(betas.nearshore, 2), round(betas.offshore.transfer, 2), round(gammas.nearshore, 2), round(R0s.offshore.transfer, 2), round(cover.offshore*100, 2))
  names(tab.offshore.transfer) = c('Susceptibility', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.near.to.off.multi = p.fit.near.to.off.multi +
    annotate(geom = "table", x = max(p.fit.near.to.off.multi$data$days.model), y = max(p.fit.near.to.off.multi$data$tissue)*0.7, label = list(tab.offshore.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.near.to.off.multi = ggplot(data = output.near.to.off.multi %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.I.fit.near.to.off.multi = ggplot(data = output.near.to.off.multi %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.D.fit.near.to.off.multi = ggplot(data = output.near.to.off.multi %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab('Day of observation period') +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  #run SIR for Midchannel based on fit from Nearshore
  site.loop = 'Midchannel'
  curr.site = 'mid'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 

  I.LS.midchannel = LS.inftiss.midchannel[1]
  S.LS.midchannel = N.LS.midchannel - I.LS.midchannel
  R.LS.midchannel = 0
  I.MS.midchannel = MS.inftiss.midchannel[1]
  S.MS.midchannel = N.MS.midchannel - I.MS.midchannel
  R.MS.midchannel = 0
  I.HS.midchannel = HS.inftiss.midchannel[1]
  S.HS.midchannel = N.HS.midchannel - I.HS.midchannel 
  R.HS.midchannel = 0
  
  P.midchannel = I.LS.midchannel + I.MS.midchannel + I.HS.midchannel
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.near.to.mid.multi = data.frame(ode(c(S.LS = S.LS.midchannel, I.LS = I.LS.midchannel, R.LS = R.LS.midchannel,
                                              S.MS = S.MS.midchannel, I.MS = I.MS.midchannel, R.MS = R.MS.midchannel,
                                              # S.HS = S.HS.midchannel, I.HS = I.HS.midchannel, R.HS = R.HS.midchannel,
                                              # P = P.midchannel),
                                              S.HS = S.HS.midchannel, I.HS = I.HS.midchannel, R.HS = R.HS.midchannel),
                                            days.model.midchannel, SIR.multi, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                                                              b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                                                              b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                                                              N.LS = N.LS.midchannel, N.MS = N.MS.midchannel, N.HS = N.HS.midchannel,
                                                                              C = cover.midchannel,
                                                                              C.LS = cover.midchannel.LS, C.MS = cover.midchannel.MS, C.HS = cover.midchannel.HS,
                                                                              l = lambda)))
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = rowSums(output.near.to.mid.multi[which(output.near.to.mid.multi$time %in% days.obs), 
                                                   which(colnames(output.near.to.mid.multi) %in% c('R.LS', 'R.MS', 'R.HS'))], 
                          na.rm = TRUE)
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops mid zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.near.to.mid.multi = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.near.to.mid.multi,
                              R_squared))
  
  output.near.to.mid.multi = output.near.to.mid.multi %>% select(-last_col())
  
  output.near.to.mid.multi = pivot_longer(output.near.to.mid.multi, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.near.to.mid.multi)[1] = 'days.model'
  colnames(output.near.to.mid.multi)[4] = 'tissue'
  
  p.fit.near.to.mid.multi = ggplot(data = output.near.to.mid.multi, aes(days.model, tissue, colour = Compartment, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment, shape = Susceptibility)) +
    scale_color_brewer(name = 'Compartment', palette = 'Set2') +
    theme_classic()
  
  beta.midchannel.transfer.LS = beta.nearshore.LS * (1 / (1 + exp(-lambda.LS * (cover.midchannel.LS))) + offset.LS)
  R0.midchannel.transfer.LS = beta.midchannel.transfer.LS / gamma.nearshore.LS
  beta.midchannel.transfer.MS = beta.nearshore.MS * (1 / (1 + exp(-lambda.MS * (cover.midchannel.MS))) + offset.MS)
  R0.midchannel.transfer.MS = beta.midchannel.transfer.MS / gamma.nearshore.MS
  beta.midchannel.transfer.HS = beta.nearshore.HS * (1 / (1 + exp(-lambda.HS * (cover.midchannel.HS))) + offset.HS)
  R0.midchannel.transfer.HS = beta.midchannel.transfer.HS / gamma.nearshore.HS
  
  betas.midchannel.transfer = c(beta.midchannel.transfer.LS, beta.midchannel.transfer.MS, beta.midchannel.transfer.HS)
  R0s.midchannel.transfer = c(R0.midchannel.transfer.LS, R0.midchannel.transfer.MS, R0.midchannel.transfer.HS)
  tab.midchannel.transfer = tibble(suscat.names, round(betas.nearshore, 2), round(betas.midchannel.transfer, 2), round(gammas.nearshore, 2), round(R0s.midchannel.transfer, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel.transfer) = c('Susceptibility', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.near.to.mid.multi = p.fit.near.to.mid.multi +
    annotate(geom = "table", x = max(p.fit.near.to.mid.multi$data$days.model), y = max(p.fit.near.to.mid.multi$data$tissue)*0.7, label = list(tab.midchannel.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.near.to.mid.multi = ggplot(data = output.near.to.mid.multi %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.I.fit.near.to.mid.multi = ggplot(data = output.near.to.mid.multi %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.D.fit.near.to.mid.multi = ggplot(data = output.near.to.mid.multi %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab('Day of observation period') +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  #run SIR for Nearshore based on fit from Offshore
  site.loop = 'Nearshore'
  curr.site = 'near'
  days.obs <- days_sites %>%
    filter(site == curr.site) %>%
    pull(days.obs) %>%
    unlist() 

  I.LS.nearshore = LS.inftiss.nearshore[1]
  S.LS.nearshore = N.LS.nearshore - I.LS.nearshore
  R.LS.nearshore = 0
  I.MS.nearshore = MS.inftiss.nearshore[1]
  S.MS.nearshore = N.MS.nearshore - I.MS.nearshore
  R.MS.nearshore = 0
  I.HS.nearshore = HS.inftiss.nearshore[1]
  S.HS.nearshore = N.HS.nearshore - I.HS.nearshore 
  R.HS.nearshore = 0
  
  P.nearshore = I.LS.nearshore + I.MS.nearshore + I.HS.nearshore
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.off.to.near.multi = data.frame(ode(c(S.LS = S.LS.nearshore, I.LS = I.LS.nearshore, R.LS = R.LS.nearshore,
                                              S.MS = S.MS.nearshore, I.MS = I.MS.nearshore, R.MS = R.MS.nearshore,
                                              # S.HS = S.HS.nearshore, I.HS = I.HS.nearshore, R.HS = R.HS.nearshore,
                                              # P = P.nearshore),
                                              S.HS = S.HS.nearshore, I.HS = I.HS.nearshore, R.HS = R.HS.nearshore),
                                            days.model.nearshore, SIR.multi, c(b.LS = beta.offshore.LS, g.LS = gamma.offshore.LS,
                                                                              b.MS = beta.offshore.MS, g.MS = gamma.offshore.MS,
                                                                              b.HS = beta.offshore.HS, g.HS = gamma.offshore.HS,
                                                                              N.LS = N.LS.nearshore, N.MS = N.MS.nearshore, N.HS = N.HS.nearshore,
                                                                              C = cover.nearshore,
                                                                              C.LS = cover.nearshore.LS, C.MS = cover.nearshore.MS, C.HS = cover.nearshore.HS,
                                                                              l = lambda)))
  
  #calculate R-squared and update error table
  # NOTE - could also fill in SSR, TSS, and observations/simulated values to error table if needed
  sim.rem.total = rowSums(output.off.to.near.multi[which(output.off.to.near.multi$time %in% days.obs), 
                                                   which(colnames(output.off.to.near.multi) %in% c('R.LS', 'R.MS', 'R.HS'))], 
                          na.rm = TRUE)
  obs.rem.total = obs.model %>%
    filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
    slice(head(row_number(), n()-DHW.modifier)) %>%
    pull(tissue)
  
  # Trim obs.rem.total from the beginning to match the length of sim.rem.total. this chops mid zero's
  if (length(obs.rem.total) > length(sim.rem.total)) {
    obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  }
  
  diff.rem.total = (sim.rem.total - obs.rem.total)
  sum_diff.total = sum(diff.rem.total^2)
  mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
  tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
  r_squared.off.to.near.multi = 1 - (sum_diff.total / tss_rem.total)
  
  error_eval <- error_eval %>%
    mutate(R_squared = ifelse(site == curr.site & 
                                host == curr.host & 
                                type == curr.type & 
                                wave == curr.wave,
                              r_squared.off.to.near.multi,
                              R_squared))
  
  output.off.to.near.multi = output.off.to.near.multi %>% select(-last_col())
  
  output.off.to.near.multi = pivot_longer(output.off.to.near.multi, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.off.to.near.multi)[1] = 'days.model'
  colnames(output.off.to.near.multi)[4] = 'tissue'
  
  p.fit.off.to.near.multi = ggplot(data = output.off.to.near.multi, aes(days.model, tissue, colour = Compartment, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment, shape = Susceptibility)) +
    scale_color_brewer(name = 'Compartment', palette = 'Set2') +
    theme_classic()
  
  beta.nearshore.transfer.LS = beta.nearshore.LS * (1 / (1 + exp(-lambda.LS * (cover.nearshore.LS))) + offset.LS)
  R0.nearshore.transfer.LS = beta.nearshore.transfer.LS / gamma.nearshore.LS
  beta.nearshore.transfer.MS = beta.nearshore.MS * (1 / (1 + exp(-lambda.MS * (cover.nearshore.MS))) + offset.MS)
  R0.nearshore.transfer.MS = beta.nearshore.transfer.MS / gamma.nearshore.MS
  beta.nearshore.transfer.HS = beta.nearshore.HS * (1 / (1 + exp(-lambda.HS * (cover.nearshore.HS))) + offset.HS)
  R0.nearshore.transfer.HS = beta.nearshore.transfer.HS / gamma.nearshore.HS
  
  betas.nearshore.transfer = c(beta.nearshore.transfer.LS, beta.nearshore.transfer.MS, beta.nearshore.transfer.HS)
  R0s.nearshore.transfer = c(R0.nearshore.transfer.LS, R0.nearshore.transfer.MS, R0.nearshore.transfer.HS)
  tab.nearshore.transfer = tibble(suscat.names, round(betas.nearshore, 2), round(betas.nearshore.transfer, 2), round(gammas.nearshore, 2), round(R0s.nearshore.transfer, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore.transfer) = c('Susceptibility', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.off.to.near.multi = p.fit.off.to.near.multi +
    annotate(geom = "table", x = max(p.fit.off.to.near.multi$data$days.model), y = max(p.fit.off.to.near.multi$data$tissue)*0.7, label = list(tab.nearshore.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.off.to.near.multi = ggplot(data = output.off.to.near.multi %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.I.fit.off.to.near.multi = ggplot(data = output.off.to.near.multi %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.D.fit.off.to.near.multi = ggplot(data = output.off.to.near.multi %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab('Day of observation period') +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  #run SIR for Midchannel based on fit from Offshore
  site.loop = 'Midchannel'
  
  I.LS.midchannel = LS.inftiss.midchannel[1]
  S.LS.midchannel = N.LS.midchannel - I.LS.midchannel
  R.LS.midchannel = 0
  I.MS.midchannel = MS.inftiss.midchannel[1]
  S.MS.midchannel = N.MS.midchannel - I.MS.midchannel
  R.MS.midchannel = 0
  I.HS.midchannel = HS.inftiss.midchannel[1]
  S.HS.midchannel = N.HS.midchannel - I.HS.midchannel 
  R.HS.midchannel = 0
  
  P.midchannel = I.LS.midchannel + I.MS.midchannel + I.HS.midchannel
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.off.to.mid.multi = data.frame(ode(c(S.LS = S.LS.midchannel, I.LS = I.LS.midchannel, R.LS = R.LS.midchannel,
                                              S.MS = S.MS.midchannel, I.MS = I.MS.midchannel, R.MS = R.MS.midchannel,
                                              # S.HS = S.HS.midchannel, I.HS = I.HS.midchannel, R.HS = R.HS.midchannel,
                                              # P = P.midchannel),
                                              S.HS = S.HS.midchannel, I.HS = I.HS.midchannel, R.HS = R.HS.midchannel),
                                           days.model.midchannel, SIR.multi, c(b.LS = beta.offshore.LS, g.LS = gamma.offshore.LS,
                                                                         b.MS = beta.offshore.MS, g.MS = gamma.offshore.MS,
                                                                         b.HS = beta.offshore.HS, g.HS = gamma.offshore.HS,
                                                                         N.LS = N.LS.midchannel, N.MS = N.MS.midchannel, N.HS = N.HS.midchannel,
                                                                         C = cover.midchannel,
                                                                         C.LS = cover.midchannel.LS, C.MS = cover.midchannel.MS, C.HS = cover.midchannel.HS,
                                                                         l = lambda)))
  
  output.off.to.mid.multi = output.off.to.mid.multi %>% select(-last_col())
  
  output.off.to.mid.multi = pivot_longer(output.off.to.mid.multi, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Susceptibility")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment)))) %>%
    mutate(Susceptibility = case_when( # NOTE - same as comment above about upstream scripts
      Susceptibility == "LS" ~ "Low",
      Susceptibility == 'MS' ~ 'Moderate',
      Susceptibility == 'HS' ~ 'High'
    ))
  
  colnames(output.off.to.mid.multi)[1] = 'days.model'
  colnames(output.off.to.mid.multi)[4] = 'tissue'
  
  p.fit.off.to.mid.multi = ggplot(data = output.off.to.mid.multi, aes(days.model, tissue, colour = Compartment, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment, shape = Susceptibility)) +
    scale_color_brewer(name = 'Compartment', palette = 'Set2') +
    theme_classic()
  
  beta.midchannel.transfer.LS = beta.offshore.LS * (1 / (1 + exp(-lambda.LS * (cover.midchannel.LS))) + offset.LS)
  R0.midchannel.transfer.LS = beta.midchannel.transfer.LS / gamma.offshore.LS
  beta.midchannel.transfer.MS = beta.offshore.MS * (1 / (1 + exp(-lambda.MS * (cover.midchannel.MS))) + offset.MS)
  R0.midchannel.transfer.MS = beta.midchannel.transfer.MS / gamma.offshore.MS
  beta.midchannel.transfer.HS = beta.offshore.HS * (1 / (1 + exp(-lambda.HS * (cover.midchannel.HS))) + offset.HS)
  R0.midchannel.transfer.HS = beta.midchannel.transfer.HS / gamma.offshore.HS
  
  betas.midchannel.transfer = c(beta.midchannel.transfer.LS, beta.midchannel.transfer.MS, beta.midchannel.transfer.HS)
  R0s.midchannel.transfer = c(R0.midchannel.transfer.LS, R0.midchannel.transfer.MS, R0.midchannel.transfer.HS)
  tab.midchannel.transfer = tibble(suscat.names, round(betas.offshore, 2), round(betas.midchannel.transfer, 2), round(gammas.offshore, 2), round(R0s.midchannel.transfer, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel.transfer) = c('Susceptibility', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.off.to.mid.multi = p.fit.off.to.mid.multi +
    annotate(geom = "table", x = max(p.fit.off.to.mid.multi$data$days.model), y = max(p.fit.off.to.mid.multi$data$tissue)*0.7, label = list(tab.midchannel.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.off.to.mid.multi = ggplot(data = output.off.to.mid.multi %>% filter(Compartment == "Susceptible"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.I.fit.off.to.mid.multi = ggplot(data = output.off.to.mid.multi %>% filter(Compartment == "Infected"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  p.D.fit.off.to.mid.multi = ggplot(data = output.off.to.mid.multi %>% filter(Compartment == "Dead"), aes(days.model, tissue, linetype = Susceptibility)) +
    xlab('Surface area of infected tissue (m2)') +
    ylab("Surface area of dead tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.multi %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue, shape = Susceptibility)) +
    theme_classic()
  
  # Combine the plots
  nearshore.to.offshore.multi = (p.fit.nearshore.multi | p.I.fit.nearshore.multi | p.D.fit.nearshore.multi) / (p.fit.offshore.multi | p.I.fit.offshore.multi | p.D.fit.offshore.multi) / (p.fit.near.to.off.multi | p.I.fit.near.to.off.multi | p.D.fit.near.to.off.multi) + plot_layout(guides = "collect",
                  axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  # nearshore.to.offshore.multi = (p.fit.nearshore.multi | p.I.fit.nearshore.multi | p.D.fit.nearshore.multi) / (p.fit.offshore.multi | p.I.fit.offshore.multi | p.D.fit.offshore.multi) / (p.fit.near.to.off.multi | p.I.fit.near.to.off.multi | p.D.fit.near.to.off.multi) + plot_layout(guides = "collect") &
  #   theme(legend.position = 'bottom') #& xlim(0, 325)
  
  nearshore.to.midchannel.multi = (p.fit.nearshore.multi | p.I.fit.nearshore.multi | p.D.fit.nearshore.multi) / (p.fit.midchannel.multi | p.I.fit.midchannel.multi | p.D.fit.midchannel.multi) / (p.fit.near.to.mid.multi | p.I.fit.near.to.mid.multi | p.D.fit.near.to.mid.multi)  + plot_layout(guides = "collect",
    axis_titles = 'collect') &
    theme(legend.position = 'bottom') # & xlim(0, 325)
  
  offshore.to.nearshore.multi = (p.fit.offshore.multi | p.I.fit.offshore.multi | p.D.fit.offshore.multi) / (p.fit.nearshore.multi | p.I.fit.nearshore.multi | p.D.fit.nearshore.multi) / (p.fit.off.to.near.multi | p.I.fit.off.to.near.multi | p.D.fit.off.to.near.multi)  + plot_layout(guides = "collect",
    axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  offshore.to.midchannel.multi = (p.fit.offshore.multi | p.I.fit.offshore.multi | p.D.fit.offshore.multi) / (p.fit.midchannel.multi | p.I.fit.midchannel.multi | p.D.fit.midchannel.multi) / (p.fit.off.to.mid.multi | p.I.fit.off.to.mid.multi | p.D.fit.off.to.mid.multi)  + plot_layout(guides = "collect",
    axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  
  
  ################################## Plots  ##################################
  
  # # Observations only [lines are observations]
  # #overlaid
  # ggplot(data = obs, aes(days, tissue, colour = Compartment, linetype = Category, shape = Site)) +
  #   xlab("Day of observation period") +
  #   ylab("Surface area of tissue (m2)") +
  #   ggtitle(paste(c("", 'All sites'), collapse="")) +
  #   geom_line() +
  #   scale_color_brewer(name = 'Susceptibility', palette = 'Set2') +
  #   theme_classic()
  
  #facet-wrapped observations
  (p.SIR.offshore.multi | p.SIR.midchannel.multi | p.SIR.nearshore.multi) + plot_layout(guides = "collect") & theme(legend.position = 'bottom') #&
  # scale_color_brewer(name = 'Susceptibility', labels = c("High", "Low", "Medium"), palette = 'Dark2')
  
  # Infection observations only [lines are observations]
  (p.I.offshore.multi | p.I.midchannel.multi | p.I.nearshore.multi) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  #whole outbreak simulations
  (p.fit.offshore.multi | p.fit.midchannel.multi | p.fit.nearshore.multi) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  #projection vs fitted
  (p.fit.offshore.multi | p.fit.near.to.off.multi)
  (p.fit.nearshore.multi | p.fit.off.to.near.multi)
  (p.fit.midchannel.multi | p.fit.near.to.mid.multi)
  
  #susceptible compartment
  (p.S.fit.offshore.multi | p.S.fit.midchannel.multi | p.S.fit.nearshore.multi) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  #infected compartment
  (p.I.fit.offshore.multi | p.I.fit.midchannel.multi | p.I.fit.nearshore.multi) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.I.fit.near.to.off.multi | p.I.fit.off.to.near.multi)
  
  #dead compartment
  (p.D.fit.offshore.multi | p.D.fit.midchannel.multi | p.D.fit.nearshore.multi) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.D.fit.near.to.off.multi | p.D.fit.off.to.near.multi | p.D.fit.near.to.mid.multi)
  
  ################################## Save workspace  ##################################
  
  # #pass workspace to downstream script
  # # save.image(file = here("output", "plots_multi_workspace_betterproj.RData"))
  # save.image(file = here("output", "plots_multi_workspace.RData"))
  # # save.image(file = here("output", "plots_multi_workspace_lower_start.RData"))
  