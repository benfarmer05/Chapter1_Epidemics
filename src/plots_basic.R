
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
  
  # Nearshore
  site.loop = 'Nearshore'
  order = 2
  N.nearshore = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.nearshore = my.SIRS.basic[[order]]
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
  
  output.nearshore = pivot_longer(output.nearshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.nearshore)[1] = 'days.model'
  colnames(output.nearshore)[3] = 'tissue'
  
  p.fit.nearshore.basic = ggplot(data = output.nearshore, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.nearshore$days.model), y = min(output.nearshore$tissue)*0.7, label = list(tab.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia')

  p.S.fit.nearshore.basic = ggplot(data = output.nearshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic = ggplot(data = output.nearshore %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.basic = ggplot(data = output.nearshore %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Midchannel
  site.loop = 'Midchannel'
  order = 1
  N.midchannel = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.midchannel = my.SIRS.basic[[order]]
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
  
  output.midchannel = pivot_longer(output.midchannel, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.midchannel)[1] = 'days.model'
  colnames(output.midchannel)[3] = 'tissue'
  
  p.fit.midchannel.basic = ggplot(data = output.midchannel, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.midchannel$days.model), y = min(output.midchannel$tissue)*0.7, label = list(tab.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel.basic = ggplot(data = output.midchannel %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.basic = ggplot(data = output.midchannel %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.basic = ggplot(data = output.midchannel %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Offshore
  site.loop = 'Offshore'
  order = 3
  N.offshore = susceptible_ref %>%
    filter(Site == site.loop) %>%
    slice(1) %>% #all values for N.site are the same between categories, so slice first row
    pull(N.site)
  
  output.offshore = my.SIRS.basic[[order]]
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
  
  output.offshore = pivot_longer(output.offshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.offshore)[1] = 'days.model'
  colnames(output.offshore)[3] = 'tissue'
  
  p.fit.offshore.basic = ggplot(data = output.offshore, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.offshore$days.model), y = min(output.offshore$tissue)*0.7, label = list(tab.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore.basic = ggplot(data = output.offshore %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.basic = ggplot(data = output.offshore %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.basic = ggplot(data = output.offshore %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  ############## TEST FOR DHW ##############
  #
  
  #Nearshore
  site.loop = 'Nearshore'
  order = 2
  
  output.nearshore.DHW = my.SIRS.basic.DHW[[order]]
  params.basic.nearshore.DHW = params.basic.DHW[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.nearshore = params.basic.nearshore.DHW[1]
  eta.nearshore = params.basic.nearshore.DHW[2]
  
  tab.nearshore = tibble(round(beta.nearshore, 2), round(beta.nearshore.adj, 2), round(gamma.nearshore, 2),
                          round(R0.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  output.nearshore.DHW = pivot_longer(output.nearshore.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.nearshore.DHW)[1] = 'days.model'
  colnames(output.nearshore.DHW)[3] = 'tissue'
  
  p.fit.nearshore.basic.DHW = ggplot(data = output.nearshore.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.nearshore.DHW$days.model), y = min(output.nearshore.DHW$tissue)*0.7, label = list(tab.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.nearshore.basic.DHW = ggplot(data = output.nearshore.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.basic.DHW = ggplot(data = output.nearshore.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.basic.DHW = ggplot(data = output.nearshore.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  #Midchannel
  site.loop = 'Midchannel'
  order = 1
  
  output.midchannel.DHW = my.SIRS.basic.DHW[[order]]
  params.basic.midchannel.DHW = params.basic.DHW[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.midchannel = params.basic.midchannel.DHW[1]
  eta.midchannel = params.basic.midchannel.DHW[2]
  
  tab.midchannel = tibble(round(beta.midchannel, 2), round(beta.midchannel.adj, 2), round(gamma.midchannel, 2),
                          round(R0.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  output.midchannel.DHW = pivot_longer(output.midchannel.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.midchannel.DHW)[1] = 'days.model'
  colnames(output.midchannel.DHW)[3] = 'tissue'
  
  p.fit.midchannel.basic.DHW = ggplot(data = output.midchannel.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.midchannel.DHW$days.model), y = min(output.midchannel.DHW$tissue)*0.7, label = list(tab.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel.basic.DHW = ggplot(data = output.midchannel.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.basic.DHW = ggplot(data = output.midchannel.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.basic.DHW = ggplot(data = output.midchannel.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  
  #Offshore
  site.loop = 'Offshore'
  order = 3
  
  output.offshore.DHW = my.SIRS.basic.DHW[[order]]
  params.basic.offshore.DHW = params.basic.DHW[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  zeta.offshore = params.basic.offshore.DHW[1]
  eta.offshore = params.basic.offshore.DHW[2]
  
  tab.offshore = tibble(round(beta.offshore, 2), round(beta.offshore.adj, 2), round(gamma.offshore, 2),
                         round(R0.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  output.offshore.DHW = pivot_longer(output.offshore.DHW, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.offshore.DHW)[1] = 'days.model'
  colnames(output.offshore.DHW)[3] = 'tissue'
  
  p.fit.offshore.basic.DHW = ggplot(data = output.offshore.DHW, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.offshore.DHW$days.model), y = min(output.offshore.DHW$tissue)*0.7, label = list(tab.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore.basic.DHW = ggplot(data = output.offshore.DHW %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.basic.DHW = ggplot(data = output.offshore.DHW %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.basic.DHW = ggplot(data = output.offshore.DHW %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
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
  
  #
  ############## TEST FOR DHW ##############
  
  
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
  
  # Fit only [lines are simulated]
  (p.fit.offshore.basic | p.fit.midchannel.basic | p.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.S.fit.offshore.basic | p.S.fit.midchannel.basic | p.S.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.I.fit.offshore.basic | p.I.fit.midchannel.basic | p.I.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.D.fit.offshore.basic | p.D.fit.midchannel.basic | p.D.fit.nearshore.basic) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  ##############################################################################################################################
  # PREDICT OUTBREAKS FROM INITIAL CONDTIONS
  
  days.model.offshore = unique(output.offshore$days.model)
  days.model.midchannel = unique(output.midchannel$days.model)
  days.model.nearshore = unique(output.nearshore$days.model)
  
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
  I.offshore = inftiss.offshore[1]
  # I.offshore = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  S.offshore = N.offshore - I.offshore
  R.offshore = 0
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.offshore.transfer = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                            days.model.offshore, SIR, c(b = beta.nearshore, g = gamma.nearshore,
                                                N = N.offshore,
                                                l = lambda.nearshore,
                                                C = cover.offshore)))

  output.offshore.transfer = pivot_longer(output.offshore.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.offshore.transfer)[1] = 'days.model'
  colnames(output.offshore.transfer)[3] = 'tissue'

  p.fit.offshore.transfer.basic = ggplot(data = output.offshore.transfer, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia')
  
  beta.offshore.transfer = beta.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.offshore))) + offset)
  # beta.offshore.transfer = (beta.nearshore * (lambda.nearshore * (1-exp(-130*(cover.offshore)))))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.offshore)))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.nearshore)))
  R0.offshore.transfer = beta.offshore.transfer / gamma.nearshore
  tab.offshore.transfer = tibble(round(beta.nearshore, 2), round(beta.offshore.transfer, 2), round(gamma.nearshore, 2), round(R0.offshore.transfer, 2), round(cover.offshore*100, 2))
  names(tab.offshore.transfer) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')

  p.fit.offshore.transfer.basic = p.fit.offshore.transfer.basic +
    annotate(geom = "table", x = min(p.fit.offshore.transfer.basic$data$days.model), y = min(p.fit.offshore.transfer.basic$data$tissue)*0.7, label = list(tab.offshore.transfer),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia')
    # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.offshore.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)

  p.S.fit.offshore.transfer.basic = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  p.I.fit.offshore.transfer.basic = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')

  p.D.fit.offshore.transfer.basic = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
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
  output.midchannel.transfer = data.frame(ode(c(S = S.midchannel, I = I.midchannel, R = R.midchannel),
                                            days.model.midchannel, SIR, c(b = beta.nearshore, g = gamma.nearshore,
                                                         N = N.midchannel,
                                                         l = lambda.nearshore,
                                                         C = cover.midchannel)))
  
  output.midchannel.transfer = pivot_longer(output.midchannel.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.midchannel.transfer)[1] = 'days.model'
  colnames(output.midchannel.transfer)[3] = 'tissue'
  
  p.fit.midchannel.transfer.basic = ggplot(data = output.midchannel.transfer, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Predicted'), collapse="")) +
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
  
  p.fit.midchannel.transfer.basic = p.fit.midchannel.transfer.basic +
    annotate(geom = "table", x = min(p.fit.midchannel.transfer.basic$data$days.model), y = min(p.fit.midchannel.transfer.basic$data$tissue)*0.7, label = list(tab.midchannel.transfer),
             family = 'Georgia', hjust = val.hjust, vjust = val.vjust)
  # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.midchannel.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.midchannel.transfer.basic = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.transfer.basic = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.transfer.basic = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
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
  output.nearshore.transfer = data.frame(ode(c(S = S.nearshore, I = I.nearshore, R = R.nearshore),
                                              days.model.nearshore, SIR, c(b = beta.offshore, g = gamma.offshore,
                                                                            N = N.nearshore,
                                                                            l = lambda.offshore,
                                                                            C = cover.nearshore)))
  
  output.nearshore.transfer = pivot_longer(output.nearshore.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.nearshore.transfer)[1] = 'days.model'
  colnames(output.nearshore.transfer)[3] = 'tissue'
  
  p.fit.nearshore.transfer.basic = ggplot(data = output.nearshore.transfer, aes(days.model, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site.loop, ' - Predicted'), collapse="")) +
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
  
  p.fit.nearshore.transfer.basic = p.fit.nearshore.transfer.basic +
    annotate(geom = "table", x = min(p.fit.nearshore.transfer.basic$data$days.model), y = min(p.fit.nearshore.transfer.basic$data$tissue)*0.7, label = list(tab.nearshore.transfer),
             family = 'Georgia', hjust = val.hjust, vjust = val.vjust)
  # geom_text(aes(x = max(days)*0.9, y = max(tissue)*0.9, label = paste("R0:", R0.nearshore.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.nearshore.transfer.basic = ggplot(data = output.nearshore.transfer %>% filter(Compartment == "Susceptible"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore.transfer.basic = ggplot(data = output.nearshore.transfer %>% filter(Compartment == "Infected"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Infected"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore.transfer.basic = ggplot(data = output.nearshore.transfer %>% filter(Compartment == "Dead"), aes(days.model, tissue)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site.loop), collapse="")) +
    geom_line() +
    geom_point(data = obs.total %>% filter(Site == site.loop, Compartment == "Dead"), aes(days.inf.site, tissue)) +
    theme_classic(base_family = 'Georgia')
  
  # Combine the plots
  nearshore.to.offshore.basic = (p.fit.nearshore.basic | p.fit.offshore.basic | p.fit.offshore.transfer.basic) + plot_layout(guides = "collect",
    axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  nearshore.to.midchannel.basic = (p.fit.nearshore.basic | p.fit.midchannel.basic | p.fit.midchannel.transfer.basic) + plot_layout(guides = "collect",
    axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  nearshore.to.offshore.basic = (p.fit.offshore.basic | p.fit.nearshore.basic | p.fit.nearshore.transfer.basic) + plot_layout(guides = "collect",
    axis_titles = 'collect') &
    theme(legend.position = 'bottom') #& xlim(0, 325)
  
  #pass workspace to downstream script
  save.image(file = here("output", "plots_basic_workspace.RData"))
  