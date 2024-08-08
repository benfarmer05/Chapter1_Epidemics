  
  rm(list=ls())
  
  load("output_multi_COVER_0-0-0_workspace.RData")
  
  library(tidyverse)
  library(patchwork)
  library(ggpmisc)
  library(ggplot2)
  library(deSolve)
  
  #ensure all 'geom_text' portions of plots are in the same font as theme_classic ('Arial' in this case)
  theme_set(theme_classic(base_family = "Arial"))
  update_geom_defaults("text", list(colour = "black", family = theme_get()$text$family))
  timepoint.midchannel = 'T7'; timepoint.offshore = 'T5'; timepoint.nearshore = 'T11'
  # display.brewer.all(colorblindFriendly = TRUE)
  
  # Nearshore
  site = 'Nearshore'
  order = 2
  
  output.nearshore = my.SIRS[[order]]
  params.nearshore = params[[order]]
  
  # order: c(min.beta.LS.tiss, min.beta.LS.tiss.adj, min.gamma.LS.tiss,
  #                 min.beta.MS.tiss, min.beta.MS.tiss.adj, min.gamma.HS.tiss,
  #                 min.beta.HS.tiss, min.beta.MS.tiss.adj, min.gamma.HS.tiss,
  #                 R0.LS, R0.MS, R0.HS,
  #                 cover.site,
  #                 cover.LS.site, cover.MS.site, cover.HS.site)
  
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
  
  suscat.names = c('LS', 'MS', 'HS')
  tab.nearshore = tibble(suscat.names, round(betas.nearshore, 2), round(betas.nearshore.adj, 2), round(gammas.nearshore, 2), round(R0s.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
  
  output.nearshore = output.nearshore %>% select(-last_col())
  
  output.nearshore = pivot_longer(output.nearshore, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.nearshore)[1] = 'days'
  colnames(output.nearshore)[4] = 'prop'
  
  
  # note - swapped colour and linetype here
  p.fit.nearshore = ggplot(data = output.nearshore, aes(days, prop, colour = Category, linetype = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Category, shape = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = max(output.nearshore$days), y = max(output.nearshore$prop)*0.7, label = list(tab.nearshore),
             vjust = 1, hjust = 1) +
    theme_classic()
  
  p.S.fit.nearshore = ggplot(data = output.nearshore %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of live tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  # note - experimented with color mapping here
  p.I.fit.nearshore = ggplot(data = output.nearshore %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line(color = '#fc8e62') +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category), color = '#fc8e62') +
    theme_classic() + 
    guides(color = 'none')
  
  p.D.fit.nearshore = ggplot(data = output.nearshore %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line(color = '#66c2a4') +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category), color = '#66c2a4') +
    theme_classic() +
    guides(color = 'none')
  
  # p.fit.nearshore;p.I.fit.nearshore;p.D.fit.nearshore
  
  # Midchannel
  site = 'Midchannel'
  order = 1
  
  output.midchannel = my.SIRS[[order]]
  params.midchannel = params[[order]]
  
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
  
  suscat.names = c('LS', 'MS', 'HS')
  tab.midchannel = tibble(suscat.names, round(betas.midchannel, 2), round(betas.midchannel.adj, 2), round(gammas.midchannel, 2), round(R0s.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
  
  output.midchannel = output.midchannel %>% select(-last_col())
  
  output.midchannel = pivot_longer(output.midchannel, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.midchannel)[1] = 'days'
  colnames(output.midchannel)[4] = 'prop'
  
  p.fit.midchannel = ggplot(data = output.midchannel, aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Compartment, shape = Category)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = max(output.midchannel$days), y = max(output.midchannel$prop)*0.7, label = list(tab.midchannel),
             vjust = 1, hjust = 1) +
    theme_classic()
  
  p.S.fit.midchannel = ggplot(data = output.midchannel %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of live tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.I.fit.midchannel = ggplot(data = output.midchannel %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line(color = '#fc8e62') +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category), color = '#fc8e62') +
    theme_classic()
  
  p.D.fit.midchannel = ggplot(data = output.midchannel %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line(color = '#66c2a4') +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category), color = '#66c2a4') +
    theme_classic()
  
  p.fit.midchannel;p.I.fit.midchannel;p.D.fit.midchannel
  
  
  # Offshore
  site = 'Offshore'
  order = 3
  
  output.offshore = my.SIRS[[order]]
  params.offshore = params[[order]]
  
  # orderr: c(min.beta.LS.tiss, min.beta.LS.tiss.adj, min.gamma.LS.tiss,
  #                 min.beta.MS.tiss, min.beta.MS.tiss.adj, min.gamma.HS.tiss,
  #                 min.beta.HS.tiss, min.beta.MS.tiss.adj, min.gamma.HS.tiss,
  #                 R0.LS, R0.MS, R0.HS,
  #                 cover.site)
  
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
  
  suscat.names = c('LS', 'MS', 'HS')
  tab.offshore = tibble(suscat.names, round(betas.offshore, 2), round(betas.offshore.adj, 2), round(gammas.offshore, 2), round(R0s.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
  
  output.offshore = output.offshore %>% select(-last_col())
  
  output.offshore = pivot_longer(output.offshore, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.offshore)[1] = 'days'
  colnames(output.offshore)[4] = 'prop'
  
  p.fit.offshore = ggplot(data = output.offshore, aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Compartment, shape = Category)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = max(output.offshore$days), y = max(output.offshore$prop)*0.7, label = list(tab.offshore),
             vjust = 1, hjust = 1) +
    theme_classic() + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore = ggplot(data = output.offshore %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of live tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic() + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.I.fit.offshore = ggplot(data = output.offshore %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line(color = '#fc8e62') +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category), color = '#fc8e62') +
    theme_classic() + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.D.fit.offshore = ggplot(data = output.offshore %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line(color = '#66c2a4') +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category), color = '#66c2a4') +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray90"))
  
  # p.fit.offshore;p.I.fit.offshore;p.D.fit.offshore
  
  # Observations only [lines are observations]
  #overlaid
  ggplot(data = obs, aes(days, prop, colour = Compartment, linetype = Category, shape = Site)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue") +
    ggtitle(paste(c("", 'All sites'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  #facet-wrapped
  (p.SIR.offshore | p.SIR.midchannel | p.SIR.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom') #&
  # scale_color_brewer(name = 'Disease compartment', labels = c("High", "Low", "Medium"), palette = 'Dark2')
  
  # Observations only [lines are observations]
  (p.I.offshore | p.I.midchannel | p.I.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  # Fit only [lines are simulated]
  (p.fit.offshore | p.fit.midchannel | p.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.S.fit.offshore | p.S.fit.midchannel | p.S.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.I.fit.offshore | p.I.fit.midchannel | p.I.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  (p.D.fit.offshore | p.D.fit.midchannel | p.D.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  ##################################################################################################################################
  
  #run SIR for Midchannel based on fit from Nearshore
  site = 'Midchannel'
  prev.timepoint = timepoint.midchannel
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  N.site = S0.snapshot$tot.sustiss
  N.LS.site = S0.snapshot$LS.sustiss
  N.MS.site = S0.snapshot$MS.sustiss
  N.HS.site = S0.snapshot$HS.sustiss
  
  LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
  MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
  HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
  
  polyp_SA = min(HS.inftiss[1:5])/5
  
  I.LS.midchannel = 0
  S.LS.midchannel = N.LS.site - I.LS.midchannel
  R.LS.midchannel = 0
  I.MS.midchannel = 0
  S.MS.midchannel = N.MS.site - I.MS.midchannel
  R.MS.midchannel = 0
  I.HS.midchannel = polyp_SA
  S.HS.midchannel = N.HS.site - I.HS.midchannel 
  R.HS.midchannel = 0
  
  P.midchannel = I.LS.midchannel + I.MS.midchannel + I.HS.midchannel
  
  output.midchannel.transfer = data.frame(ode(c(S.LS = S.LS.midchannel, I.LS = I.LS.midchannel, R.LS = R.LS.midchannel,
                                  S.MS = S.MS.midchannel, I.MS = I.MS.midchannel, R.MS = R.MS.midchannel,
                                  S.HS = S.HS.midchannel, I.HS = I.HS.midchannel, R.HS = R.HS.midchannel,
                                  P = P.midchannel),
                                time, SIR, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                             b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                             b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                             N.LS = N.LS.site, N.MS = N.MS.site, N.HS = N.HS.site,
                                             C = cover.midchannel,
                                             C.LS = cover.midchannel.LS, C.MS = cover.midchannel.MS, C.HS = cover.midchannel.HS,
                                             l = lambda)))
  
  output.midchannel.transfer = output.midchannel.transfer %>% select(-last_col())
  
  output.midchannel.transfer = pivot_longer(output.midchannel.transfer, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.midchannel.transfer)[1] = 'days'
  colnames(output.midchannel.transfer)[4] = 'prop'
  
  p.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer, aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Compartment, shape = Category)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
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
  names(tab.midchannel.transfer) = c('Category', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.midchannel.transfer = p.fit.midchannel.transfer +
    annotate(geom = "table", x = max(p.fit.midchannel.transfer$data$days), y = max(p.fit.midchannel.transfer$data$prop)*0.7, label = list(tab.midchannel.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.I.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.D.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab('Surface area of infected tissue (m2)') +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  nearshore.to.midchannel = (p.fit.nearshore | p.I.fit.nearshore | p.D.fit.nearshore) / (p.fit.midchannel | p.I.fit.midchannel | p.D.fit.midchannel) / (p.fit.midchannel.transfer | p.I.fit.midchannel.transfer | p.D.fit.midchannel.transfer)
  
  #run SIR for Offshore based on fit from Nearshore
  site = 'Offshore'
  prev.timepoint = timepoint.offshore
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  N.site = S0.snapshot$tot.sustiss
  N.LS.site = S0.snapshot$LS.sustiss
  N.MS.site = S0.snapshot$MS.sustiss
  N.HS.site = S0.snapshot$HS.sustiss
  
  LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
  MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
  HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
  
  polyp_SA = min(HS.inftiss[1:5])/5
  
  I.LS.offshore = 0
  S.LS.offshore = N.LS.site - I.LS.offshore
  R.LS.offshore = 0
  I.MS.offshore = 0
  S.MS.offshore = N.MS.site - I.MS.offshore
  R.MS.offshore = 0
  I.HS.offshore = polyp_SA
  S.HS.offshore = N.HS.site - I.HS.offshore 
  R.HS.offshore = 0
  
  P.offshore = I.LS.offshore + I.MS.offshore + I.HS.offshore
  
  output.offshore.transfer = data.frame(ode(c(S.LS = S.LS.offshore, I.LS = I.LS.offshore, R.LS = R.LS.offshore,
                                                S.MS = S.MS.offshore, I.MS = I.MS.offshore, R.MS = R.MS.offshore,
                                                S.HS = S.HS.offshore, I.HS = I.HS.offshore, R.HS = R.HS.offshore,
                                                P = P.offshore),
                                              time, SIR, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                                           b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                                           b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                                           N.LS = N.LS.site, N.MS = N.MS.site, N.HS = N.HS.site,
                                                           C = cover.offshore,
                                                           C.LS = cover.offshore.LS, C.MS = cover.offshore.MS, C.HS = cover.offshore.HS,
                                                           l = lambda)))
  
  output.offshore.transfer = output.offshore.transfer %>% select(-last_col())
  
  output.offshore.transfer = pivot_longer(output.offshore.transfer, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.offshore.transfer)[1] = 'days'
  colnames(output.offshore.transfer)[4] = 'prop'
  
  p.fit.offshore.transfer = ggplot(data = output.offshore.transfer, aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Compartment, shape = Category)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
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
  names(tab.offshore.transfer) = c('Category', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.offshore.transfer = p.fit.offshore.transfer +
    annotate(geom = "table", x = max(p.fit.offshore.transfer$data$days), y = max(p.fit.offshore.transfer$data$prop)*0.7, label = list(tab.offshore.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.offshore.transfer = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.I.fit.offshore.transfer = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.D.fit.offshore.transfer = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab('Surface area of infected tissue (m2)') +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  nearshore.to.offshore = (p.fit.nearshore | p.I.fit.nearshore | p.D.fit.nearshore) / (p.fit.offshore | p.I.fit.offshore | p.D.fit.offshore) / (p.fit.offshore.transfer | p.I.fit.offshore.transfer | p.D.fit.offshore.transfer)  + plot_layout(guides = "collect",
                                                                                                                                                                                                                                                axis_titles = 'collect') &
    theme(legend.position = 'bottom') &
    xlim(0, 325)
  
  #run SIR for Nearshore based on fit from Offshore
  site = 'Nearshore'
  prev.timepoint = timepoint.nearshore
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  N.site = S0.snapshot$tot.sustiss
  N.LS.site = S0.snapshot$LS.sustiss
  N.MS.site = S0.snapshot$MS.sustiss
  N.HS.site = S0.snapshot$HS.sustiss
  
  LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
  MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
  HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
  
  polyp_SA = min(HS.inftiss[1:5])/5
  
  I.LS.nearshore = 0
  S.LS.nearshore = N.LS.site - I.LS.nearshore
  R.LS.nearshore = 0
  I.MS.nearshore = 0
  S.MS.nearshore = N.MS.site - I.MS.nearshore
  R.MS.nearshore = 0
  I.HS.nearshore = polyp_SA
  S.HS.nearshore = N.HS.site - I.HS.nearshore 
  R.HS.nearshore = 0
  
  P.nearshore = I.LS.nearshore + I.MS.nearshore + I.HS.nearshore
  
  output.nearshore.transfer = data.frame(ode(c(S.LS = S.LS.nearshore, I.LS = I.LS.nearshore, R.LS = R.LS.nearshore,
                                              S.MS = S.MS.nearshore, I.MS = I.MS.nearshore, R.MS = R.MS.nearshore,
                                              S.HS = S.HS.nearshore, I.HS = I.HS.nearshore, R.HS = R.HS.nearshore,
                                              P = P.nearshore),
                                            time, SIR, c(b.LS = beta.offshore.LS, g.LS = gamma.offshore.LS,
                                                         b.MS = beta.offshore.MS, g.MS = gamma.offshore.MS,
                                                         b.HS = beta.offshore.HS, g.HS = gamma.offshore.HS,
                                                         N.LS = N.LS.site, N.MS = N.MS.site, N.HS = N.HS.site,
                                                         C = cover.nearshore,
                                                         C.LS = cover.nearshore.LS, C.MS = cover.nearshore.MS, C.HS = cover.nearshore.HS,
                                                         l = lambda)))
  
  output.nearshore.transfer = output.nearshore.transfer %>% select(-last_col())
  
  output.nearshore.transfer = pivot_longer(output.nearshore.transfer, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.nearshore.transfer)[1] = 'days'
  colnames(output.nearshore.transfer)[4] = 'prop'
  
  p.fit.nearshore.transfer = ggplot(data = output.nearshore.transfer, aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Compartment, shape = Category)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  beta.nearshore.transfer.LS = beta.offshore.LS * (1 / (1 + exp(-lambda.LS * (cover.nearshore.LS))) + offset.LS)
  R0.nearshore.transfer.LS = beta.nearshore.transfer.LS / gamma.offshore.LS
  beta.nearshore.transfer.MS = beta.offshore.MS * (1 / (1 + exp(-lambda.MS * (cover.nearshore.MS))) + offset.MS)
  R0.nearshore.transfer.MS = beta.nearshore.transfer.MS / gamma.offshore.MS
  beta.nearshore.transfer.HS = beta.offshore.HS * (1 / (1 + exp(-lambda.HS * (cover.nearshore.HS))) + offset.HS)
  R0.nearshore.transfer.HS = beta.nearshore.transfer.HS / gamma.offshore.HS
  
  betas.nearshore.transfer = c(beta.nearshore.transfer.LS, beta.nearshore.transfer.MS, beta.nearshore.transfer.HS)
  R0s.nearshore.transfer = c(R0.nearshore.transfer.LS, R0.nearshore.transfer.MS, R0.nearshore.transfer.HS)
  tab.nearshore.transfer = tibble(suscat.names, round(betas.offshore, 2), round(betas.nearshore.transfer, 2), round(gammas.offshore, 2), round(R0s.nearshore.transfer, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore.transfer) = c('Category', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.nearshore.transfer = p.fit.nearshore.transfer +
    annotate(geom = "table", x = max(p.fit.nearshore.transfer$data$days), y = max(p.fit.nearshore.transfer$data$prop)*0.7, label = list(tab.nearshore.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.nearshore.transfer = ggplot(data = output.nearshore.transfer %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.I.fit.nearshore.transfer = ggplot(data = output.nearshore.transfer %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.D.fit.nearshore.transfer = ggplot(data = output.nearshore.transfer %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab('Surface area of infected tissue (m2)') +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  offshore.to.nearshore = (p.fit.offshore | p.I.fit.offshore | p.D.fit.offshore) / (p.fit.nearshore | p.I.fit.nearshore | p.D.fit.nearshore) / (p.fit.nearshore.transfer | p.I.fit.nearshore.transfer | p.D.fit.nearshore.transfer)
  
  #run SIR for Midchannel based on fit from Offshore
  site = 'Midchannel'
  prev.timepoint = timepoint.midchannel
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  N.site = S0.snapshot$tot.sustiss
  N.LS.site = S0.snapshot$LS.sustiss
  N.MS.site = S0.snapshot$MS.sustiss
  N.HS.site = S0.snapshot$HS.sustiss
  
  LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
  MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
  HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
  
  polyp_SA = min(HS.inftiss[1:5])/5
  
  I.LS.midchannel = 0
  S.LS.midchannel = N.LS.site - I.LS.midchannel
  R.LS.midchannel = 0
  I.MS.midchannel = 0
  S.MS.midchannel = N.MS.site - I.MS.midchannel
  R.MS.midchannel = 0
  I.HS.midchannel = polyp_SA
  S.HS.midchannel = N.HS.site - I.HS.midchannel 
  R.HS.midchannel = 0
  
  P.midchannel = I.LS.midchannel + I.MS.midchannel + I.HS.midchannel
  
  output.midchannel.transfer = data.frame(ode(c(S.LS = S.LS.midchannel, I.LS = I.LS.midchannel, R.LS = R.LS.midchannel,
                                               S.MS = S.MS.midchannel, I.MS = I.MS.midchannel, R.MS = R.MS.midchannel,
                                               S.HS = S.HS.midchannel, I.HS = I.HS.midchannel, R.HS = R.HS.midchannel,
                                               P = P.midchannel),
                                             time, SIR, c(b.LS = beta.offshore.LS, g.LS = gamma.offshore.LS,
                                                          b.MS = beta.offshore.MS, g.MS = gamma.offshore.MS,
                                                          b.HS = beta.offshore.HS, g.HS = gamma.offshore.HS,
                                                          N.LS = N.LS.site, N.MS = N.MS.site, N.HS = N.HS.site,
                                                          C = cover.midchannel,
                                                          C.LS = cover.midchannel.LS, C.MS = cover.midchannel.MS, C.HS = cover.midchannel.HS,
                                                          l = lambda)))
  
  output.midchannel.transfer = output.midchannel.transfer %>% select(-last_col())
  
  output.midchannel.transfer = pivot_longer(output.midchannel.transfer, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>% 
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible', 
                                ifelse(Compartment == 'I.', 'Infected',
                                       ifelse(Compartment == 'R.', 'Dead', Compartment))))
  colnames(output.midchannel.transfer)[1] = 'days'
  colnames(output.midchannel.transfer)[4] = 'prop'
  
  p.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer, aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site), aes(days, prop, colour = Compartment, shape = Category)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
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
  names(tab.midchannel.transfer) = c('Category', 'beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover')
  
  p.fit.midchannel.transfer = p.fit.midchannel.transfer +
    annotate(geom = "table", x = max(p.fit.midchannel.transfer$data$days), y = max(p.fit.midchannel.transfer$data$prop)*0.7, label = list(tab.midchannel.transfer),
             vjust = 1, hjust = 1)
  
  p.S.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.I.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab('Surface area of infected tissue (m2)') +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Infected"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  p.D.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab('Surface area of infected tissue (m2)') +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs %>% filter(Site == site, Compartment == "Dead"), aes(days, prop, shape = Category)) +
    theme_classic()
  
  # 
  # t.1 = (p.fit.offshore + ggtitle(NULL) | p.I.fit.offshore + ggtitle(NULL) | p.D.fit.offshore + ggtitle(NULL)) + plot_layout(guides = "collect", axis_titles = 'collect')
  # t.1 = t.1 & plot_annotation(title = "Offshore - Fitted")
  # t.2 = (p.fit.midchannel + ggtitle(NULL) | p.I.fit.midchannel + ggtitle(NULL) | p.D.fit.midchannel + ggtitle(NULL)) & plot_annotation(title = "Midchannel - Fitted")
  # t.3 = (p.fit.midchannel.transfer + ggtitle(NULL) | p.I.fit.midchannel.transfer + ggtitle(NULL) | p.D.fit.midchannel.transfer + ggtitle(NULL)) & plot_annotation(title = "Midchannel - Predicted")
  # 
  # offshore.to.midchannel = (t.1 / t.2 / t.3) + plot_layout(guides = "collect", axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') &
  #   xlim(0, 325)
  # 
  # nearshore.to.offshore = (p.fit.nearshore | p.fit.offshore | p.fit.offshore.transfer) + plot_layout(guides = "collect",
  #                                                                                                    axis_titles = 'collect') &
  #   theme(legend.position = 'bottom') &
  #   xlim(0, 325)
  # 
  # nearshore.to.midchannel
  # nearshore.to.offshore
  # offshore.to.nearshore
  # offshore.to.midchannel

  # # Save workspace
  # save.image(file = "plots_multi_COVER_0-0-0_workspace.RData")
  