
  #to-do:
# - optimize the prediction of nearshore parameters onto offshore. and vice versa probably
# - keep testing effect of simply not fitting to infections at all in the first optimization step
# - also, thinking about the Dobbelaere approach. I could try a fit that accounts for all sites summed together (or averaged)
#     - what would the prediction look like? similar to the 2020 paper? there are important differences in how they accounted for area
#     - could allow beta to do whatever it wants, and then calculate a gamma that a constant R0 dictates...does this make sense?
#         - might allow the actual shape of the curve to change between sites
#         - could also always test that more with the multi-scale models or a UVI-specific wheels study...eventually
#         - also though, my theory is really that only transmission (beta) itself should change with site density. so I like what I did, too
#
#     - our model is tissue-based too - which is a huge plus. while its predictions aren't perfect, it's probably a much closer reality
#         than just assuming a bunch of whole coral colonies die
#         - and ours is multi-species!

  rm(list=ls())
  
  load("output_basic_COVER_1.0_workspace.RData")
  
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
  
  #for plotting parameter value output in a tibble table format
  val.hjust = 0
  val.vjust = -4
  
  # Nearshore
  site = 'Nearshore'
  order = 2
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == timepoint.nearshore)
  S0.tot.sustiss = S0.snapshot$tot.sustiss
  
  output.nearshore = my.SIRS[[order]]
  params.nearshore = params[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.nearshore = params.nearshore[1]
  beta.nearshore.adj = params.nearshore[2]
  gamma.nearshore = params.nearshore[3]
  lambda.nearshore = params.nearshore[4]
  R0.nearshore = params.nearshore[5]
  cover.nearshore = params.nearshore[6]
  
  tab.nearshore = tibble(round(beta.nearshore, 2), round(beta.nearshore.adj, 2), round(gamma.nearshore, 2),
                          round(R0.nearshore, 2), round(cover.nearshore*100, 2))
  names(tab.nearshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  

  output.nearshore = pivot_longer(output.nearshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.nearshore)[1] = 'days'
  colnames(output.nearshore)[3] = 'prop'
  
  p.fit.nearshore = ggplot(data = output.nearshore, aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site, " - Fitted"), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site), aes(days, prop, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.nearshore$days), y = min(output.nearshore$prop)*0.7, label = list(tab.nearshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia')

  p.S.fit.nearshore = ggplot(data = output.nearshore %>% filter(Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.nearshore = ggplot(data = output.nearshore %>% filter(Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Infected"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.nearshore = ggplot(data = output.nearshore %>% filter(Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Dead"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  # Midchannel
  site = 'Midchannel'
  order = 1
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == timepoint.midchannel)
  S0.tot.sustiss = S0.snapshot$tot.sustiss
  
  output.midchannel = my.SIRS[[order]]
  params.midchannel = params[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.midchannel = params.midchannel[1]
  beta.midchannel.adj = params.midchannel[2]
  gamma.midchannel = params.midchannel[3]
  lambda.midchannel = params.midchannel[4]
  R0.midchannel = params.midchannel[5]
  cover.midchannel = params.midchannel[6]
  
  tab.midchannel = tibble(round(beta.midchannel, 2), round(beta.midchannel.adj, 2), round(gamma.midchannel, 2),
                          round(R0.midchannel, 2), round(cover.midchannel*100, 2))
  names(tab.midchannel) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  output.midchannel = pivot_longer(output.midchannel, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.midchannel)[1] = 'days'
  colnames(output.midchannel)[3] = 'prop'
  
  p.fit.midchannel = ggplot(data = output.midchannel, aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site), aes(days, prop, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.midchannel$days), y = min(output.midchannel$prop)*0.7, label = list(tab.midchannel),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') +
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.midchannel = ggplot(data = output.midchannel %>% filter(Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel = ggplot(data = output.midchannel %>% filter(Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Infected"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel = ggplot(data = output.midchannel %>% filter(Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Dead"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  # Offshore
  site = 'Offshore'
  order = 3
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == timepoint.offshore)
  S0.tot.sustiss = S0.snapshot$tot.sustiss
  
  output.offshore = my.SIRS[[order]]
  params.offshore = params[[order]] #beta, cover-adjusted beta, gamma, lambda, R0, cover
  beta.offshore = params.offshore[1]
  beta.offshore.adj = params.offshore[2]
  gamma.offshore = params.offshore[3]
  lambda.offshore = params.offshore[4]
  R0.offshore = params.offshore[5]
  cover.offshore = params.offshore[6]
  
  tab.offshore = tibble(round(beta.offshore, 2), round(beta.offshore.adj, 2), round(gamma.offshore, 2),
                          round(R0.offshore, 2), round(cover.offshore*100, 2))
  names(tab.offshore) = c('beta', 'Adj. beta', 'gamma', 'R0', 'Cover (%)')
  
  output.offshore = pivot_longer(output.offshore, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  
  colnames(output.offshore)[1] = 'days'
  colnames(output.offshore)[3] = 'prop'
  
  p.fit.offshore = ggplot(data = output.offshore, aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site, ' - Fitted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site), aes(days, prop, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    annotate(geom = "table", x = min(output.offshore$days), y = min(output.offshore$prop)*0.7, label = list(tab.offshore),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
    theme_classic(base_family = 'Georgia') + 
    theme(panel.background = element_rect(fill = "gray90"))
  
  p.S.fit.offshore = ggplot(data = output.offshore %>% filter(Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of susceptible tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore = ggplot(data = output.offshore %>% filter(Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Infected"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore = ggplot(data = output.offshore %>% filter(Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of dead tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Dead"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  
  # # Observations only [lines are observations]
  # #overlaid
  # ggplot(data = obs, aes(days, prop, colour = Compartment, linetype = Category, shape = Site)) +
  #   xlab("Day of observation period") +
  #   ylab("Proportion of tissue") +
  #   ggtitle(paste(c("", 'All sites'), collapse="")) +
  #   geom_line() +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   theme_classic(base_family = 'Georgia')
  # 
  # #facet-wrapped
  # (p.SIR.offshore | p.SIR.midchannel | p.SIR.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom') #&
  # # scale_color_brewer(name = 'Disease compartment', labels = c("High", "Low", "Medium"), palette = 'Dark2')
  # 
  # # Observations only [lines are observations]
  # (p.I.offshore | p.I.midchannel | p.I.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  # # Fit only [lines are simulated]
  # (p.fit.offshore | p.fit.midchannel | p.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  # (p.S.fit.offshore | p.S.fit.midchannel | p.S.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  # (p.I.fit.offshore | p.I.fit.midchannel | p.I.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  # (p.D.fit.offshore | p.D.fit.midchannel | p.D.fit.nearshore) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  ##############################################################################################################################

  #run SIR for Offshore based on fit from Nearshore
  # site = 'Nearshore'
  # prev.timepoint = 'T11'
  # S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  # S0.tot.sustiss = S0.snapshot$tot.sustiss

  #Offshore
  site = 'Offshore'
  prev.timepoint = 'T5'
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  S0.tot.sustiss = S0.snapshot$tot.sustiss

  LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
  MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
  HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
  inftiss = LS.inftiss+MS.inftiss+HS.inftiss

  # polyp_SA = min(inftiss[1:5])/5
  polyp_SA = inftiss[1]/5

  # I.offshore = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  I.offshore = polyp_SA
  S.offshore = S0.tot.sustiss - I.offshore
  R.offshore = 0

  #simulation using initial state variables from naive site and parameters from fitted site
  output.offshore.transfer = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                   time, SIR, c(b = beta.nearshore, g = gamma.nearshore,
                                                N = S0.tot.sustiss,
                                                l = lambda.nearshore,
                                                C = cover.offshore)))
  
  output.offshore.transfer = pivot_longer(output.offshore.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.offshore.transfer)[1] = 'days'
  colnames(output.offshore.transfer)[3] = 'prop'

  p.fit.offshore.transfer = ggplot(data = output.offshore.transfer, aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site), aes(days, prop, colour = Compartment)) +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic(base_family = 'Georgia')
  
  beta.offshore.transfer = beta.nearshore * (1 / (1 + exp(-lambda.modifier * (cover.offshore))) + offset)
  # beta.offshore.transfer = (beta.nearshore * (lambda.nearshore * (1-exp(-130*(cover.offshore)))))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.offshore)))
  # (beta.nearshore * (1 + lambda.nearshore * sqrt(cover.nearshore)))
  R0.offshore.transfer = beta.offshore.transfer / gamma.nearshore
  tab.offshore.transfer = tibble(round(beta.nearshore, 2), round(beta.offshore.transfer, 2), round(gamma.nearshore, 2), round(R0.offshore.transfer, 2), round(cover.offshore*100, 2))
  names(tab.offshore.transfer) = c('beta', 'Adj. beta', 'gamma', 'Adj. R0', 'Cover (%)')
  
  p.fit.offshore.transfer = p.fit.offshore.transfer +
    annotate(geom = "table", x = min(p.fit.offshore.transfer$data$days), y = min(p.fit.offshore.transfer$data$prop)*0.7, label = list(tab.offshore.transfer),
             vjust = val.vjust, hjust = val.hjust, family = 'Georgia')
    # geom_text(aes(x = max(days)*0.9, y = max(prop)*0.9, label = paste("R0:", R0.offshore.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.offshore.transfer = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.offshore.transfer = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site, ' - Projected'), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Infected"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.offshore.transfer = ggplot(data = output.offshore.transfer %>% filter(Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Dead"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  # p.fit.nearshore
  # p.fit.offshore
  # p.fit.offshore.transfer
  # p.S.fit.nearshore
  # p.S.fit.offshore
  # p.S.fit.offshore.transfer
  # p.I.fit.nearshore
  # p.I.fit.offshore.transfer
  # p.D.fit.offshore.transfer
  # p.fit.offshore.transfer;p.S.fit.offshore.transfer;p.I.fit.offshore.transfer;p.D.fit.offshore.transfer
  # (p.fit.nearshore | p.fit.offshore | p.fit.offshore.transfer)
  
  # (p.fit.nearshore | p.I.fit.nearshore | p.D.fit.nearshore) / (p.fit.offshore.transfer | p.I.fit.offshore.transfer | p.D.fit.offshore.transfer)
  
  
  #run SIR for Midchannel based on fit from Nearshore
  # site = 'Nearshore'
  # prev.timepoint = 'T11'
  # S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  # S0.tot.sustiss = S0.snapshot$tot.sustiss
  
  #Midchannel
  site = 'Midchannel'
  prev.timepoint = 'T5'
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  S0.tot.sustiss = S0.snapshot$tot.sustiss
  
  LS.inftiss = obs %>% filter(Site == site, Category == "LS", Compartment == "Infected") %>% pull(prop)
  MS.inftiss = obs %>% filter(Site == site, Category == "MS", Compartment == "Infected") %>% pull(prop)
  HS.inftiss = obs %>% filter(Site == site, Category == "HS", Compartment == "Infected") %>% pull(prop)
  inftiss = LS.inftiss+MS.inftiss+HS.inftiss
  
  # polyp_SA = min(inftiss[1:5])/5
  polyp_SA = inftiss[1]/5
  
  # I.midchannel = 1e-4 #m2 - equivalent to 100 mm2, which is a rough approximation of a fully infected medium-sized coral polyp
  I.midchannel = polyp_SA
  S.midchannel = S0.tot.sustiss - I.midchannel
  R.midchannel = 0
  
  #simulation using initial state variables from naive site and parameters from fitted site
  output.midchannel.transfer = data.frame(ode(c(S = S.midchannel, I = I.midchannel, R = R.midchannel),
                                            time, SIR, c(b = beta.nearshore, g = gamma.nearshore,
                                                         N = S0.tot.sustiss,
                                                         l = lambda.nearshore,
                                                         C = cover.midchannel)))
  
  output.midchannel.transfer = pivot_longer(output.midchannel.transfer, cols = -1, names_to = c("Compartment")) %>%
    mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
    mutate(Compartment = ifelse(Compartment == 'S', 'Susceptible',
                                ifelse(Compartment == 'I', 'Infected',
                                       ifelse(Compartment == 'R', 'Dead', Compartment))))
  colnames(output.midchannel.transfer)[1] = 'days'
  colnames(output.midchannel.transfer)[3] = 'prop'
  
  p.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer, aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Surface area of tissue (m2)") +
    ggtitle(paste(c("", site, ' - Predicted'), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site), aes(days, prop, colour = Compartment)) +
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
  
  p.fit.midchannel.transfer = p.fit.midchannel.transfer +
    annotate(geom = "table", x = min(p.fit.midchannel.transfer$data$days), y = min(p.fit.midchannel.transfer$data$prop)*0.7, label = list(tab.midchannel.transfer),
             family = 'Georgia', hjust = val.hjust, vjust = val.vjust)
  # geom_text(aes(x = max(days)*0.9, y = max(prop)*0.9, label = paste("R0:", R0.midchannel.transfer)), color = "black", hjust = 1, family = 'Georgia', vjust = 1, size = 4)
  
  p.S.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Susceptible"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.I.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Infected"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  p.D.fit.midchannel.transfer = ggplot(data = output.midchannel.transfer %>% filter(Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Surface area of infected tissue") +
    ggtitle(paste(c("", site), collapse="")) +
    geom_line() +
    geom_point(data = obs.basic %>% filter(Site == site, Compartment == "Dead"), aes(days, prop)) +
    theme_classic(base_family = 'Georgia')
  
  # p.fit.nearshore
  # p.fit.midchannel
  # p.fit.midchannel.transfer
  # p.S.fit.nearshore
  # p.S.fit.midchannel.transfer
  # p.I.fit.nearshore
  # p.I.fit.midchannel.transfer
  # p.D.fit.midchannel.transfer
  # p.fit.midchannel.transfer;p.S.fit.midchannel.transfer;p.I.fit.midchannel.transfer;p.D.fit.midchannel.transfer
  
  # (p.fit.nearshore | p.I.fit.nearshore | p.D.fit.nearshore) / (p.fit.midchannel.transfer | p.I.fit.midchannel.transfer | p.D.fit.midchannel.transfer)
  
  # Remove legend from individual plots
  # p.fit.nearshore <- p.fit.nearshore + theme(legend.position = 'none')
  # p.fit.midchannel <- p.fit.midchannel + theme(legend.position = 'none')
  # p.fit.midchannel.transfer <- p.fit.midchannel.transfer + theme(legend.position = 'none')
  
  # Combine the plots
  nearshore.to.offshore = (p.fit.nearshore | p.fit.offshore | p.fit.offshore.transfer) + plot_layout(guides = "collect",
                                                                                                            axis_titles = 'collect') &
    theme(legend.position = 'bottom') &
    xlim(0, 325)
  
  nearshore.to.midchannel = (p.fit.nearshore | p.fit.midchannel | p.fit.midchannel.transfer) + plot_layout(guides = "collect",
                                                                                                  axis_titles = 'collect') &
    theme(legend.position = 'bottom') &
    xlim(0, 325)
  
  # # Save workspace
  # save.image(file = "plots_basic_COVER_1.0_workspace.RData")
  