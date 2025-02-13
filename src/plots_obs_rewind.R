  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(ggplot2)
  library(dplyr)
  
  #import workspace from FLKEYS_data_processing.R
  load(here("output", "data_processing_workspace_rewind.RData"))
  
  obs.total = obs %>%
    filter(Category == 'Total')
  obs.multi = obs %>%
    filter(Category != "Total")
  
  # MULTI-GROUP SIR
  # scaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Offshore"), aes(days.inf.site, tissue.scaled, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  # #attempt at plotting SST
  # # Adjust the scaling_factor for plotting SST with epidemics
  # scaling_factor <- max(obs.multi$tissue) / max(obs.multi$SST_90th_HS, na.rm = TRUE)
  # scaling_factor.scaled <- max(obs.multi$tissue.scaled) / max(obs.multi$SST_90th_HS, na.rm = TRUE)
  # 
  # p.SIR.offshore.scale.multi2 = ggplot(data = obs.multi %>% filter(Site == "Offshore"), aes(x = days.inf.site)) +
  #   xlab("Day of observation period") +
  #   ylab("Proportion of tissue.scaled") +
  #   ggtitle('Offshore') +
  # 
  #   # Epidemic outbreak lines
  #   geom_line(aes(y = tissue.scaled, colour = Compartment, linetype = Category)) +
  # 
  #   # SST line, transformed by the scaling factor
  #   geom_line(aes(y = SST_90th_HS * scaling_factor), color = "red", linetype = "dashed", size = 1.2, show.legend = TRUE) +
  # 
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  # 
  #   # Theme adjustments
  #   theme_classic() +
  # 
  #   # Adjust the primary y-axis without forcing it to match the secondary axis
  #   scale_y_continuous(
  #     name = "Proportion of tissue.scaled",
  #     sec.axis = sec_axis(~ . / scaling_factor, name = "Sea Surface Temperature (SST)")
  #   )
  
  p.I.offshore.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel"), aes(days.inf.site, tissue.scaled, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.scale.multi = ggplot(data = obs.multi %>%
                                   filter(Site == "Nearshore"),
                                   # filter(Site == "Nearshore", Category == 'High'),
                                 aes(days.inf.site, tissue.scaled, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.nearshore.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.scale.multi = ggplot(data = obs.multi %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days.inf.site, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # unscaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.multi = ggplot(data = obs.multi %>%
                            filter(Site == "Offshore"),
                            # filter(Site == "Offshore", Category == 'High'),
                          aes(days.inf.site, tissue, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.offshore.multi = ggplot(data = obs.multi %>% filter(Site == "Offshore", Compartment == "Susceptible"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.offshore.multi = ggplot(data = obs.multi %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.multi = ggplot(data = obs.multi %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel"), aes(days.inf.site, tissue, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.multi = ggplot(data = obs.multi %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.multi = ggplot(data = obs.multi %>%
                             filter(Site == "Nearshore"),
                             # filter(Site == "Nearshore", Category == 'High'),
                           aes(days.inf.site, tissue, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.nearshore.multi = ggplot(data = obs.multi %>% filter(Site == "Nearshore", Compartment == "Susceptible"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.nearshore.multi = ggplot(data = obs.multi %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.multi = ggplot(data = obs.multi %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days.inf.site, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # BASIC SIR
  # scaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Offshore"), aes(days.inf.site, tissue.scaled, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.offshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days.inf.site, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days.inf.site, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.scale.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel"), aes(days.inf.site, tissue.scaled, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.midchannel.scale.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days.inf.site, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.scale.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days.inf.site, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore"), aes(days.inf.site, tissue.scaled, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.nearshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days.inf.site, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days.inf.site, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # unscaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore"), aes(days.inf.site, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel"), aes(days.inf.site, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # #attempt at plotting SST
  # scaling_factor <- max(obs.multi %>% filter(Site == "Midchannel", Compartment == "Infected") %>% pull(tissue)) /
  #   max(obs.multi$SST_90th_HS, na.rm = TRUE)
  # 
  # p.I.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days.inf.site, tissue)) +
  #   xlab("Day of observation period") +
  #   ylab("Tissue Surface Area (m2)") +
  #   ggtitle(paste(c("", 'Midchannel'), collapse="")) +
  #   geom_line(aes(y = tissue, colour = Compartment, linetype = Category)) +
  #   geom_line(aes(y = SST_90th_HS * scaling_factor), color = "red", linetype = "dashed", linewidth = 1.2, show.legend = TRUE) +
  #   geom_hline(yintercept = 30 * scaling_factor, linetype = "dashed", color = "red", linewidth = 0.5, alpha = 0.5) +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   theme_classic() +
  #   scale_y_continuous(
  #     name = "Proportion of tissue.scaled",
  #     sec.axis = sec_axis(~ . / scaling_factor, name = "Sea Surface Temperature (SST)")
  #   )
  
  p.I.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore"), aes(days.inf.site, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Susceptible"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days.inf.site, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #save workspace for returning to plots
  save.image(file = here("output", "plots_obs_workspace_rewind.RData"))
  