  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())

  library(here)
  library(ggplot2)
  library(dplyr)
  
  #import workspace from FLKEYS_data_processing.R
  load(here("output", "data_processing_workspace.RData"))
  
  obs.total = obs %>%
    filter(Category == 'Total')
  obs = obs %>%
    filter(Category != "Total")
  
  # MULTI-GROUP SIR
  # scaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.scale = ggplot(data = obs %>% filter(Site == "Offshore"), aes(days, tissue.scaled, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.offshore.scale = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.scale = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.scale = ggplot(data = obs %>% filter(Site == "Midchannel"), aes(days, tissue.scaled, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.scale = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel.scale = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.scale = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.scale = ggplot(data = obs %>%
                                   filter(Site == "Nearshore"),
                                   # filter(Site == "Nearshore", Category == 'High'),
                                 aes(days, tissue.scaled, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.nearshore.scale = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.scale = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, tissue.scaled, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue.scaled") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # unscaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore = ggplot(data = obs %>%
                            # filter(Site == "Offshore"),
                            filter(Site == "Offshore", Category == 'High'),
                          aes(days, tissue, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.offshore = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Susceptible"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.offshore = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel"), aes(days, tissue, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore = ggplot(data = obs %>%
                             filter(Site == "Nearshore"),
                             # filter(Site == "Nearshore", Category == 'High'),
                           aes(days, tissue, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Susceptible"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, tissue, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # BASIC SIR
  # scaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Offshore"), aes(days, tissue.scaled, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.offshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.scale.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel"), aes(days, tissue.scaled, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.midchannel.scale.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.scale.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore"), aes(days, tissue.scaled, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.nearshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.scale.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, tissue.scaled)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # unscaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore"), aes(days, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Susceptible"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.basic = ggplot(data = obs.total %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel"), aes(days, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.basic = ggplot(data = obs.total %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore"), aes(days, tissue, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Susceptible"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.basic = ggplot(data = obs.total %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, tissue)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # # Save workspace
  # save.image(file = "FLKEYS_workspace.RData")
  
  
  
  
  
  
  
  
  # # development
  # #
  # # Create a long format of the summary data for SIR metrics
  # summary_long <- summary %>%
  #   select(site, days.survey, 
  #          low.susnum, moderate.susnum, high.susnum,
  #          low.infnum, moderate.infnum, high.infnum,
  #          low.deadnum, moderate.deadnum, high.deadnum) %>%
  #   pivot_longer(
  #     cols = -c(site, days.survey), 
  #     names_to = c("susceptibility", "state"),
  #     names_pattern = "(.*)\\.(.*)"
  #   ) %>%
  #   # Replace state labels with user-friendly names
  #   mutate(state = recode(state,
  #                         susnum = "Susceptible",
  #                         infnum = "Infected",
  #                         deadnum = "Dead"
  #   )) %>%
  #   # Format susceptibility labels with capital letters
  #   mutate(susceptibility = recode(susceptibility,
  #                                  low = "Low",
  #                                  moderate = "Moderate",
  #                                  high = "High"
  #   ))
  # 
  # # Check the structure of summary_long to ensure correct reshaping
  # str(summary_long)
  # 
  # # Plot using ggplot
  # ggplot(summary_long, aes(x = days.survey, y = value, color = state, linetype = susceptibility)) +
  #   geom_line() +
  #   facet_wrap(~ site, scales = "free_y") +
  #   labs(title = "SIR Metrics through Time by Site",
  #        x = "Days since First Infection",
  #        y = "Count",
  #        color = "State",
  #        linetype = "Susceptibility") +
  #   theme_minimal()  
  # 
  # # Update the summary_long to include total counts
  # tot.summary_long <- summary %>%
  #   pivot_longer(cols = c(tot.susnum, tot.infnum, tot.deadnum),
  #                names_to = "state",
  #                values_to = "value") %>%
  #   mutate(
  #     state = case_when(
  #       state == "tot.susnum" ~ "Susceptible",
  #       state == "tot.infnum" ~ "Infected",
  #       state == "tot.deadnum" ~ "Dead"
  #     )
  #   ) %>%
  #   select(site, days.survey, state, value)
  # 
  # # Combined plot with coloring by both State and Site
  # combined_plot <- ggplot(tot.summary_long, aes(x = days.survey, y = value, color = state, group = interaction(state, site))) +
  #   geom_line(aes(linetype = site)) +
  #   labs(title = "Combined SIR Metrics for All Sites",
  #        color = "State",
  #        linetype = "Site") +
  #   theme_minimal()
