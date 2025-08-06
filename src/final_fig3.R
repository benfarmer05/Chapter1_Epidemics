 
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(patchwork)
  
  ################################## Set-up ##################################
  
  #import workspace from upstream script
  load(here("output/tables_figures_workspace.RData"))
  
  ################################## create new Fig 3e ##################################
  
  # Load both datasets
  panel_e_ws1 <- readRDS(here("output", "panel_e_data_workspace_1.rds"))
  panel_e_ws2 <- readRDS(here("output", "panel_e_data_workspace_2.rds"))
  
  # Load observed data (assuming it's the same for both)
  obs_e_ws1 <- readRDS(here("output", "panel_e_obs_workspace_1.rds"))
  
  # # Change "Recovered" to "Removed" in all datasets
  # panel_e_ws1$Compartment <- ifelse(panel_e_ws1$Compartment == "Recovered", "Removed", panel_e_ws1$Compartment)
  # panel_e_ws2$Compartment <- ifelse(panel_e_ws2$Compartment == "Recovered", "Removed", panel_e_ws2$Compartment)
  # 
  # # Handle factor conversion for obs_e_ws1
  # obs_e_ws1$Compartment <- as.character(obs_e_ws1$Compartment)
  # obs_e_ws1$Compartment <- ifelse(obs_e_ws1$Compartment == "Recovered", "Removed", obs_e_ws1$Compartment)
  # obs_e_ws1$Compartment <- factor(obs_e_ws1$Compartment, levels = c("Susceptible", "Infected", "Removed"))
  
  # Filter out 'Infected' compartment from workspace_2 (NL Dens.)
  panel_e_ws2_filtered <- panel_e_ws2 %>%
    filter(Compartment != "Infected")
  
  # Rename workspace labels - only NL Dens. will show in legend
  panel_e_ws1$workspace <- "Freq."
  panel_e_ws2_filtered$workspace <- "NL Dens."
  
  # Combine the datasets
  combined_panel_e <- rbind(panel_e_ws1, panel_e_ws2_filtered)
  
  # Calculate consistent y-axis limits
  # For panels A&B (nearshore data) - get max from both single and multi-host nearshore data
  max_nearshore <- max(c(
    max(data_fig3 %>% filter(Site == "Nearshore") %>% pull(tissue), na.rm = TRUE),
    max(combined_panel_e$tissue, na.rm = TRUE)  # in case panel E data is higher
  ), na.rm = TRUE)
  
  # For panels C/D/E/F (offshore data) - get max from offshore and panel E data
  max_offshore <- max(c(
    max(data_fig3 %>% filter(Site == "Offshore") %>% pull(tissue), na.rm = TRUE),
    max(combined_panel_e$tissue, na.rm = TRUE)
  ), na.rm = TRUE)
  
  # Create the combined plot with only NL Dens. in the legend
  p_combined_panel_e <- ggplot(data = combined_panel_e, 
                               aes(x = days.model, y = tissue, 
                                   color = workspace, linetype = workspace)) +
    geom_line(aes(group = interaction(Compartment, workspace)), linewidth = 0.75) +
    geom_point(data = obs_e_ws1, 
               aes(x = days.inf.site, y = tissue, shape = Compartment), 
               color = "black", size = 1.3, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0, max(combined_panel_e$days.model, na.rm = TRUE))) +
    scale_y_continuous(limits = c(0, max_offshore)) +  # Set consistent y-axis for offshore panels
    scale_shape_manual(values = c("Susceptible" = 16, "Infected" = 17, "Removed" = 15),
                       name = "") +  # Remove "Compartment" title to match your other plots
    scale_color_manual(values = c("Freq." = "black", "NL Dens." = "#228B22"),  # Black and forest green
                       name = "Model",
                       breaks = c("NL Dens.")) +  # Only show NL Dens. in legend
    scale_linetype_manual(values = c("Freq." = "solid", "NL Dens." = "dashed"),
                          guide = "none") +  # Hide linetype legend since both are solid
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.key.height = unit(0, "cm"))  # Match your other plot styling
  
  p5.fig3.final = p_combined_panel_e
  
  # You'll also need to add y-axis limits to your other plots
  # Add these lines to your existing plot code before creating fig3_final:
  
  # For panels A&B (p1.fig3 and p2.fig3), add this:
  p1.fig3 <- p1.fig3 + scale_y_continuous(limits = c(0, max_nearshore))
  p2.fig3 <- p2.fig3 + scale_y_continuous(limits = c(0, max_nearshore))
  
  # For panels C, D, F (p3.fig3, p4.fig3, p6.fig3), add this:
  p3.fig3 <- p3.fig3 + scale_y_continuous(limits = c(0, max_offshore))
  p4.fig3 <- p4.fig3 + scale_y_continuous(limits = c(0, max_offshore))
  p6.fig3 <- p6.fig3 + scale_y_continuous(limits = c(0, max_offshore))
  
  fig3_final =
    p1.fig3 + p2.fig3 + p3.fig3 + p4.fig3 + p5.fig3.final + p6.fig3 +
    plot_layout(axis_titles = "collect",
                guide = 'collect',
                design = "AB\nCD\nEF") &
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize, color = 'black'),
          axis.text = element_text(size = textsize, color = 'black'),
          # legend.box = "vertical",  # Stack legends vertically
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black")
          # legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
          # plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  # add annotation separately (must be done to avoid breaking axes collect)
  fig3_final <- fig3_final + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)

  fig3_final

  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "fig3_final.pdf"), type = "pdf")

  #ggplot-export to image
  ggsave(filename = here("output", "fig3_final.png"), device = "png", width = 7.087, height = 5, dpi = 1200)

  # Close the Quartz device
  dev.off()
  
  ################################## create new Fig S4 ##################################
  
  # Load both datasets
  panel_e_ws1 <- readRDS(here("output", "panel_e_data_workspace_1.rds"))
  panel_e_ws2 <- readRDS(here("output", "panel_e_data_workspace_2.rds"))
  
  # Load observed data (assuming it's the same for both)
  obs_e_ws1 <- readRDS(here("output", "panel_e_obs_workspace_1.rds"))
  
  # Filter out 'Infected' compartment
  obs_e_ws1_filtered <- obs_e_ws1 %>%
    filter(Compartment == "Infected")
  panel_e_ws1_filtered <- panel_e_ws1 %>%
    filter(Compartment == "Infected")
  panel_e_ws2_filtered <- panel_e_ws2 %>%
    filter(Compartment == "Infected")
  
  # Rename workspace labels - only NL Dens. will show in legend
  panel_e_ws1_filtered$workspace <- "Freq."
  panel_e_ws2_filtered$workspace <- "NL Dens."
  
  # Combine the datasets
  combined_panel_e <- rbind(panel_e_ws1_filtered, panel_e_ws2_filtered)
  
  # Calculate consistent y-axis limits
  # For panels A&B (nearshore data) - get max from both single and multi-host nearshore data
  max_nearshore <- max(c(
    max(data_fig3 %>% filter(Site == "Nearshore", Compartment == 'Infected') %>% pull(tissue), na.rm = TRUE),
    max(combined_panel_e$tissue, na.rm = TRUE)  # in case panel E data is higher
  ), na.rm = TRUE)
  
  # For panels C/D/E/F (offshore data) - get max from offshore and panel E data
  max_offshore <- max(c(
    max(data_fig3 %>% filter(Site == "Offshore", Compartment == 'Infected') %>% pull(tissue), na.rm = TRUE),
    max(combined_panel_e$tissue, na.rm = TRUE)
  ), na.rm = TRUE)
  
  # Create the combined plot with only NL Dens. in the legend
  p_combined_panel_e <- ggplot(data = combined_panel_e, 
                               aes(x = days.model, y = tissue, 
                                   color = workspace, linetype = workspace)) +
    geom_line(aes(group = interaction(Compartment, workspace)), linewidth = 0.75) +
    geom_point(data = obs_e_ws1_filtered, 
               aes(x = days.inf.site, y = tissue), 
               color = "black", size = 1.3, shape = 17) +  # Removed shape aesthetic and set fixed shape
    scale_x_continuous(limits = c(0, max(combined_panel_e$days.model, na.rm = TRUE))) +
    scale_y_continuous(limits = c(0, max_offshore)) +  # Set consistent y-axis for offshore panels
    # Removed scale_shape_manual entirely since we don't want it in the legend
    scale_color_manual(values = c("Freq." = "black", "NL Dens." = "#228B22"),  # Black and forest green
                       name = "Model",
                       breaks = c("NL Dens.")) +  # Only show NL Dens. in legend
    scale_linetype_manual(values = c("Freq." = "solid", "NL Dens." = "solid"),
                          guide = "none") +  # Hide linetype legend since both are solid
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.key.height = unit(0, "cm"))  # Match your other plot styling
  
  p5.figS4.final = p_combined_panel_e
  
  # You'll also need to add y-axis limits to your other plots
  # Add these lines to your existing plot code before creating fig3_final:
  
  # For panels A&B (p1.figS4 and p2.figS4), add this:
  p1.figS4 <- p1.figS4 + scale_y_continuous(limits = c(0, max_nearshore))
  p2.figS4 <- p2.figS4 + scale_y_continuous(limits = c(0, max_nearshore))
  
  # For panels C, D, F (p3.figS4, p4.figS4, p6.figS4), add this:
  p3.figS4 <- p3.figS4 + scale_y_continuous(limits = c(0, max_offshore))
  p4.figS4 <- p4.figS4 + scale_y_continuous(limits = c(0, max_offshore))
  p6.figS4 <- p6.figS4 + scale_y_continuous(limits = c(0, max_offshore))
  
  figS4_final =
    p1.figS4 + p2.figS4 + p3.figS4 + p4.figS4 + p5.figS4.final + p6.figS4 +
    plot_layout(axis_titles = "collect",
                guide = 'collect',
                design = "AB\nCD\nEF") &
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize, color = 'black'),
          axis.text = element_text(size = textsize, color = 'black'),
          # legend.box = "vertical",  # Stack legends vertically
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black")
          # legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
          # plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  # add annotation separately (must be done to avoid breaking axes collect)
  figS4_final <- figS4_final + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)
  
  figS4_final
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "figS4_final.pdf"), type = "pdf")
  
  #ggplot-export to image
  ggsave(filename = here("output", "figS4_final.png"), device = "png", width = 7.087, height = 5, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
  
  
  ################################## create new Fig S5 ##################################
  
  # Load both datasets
  panel_e_ws1 <- readRDS(here("output", "panel_e_data_workspace_1.rds"))
  panel_e_ws2 <- readRDS(here("output", "panel_e_data_workspace_2.rds"))
  
  # Load observed data (assuming it's the same for both)
  obs_e_ws1 <- readRDS(here("output", "panel_e_obs_workspace_1.rds"))
  
  # Filter out 'Recovered' compartment
  obs_e_ws1_filtered <- obs_e_ws1 %>%
    filter(Compartment == "Recovered")
  panel_e_ws1_filtered <- panel_e_ws1 %>%
    filter(Compartment == "Recovered")
  panel_e_ws2_filtered <- panel_e_ws2 %>%
    filter(Compartment == "Recovered")
  
  # Rename workspace labels - only NL Dens. will show in legend
  panel_e_ws1_filtered$workspace <- "Freq."
  panel_e_ws2_filtered$workspace <- "NL Dens."
  
  # Combine the datasets
  combined_panel_e <- rbind(panel_e_ws1_filtered, panel_e_ws2_filtered)
  
  # Calculate consistent y-axis limits
  # For panels A&B (nearshore data) - get max from both single and multi-host nearshore data
  max_nearshore <- max(c(
    max(data_fig3 %>% filter(Site == "Nearshore", Compartment == 'Recovered') %>% pull(tissue), na.rm = TRUE),
    max(combined_panel_e$tissue, na.rm = TRUE)  # in case panel E data is higher
  ), na.rm = TRUE)
  
  # For panels C/D/E/F (offshore data) - get max from offshore and panel E data
  max_offshore <- max(c(
    max(data_fig3 %>% filter(Site == "Offshore", Compartment == 'Recovered') %>% pull(tissue), na.rm = TRUE),
    max(combined_panel_e$tissue, na.rm = TRUE)
  ), na.rm = TRUE)
  
  # Create the combined plot with only NL Dens. in the legend
  p_combined_panel_e <- ggplot(data = combined_panel_e, 
                               aes(x = days.model, y = tissue, 
                                   color = workspace, linetype = workspace)) +
    geom_line(aes(group = interaction(Compartment, workspace)), linewidth = 0.75) +
    geom_point(data = obs_e_ws1_filtered, 
               aes(x = days.inf.site, y = tissue), 
               color = "black", size = 1.3, shape = 17) +  # Removed shape aesthetic and set fixed shape
    scale_x_continuous(limits = c(0, max(combined_panel_e$days.model, na.rm = TRUE))) +
    scale_y_continuous(limits = c(0, max_offshore)) +  # Set consistent y-axis for offshore panels
    # Removed scale_shape_manual entirely since we don't want it in the legend
    scale_color_manual(values = c("Freq." = "black", "NL Dens." = "#228B22"),  # Black and forest green
                       name = "Model",
                       breaks = c("NL Dens.")) +  # Only show NL Dens. in legend
    scale_linetype_manual(values = c("Freq." = "solid", "NL Dens." = "solid"),
                          guide = "none") +  # Hide linetype legend since both are solid
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.key.height = unit(0, "cm"))  # Match your other plot styling
  
  p5.figS5.final = p_combined_panel_e
  
  # You'll also need to add y-axis limits to your other plots
  # Add these lines to your existing plot code before creating fig3_final:
  
  # For panels A&B (p1.figS5 and p2.figS5), add this:
  p1.figS5 <- p1.figS5 + scale_y_continuous(limits = c(0, max_nearshore))
  p2.figS5 <- p2.figS5 + scale_y_continuous(limits = c(0, max_nearshore))
  
  # For panels C, D, F (p3.figS5, p4.figS5, p6.figS5), add this:
  p3.figS5 <- p3.figS5 + scale_y_continuous(limits = c(0, max_offshore))
  p4.figS5 <- p4.figS5 + scale_y_continuous(limits = c(0, max_offshore))
  p6.figS5 <- p6.figS5 + scale_y_continuous(limits = c(0, max_offshore))
  
  figS5_final =
    p1.figS5 + p2.figS5 + p3.figS5 + p4.figS5 + p5.figS5.final + p6.figS5 +
    plot_layout(axis_titles = "collect",
                guide = 'collect',
                design = "AB\nCD\nEF") &
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize, color = 'black'),
          axis.text = element_text(size = textsize, color = 'black'),
          # legend.box = "vertical",  # Stack legends vertically
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black")
          # legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
          # plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  # add annotation separately (must be done to avoid breaking axes collect)
  figS5_final <- figS5_final + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)
  
  figS5_final
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "figS5_final.pdf"), type = "pdf")
  
  #ggplot-export to image
  ggsave(filename = here("output", "figS5_final.png"), device = "png", width = 7.087, height = 5, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
  
  