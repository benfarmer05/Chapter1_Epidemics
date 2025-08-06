  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(patchwork)
  
  ################################## Set-up ##################################
  
  load(here("output/tables_figures_workspace.RData"))
  
  sst_min <- 23
  sst_max <- 32.5
  
  #import figure panels
  fig4 = readRDS(here('output', 'fig4.rds'))
  
  ################################## create new Fig 2 ##################################
  
  # Option 1: Simple side-by-side arrangement
  combined_fig <- fig2a + fig2b
  
  # Option 2: With panel labels
  combined_fig <- fig2a + fig2b + 
    plot_annotation(tag_levels = 'A')
  
  # # Option 3: More control over layout and labels
  # combined_fig <- fig2a + fig2b + 
  #   plot_annotation(tag_levels = 'A') +
  #   plot_layout(ncol = 2)
  
  # # Option 4: Vertical arrangement
  # combined_fig <- fig2a / fig2b + 
  #   plot_annotation(tag_levels = 'A')
  
  # Display the combined figure
  combined_fig
  
  ################################## export to PNG / PDF ##################################
  
  #max dimensions are 7.087 in. wide by 9.45 in. tall (3.35 inches preferred)
  quartz(h = 3, w = 7.087)
  
  combined_fig
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "fig2.pdf"), type = "pdf")
  
  #ggplot-export to image
  ggsave(filename = here("output", "fig2.png"), device = "png", width = 7.087, height = 3, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
