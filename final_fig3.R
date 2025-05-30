# =============================================================================
# PART 2: AFTER RUNNING BOTH WORKSPACES, USE THIS TO COMBINE AND PLOT
# =============================================================================

# Load both datasets
panel_e_ws1 <- readRDS(here("output", "panel_e_data_workspace_1.rds"))
panel_e_ws2 <- readRDS(here("output", "panel_e_data_workspace_2.rds"))

# Load observed data (assuming it's the same for both)
obs_e_ws1 <- readRDS(here("output", "panel_e_obs_workspace_1.rds"))
# obs_e_ws2 <- readRDS(here("output", "panel_e_obs_workspace_2.rds"))  # if different

# Combine the datasets
combined_panel_e <- rbind(panel_e_ws1, panel_e_ws2)

# Create the combined plot
p_combined_panel_e <- ggplot(data = combined_panel_e, 
                             aes(x = days.model, y = tissue, 
                                 color = workspace, linetype = workspace)) +
  geom_line(aes(group = interaction(Compartment, workspace)), linewidth = 0.75) +
  geom_point(data = obs_e_ws1, 
             aes(x = days.inf.site, y = tissue, shape = Compartment), 
             color = "black", size = 1.3, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, max(combined_panel_e$days.model, na.rm = TRUE))) +
  scale_shape_manual(values = c("Susceptible" = 16, "Infected" = 17, "Recovered" = 15),
                     name = "Compartment") +
  scale_color_manual(values = c("workspace_1" = "blue", "workspace_2" = "red"),
                     name = "Workspace") +
  scale_linetype_manual(values = c("workspace_1" = "solid", "workspace_2" = "dashed"),
                        name = "Workspace") +
  xlab("Day of outbreak") +
  ylab("Surface area of tissue (m²)") +
  theme_classic(base_family = "Georgia") +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))

# Display the plot
print(p_combined_panel_e)

# Save the combined plot
ggsave(filename = here("output", "panel_e_combined_workspaces.png"), 
       plot = p_combined_panel_e,
       device = "png", width = 6, height = 4, dpi = 300)