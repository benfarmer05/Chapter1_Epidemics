  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(patchwork)
  library(extrafont)
  library(ggthemes)
  library(scales)
  
  ################################## Set-up ##################################
  
  # #this was required for me to run at least once on an M1 Macbook (last updated December 2024)
  # font_import() 
  # loadfonts(device = "pdf")
  
  #very helpful:
  # https://www.rforecology.com/post/exporting-plots-in-r/
  
  #import workspace from upstream script
  load(here("output/error_eval_workspace.RData"))

  ################################## Table 1 ##################################
  
  # NOTE - the coral cover conversions done here are not ratios based on particular taxa. this could be improved by re-analyzing CPCe output,
  #         but may not be worth it for this study. the main thing it would affect is presentation of this table, and some interpretation
  #         of the community model
  
  species_names <- c(
  "AAGA" = "Agaricia agaricites",
  "ACER" = "Acropora cervicornis",
  "CNAT" = "Colpophyllia natans",
  "DLAB" = "Diploria labyrinthiformis",
  "DSTO" = "Dichocoenia stokesi",
  "MCAV" = "Montastraea cavernosa",
  "MMEA" = "Meandrina meandrites",
  "MYCE" = "Mycetophyllia spp.",
  "OANN" = "Orbicella annularis",
  "OCUL" = "Oculina spp.",
  "ODIF" = "Oculina diffusa",
  "OFAV" = "Orbicella faveolata",
  "PAST" = "Porites astreoides",
  "PCLI" = "Pseudodiploria clivosa",
  "PDIV" = "Porites divaricata",
  "PPOR" = "Porites porites",
  "PSTR" = "Pseudodiploria strigosa",
  "SBOU" = "Solenastrea bournoni",
  "SINT" = "Stephanocoenia intersepta",
  "SRAD" = "Siderastrea radians",
  "SSID" = "Siderastrea siderea"
  )
  
  # Aggregate data by species and summarize tissue as cover values for TP = 01 only
  summary_species_cover <- survey_tissue %>%
    mutate(TP = sprintf("%02d", dense_rank(date))) %>% # Format time points as "01", "02", etc.
    filter(TP == "01") %>% # Only include data for TP = 01
    group_by(site, date, TP, spp, susc) %>% # Group by species and TP
    mutate(
      is_susceptible = ifelse(is.na(inftiss) & (is.na(cum_percloss) | cum_percloss != 100), 1, 0),  
      is_infected = ifelse(inftiss > 0, 1, 0),  
      is_dead = ifelse(cum_percloss == 100, 1, 0)  
  ) %>%
  summarise(
    tottiss = sum(remaintiss, inftiss, cum_tissloss, na.rm = TRUE), # Total tissue
    sustiss = sum(remaintiss, na.rm = TRUE),  # Susceptible tissue
    inftiss = sum(inftiss, na.rm = TRUE),     # Infected tissue
    deadtiss = sum(cum_tissloss, na.rm = TRUE), # Dead tissue
    totnum = n(),  # Count of all entries for each species
    susnum = sum(is_susceptible, na.rm = TRUE),  
    infnum = sum(is_infected, na.rm = TRUE),  
    deadnum = sum(is_dead, na.rm = TRUE),  
    .groups = 'drop'  # Ungroup after summarizing
  ) %>%
  complete(site, spp, TP, fill = list(
    tottiss = 0, sustiss = 0, inftiss = 0, deadtiss = 0,
    totnum = 0, susnum = 0, infnum = 0, deadnum = 0
  )) %>% # Ensure all combinations of site, spp, and TP exist, filling missing values with 0
  ungroup()
  
  # #checking it's working right
  # summation_cover_moderate <- summary_species_cover %>%
  #   filter(susc == "high", site == "near") %>%
  #   summarise(tot_tissue = sum(tottiss, na.rm = TRUE))
  
  # Ensure unique SA.cover.ratio per Site in susceptible_ref
  susceptible_ref_unique <- susceptible_ref %>%
    select(Site, SA.cover.ratio) %>%
    distinct(Site, .keep_all = TRUE)  # Remove duplicates, keeping the first occurrence

  # Calculate cover and add species full names
  final_table <- summary_species_cover %>%
    mutate(site = case_when(
      site == "mid"  ~ "Midchannel",
      site == "near" ~ "Nearshore",
      site == "off"  ~ "Offshore",
      TRUE ~ site  # Keep unchanged if not one of the above
    )) %>%
    left_join(susceptible_ref_unique, by = c("site" = "Site")) %>%
    mutate(
      cover = tottiss / SA.cover.ratio,
      species_name = recode(spp, !!!species_names)
    ) %>%
    group_by(species_name) %>%
    mutate(susc = ifelse(is.na(susc), first(na.omit(susc)), susc)) %>%
    ungroup() %>%
    select(site, susc, species_name, tottiss, totnum, cover)
  
  # #checking it's working right
  # summation_cover_moderate <- final_table %>%
  #   filter(susc == "moderate") %>%
  #   summarise(total_cover = sum(cover, na.rm = TRUE))
  
  # Summarize and assign the cover in the desired format (X, X, X) for each species
  final_table_summary <- final_table %>%
    group_by(susc, species_name, site) %>%
    summarise(
      total_SA = sum(tottiss, na.rm = TRUE),
      total_cover = sum(cover, na.rm = TRUE),
      total_num = sum(totnum, na.rm = TRUE),
      .groups = "drop"
    )
  
  # # Step 2: Format the cover and SA in (X, X, X) format for each species and site.
  # final_table_formatted <- final_table_summary %>%
  #   group_by(susc, species_name) %>%
  #   reframe(
  #     cover = paste(
  #       ifelse(round(total_cover[site == "Offshore"], 2) == 0, "<0.01", round(total_cover[site == "Offshore"], 2)), 
  #       ifelse(round(total_cover[site == "Midchannel"], 2) == 0, "<0.01", round(total_cover[site == "Midchannel"], 2)), 
  #       ifelse(round(total_cover[site == "Nearshore"], 2) == 0, "<0.01", round(total_cover[site == "Nearshore"], 2)), 
  #       sep = ", "
  #     ),
  #     number = paste(
  #       ifelse(round(total_num[site == "Offshore"], 2) == 0, "<0.01", round(total_num[site == "Offshore"], 2)), 
  #       ifelse(round(total_num[site == "Midchannel"], 2) == 0, "<0.01", round(total_num[site == "Midchannel"], 2)), 
  #       ifelse(round(total_num[site == "Nearshore"], 2) == 0, "<0.01", round(total_num[site == "Nearshore"], 2)), 
  #       sep = ", "
  #     ),
  #     SA = paste(
  #       ifelse(round(total_SA[site == "Offshore"], 2) == 0, "<0.01", round(total_SA[site == "Offshore"], 2)), 
  #       ifelse(round(total_SA[site == "Midchannel"], 2) == 0, "<0.01", round(total_SA[site == "Midchannel"], 2)), 
  #       ifelse(round(total_SA[site == "Nearshore"], 2) == 0, "<0.01", round(total_SA[site == "Nearshore"], 2)), 
  #       sep = ", "
  #     )
  #   ) %>%
  #   arrange(susc, species_name)
  
  # Step 2: Format the cover and SA in (X, X, X) format for each species and site.
  final_table_formatted <- final_table_summary %>%
    group_by(susc, species_name) %>%
    reframe(
      cover = paste(
        ifelse(total_cover[site == "Offshore"] == 0, "0.00", 
               ifelse(round(total_cover[site == "Offshore"], 2) == 0, "<0.01", round(total_cover[site == "Offshore"], 2))), 
        ifelse(total_cover[site == "Midchannel"] == 0, "0.00", 
               ifelse(round(total_cover[site == "Midchannel"], 2) == 0, "<0.01", round(total_cover[site == "Midchannel"], 2))), 
        ifelse(total_cover[site == "Nearshore"] == 0, "0.00", 
               ifelse(round(total_cover[site == "Nearshore"], 2) == 0, "<0.01", round(total_cover[site == "Nearshore"], 2))), 
        sep = ", "
      ),
      number = paste(
        ifelse(total_num[site == "Offshore"] == 0, "0", 
               ifelse(round(total_num[site == "Offshore"], 2) == 0, "<0.01", round(total_num[site == "Offshore"], 2))), 
        ifelse(total_num[site == "Midchannel"] == 0, "0", 
               ifelse(round(total_num[site == "Midchannel"], 2) == 0, "<0.01", round(total_num[site == "Midchannel"], 2))), 
        ifelse(total_num[site == "Nearshore"] == 0, "0", 
               ifelse(round(total_num[site == "Nearshore"], 2) == 0, "<0.01", round(total_num[site == "Nearshore"], 2))), 
        sep = ", "
      ),
      SA = paste(
        ifelse(total_SA[site == "Offshore"] == 0, "0.00", 
               ifelse(round(total_SA[site == "Offshore"], 2) == 0, "<0.01", round(total_SA[site == "Offshore"], 2))), 
        ifelse(total_SA[site == "Midchannel"] == 0, "0.00", 
               ifelse(round(total_SA[site == "Midchannel"], 2) == 0, "<0.01", round(total_SA[site == "Midchannel"], 2))), 
        ifelse(total_SA[site == "Nearshore"] == 0, "0.00", 
               ifelse(round(total_SA[site == "Nearshore"], 2) == 0, "<0.01", round(total_SA[site == "Nearshore"], 2))), 
        sep = ", "
      )
    ) %>%
    arrange(susc, species_name)
  
  
  # Step 3: Add in total and grand total
  final_table_output <- final_table_formatted %>%
    mutate(across(where(is.numeric), ~sprintf("%.2f", .))) %>%
    # Add the "Total" row
    bind_rows(
      final_table_summary %>%
        group_by(susc) %>%
        summarise(
          species_name = "Total",
          cover = paste(
            ifelse(round(sum(total_cover[site == "Offshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_cover[site == "Offshore"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_cover[site == "Midchannel"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_cover[site == "Midchannel"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_cover[site == "Nearshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_cover[site == "Nearshore"], na.rm = TRUE), 2))),
            sep = ", "
          ),
          number = paste(
            sum(total_num[site == "Offshore"], na.rm = TRUE),
            sum(total_num[site == "Midchannel"], na.rm = TRUE),
            sum(total_num[site == "Nearshore"], na.rm = TRUE),
            sep = ", "
          ),
          SA = paste(
            ifelse(round(sum(total_SA[site == "Offshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_SA[site == "Offshore"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_SA[site == "Midchannel"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_SA[site == "Midchannel"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_SA[site == "Nearshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_SA[site == "Nearshore"], na.rm = TRUE), 2))),
            sep = ", "
          ),
          .groups = "drop"
        ) %>%
        ungroup()
    ) %>%
    # Add the "Grand Total" row at the bottom
    bind_rows(
      final_table_summary %>%
        summarise(
          species_name = "Grand Total",
          cover = paste(
            ifelse(round(sum(total_cover[site == "Offshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_cover[site == "Offshore"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_cover[site == "Midchannel"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_cover[site == "Midchannel"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_cover[site == "Nearshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_cover[site == "Nearshore"], na.rm = TRUE), 2))),
            sep = ", "
          ),
          number = paste(
            sum(total_num[site == "Offshore"], na.rm = TRUE),
            sum(total_num[site == "Midchannel"], na.rm = TRUE),
            sum(total_num[site == "Nearshore"], na.rm = TRUE),
            sep = ", "
          ),
          SA = paste(
            ifelse(round(sum(total_SA[site == "Offshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_SA[site == "Offshore"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_SA[site == "Midchannel"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_SA[site == "Midchannel"], na.rm = TRUE), 2))),
            ifelse(round(sum(total_SA[site == "Nearshore"], na.rm = TRUE), 2) == 0, "<0.01", sprintf("%.2f", round(sum(total_SA[site == "Nearshore"], na.rm = TRUE), 2))),
            sep = ", "
          ),
          .groups = "drop"
        )
    ) %>%
    # Sort so that "Total" and "Grand Total" are at the bottom of each 'susc' group
    arrange(susc, species_name == "Total", species_name == "Grand Total", species_name)
  
  # #option to remove comma separation
  # final_table_output2 <- final_table_output %>%
  #   mutate(
  #     cover = gsub(", ", "\t", cover)  # Replace commas with tabs
  #   )
  
  # final_table_output3 = final_table_output %>%
  #   select(cover)
  
  # Use the "here" package to specify the file path
  output_path <- here("output", "formatted_table.csv")
  
  # Write the table to a CSV file without row names
  write.csv(final_table_output, output_path, row.names = FALSE)

  # Print the file path for confirmation
  cat("Table saved to: ", output_path, "\n")
  # Print the final table
  print(final_table_output)
  
  
  ################################## Figure prep ##################################
  # # List of all data frames
  # dfs <- list(
  #   output.basic.midchannel = output.basic.midchannel,
  #   output.basic.midchannel.DHW = output.basic.midchannel.DHW,
  #   output.near.to.mid.basic = output.basic.midchannel.transfer,
  #   output.basic.offshore = output.basic.offshore,
  #   output.basic.offshore.DHW = output.basic.offshore.DHW,
  #   output.near.to.off.basic = output.basic.offshore.transfer,
  #   output.basic.nearshore = output.basic.nearshore,
  #   output.basic.nearshore.DHW = output.basic.nearshore.DHW,
  #   output.off.to.near.basic = output.basic.nearshore.transfer,
  # 
  #   output.midchannel = output.midchannel,
  #   output.offshore = output.offshore,
  #   output.nearshore = output.nearshore,
  #   output.near.to.off.multi = output.near.to.off.multi
  # )
  # 
  # # Function to extract metadata from names
  # extract_metadata <- function(name) {
  #   site <- case_when(
  #     grepl("midchannel", name) ~ "Midchannel",
  #     grepl("offshore", name) ~ "Offshore",
  #     grepl("nearshore", name) ~ "Nearshore",
  #     TRUE ~ "Unknown"
  #   )
  #   
  #   host <- ifelse(grepl("basic", name), "Single-host", "Multi-host")
  #   
  #   type <- case_when(
  #     grepl("DHW", name) ~ "DHW",
  #     grepl("transfer", name) | grepl("\\.to\\.", name) ~ "Projected",
  #     TRUE ~ "Fitted"
  #   )
  #   
  #   tibble(Site = site, Host = host, Type = type)
  # }
  # 
  # # Concatenate all data frames with new metadata columns
  # data_fig3 <- imap_dfr(dfs, ~ mutate(.x, Site = extract_metadata(.y)$Site, 
  #                                        Host = extract_metadata(.y)$Host, 
  #                                        Type = extract_metadata(.y)$Type))
  
  dfs <- list(
    output.basic.midchannel = output.basic.midchannel,
    output.basic.midchannel.DHW = output.basic.midchannel.DHW,
    output.near.to.mid.basic = output.basic.midchannel.transfer,
    output.basic.offshore = output.basic.offshore,
    output.basic.offshore.DHW = output.basic.offshore.DHW,
    output.near.to.off.basic = output.basic.offshore.transfer,
    output.basic.nearshore = output.basic.nearshore,
    output.basic.nearshore.DHW = output.basic.nearshore.DHW,
    output.off.to.near.basic = output.basic.nearshore.transfer,
    
    output.midchannel = output.midchannel,
    output.offshore = output.offshore,
    output.nearshore = output.nearshore,
    output.near.to.off.multi = output.near.to.off.multi
  )
  
  # Function to extract metadata from names
  extract_metadata <- function(name) {
    site <- case_when(
      grepl("\\.to\\.", name) ~ case_when(
        grepl("to\\.mid", name) ~ "Midchannel",
        grepl("to\\.off", name) ~ "Offshore",
        grepl("to\\.near", name) ~ "Nearshore",
        TRUE ~ "Unknown"
      ),
      grepl("midchannel", name) ~ "Midchannel",
      grepl("offshore", name) ~ "Offshore",
      grepl("nearshore", name) ~ "Nearshore",
      TRUE ~ "Unknown"
    )
    
    host <- ifelse(grepl("basic", name), "Single-host", "Multi-host")
    
    type <- case_when(
      grepl("DHW", name) ~ "DHW",
      grepl("transfer", name) | grepl("\\.to\\.", name) ~ "Projected",
      TRUE ~ "Fitted"
    )
    
    tibble(Site = site, Host = host, Type = type)
  }
  
  # Concatenate all data frames with new metadata columns
  data_fig3 <- imap_dfr(dfs, ~ mutate(.x, 
                                      Site = extract_metadata(.y)$Site, 
                                      Host = extract_metadata(.y)$Host, 
                                      Type = extract_metadata(.y)$Type))
  
  
  
  ################################## Figure 2  ##################################
  
  # NOTE - need to make the top row of fig2 proportional (see teams chat Jan 30 2025)
  
  #version where "Total" is its own thing, not part of the legend
  obs.total.figures <- obs.total %>%
    rename(Susceptibility = Category) %>%
    mutate(Site = factor(Site, levels = c("Offshore", "Midchannel", "Nearshore"))) %>%
    mutate(Compartment = case_when(Compartment == "Dead" ~ "Recovered", TRUE ~ Compartment)) %>%
    mutate(Compartment = factor(Compartment, levels = c("Susceptible", "Infected", "Recovered")))
  
  obs.multi.figures <- obs.multi %>%
   mutate(Susceptibility = factor(Susceptibility, levels = c("Low", "Moderate", "High"))) %>%
   mutate(Site = factor(Site, levels = c("Offshore", "Midchannel", "Nearshore"))) %>%
   mutate(Compartment = case_when(Compartment == "Dead" ~ "Recovered", TRUE ~ Compartment)) %>%
   mutate(Compartment = factor(Compartment, levels = c("Susceptible", "Infected", "Recovered")))
  
  # Find the first date after the initial epidemic wave where SST exceeds threshold set during modeling
  target_date <- DHW.CRW %>%
   filter(date > as.Date(date_threshold) & SST.90th_HS > SST_threshold_value) %>%
   slice(1) %>%
   pull(date)
  
  target_days <- obs.total.figures %>%
   filter(date == target_date) %>%
   select(Site, days.inf.site)
  
  # Define max days across all sites
  #   - important for creating tidy, comparable time-based plots when necessary
  max_days_all_sites <- max(obs.total.figures$days.inf.site, na.rm = TRUE)
  
  # Absolute Values Plot (Free Scales)
  absolute_plot <- ggplot() +
   geom_ribbon(data = obs.total.figures %>% filter(Compartment == "Infected"),
               aes(x = days.inf.site, ymin = 0, ymax = tissue), 
               fill = "gray80", alpha = 0.8) +
   geom_line(data = obs.multi.figures %>% filter(Compartment == "Infected"),
             aes(x = days.inf.site, y = tissue, color = Susceptibility),
             linewidth = 0.60) +
   scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                 "Moderate" = "#FFD700", # Yellow
                                 "High" = "#FF1493")) +  # Deep Pink
   facet_wrap(~ Site, scales = "free") +
   # scale_x_continuous(expand = c(0, 0)) +
   xlab("Day of outbreak") +
   ylab("Tissue surface area (m²)") +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +
   scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
   scale_linetype_manual(values = c("dotdash", "longdash", "solid")) +
   theme_classic(base_family = "Georgia") +
   theme(legend.position = "bottom",
         axis.title = element_text(size = 9),
         axis.text = element_text(size = 7),
         strip.text = element_text(size = 8),
         legend.text = element_text(size = 7),
         legend.title = element_text(size = 9),
         legend.key.height = unit(0, "cm"))
  
  max_value.fig2 <- max(obs.total.figures %>% filter(Compartment == "Infected") %>% pull(tissue), na.rm = TRUE)
  
  # Relative Values Plot (Fixed Scales)
  relative_plot <- ggplot() +
   geom_ribbon(data = obs.total.figures %>% filter(Compartment == "Infected"),
               aes(x = days.inf.site, ymin = 0, ymax = tissue), 
               fill = "gray80", alpha = 0.8) +
   geom_line(data = obs.multi.figures %>% filter(Compartment == "Infected"),
             aes(x = days.inf.site, y = tissue, color = Susceptibility),
             linewidth = 0.60) +
   scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                 "Moderate" = "#FFD700", # Yellow
                                 "High" = "#FF1493")) +  # Deep Pink
   # facet_wrap(~ Site, scales = "free_x") +  # Fixed scales
   facet_wrap(~ Site, scales = "free") +  # Fixed scales
   xlab("Day of outbreak") +
   ylab("Tissue surface area (m²)") +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +  # Apply global max
   scale_y_continuous(limits = c(0, max_value.fig2), labels = scales::label_number(accuracy = 0.001)) +
   scale_linetype_manual(values = c("dotdash", "longdash", "solid")) +
   theme_classic(base_family = "Georgia") +
   theme(legend.position = "bottom",
         axis.title = element_text(size = 9),
         axis.text = element_text(size = 7),
         # strip.text = element_text(size = 8),
         strip.text = element_blank(),
         legend.text = element_text(size = 7),
         legend.title = element_text(size = 9),
         legend.key.height = unit(0, "cm"))
  
  fig2 <- (absolute_plot / relative_plot) +
   plot_layout(guides = "collect", axes = "collect") &  # Collect the legends
   theme(legend.position = "bottom",
         legend.box.spacing = unit(0, "cm"),
         plot.margin = margin(t = 5, r = 5, b = 0, l = 5))
  # legend.margin = margin(t = -5, r = 0, b = -5, l = 0)
  # theme(legend.position = "bottom", strip.placement = "outside", strip.text = element_text(size = 15))   # Set legend position
  
  #plot_layout(heights = c(1, 1.2))
  
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  # quartz(h = 6, w = 7.087)
  quartz(h = 3, w = 5)
  
  fig2
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig2.pdf"), type = "pdf")
  
  # Close the Quartz device
  dev.off()
  
  # # Save the plot with specified width and height, using here for the file path
  # # ggsave(here("output", "fig2.pdf"), plot = fig2)
  # ggsave(here("output", "fig2.pdf"), plot = combined_plot, width = 5, height = 3, units = "in")   #18 inches max width
  
  
  ################################## Figure 3 ##################################
  p.fit.nearshore.basic / p.fit.offshore.basic / p.fit.near.to.off.basic
  p.fit.nearshore.multi / p.fit.offshore.multi / p.fit.near.to.off.multi
  p.fit.nearshore.basic.DHW
  
  #this will be a 9-panel figure
  #    leftmost column is nearshore-fitted, offshore-fitted, and offshore-projected for single-host model
  #    middle column is the same as the leftmost, but for the multi-host model
  #    rightmost column is the same as leftmost, but for the single-host model with inclusion of thermal stress
  
  
  # NOTE - could choose to panel columns / rows with facet_wrap; would provide some flexibility on a fixed y-axis. then consider inserts with higher detail ?
  
  
  # # COLUMN 1
  # p.fit.nearshore.basic / p.fit.offshore.basic / p.fit.near.to.off.basic
  # p.fit.nearshore.basic + p.fit.offshore.basic + p.fit.near.to.off.basic
  
  data_fig3 = data_fig3 %>%
    mutate(Compartment = case_when(
      Compartment == "Dead" ~ "Recovered",
      TRUE ~ Compartment) 
    )
  
  max_value.fig3 <- max(data_fig3 %>% 
                     select(tissue) %>%
                     pull(tissue), na.rm = TRUE)
  
  site.loop = 'Nearshore'
  p.fit.nearshore.single.figures =
  ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Fitted"),
         aes(x = days.model, y = tissue, group = Compartment)) +
   xlab("Day of outbreak") +
   ylab("Surface area of tissue (m2)") +
   # ggtitle(paste0(site.loop, " - Fitted")) +
   geom_line(aes(group = Compartment), color = "black", linewidth = 0.4) +
   geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = 1.3) +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +
   # scale_y_continuous(limits = c(0, max_value.fig3)) +
   scale_shape_manual(values = c("Susceptible" = 16,
                                 "Infected" = 17,
                                 "Recovered" = 15),
                      name = "") + # Removes "Compartment" from the legend title
    # scale_shape_manual(values = c("Susceptible" = 16,
    #                               "Infected" = 1,
    #                               "Recovered" = 13),
    #                    name = "") + # Removes "Compartment" from the legend title
   theme_classic(base_family = 'Georgia')
  
  site.loop = 'Offshore'
  p.fit.offshore.single.figures =
  ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Fitted"), aes(x = days.model, y = tissue, group = Compartment)) +
   xlab("Day of outbreak") +
   ylab("Surface area of tissue (m2)") +
   # ggtitle(paste0(site.loop, " - Fitted")) +
   geom_line(aes(group = Compartment), color = "black", linewidth = 0.4) +
   geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = 1.3) +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +
   # scale_y_continuous(limits = c(0, max_value.fig3)) +
   scale_shape_manual(values = c("Susceptible" = 16,
                                 "Infected" = 17,
                                 "Recovered" = 15),
                      name = "") + # Removes "Compartment" from the legend title
   theme_classic(base_family = 'Georgia')
  
  p.fit.near.to.off.single.figures =
  ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Projected"), aes(x = days.model, y = tissue, group = Compartment)) +
   xlab("Day of outbreak") +
   ylab("Surface area of tissue (m2)") +
   # ggtitle(paste0(site.loop, " - Projected")) +
   geom_line(aes(group = Compartment), color = "black", linewidth = 0.4) +
   geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = 1.3) +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +
   # scale_y_continuous(limits = c(0, max_value.fig3)) +
   scale_shape_manual(values = c("Susceptible" = 16,
                                 "Infected" = 17,
                                 "Recovered" = 15),
                      name = "") + # Removes "Compartment" from the legend title
   theme_classic(base_family = 'Georgia')
  
  fig3_col1 = (p.fit.nearshore.single.figures / p.fit.offshore.single.figures / p.fit.near.to.off.single.figures) + 
   plot_layout(guides = "collect", axes = "collect") &  # Collect the legends
   theme(legend.position = "bottom",
         legend.box.spacing = unit(0, "cm"),
         legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
         legend.text = element_text(size = 8),  # Reduce legend text size
         axis.title.x = element_blank(),  # Remove x-axis label
         # legend.spacing.y = unit(0, "cm"),    # Reduce space between legend items (doesn't work)
         plot.margin = margin(t = 5, r = 5, b = 0, l = 5) #line = element_line(linewidth = 1)  # Change the universal linewidth
   )
   # ) &
   #  labs(x = "Day of outbreak")  # Apply x-axis label at the composite level
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  # quartz(h = 6, w = 7.087)
  quartz(h = 5, w = 3)
  
  fig3_col1
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig1.pdf"), type = "pdf")
  
  # Close the Quartz device
  dev.off()
  
  
  
  
  #COLUMN 2
  # p.fit.nearshore.multi / p.fit.offshore.multi / p.fit.near.to.off.multis
  
  # NOTE - could consider a shaded ribbon to display the total accumulated recovered compartment
  #         (would accentuate how high-susceptibility corals dominate the mortality)
  
  site.loop = 'Nearshore'
  p.fit.nearshore.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility, shape = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m2)") +
    # ggtitle(paste0(site.loop, " - Fitted")) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility), size = 1.3) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = 'Georgia')
  
  site.loop = 'Offshore'
  p.fit.offshore.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility, shape = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m2)") +
    # ggtitle(paste0(site.loop, " - Fitted")) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility), size = 1.3) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = 'Georgia')
  
  p.fit.near.to.off.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Multi-host", Type == "Projected"), aes(x = days.model, y = tissue, color = Susceptibility, shape = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m2)") +
    # ggtitle(paste0(site.loop, " - Projected")) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility), size = 1.3) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    # scale_y_continuous(limits = c(0, max_value.fig3)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = 'Georgia')
  
  fig3_col2 = (p.fit.nearshore.multi.figures / p.fit.offshore.multi.figures / p.fit.near.to.off.multi.figures) + 
    plot_layout(guides = "collect", axes = "collect") &  # Collect the legends
    theme(legend.position = "bottom",
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
          legend.text = element_text(size = 8),  # Reduce legend text size
          axis.title.y = element_blank(),  # Remove y-axis label
          # axis.title.x = element_blank(),  # Remove x-axis label
          # legend.spacing.y = unit(0, "cm"),    # Reduce space between legend items (doesn't work)
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)#, line = element_line(linewidth = 1)  # Change the universal linewidth
    )
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  # quartz(h = 6, w = 7.087)
  quartz(h = 5, w = 3)
  
  fig3_col2
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig1.pdf"), type = "pdf")
  
  # Close the Quartz device
  dev.off()
  
  
  
  
  
  
  
  #COLUMN 3
  # p.fit.nearshore.single.DHW
  
  
  # p.fit.nearshore.single.DHW = ggplot(data = output.single.nearshore.DHW, aes(days.model, tissue, colour = Compartment)) +
  #   xlab("Day of outbreak") +
  #   ylab("Surface area of tissue (m2)") +
  #   ggtitle(paste(c("", site.loop, ' - Fitted'), collapse="")) +
  #   geom_line() +
  #   geom_point(data = obs.total %>% filter(Site == site.loop), aes(days.inf.site, tissue, colour = Compartment)) +
  #   scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  #   annotate(geom = "table", x = min(output.single.nearshore.DHW$days.model), y = min(output.single.nearshore.DHW$tissue)*0.7, label = list(tab.nearshore),
  #            vjust = val.vjust, hjust = val.hjust, family = 'Georgia') +
  #   theme_classic(base_family = 'Georgia') +
  #   theme(panel.background = element_rect(fill = "gray90"))
  
  
  site.loop = 'Nearshore'
  p.fit.nearshore.single.DHW.figures = 
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "DHW"),
           aes(x = days.model, y = tissue, group = Compartment)) +
    # xlab("Day of outbreak") +
    ylab("Surface area of tissue (m2)") +
    # ggtitle(paste0(site.loop, " - Fitted")) +
    geom_line(aes(group = Compartment), color = "black", linewidth = 0.4) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = 1.3) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    # scale_y_continuous(limits = c(0, max_value.fig3)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    # scale_shape_manual(values = c("Susceptible" = 16,
    #                               "Infected" = 1,
    #                               "Recovered" = 13),
    #                    name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = 'Georgia')
  
  site.loop = 'Offshore'
  p.fit.offshore.single.DHW.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "DHW"), aes(x = days.model, y = tissue, group = Compartment)) +
    # xlab("Day of outbreak") +
    ylab("Surface area of tissue (m2)") +
    # ggtitle(paste0(site.loop, " - Fitted")) +
    geom_line(aes(group = Compartment), color = "black", linewidth = 0.4) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = 1.3) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    # scale_y_continuous(limits = c(0, max_value.fig3)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = 'Georgia')
  
  # NOTE - have not attempted to project an epidemic with DHW information ... I guess I could try ?
  p.fit.near.to.off.single.DHW.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Projected"), aes(x = days.model, y = tissue, group = Compartment)) +
    # xlab("Day of outbreak") +
    ylab("Surface area of tissue (m2)") +
    # ggtitle(paste0(site.loop, " - Projected")) +
    # geom_line(aes(group = Compartment), color = "black", linewidth = 0.4) +
    # geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = 1.3) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    # scale_y_continuous(limits = c(0, max_value.fig3)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = 'Georgia')
  
  
  # # Create a blank plot for p.fit.near.to.off.single.figures
  # p.fit.near.to.off.single.figures <- ggplot() + 
  #   theme_void()  # Ensures it appears blank with no axes or labels
  
  
  fig3_col3 = (p.fit.nearshore.single.DHW.figures / p.fit.offshore.single.DHW.figures / p.fit.near.to.off.single.DHW.figures) +
  # fig3_col3 = (p.fit.nearshore.single.figures / p.fit.offshore.single.figures) + 
    plot_layout(guides = "collect", axes = "collect") &  # Collect the legends
    theme(legend.position = "bottom",
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
          legend.text = element_text(size = 8),  # Reduce legend text size
          axis.title.y = element_blank(),  # Remove y-axis label
          axis.title.x = element_blank(),  # Remove x-axis label
          # legend.spacing.y = unit(0, "cm"),    # Reduce space between legend items (doesn't work)
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)#, line = element_line(linewidth = 1)  # Change the universal linewidth
    )
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  # quartz(h = 6, w = 7.087)
  quartz(h = 5, w = 3)
  
  fig3_col3
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig1.pdf"), type = "pdf")
  
  # Close the Quartz device
  dev.off()
  
  
  
  
  
  
  
  #COMPOSITE
  fig3_col1 = fig3_col1 & theme(legend.position = "none")  # Remove legend from fig3_col1
  fig3_col3 = fig3_col3 & theme(legend.position = "none")  # Remove legend from fig3_col1
  
  # #failed attempt at better formatting the shared legend
  # fig3_col2 = fig3_col2 & 
  #   theme(legend.position = "bottom",
  #         legend.direction = "vertical",  # Stack legends vertically
  #         legend.box = "vertical",        # Ensure they are treated as separate boxes
  #         legend.spacing.y = unit(0.3, "cm"))  # Adjust vertical spacing
  
  fig3 = (fig3_col1 | fig3_col2) +
  # fig3 = (fig3_col1 | fig3_col2 | fig3_col3) +
  # fig3 = (fig3_col1 | fig3_col3) +
    plot_layout(guides = "collect", axes = "collect")  #&
    # ylab("Surface area of tissue (m²)") &
    # xlab("Day of ...  ")
    # labs(x = "Day of outbreak",
    #      y = "Surface area of tissue (m²)")
  
  # add annotation separately (must be done to avoid breaking axes collect)
  fig3 <- fig3 + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  # quartz(h = 6, w = 7.087)
  quartz(h = 4, w = 6)
  
  fig3
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig1.pdf"), type = "pdf")
  
  # # Close the Quartz device
  # dev.off()
  
  
  ################################## Save output ##################################

  # save.image(file = here("output", "tables_figures_workspace.RData"))
  