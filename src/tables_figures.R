  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(patchwork)
  library(extrafont)
  library(ggthemes)
  library(scales)
  library(kableExtra)
  library(cowplot)
  
  ################################## Set-up ##################################
  
  # #this was required for me to run at least once on an M1 Macbook (last updated December 2024)
  # font_import() 
  # loadfonts(device = "pdf")
  
  #very helpful:
  # https://www.rforecology.com/post/exporting-plots-in-r/
  
  #import workspace from upstream script
  load(here("output/error_eval_workspace.RData"))
  
  # ################################## Proportions ##################################
  # 
  # 
  # # SINGLE-HOST
  # #comparison of peak infection prevalence
  # max(output.basic.nearshore.full.raw$I) / max(output.basic.nearshore.full.raw$S) * 100
  # max(output.basic.offshore.transfer.raw$I) / max(output.basic.offshore.transfer.raw$S) * 100
  # max(output.basic.offshore.full.raw$I) / max(output.basic.offshore.full.raw$S) * 100
  # 
  # #comparison of proportions of tissue lost
  # max(output.basic.nearshore.full.raw$R) / max(output.basic.nearshore.full.raw$S) * 100
  # max(output.basic.offshore.transfer.raw$R) / max(output.basic.offshore.transfer.raw$S) * 100
  # max(output.basic.offshore.full.raw$R) / max(output.basic.offshore.full.raw$S) * 100
  # 
  # 
  # 
  # # MULTI-HOST
  # 
  # #test of proportions - confirms that though shapes of projected outbreaks look identical to source fit, they aren't exactly
  # max(output.raw.nearshore$I.LS) / max(output.raw.nearshore$I.HS)
  # max(output.raw.near.to.off.multi$I.LS) / max(output.raw.near.to.off.multi$I.HS)
  # #
  # #comparison of peak infection prevalence
  # (max(output.raw.nearshore$I.LS) + max(output.raw.nearshore$I.MS) + max(output.raw.nearshore$I.HS)) /
  #   (max(output.raw.nearshore$S.LS) + max(output.raw.nearshore$S.MS) + max(output.raw.nearshore$S.HS)) *
  #   100
  # (max(output.raw.near.to.off.multi$I.LS) + max(output.raw.near.to.off.multi$I.MS) + max(output.raw.near.to.off.multi$I.HS)) /
  #   (max(output.raw.near.to.off.multi$S.LS) + max(output.raw.near.to.off.multi$S.MS) + max(output.raw.near.to.off.multi$S.HS)) *
  #   100
  # #comparison of proportions of tissue lost
  # max(output.raw.nearshore$R.LS) / max(output.raw.nearshore$S.LS)
  # max(output.raw.near.to.off.multi$R.LS) / max(output.raw.near.to.off.multi$S.LS)
  # max(output.raw.nearshore$R.MS) / max(output.raw.nearshore$S.MS)
  # max(output.raw.near.to.off.multi$R.MS) / max(output.raw.near.to.off.multi$S.MS)
  # max(output.raw.nearshore$R.HS) / max(output.raw.nearshore$S.HS)
  # max(output.raw.near.to.off.multi$R.HS) / max(output.raw.near.to.off.multi$S.HS)
  # 
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
    output.basic.midchannel = output.basic.midchannel.full,
    # output.basic.midchannel.DHW = output.basic.midchannel.DHW,
    output.near.to.mid.basic = output.basic.midchannel.transfer,
    output.basic.offshore = output.basic.offshore.full,
    # output.basic.offshore.DHW = output.basic.offshore.DHW,
    output.near.to.off.basic = output.basic.offshore.transfer,
    output.basic.nearshore = output.basic.nearshore.full,
    # output.basic.nearshore.DHW = output.basic.nearshore.DHW,
    output.off.to.near.basic = output.basic.nearshore.transfer,
    
    output.midchannel = output.midchannel,
    output.offshore = output.offshore,
    output.nearshore = output.nearshore,
    output.near.to.off.multi = output.near.to.off.multi,
    output.near.to.mid.multi = output.near.to.off.multi
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
  
  
  
  ################################## Figure S1 ##################################
  
  linewidths =  1.3 #0.4 #0.75 is roughly 1 pt. ggplot measures these in mm, not points
  symbsizes = 0.4
  titlesize = 10 #9  #text sizes in ggplot are actually in units of points, when specified using element_text
  textsize = 9
  palette = 'viridis'
  symbalpha = 1
  tensoralpha = 0.7
  GAMalpha = 1
  
  figS1.main = ggplot(SA_predictions, aes(x = x, y = y)) +
    geom_point(color = "black", size = symbsizes, alpha = symbalpha, show.legend = FALSE) + #alpha = 0.7
    # geom_line(aes(y = pred_linear, color = "Linear"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam, color = "GAM"), linewidth = linewidths, alpha = GAMalpha) +
    geom_line(aes(y = pred_gam_gamma_tp, color = "GAM with Tensor Product"), linewidth = linewidths, alpha = tensoralpha) +
    # geom_line(aes(y = pred_gam_gamma_ps, color = "GAM with P-splines"), linewidth = linewidths, alpha = 0.7) +
    labs(
      x = "Maximum colony diameter (cm)",
      y = "Colony surface area (m²)",
      color = "Model"
    ) +
    scale_color_manual(
      values = c(
        "GAM" = "#00BFC4", # "#1B9E77",    # red        
        "GAM with Tensor Product" = "#E7298A" #darkorange
      )
    ) +
    theme_classic(base_family = "Georgia") +
    theme(
      axis.title = element_text(size = titlesize, color = 'black'),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.ticks = element_line(color = "black"),
      legend.text = element_text(size = textsize),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box.spacing = unit(0, "cm"),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  

  
  
  # Create the inset (zoomed-in) plot
  # Note: The inset plot won't have its own legend to avoid duplication
  figS1.inset = ggplot(SA_predictions, aes(x = x, y = y)) +
    geom_point(color = "black", size = symbsizes, alpha = symbalpha) +
    geom_line(aes(y = pred_gam, color = "GAM"), linewidth = linewidths, alpha = GAMalpha) +
    geom_line(aes(y = pred_gam_gamma_tp, color = "GAM with Tensor Product"), linewidth = linewidths, alpha = tensoralpha) +
    labs(
      x = "Maximum colony diameter (cm)",
      y = "Colony surface area (m²)"
    ) +
    scale_color_manual(
      values = c(
        "GAM" = "#00BFC4", # "#1B9E77",            #red
        "GAM with Tensor Product" = "#E7298A"  #darkorange
      )
    ) +
    xlim(0, 20) +
    ylim(0, 0.10) +
    theme_classic(base_family = "Georgia") +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = textsize - 2, color = 'black'),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"  # No legend for the inset
    )
  
  # Combine the plots using cowplot's ggdraw and draw_plot functions
  figS1 = ggdraw() +
    draw_plot(figS1.main) +
    # Add the inset plot - adjust these parameters as needed:
    # x, y: position of the bottom-left corner of the inset (0,0 is bottom-left, 1,1 is top-right)
    # width, height: size of the inset as a proportion of the main plot
    draw_plot(figS1.inset, x = 0.15, y = 0.55, width = 0.45, height = 0.40) +
    # Optional: Add a border around the inset
    draw_label("", x = 0.15, y = 0.55, hjust = 0, vjust = 0, size = 12)
  
  

  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  quartz(h = 3, w = 3.35)
  # quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)
  
  figS1
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "figS1.pdf"), type = "pdf")
  # 
  # #ggplot-export to image
  # ggsave(filename = here("output", "figS1.png"), device = "png", width = 3.35, height = 3, dpi = 1200)
  
  dev.off()
  
  ################################## Figure 2  ##################################
  
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
  
  #calculate proportion of tissue in each SIR stage by timepoint
  site_mapping <- data.frame(
    Site = c("Midchannel", "Nearshore", "Offshore"),
    Site_short = c("mid", "near", "off"),
    stringsAsFactors = FALSE
  )
  obs_total_with_n <- obs.total.figures %>%
    # filter(Compartment == "Infected") %>%
    left_join(site_mapping, by = "Site") %>%
    left_join(site_ref, by = c("Site_short" = "Site")) %>%
    mutate(Site = factor(Site, levels = c("Offshore", "Midchannel", "Nearshore"))) %>%
    group_by(Site) %>%
    mutate(tissue_normalized = tissue / N.site.y * 100)
  
  obs_with_n <- obs.multi.figures %>% 
    # filter(Compartment == "Infected") %>%
    left_join(site_mapping, by = "Site") %>%  # Add the short names
    left_join(site_ref, by = c("Site_short" = "Site")) %>%  # Join with site_ref using short names
    mutate(Site = factor(Site, levels = c("Offshore", "Midchannel", "Nearshore"))) %>%
    mutate(tissue_normalized = tissue / N.site.y * 100)
  
  #calculate proportions of *colonies* in each SIR stage by timepoint
  obs_total_with_n = obs_total_with_n %>%
    group_by(Site) %>%
    mutate(max_susceptible_count = max(count)) %>%
    # filter(Compartment == "Infected") %>%
    mutate(prevalence = count / max_susceptible_count * 100) %>%
    ungroup()
  
  obs_with_n = obs_with_n %>%
    group_by(Site) %>%
    mutate(max_susceptible_count = max(count)) %>%
    # filter(Compartment == "Infected") %>%
    mutate(prevalence = count / max_susceptible_count * 100) %>%
    ungroup()
  
  #pull maximum proportion infected and prevalence of infected colonies per site, for reference
  max_prop_prev = obs_total_with_n %>%
    filter(Compartment != 'Susceptible') %>%
    group_by(Site, Compartment) %>%
    summarize(max_tissue_normalized = max(tissue_normalized, na.rm = TRUE),
              max_prevalence = max(prevalence, na.rm = TRUE))
  
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
  proportion_plot <- ggplot() +
    geom_ribbon(data = obs_total_with_n %>% filter(Compartment == "Infected"),
                aes(x = days.inf.site, ymin = 0, ymax = tissue_normalized), 
                fill = "gray80", alpha = 0.8) +
    geom_line(data = obs_with_n %>% filter(Compartment == "Infected"),
              aes(x = days.inf.site, y = tissue_normalized, color = Susceptibility),
              linewidth = 1) + #0.6
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    facet_wrap(~ Site, scales = "free") +
    # scale_x_continuous(expand = c(0, 0)) +
    xlab("Day of outbreak") +
    ylab("% tissue infected") +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_linetype_manual(values = c("dotdash", "longdash", "solid")) +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 7),
          # Remove x-axis text but keep the line
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),  # This removes the tick marks
          strip.text = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          legend.key.height = unit(0, "cm"))
  
  # Relative Values Plot (Fixed Scales)
  max_value.fig2.row2 <- max(obs.total.figures %>% filter(Compartment == "Infected") %>% pull(tissue), na.rm = TRUE)
  
  relative_plot <- ggplot() +
   geom_ribbon(data = obs.total.figures %>% filter(Compartment == "Infected"),
               aes(x = days.inf.site, ymin = 0, ymax = tissue), 
               fill = "gray80", alpha = 0.8) +
   geom_line(data = obs.multi.figures %>% filter(Compartment == "Infected"),
             aes(x = days.inf.site, y = tissue, color = Susceptibility),
             linewidth = 1) + #0.6
   scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                 "Moderate" = "#FFD700", # Yellow
                                 "High" = "#FF1493")) +  # Deep Pink
   # facet_wrap(~ Site, scales = "free_x") +  # Fixed scales
   facet_wrap(~ Site, scales = "free") +  # Fixed scales
   xlab("Day of outbreak") +
   ylab("SA infected (m²)") +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +  # Apply global max
   scale_y_continuous(limits = c(0, max_value.fig2.row2), labels = scales::label_number(accuracy = 0.1)) +
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
  
  fig2 <- (proportion_plot / relative_plot) +
   plot_layout(guides = "collect") &
   # plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom",
         legend.box.spacing = unit(0, "cm"),
         plot.margin = margin(t = 5, r = 5, b = 0, l = 5))
  # legend.margin = margin(t = -5, r = 0, b = -5, l = 0)
  # theme(legend.position = "bottom", strip.placement = "outside", strip.text = element_text(size = 15))   # Set legend position
  
  #plot_layout(heights = c(1, 1.2))
  
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 2, w = 3.35)
  quartz(h = 3.5, w = 7.087)
  # quartz(h = 3, w = 5)
  
  fig2
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig2.pdf"), type = "pdf")
  # 
  # #ggplot-export to image
  # ggsave(filename = here("output", "fig2.png"), device = "png", width = 7.087, height = 3, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
  
  # # Save the plot with specified width and height, using here for the file path
  # # ggsave(here("output", "fig2.pdf"), plot = fig2)
  # ggsave(here("output", "fig2.pdf"), plot = fig2, width = 5, height = 3, units = "in")   #18 inches max width
  # ggsave(here("output", "fig2.png"), plot = fig2, width = 5, height = 3, units = "in")   #18 inches max width
  
  ################################## Figure 3 ##################################
  
  #this is a 6-panel figure
  #    left column is nearshore-fitted, offshore-fitted, and offshore-projected for single-host model
  #    right column is the same as the leftmost, but for the multi-host model

  # NOTE - could choose to panel columns / rows with facet_wrap; would provide some flexibility on a fixed y-axis. then consider inserts with higher detail ?
  
  linewidths =  0.75 #0.4 #0.75 is roughly 1 pt. ggplot measures these in mm, not points
  symbsizes = 1.3
  titlesize = 10 #9  #text sizes in ggplot are actually in units of points, when specified using element_text
  textsize = 9
  
  # COLUMN 1
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
   ylab("Surface area of tissue (m²)") +
   geom_line(aes(group = Compartment), color = "black", linewidth = linewidths) +
   geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = symbsizes) +
   scale_x_continuous(limits = c(0, max_days_all_sites)) +
   scale_shape_manual(values = c("Susceptible" = 16,
                                 "Infected" = 17,
                                 "Recovered" = 15),
                      name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          # axis.text.x = element_blank(), #remove x-axis text but keep the line
          # axis.ticks.x = element_blank(), #remove the tick marks
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
    
    # site.loop = 'Midchannel'
    # p.fit.midchannel.single.figures =
    #   ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Fitted"), aes(x = days.model, y = tissue, group = Compartment)) +
    #   xlab("Day of outbreak") +
    #   ylab("SA of tissue (m2)") +
    #   # ggtitle(paste0(site.loop, " - Fitted")) +
    #   geom_line(aes(group = Compartment), color = "black", linewidth = linewidths) +
    #   geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = symbsizes) +
    #   scale_x_continuous(limits = c(0, max_days_all_sites)) +
    #   # scale_y_continuous(limits = c(0, max_value.fig3)) +
    #   scale_shape_manual(values = c("Susceptible" = 16,
    #                                 "Infected" = 17,
    #                                 "Recovered" = 15),
    #                      name = "") + # Removes "Compartment" from the legend title
    #   theme_classic(base_family = 'Georgia')
    # 
    # p.fit.near.to.mid.single.figures =
    #   ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Projected"), aes(x = days.model, y = tissue, group = Compartment)) +
    #   xlab("Day of outbreak") +
    #   ylab("Surface area of tissue (m2)") +
    #   # ggtitle(paste0(site.loop, " - Projected")) +
    #   geom_line(aes(group = Compartment), color = "black", linewidth = linewidths) +
    #   geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = symbsizes) +
    #   scale_x_continuous(limits = c(0, max_days_all_sites)) +
    #   # scale_y_continuous(limits = c(0, max_value.fig3)) +
    #   scale_shape_manual(values = c("Susceptible" = 16,
    #                                 "Infected" = 17,
    #                                 "Recovered" = 15),
    #                      name = "") + # Removes "Compartment" from the legend title
    #   theme_classic(base_family = 'Georgia')
    
  site.loop = 'Offshore'
  p.fit.offshore.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Fitted"),
           aes(x = days.model, y = tissue, group = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(aes(group = Compartment), color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          # axis.text.x = element_blank(), #remove x-axis text but keep the line
          # axis.ticks.x = element_blank(), #remove the tick marks
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
    site.loop = 'Offshore'
    p.fit.near.to.off.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Single-host", Type == "Projected"), aes(x = days.model, y = tissue, group = Compartment)) +
     xlab("Day of outbreak") +
     ylab("Surface area of tissue (m²)") +
     geom_line(aes(group = Compartment), color = "black", linewidth = linewidths) +
     geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = symbsizes) +
     scale_x_continuous(limits = c(0, max_days_all_sites)) +
     # scale_y_continuous(limits = c(0, max_value.fig3)) +
     scale_shape_manual(values = c("Susceptible" = 16,
                                   "Infected" = 17,
                                   "Recovered" = 15),
                        name = "") + # Removes "Compartment" from the legend title
      theme_classic(base_family = "Georgia") +
      theme(legend.position = "bottom",
            axis.title = element_text(size = titlesize),
            axis.text = element_text(size = textsize),
            # axis.text.x = element_blank(), #remove x-axis text but keep the line
            # axis.ticks.x = element_blank(), #remove the tick marks
            # axis.title.x = element_blank(),
            # axis.title.y = element_blank(),
            strip.text = element_text(size = 8),
            legend.text = element_text(size = textsize),
            legend.title = element_text(size = titlesize),
            legend.key.height = unit(0, "cm"))
    
    fig3_col1 = (p.fit.nearshore.single.figures / p.fit.offshore.single.figures / p.fit.near.to.off.single.figures) + 
      plot_layout(guides = "collect", axis_titles = 'collect') &
      theme(legend.position = "bottom",
            legend.box.spacing = unit(0, "cm"),
            legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
            legend.text = element_text(size = 8),  # Reduce legend text size
            # axis.title.x = element_blank(),  # Remove x-axis label
            # legend.spacing.y = unit(0, "cm"),    # Reduce space between legend items (doesn't work)
            plot.margin = margin(t = 5, r = 5, b = 0, l = 5)) &
      labs(x = "Day of outbreak")  # Apply x-axis label at the composite level

  # # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # # quartz(h = 5, w = 3.35)
  # # quartz(h = 6, w = 7.087)
  # quartz(h = 5, w = 3)
  # 
  # fig3_col1
  
  # #for Benthics
  # quartz(h = 4, w = 4)
  # fig3_col1
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig1.pdf"), type = "pdf")
  
  # #ggplot-export to image
  # ggsave(filename = here("output", "fig3_col1.png"), device = "png", width = 4, height = 4, dpi = 1200)

  # # Close the Quartz device
  # dev.off()
  
  #COLUMN 2
  # p.fit.nearshore.multi / p.fit.offshore.multi / p.fit.near.to.off.multis
  
  # NOTE - could consider a shaded ribbon to display the total accumulated recovered compartment
  #         (would accentuate how high-susceptibility corals dominate the mortality)
  
  site.loop = 'Nearshore'
  p.fit.nearshore.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility, shape = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility), size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          # axis.text.x = element_blank(), #remove x-axis text but keep the line
          # axis.ticks.x = element_blank(), #remove the tick marks
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.offshore.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility, shape = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility), size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          # axis.text.x = element_blank(), #remove x-axis text but keep the line
          # axis.ticks.x = element_blank(), #remove the tick marks
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.near.to.off.multi.figures =
    ggplot(data = data_fig3 %>%
             # filter(Site == site.loop, Host == "Multi-host", Type == "Projected"),
             # aes(x = days.model, y = tissue, color = Susceptibility, shape = Compartment)) +
           filter(Site == site.loop, Host == "Multi-host", Type == "Projected") %>%
             group_by(days.model, Compartment) %>%
             summarize(tissue = sum(tissue), .groups = "drop"),
           aes(x = days.model, y = tissue, shape = Compartment)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    # geom_point(data = obs.multi.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility), size = symbsizes) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop), aes(x = days.inf.site, y = tissue, shape = Compartment), color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          # axis.text.x = element_blank(), #remove x-axis text but keep the line
          # axis.ticks.x = element_blank(), #remove the tick marks
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  fig3_col2 = (p.fit.nearshore.multi.figures / p.fit.offshore.multi.figures / p.fit.near.to.off.multi.figures) + 
    plot_layout(axis_titles = 'collect', guides = "collect") &
    # plot_layout(guides = "collect", axes = "collect") &  # Collect the legends
    # theme(legend.position = "bottom",
    #       legend.box.spacing = unit(0, "cm"),
    #       legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
    #       legend.text = element_text(size = 8),  # Reduce legend text size
    #       # axis.title.y = element_blank(),  # Remove y-axis label
    #       # axis.title.x = element_blank(),  # Remove x-axis label
    #       # legend.spacing.y = unit(0, "cm"),    # Reduce space between legend items (doesn't work)
    #       legend.title = element_blank(),  # Add this line to remove legend titles
    #       plot.margin = margin(t = 5, r = 5, b = 0, l = 5)#, line = element_line(linewidth = 1)  # Change the universal linewidth
    # )
    theme(legend.position = "bottom",
        legend.box = "vertical",  # Stack legends vertically
        legend.box.spacing = unit(0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  )
  
  # # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # # quartz(h = 5, w = 3.35)
  # # quartz(h = 6, w = 7.087)
  # quartz(h = 4, w = 4)
  # 
  # fig3_col2
  # 
  # #ggplot-export to image
  # ggsave(filename = here("output", "fig3_col2_blacklines.png"), device = "png", width = 4, height = 4, dpi = 1200)
  
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig1.pdf"), type = "pdf")
  
  # # Close the Quartz device
  # dev.off()
  
  #COMPOSITE
  # fig3_col1 = fig3_col1 & theme(legend.position = "none")  # Remove legend from fig3_col1
  
  fig3 = (fig3_col1 | fig3_col2) +
    plot_layout(guides = "collect", axis_titles = "collect")  &
    ylab("Surface area of tissue (m²)") &
    xlab("Day of outbreak") &
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          legend.box = "vertical",  # Stack legends vertically
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  p1 = p.fit.nearshore.single.figures
  p2 = p.fit.nearshore.multi.figures
  p3 = p.fit.offshore.single.figures
  p4 = p.fit.offshore.multi.figures
  p5 = p.fit.near.to.off.single.figures
  p6 = p.fit.near.to.off.multi.figures
  
  fig3 =
  p1 + p2 + p3 + p4 + p5 + p6 +
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
  fig3 <- fig3 + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)
  
  fig3
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "fig3.pdf"), type = "pdf")

  #ggplot-export to image
  ggsave(filename = here("output", "fig3.png"), device = "png", width = 7.087, height = 5, dpi = 1200)

  # Close the Quartz device
  dev.off()
  
  
  ################################## Figure S5 ##################################
  
  #figure 3, but just infections
  
  linewidths =  0.75 #0.4 #0.75 is roughly 1 pt. ggplot measures these in mm, not points
  symbsizes = 1.3
  titlesize = 10 #9  #text sizes in ggplot are actually in units of points, when specified using element_text
  textsize = 9
  
  # COLUMN 1
  
  site.loop = 'Nearshore'
  p.fit.nearshore.inf.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Infected', Host == "Single-host", Type == "Fitted"),
           aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Infected'), aes(x = days.inf.site, y = tissue),
               shape = 17, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.offshore.inf.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Infected', Host == "Single-host", Type == "Fitted"),
           aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Infected'), aes(x = days.inf.site, y = tissue),
               shape = 17, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    # ylim(0, 0.018) +
    ylim(0, 0.105) +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.near.to.off.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Infected', Host == "Single-host", Type == "Projected"), aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Infected'), aes(x = days.inf.site, y = tissue),
               shape = 17, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    # ylim(0, 0.018) +
    ylim(0, 0.105) +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  figS5_col1 = (p.fit.nearshore.inf.single.figures / p.fit.offshore.inf.single.figures / p.fit.near.to.off.single.figures) + 
    plot_layout(guides = "collect", axis_titles = 'collect') &
    theme(legend.position = "bottom",
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
          legend.text = element_text(size = 8),  # Reduce legend text size
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)) &
    labs(x = "Day of outbreak")  # Apply x-axis label at the composite level
  
  site.loop = 'Nearshore'
  p.fit.nearshore.inf.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Infected', Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop, Compartment == 'Infected'), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility),
               shape = 17, size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.offshore.inf.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Infected', Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop, Compartment == 'Infected'), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility),
               shape = 17, size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.near.to.off.multi.figures =
    ggplot(data = data_fig3 %>%
             filter(Site == site.loop, Compartment == 'Infected', Host == "Multi-host", Type == "Projected") %>%
             group_by(days.model, Compartment) %>%
             summarize(tissue = sum(tissue), .groups = "drop"),
           aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Infected'), aes(x = days.inf.site, y = tissue),
               shape = 17, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    # ylim(0, 0.018) +
    ylim(0, 0.105) +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  figS5_col2 = (p.fit.nearshore.inf.multi.figures / p.fit.offshore.inf.multi.figures / p.fit.near.to.off.multi.figures) + 
    plot_layout(axis_titles = 'collect', guides = "collect") &
    theme(legend.position = "bottom",
          legend.box = "vertical",  # Stack legends vertically
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  p1 = p.fit.nearshore.inf.single.figures
  p2 = p.fit.nearshore.inf.multi.figures
  p3 = p.fit.offshore.inf.single.figures
  p4 = p.fit.offshore.inf.multi.figures
  p5 = p.fit.near.to.off.single.figures
  p6 = p.fit.near.to.off.multi.figures
  
  figS5 =
    p1 + p2 + p3 + p4 + p5 + p6 +
    plot_layout(axis_titles = "collect",
                guide = 'collect',
                design = "AB\nCD\nEF") &
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize, color = 'black'),
          axis.text = element_text(size = textsize, color = 'black'),
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black")
    )
  
  # add annotation separately (must be done to avoid breaking axes collect)
  figS5 <- figS5 + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)
  
  figS5
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "figS5.pdf"), type = "pdf")
  
  #ggplot-export to image
  ggsave(filename = here("output", "figS5.png"), device = "png", width = 7.087, height = 5, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
  
  
  ################################## Figure S6 ##################################
  
  #figure 3, but just removal
  
  linewidths =  0.75 #0.4 #0.75 is roughly 1 pt. ggplot measures these in mm, not points
  symbsizes = 1.3
  titlesize = 10 #9  #text sizes in ggplot are actually in units of points, when specified using element_text
  textsize = 9
  
  # COLUMN 1
  
  site.loop = 'Nearshore'
  p.fit.nearshore.rem.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Recovered', Host == "Single-host", Type == "Fitted"),
           aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Recovered'), aes(x = days.inf.site, y = tissue),
               shape = 15, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.offshore.rem.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Recovered', Host == "Single-host", Type == "Fitted"),
           aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Recovered'), aes(x = days.inf.site, y = tissue),
               shape = 15, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    # ylim(0, 1.5) +
    ylim(0, 3.5) +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.near.to.off.single.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Recovered', Host == "Single-host", Type == "Projected"), aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(color = "black", linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Recovered'), aes(x = days.inf.site, y = tissue),
               shape = 15, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    # ylim(0, 1.5) +
    ylim(0, 3.5) +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  figS6_col1 = (p.fit.nearshore.rem.single.figures / p.fit.offshore.rem.single.figures / p.fit.near.to.off.single.figures) + 
    plot_layout(guides = "collect", axis_titles = 'collect') &
    theme(legend.position = "bottom",
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),  # Shrink legend key size (symbols)
          legend.text = element_text(size = 8),  # Reduce legend text size
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)) &
    labs(x = "Day of outbreak")  # Apply x-axis label at the composite level
  
  site.loop = 'Nearshore'
  p.fit.nearshore.rem.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Recovered', Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop, Compartment == 'Recovered'), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility),
               shape = 15, size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.offshore.rem.multi.figures =
    ggplot(data = data_fig3 %>% filter(Site == site.loop, Compartment == 'Recovered', Host == "Multi-host", Type == "Fitted"), aes(x = days.model, y = tissue, color = Susceptibility)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    scale_color_manual(values = c("Low" = "#1E90FF",   # Blue
                                  "Moderate" = "#FFD700", # Yellow
                                  "High" = "#FF1493")) +  # Deep Pink
    geom_point(data = obs.multi.figures %>% filter(Site == site.loop, Compartment == 'Recovered'), aes(x = days.inf.site, y = tissue, shape = Compartment, fill = Susceptibility),
               shape = 15, size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    theme_classic(base_family = "Georgia") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  site.loop = 'Offshore'
  p.fit.near.to.off.multi.figures =
    ggplot(data = data_fig3 %>%
             filter(Site == site.loop, Compartment == 'Recovered', Host == "Multi-host", Type == "Projected") %>%
             group_by(days.model, Compartment) %>%
             summarize(tissue = sum(tissue), .groups = "drop"),
           aes(x = days.model, y = tissue)) +
    xlab("Day of outbreak") +
    ylab("Surface area of tissue (m²)") +
    geom_line(linewidth = linewidths) +
    geom_point(data = obs.total.figures %>% filter(Site == site.loop, Compartment == 'Recovered'), aes(x = days.inf.site, y = tissue),
               shape = 15, color = "black", size = symbsizes) +
    scale_x_continuous(limits = c(0, max_days_all_sites)) +
    scale_shape_manual(values = c("Susceptible" = 16,
                                  "Infected" = 17,
                                  "Recovered" = 15),
                       name = "") + # Removes "Compartment" from the legend title
    theme_classic(base_family = "Georgia") +
    # ylim(0, 1.5) +
    ylim(0, 3.5) +
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize),
          axis.text = element_text(size = textsize),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize),
          legend.key.height = unit(0, "cm"))
  
  figS6_col2 = (p.fit.nearshore.rem.multi.figures / p.fit.offshore.rem.multi.figures / p.fit.near.to.off.multi.figures) + 
    plot_layout(axis_titles = 'collect', guides = "collect") &
    theme(legend.position = "bottom",
          legend.box = "vertical",  # Stack legends vertically
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          legend.margin = margin(5, 0, 0, 0),  # Reduce margin around each legend
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  p1 = p.fit.nearshore.rem.single.figures
  p2 = p.fit.nearshore.rem.multi.figures
  p3 = p.fit.offshore.rem.single.figures
  p4 = p.fit.offshore.rem.multi.figures
  p5 = p.fit.near.to.off.single.figures
  p6 = p.fit.near.to.off.multi.figures
  
  figS6 =
    p1 + p2 + p3 + p4 + p5 + p6 +
    plot_layout(axis_titles = "collect",
                guide = 'collect',
                design = "AB\nCD\nEF") &
    theme(legend.position = "bottom",
          axis.title = element_text(size = titlesize, color = 'black'),
          axis.text = element_text(size = textsize, color = 'black'),
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black")
    )
  
  # add annotation separately (must be done to avoid breaking axes collect)
  figS6 <- figS6 + plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(1, 1))  # Move labels to top right
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # quartz(h = 5, w = 3.35)
  quartz(h = 5, w = 7.087)
  # quartz(h = 6, w = 5)
  
  figS6
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", "figS6.pdf"), type = "pdf")
  
  #ggplot-export to image
  ggsave(filename = here("output", "figS6.png"), device = "png", width = 7.087, height = 5, dpi = 1200)
  
  # Close the Quartz device
  dev.off()
  
  
  
  ################################## Table X ##################################
  
  # First reshape the data to get the format you need
  reformatted_table_format <- max_prop_prev %>%
    # Convert from long to wide format
    pivot_wider(
      id_cols = Site,
      names_from = Compartment,
      values_from = c(max_tissue_normalized, max_prevalence)
    ) %>%
    # Select and rename columns in the correct order
    transmute(
      Site = Site,
      "Tissue prev." = round(max_tissue_normalized_Infected, 2),
      "Colony prev." = round(max_prevalence_Infected, 2),
      "% Tissue lost" = round(max_tissue_normalized_Recovered, 2),
      "% Colonies lost" = round(max_prevalence_Recovered, 2)
    )
  
  mean_row <- tibble(
    Site = "Mean",
    "Tissue prev." = round(mean(reformatted_table_format$`Tissue prev.`, na.rm = TRUE), 2),
    "Colony prev." = round(mean(reformatted_table_format$`Colony prev.`, na.rm = TRUE), 2),
    "% Tissue lost" = round(mean(reformatted_table_format$`% Tissue lost`, na.rm = TRUE), 2),
    "% Colonies lost" = round(mean(reformatted_table_format$`% Colonies lost`, na.rm = TRUE), 2)
  )
  
  # Combine data with mean row
  final_table2 <- bind_rows(reformatted_table_format, mean_row)
  
  # #optional html table for display in R
  # formatted_table <- final_table2 %>%
  #   kable("html") %>%
  #   kable_styling(bootstrap_options = c("striped", "hover"),
  #                 full_width = FALSE)
  
  ################################## Table X2 - Freq. ##################################
  
  # NOTE - 21 April 2025
  #   - ran the below using cover model workspace (i.e., single-host model w/ nonlinear density dependence)
  #   - the below also specifically notates objects using 'freq' to indicates frequency-dependence - but can be
  #       used to summarize output metrics from density-dependent runs too. ideally could be re-coded accordingly

  # error_eval stuff here
  #
  table_x2 = error_metrics %>%
    select(-NRMSE_range, -NRMSE_mean, -sMAPE) %>%
    filter(!(type == "DHW" | wave == "Pre-heat")) %>%
    select(-wave)



  # param stuff here
  #
  #frequency-dependent
  #basic
  beta.offshore.single.freq = params.basic.offshore.full[1]
  # min.beta.tiss.adj = min.beta.tiss * (1 / (1 + exp(-lambda.modifier * (cover.site))) + offset)
  # beta.offshore.adj.full = params.basic.offshore.full[2]
  gamma.offshore.single.freq = params.basic.offshore.full[3]
  R0.offshore.single.freq = beta.offshore.single.freq / gamma.offshore.single.freq
  #
  beta.midchannel.single.freq = params.basic.midchannel.full[1]
  gamma.midchannel.single.freq = params.basic.midchannel.full[3]
  R0.midchannel.single.freq = beta.midchannel.single.freq / gamma.midchannel.single.freq
  #
  beta.nearshore.single.freq = params.basic.nearshore.full[1]
  gamma.nearshore.single.freq = params.basic.nearshore.full[3]
  R0.nearshore.single.freq = beta.nearshore.single.freq / gamma.nearshore.single.freq
  #
  #
  #
  #mixed dependency
  #multi
  beta.LS.offshore.multi = beta.offshore.LS
  beta.MS.offshore.multi = beta.offshore.MS
  beta.HS.offshore.multi = beta.offshore.HS
  # betas.offshore.adj = c(beta.offshore.LS.adj, beta.offshore.MS.adj, beta.offshore.HS.adj)
  gamma.LS.offshore.multi = gamma.offshore.LS
  gamma.MS.offshore.multi = gamma.offshore.MS
  gamma.HS.offshore.multi = gamma.offshore.HS
  R0s.offshore.multi = c(beta.LS.offshore.multi, beta.MS.offshore.multi, beta.HS.offshore.multi) / c(gamma.LS.offshore.multi, gamma.MS.offshore.multi, gamma.HS.offshore.multi)
  #
  beta.LS.midchannel.multi = beta.midchannel.LS
  beta.MS.midchannel.multi = beta.midchannel.MS
  beta.HS.midchannel.multi = beta.midchannel.HS
  gamma.LS.midchannel.multi = gamma.midchannel.LS
  gamma.MS.midchannel.multi = gamma.midchannel.MS
  gamma.HS.midchannel.multi = gamma.midchannel.HS
  R0s.midchannel.multi = c(beta.LS.midchannel.multi, beta.MS.midchannel.multi, beta.HS.midchannel.multi) / c(gamma.LS.midchannel.multi, gamma.MS.midchannel.multi, gamma.HS.midchannel.multi)
  #
  beta.LS.nearshore.multi = beta.nearshore.LS
  beta.MS.nearshore.multi = beta.nearshore.MS
  beta.HS.nearshore.multi = beta.nearshore.HS
  gamma.LS.nearshore.multi = gamma.nearshore.LS
  gamma.MS.nearshore.multi = gamma.nearshore.MS
  gamma.HS.nearshore.multi = gamma.nearshore.HS
  R0s.nearshore.multi = c(beta.LS.nearshore.multi, beta.MS.nearshore.multi, beta.HS.nearshore.multi) / c(gamma.LS.nearshore.multi, gamma.MS.nearshore.multi, gamma.HS.nearshore.multi)

  # Create a new dataframe with the updated structure - now with separate Beta and Gamma columns
  reformatted_table <- data.frame(
    Simulation = character(),
    Host = character(),
    Dependence = character(),
    Betas = character(),      # New column for Beta values
    Gammas = character(),     # New column for Gamma values
    Effective_R0 = character(),
    R2 = numeric(),
    RMSE = numeric(),
    stringsAsFactors = FALSE
  )

  # Map site/type/host combinations to the desired output format
  for (i in 1:nrow(table_x2)) {
    row <- table_x2[i,]

    # Determine simulation name based on site and type
    simulation <- ""
    if (row$site == "near") {
      if (row$type == "Fitted") {
        simulation <- "Nearshore"
      } else {
        simulation <- "Off -> Near"
      }
    } else if (row$site == "mid") {
      simulation <- "Midchannel"
    } else if (row$site == "off") {
      if (row$type == "Fitted") {
        simulation <- "Offshore"
      } else {
        simulation <- "Near -> Off"
      }
    }

    # Determine dependence type
    dependence <- ifelse(row$host == "Single", "Frequency", "Mixed")

    # Initialize empty beta, gamma and R0 values (will be filled later)
    betas <- ""
    gammas <- ""
    effective_r0 <- ""

    # Add to the new dataframe
    reformatted_table <- rbind(reformatted_table, data.frame(
      Simulation = simulation,
      Host = row$host,
      Dependence = dependence,
      Betas = betas,
      Gammas = gammas,
      Effective_R0 = effective_r0,
      R2 = round(row$R_squared, 2),
      RMSE = round(row$RMSE, 2)
    ))
  }

  # Define the order of simulations and hosts for proper sorting
  order_simulations <- c("Offshore", "Midchannel", "Nearshore", "Near -> Off", "Off -> Near")
  order_hosts <- c("Single", "Multi")

  # Create sorting keys
  simulation_order <- match(reformatted_table$Simulation, order_simulations)
  host_order <- match(reformatted_table$Host, order_hosts)

  # Sort the table
  reformatted_table <- reformatted_table[order(host_order, simulation_order), ]

  # Now fill in the parameter values and R0s from the provided variable names

  # Offshore Single Frequency
  idx <- which(reformatted_table$Simulation == "Offshore" & reformatted_table$Host == "Single")
  reformatted_table$Betas[idx] <- as.character(round(beta.offshore.single.freq, 2))
  reformatted_table$Gammas[idx] <- as.character(round(gamma.offshore.single.freq, 2))
  reformatted_table$Effective_R0[idx] <- as.character(round(R0.offshore.single.freq, 2))

  # Midchannel Single Frequency
  idx <- which(reformatted_table$Simulation == "Midchannel" & reformatted_table$Host == "Single")
  reformatted_table$Betas[idx] <- as.character(round(beta.midchannel.single.freq, 2))
  reformatted_table$Gammas[idx] <- as.character(round(gamma.midchannel.single.freq, 2))
  reformatted_table$Effective_R0[idx] <- as.character(round(R0.midchannel.single.freq, 2))

  # Nearshore Single Frequency
  idx <- which(reformatted_table$Simulation == "Nearshore" & reformatted_table$Host == "Single")
  reformatted_table$Betas[idx] <- as.character(round(beta.nearshore.single.freq, 2))
  reformatted_table$Gammas[idx] <- as.character(round(gamma.nearshore.single.freq, 2))
  reformatted_table$Effective_R0[idx] <- as.character(round(R0.nearshore.single.freq, 2))

  # Near -> Off Single (uses Nearshore parameters)
  idx <- which(reformatted_table$Simulation == "Near -> Off" & reformatted_table$Host == "Single")
  reformatted_table$Betas[idx] <- as.character(round(beta.nearshore.single.freq, 2))
  reformatted_table$Gammas[idx] <- as.character(round(gamma.nearshore.single.freq, 2))
  reformatted_table$Effective_R0[idx] <- as.character(round(R0.nearshore.single.freq, 2))

  # Off -> Near Single (uses Offshore parameters)
  idx <- which(reformatted_table$Simulation == "Off -> Near" & reformatted_table$Host == "Single")
  reformatted_table$Betas[idx] <- as.character(round(beta.offshore.single.freq, 2))
  reformatted_table$Gammas[idx] <- as.character(round(gamma.offshore.single.freq, 2))
  reformatted_table$Effective_R0[idx] <- as.character(round(R0.offshore.single.freq, 2))

  # Offshore Multi Mixed
  idx <- which(reformatted_table$Simulation == "Offshore" & reformatted_table$Host == "Multi")
  reformatted_table$Betas[idx] <- paste(round(beta.LS.offshore.multi, 2), round(beta.MS.offshore.multi, 2), round(beta.HS.offshore.multi, 2), sep = ", ")
  reformatted_table$Gammas[idx] <- paste(round(gamma.LS.offshore.multi, 2), round(gamma.MS.offshore.multi, 2), round(gamma.HS.offshore.multi, 2), sep = ", ")
  reformatted_table$Effective_R0[idx] <- paste(round(R0s.offshore.multi, 2), collapse = ", ")

  # Midchannel Multi Mixed
  idx <- which(reformatted_table$Simulation == "Midchannel" & reformatted_table$Host == "Multi")
  reformatted_table$Betas[idx] <- paste(round(beta.LS.midchannel.multi, 2), round(beta.MS.midchannel.multi, 2), round(beta.HS.midchannel.multi, 2), sep = ", ")
  reformatted_table$Gammas[idx] <- paste(round(gamma.LS.midchannel.multi, 2), round(gamma.MS.midchannel.multi, 2), round(gamma.HS.midchannel.multi, 2), sep = ", ")
  reformatted_table$Effective_R0[idx] <- paste(round(R0s.midchannel.multi, 2), collapse = ", ")

  # Nearshore Multi Mixed
  idx <- which(reformatted_table$Simulation == "Nearshore" & reformatted_table$Host == "Multi")
  reformatted_table$Betas[idx] <- paste(round(beta.LS.nearshore.multi, 2), round(beta.MS.nearshore.multi, 2), round(beta.HS.nearshore.multi, 2), sep = ", ")
  reformatted_table$Gammas[idx] <- paste(round(gamma.LS.nearshore.multi, 2), round(gamma.MS.nearshore.multi, 2), round(gamma.HS.nearshore.multi, 2), sep = ", ")
  reformatted_table$Effective_R0[idx] <- paste(round(R0s.nearshore.multi, 2), collapse = ", ")

  # Near -> Off Multi (uses Nearshore parameters)
  idx <- which(reformatted_table$Simulation == "Near -> Off" & reformatted_table$Host == "Multi")
  reformatted_table$Betas[idx] <- paste(round(beta.LS.nearshore.multi, 2), round(beta.MS.nearshore.multi, 2), round(beta.HS.nearshore.multi, 2), sep = ", ")
  reformatted_table$Gammas[idx] <- paste(round(gamma.LS.nearshore.multi, 2), round(gamma.MS.nearshore.multi, 2), round(gamma.HS.nearshore.multi, 2), sep = ", ")
  reformatted_table$Effective_R0[idx] <- paste(round(R0s.nearshore.multi, 2), collapse = ", ")

  # Off -> Near Multi (uses Offshore parameters)
  idx <- which(reformatted_table$Simulation == "Off -> Near" & reformatted_table$Host == "Multi")
  reformatted_table$Betas[idx] <- paste(round(beta.LS.offshore.multi, 2), round(beta.MS.offshore.multi, 2), round(beta.HS.offshore.multi, 2), sep = ", ")
  reformatted_table$Gammas[idx] <- paste(round(gamma.LS.offshore.multi, 2), round(gamma.MS.offshore.multi, 2), round(gamma.HS.offshore.multi, 2), sep = ", ")
  reformatted_table$Effective_R0[idx] <- paste(round(R0s.offshore.multi, 2), collapse = ", ")

  # Print the final table
  print(reformatted_table)

  #write to CSV for easy dumping of data into Word table in manuscript
  output_path <- here("output", "tablex2_freqency.csv")

  # Write the table to a CSV file without row names
  write.csv(reformatted_table, output_path, row.names = FALSE)

  # Print the file path for confirmation
  cat("Table saved to: ", output_path, "\n")
  # Print the final table
  print(reformatted_table)

  
  ################################## Force of infection / R0 ##################################
  
  # TEST - stopping point - note - 22 April 2025
  # - working out Jacobian matrix for group-wide R0
  # Define parameters
  
  # OFFSHORE
  betas.offshore.jacobian.freq <- c(beta.LS.offshore.multi, beta.MS.offshore.multi, beta.HS.offshore.multi)
  gammas.offshore.jacobian.freq <- c(gamma.LS.offshore.multi, gamma.MS.offshore.multi, gamma.HS.offshore.multi)
  Ns.offshore <- c(N.LS.offshore, N.MS.offshore, N.HS.offshore)
  
  # Define group names for clear labeling
  group_names <- c("LS", "MS", "HS")
  
  # Construct F matrix
  F.offshore.freq <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      F.offshore.freq[i, j] <- betas.offshore.jacobian.freq[j] * Ns.offshore[i] / Ns.offshore[j]
    }
  }
  
  # Add row and column names to F matrix
  rownames(F.offshore.freq) <- paste0("To_", group_names)
  colnames(F.offshore.freq) <- paste0("From_", group_names)
  
  # Construct V matrix (diagonal with gammas)
  V.offshore.freq <- diag(gammas.offshore.jacobian.freq)
  
  # Add row and column names to V matrix
  rownames(V.offshore.freq) <- group_names
  colnames(V.offshore.freq) <- group_names
  
  # Next-generation matrix
  K.offshore.freq <- F.offshore.freq %*% solve(V.offshore.freq)
  
  # Add row and column names to K matrix
  rownames(K.offshore.freq) <- paste0("To_", group_names)
  colnames(K.offshore.freq) <- paste0("From_", group_names)
  
  # R0 is the dominant eigenvalue
  eigenvalues.offshore.freq <- eigen(K.offshore.freq)$values
  R0.offshore.jacobian.freq <- max(Re(eigenvalues.offshore.freq))
  
  
  
  
  # MIDCHANNEL
  betas.midchannel.jacobian.freq <- c(beta.LS.midchannel.multi, beta.MS.midchannel.multi, beta.HS.midchannel.multi)
  gammas.midchannel.jacobian.freq <- c(gamma.LS.midchannel.multi, gamma.MS.midchannel.multi, gamma.HS.midchannel.multi)
  Ns.midchannel <- c(N.LS.midchannel, N.MS.midchannel, N.HS.midchannel)
  
  # Construct F matrix
  F.midchannel.freq <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      F.midchannel.freq[i, j] <- betas.midchannel.jacobian.freq[j] * Ns.midchannel[i] / Ns.midchannel[j]
    }
  }
  
  # Add row and column names to F matrix
  rownames(F.midchannel.freq) <- paste0("To_", group_names)
  colnames(F.midchannel.freq) <- paste0("From_", group_names)
  
  # Construct V matrix (diagonal with gammas)
  V.midchannel.freq <- diag(gammas.midchannel.jacobian.freq)
  
  # Add row and column names to V matrix
  rownames(V.midchannel.freq) <- group_names
  colnames(V.midchannel.freq) <- group_names
  
  # Next-generation matrix
  K.midchannel.freq <- F.midchannel.freq %*% solve(V.midchannel.freq)
  
  # Add row and column names to K matrix
  rownames(K.midchannel.freq) <- paste0("To_", group_names)
  colnames(K.midchannel.freq) <- paste0("From_", group_names)
  
  # R0 is the dominant eigenvalue
  eigenvalues.midchannel.freq <- eigen(K.midchannel.freq)$values
  R0.midchannel.jacobian.freq <- max(Re(eigenvalues.midchannel.freq))
  
  
  
  
  # NEARSHORE
  betas.nearshore.jacobian.freq <- c(beta.LS.nearshore.multi, beta.MS.nearshore.multi, beta.HS.nearshore.multi)
  gammas.nearshore.jacobian.freq <- c(gamma.LS.nearshore.multi, gamma.MS.nearshore.multi, gamma.HS.nearshore.multi)
  Ns.nearshore <- c(N.LS.nearshore, N.MS.nearshore, N.HS.nearshore)
  
  # Construct F matrix
  F.nearshore.freq <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      F.nearshore.freq[i, j] <- betas.nearshore.jacobian.freq[j] * Ns.nearshore[i] / Ns.nearshore[j]
    }
  }
  
  # Add row and column names to F matrix
  rownames(F.nearshore.freq) <- paste0("To_", group_names)
  colnames(F.nearshore.freq) <- paste0("From_", group_names)
  
  # Construct V matrix (diagonal with gammas)
  V.nearshore.freq <- diag(gammas.nearshore.jacobian.freq)
  
  # Add row and column names to V matrix
  rownames(V.nearshore.freq) <- group_names
  colnames(V.nearshore.freq) <- group_names
  
  # Next-generation matrix
  K.nearshore.freq <- F.nearshore.freq %*% solve(V.nearshore.freq)
  
  # Add row and column names to K matrix
  rownames(K.nearshore.freq) <- paste0("To_", group_names)
  colnames(K.nearshore.freq) <- paste0("From_", group_names)
  
  # R0 is the dominant eigenvalue
  eigenvalues.nearshore.freq <- eigen(K.nearshore.freq)$values
  R0.nearshore.jacobian.freq <- max(Re(eigenvalues.nearshore.freq))
  
  
  
  # NEAR TO OFF
  betas.nearshore.jacobian.freq
  gammas.nearshore.jacobian.freq
  Ns.offshore
  
  # Construct F matrix
  F.near.to.off.freq <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      F.near.to.off.freq[i, j] <- betas.nearshore.jacobian.freq[j] * Ns.offshore[i] / Ns.offshore[j]
    }
  }
  
  # Add row and column names to F matrix
  rownames(F.near.to.off.freq) <- paste0("To_", group_names)
  colnames(F.near.to.off.freq) <- paste0("From_", group_names)
  
  # Construct V matrix (diagonal with gammas)
  V.near.to.off.freq <- diag(gammas.nearshore.jacobian.freq)
  
  # Add row and column names to V matrix
  rownames(V.near.to.off.freq) <- group_names
  colnames(V.near.to.off.freq) <- group_names
  
  # Next-generation matrix
  K.near.to.off.freq <- F.near.to.off.freq %*% solve(V.near.to.off.freq)
  
  # Add row and column names to K matrix
  rownames(K.near.to.off.freq) <- paste0("To_", group_names)
  colnames(K.near.to.off.freq) <- paste0("From_", group_names)
  
  # R0 is the dominant eigenvalue
  eigenvalues.near.to.off.freq <- eigen(K.near.to.off.freq)$values
  R0.near.to.off.jacobian.freq <- max(Re(eigenvalues.near.to.off.freq))
  
  
  
  # OFF TO NEAR
  betas.offshore.jacobian.freq
  gammas.offshore.jacobian.freq
  Ns.nearshore
  
  # Construct F matrix
  F.off.to.near.freq <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      F.off.to.near.freq[i, j] <- betas.offshore.jacobian.freq[j] * Ns.nearshore[i] / Ns.nearshore[j]
    }
  }
  
  # Add row and column names to F matrix
  rownames(F.off.to.near.freq) <- paste0("To_", group_names)
  colnames(F.off.to.near.freq) <- paste0("From_", group_names)
  
  # Construct V matrix (diagonal with gammas)
  V.off.to.near.freq <- diag(gammas.offshore.jacobian.freq)
  
  # Add row and column names to V matrix
  rownames(V.off.to.near.freq) <- group_names
  colnames(V.off.to.near.freq) <- group_names
  
  # Next-generation matrix
  K.off.to.near.freq <- F.off.to.near.freq %*% solve(V.off.to.near.freq)
  
  # Add row and column names to K matrix
  rownames(K.off.to.near.freq) <- paste0("To_", group_names)
  colnames(K.off.to.near.freq) <- paste0("From_", group_names)
  
  # R0 is the dominant eigenvalue
  eigenvalues.off.to.near.freq <- eigen(K.off.to.near.freq)$values
  R0.off.to.near.jacobian.freq <- max(Re(eigenvalues.off.to.near.freq))
  
  
  
  
  # RESULTS
  # # Print F matrix with interpretation
  # # cat("Offshore F Matrix (new infections):\n")
  # print(F.offshore.freq)
  # # cat("\nInterpretation: F[i,j] represents infection rate from group j to group i\n")
  # # cat("For example, F['To_LS', 'From_MS'] is the rate at which MS group infects LS group\n\n")
  # 
  # # Print V matrix with interpretation
  # # cat("Offshore V Matrix (transition rates):\n")
  # print(V.offshore.freq)
  # # cat("\nInterpretation: V contains the recovery rates for each group on the diagonal\n\n")
  # 
  # # Print K matrix with interpretation
  # # cat("Offshore K Matrix (next-generation matrix):\n")
  # print(K.offshore.freq)
  # # cat("\nInterpretation: K[i,j] is the expected number of new infections in group i\n")
  # # cat("produced by a single infected individual in group j during its infectious period\n\n")
  
  # cat("Offshore R0 value:", R0.offshore.jacobian.freq, "\n")
  # 
  # print(F.midchannel.freq)
  # print(V.midchannel.freq)
  # print(K.midchannel.freq)
  # cat("Midchannel R0 value:", R0.midchannel.jacobian.freq, "\n")
  # 
  # print(F.nearshore.freq)
  # print(V.nearshore.freq)
  # print(K.nearshore.freq)
  # cat("Nearshore R0 value:", R0.nearshore.jacobian.freq, "\n")
  
  #R0s
  cat("Offshore R0 value:", R0.offshore.jacobian.freq, "\n")
  cat("Midchannel R0 value:", R0.midchannel.jacobian.freq, "\n")
  cat("Nearshore R0 value:", R0.nearshore.jacobian.freq, "\n")
  cat("Near -> Off R0 value:", R0.near.to.off.jacobian.freq, "\n")
  cat("Off -> Near R0 value:", R0.off.to.near.jacobian.freq, "\n")
  
  #Fs
  print(F.offshore.freq)
  print(F.midchannel.freq)
  print(F.nearshore.freq)
  print(F.near.to.off.freq)
  print(F.off.to.near.freq)
  
  #Ks
  print(K.offshore.freq)
  print(K.midchannel.freq)
  print(K.nearshore.freq)
  print(K.near.to.off.freq)
  print(K.off.to.near.freq)
  
  #forces of infection exerted on species j by the rest of the community (Dobson 2004)
  print(rowSums(K.offshore.freq))
  print(rowSums(K.midchannel.freq))
  print(rowSums(K.nearshore.freq))
  print(rowSums(K.near.to.off.freq))
  print(rowSums(K.off.to.near.freq))
  
  #forces of infection placed on other species due to presence of infections in species j (Dobson 2004):
  print(colSums(K.offshore.freq))
  print(colSums(K.midchannel.freq))
  print(colSums(K.nearshore.freq))
  print(colSums(K.near.to.off.freq))
  print(colSums(K.off.to.near.freq))
  
  #re-write table with community-wide R0 instead of within-group
  idx <- which(reformatted_table$Simulation == "Offshore" & reformatted_table$Host == "Multi")
  reformatted_table$Effective_R0[idx] <- paste(round(R0.offshore.jacobian.freq, 2), collapse = ", ")
  
  idx <- which(reformatted_table$Simulation == "Midchannel" & reformatted_table$Host == "Multi")
  reformatted_table$Effective_R0[idx] <- paste(round(R0.midchannel.jacobian.freq, 2), collapse = ", ")
  
  idx <- which(reformatted_table$Simulation == "Nearshore" & reformatted_table$Host == "Multi")
  reformatted_table$Effective_R0[idx] <- paste(round(R0.nearshore.jacobian.freq, 2), collapse = ", ")
  
  idx <- which(reformatted_table$Simulation == "Near -> Off" & reformatted_table$Host == "Multi")
  reformatted_table$Effective_R0[idx] <- paste(round(R0.near.to.off.jacobian.freq, 2), collapse = ", ")
  
  idx <- which(reformatted_table$Simulation == "Off -> Near" & reformatted_table$Host == "Multi")
  reformatted_table$Effective_R0[idx] <- paste(round(R0.off.to.near.jacobian.freq, 2), collapse = ", ")
  
  output_path <- here("output", "tablex2_freqency.csv")
  
  # Write the table to a CSV file without row names
  write.csv(reformatted_table, output_path, row.names = FALSE)
  
  # Print the file path for confirmation
  cat("Table saved to: ", output_path, "\n")
  # Print the final table
  print(reformatted_table)
  
  
  # ################################## Table X2 - Dens. ##################################
  # 
  # # error_eval stuff here
  # #
  # table_x2 = error_metrics %>%
  #   select(-NRMSE_range, -NRMSE_mean, -sMAPE) %>%
  #   filter(!(type == "DHW" | wave == "Pre-heat")) %>%
  #   select(-wave)
  # 
  # # param stuff here
  # #
  # #density-dependent
  # #basic
  # beta.offshore.single.dens = params.basic.offshore.full[1] * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.offshore.full)) / (1 - exp(-k_val))))
  # gamma.offshore.single.dens = params.basic.offshore.full[3]
  # R0.offshore.single.dens = beta.offshore.single.dens / gamma.offshore.single.dens
  # #
  # beta.midchannel.single.dens = params.basic.midchannel.full[1] * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.midchannel.full)) / (1 - exp(-k_val))))
  # gamma.midchannel.single.dens = params.basic.midchannel.full[3]
  # R0.midchannel.single.dens = beta.midchannel.single.dens / gamma.midchannel.single.dens
  # #
  # beta.nearshore.single.dens = params.basic.nearshore.full[1] * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.nearshore.full)) / (1 - exp(-k_val))))
  # gamma.nearshore.single.dens = params.basic.nearshore.full[3]
  # R0.nearshore.single.dens = beta.nearshore.single.dens / gamma.nearshore.single.dens
  # #
  # beta.near.to.off.single.dens = params.basic.nearshore.full[1] * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.offshore.full)) / (1 - exp(-k_val))))
  # R0.near.to.off.single.dens = beta.near.to.off.single.dens / gamma.nearshore.single.dens
  # #
  # beta.off.to.near.single.dens = params.basic.offshore.full[1] * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.nearshore.full)) / (1 - exp(-k_val))))
  # R0.off.to.near.single.dens = beta.off.to.near.single.dens / gamma.offshore.single.dens
  # #
  # #
  # #
  # #mixed dependency
  # #multi
  # beta.LS.offshore.multi = beta.offshore.LS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.offshore.LS)) / (1 - exp(-k_val))))
  # beta.MS.offshore.multi = beta.offshore.MS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.offshore.MS)) / (1 - exp(-k_val))))
  # beta.HS.offshore.multi = beta.offshore.HS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.offshore.HS)) / (1 - exp(-k_val))))
  # gamma.LS.offshore.multi = gamma.offshore.LS
  # gamma.MS.offshore.multi = gamma.offshore.MS
  # gamma.HS.offshore.multi = gamma.offshore.HS
  # R0s.offshore.multi = c(beta.LS.offshore.multi, beta.MS.offshore.multi, beta.HS.offshore.multi) / c(gamma.LS.offshore.multi, gamma.MS.offshore.multi, gamma.HS.offshore.multi)
  # #
  # beta.LS.midchannel.multi = beta.midchannel.LS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.midchannel.LS)) / (1 - exp(-k_val))))
  # beta.MS.midchannel.multi = beta.midchannel.MS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.midchannel.MS)) / (1 - exp(-k_val))))
  # beta.HS.midchannel.multi = beta.midchannel.HS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.midchannel.HS)) / (1 - exp(-k_val))))
  # gamma.LS.midchannel.multi = gamma.midchannel.LS
  # gamma.MS.midchannel.multi = gamma.midchannel.MS
  # gamma.HS.midchannel.multi = gamma.midchannel.HS
  # R0s.midchannel.multi = c(beta.LS.midchannel.multi, beta.MS.midchannel.multi, beta.HS.midchannel.multi) / c(gamma.LS.midchannel.multi, gamma.MS.midchannel.multi, gamma.HS.midchannel.multi)
  # #
  # beta.LS.nearshore.multi = beta.nearshore.LS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.nearshore.LS)) / (1 - exp(-k_val))))
  # beta.MS.nearshore.multi = beta.nearshore.MS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.nearshore.MS)) / (1 - exp(-k_val))))
  # beta.HS.nearshore.multi = beta.nearshore.HS * ((1 - alpha_val) + alpha_val*((1 - exp(-k_val*cover.nearshore.HS)) / (1 - exp(-k_val))))
  # gamma.LS.nearshore.multi = gamma.nearshore.LS
  # gamma.MS.nearshore.multi = gamma.nearshore.MS
  # gamma.HS.nearshore.multi = gamma.nearshore.HS
  # R0s.nearshore.multi = c(beta.LS.nearshore.multi, beta.MS.nearshore.multi, beta.HS.nearshore.multi) / c(gamma.LS.nearshore.multi, gamma.MS.nearshore.multi, gamma.HS.nearshore.multi)
  # 
  # # Create a new dataframe with the updated structure - now with separate Beta and Gamma columns
  # reformatted_table <- data.frame(
  #   Simulation = character(),
  #   Host = character(),
  #   Dependence = character(),
  #   Betas = character(),      # New column for Beta values
  #   Gammas = character(),     # New column for Gamma values
  #   Effective_R0 = character(),
  #   R2 = numeric(),
  #   RMSE = numeric(),
  #   stringsAsFactors = FALSE
  # )
  # 
  # # Map site/type/host combinations to the desired output format
  # for (i in 1:nrow(table_x2)) {
  #   row <- table_x2[i,]
  #   
  #   # Determine simulation name based on site and type
  #   simulation <- ""
  #   if (row$site == "near") {
  #     if (row$type == "Fitted") {
  #       simulation <- "Nearshore"
  #     } else {
  #       simulation <- "Off -> Near"
  #     }
  #   } else if (row$site == "mid") {
  #     simulation <- "Midchannel"
  #   } else if (row$site == "off") {
  #     if (row$type == "Fitted") {
  #       simulation <- "Offshore"
  #     } else {
  #       simulation <- "Near -> Off"
  #     }
  #   }
  #   
  #   # Determine dependence type
  #   dependence <- ifelse(row$host == "Single", "Density", "Mixed")
  #   
  #   # Initialize empty beta, gamma and R0 values (will be filled later)
  #   betas <- ""
  #   gammas <- ""
  #   effective_r0 <- ""
  #   
  #   # Add to the new dataframe
  #   reformatted_table <- rbind(reformatted_table, data.frame(
  #     Simulation = simulation,
  #     Host = row$host,
  #     Dependence = dependence,
  #     Betas = betas,
  #     Gammas = gammas,
  #     Effective_R0 = effective_r0,
  #     R2 = round(row$R_squared, 2),
  #     RMSE = round(row$RMSE, 2)
  #   ))
  # }
  # 
  # # Define the order of simulations and hosts for proper sorting
  # order_simulations <- c("Offshore", "Midchannel", "Nearshore", "Near -> Off", "Off -> Near")
  # order_hosts <- c("Single", "Multi")
  # 
  # # Create sorting keys
  # simulation_order <- match(reformatted_table$Simulation, order_simulations)
  # host_order <- match(reformatted_table$Host, order_hosts)
  # 
  # # Sort the table
  # reformatted_table <- reformatted_table[order(host_order, simulation_order), ]
  # 
  # # Now fill in the parameter values and R0s from the provided variable names
  # 
  # # Offshore Single Frequency
  # idx <- which(reformatted_table$Simulation == "Offshore" & reformatted_table$Host == "Single")
  # reformatted_table$Betas[idx] <- as.character(round(beta.offshore.single.dens, 2))
  # reformatted_table$Gammas[idx] <- as.character(round(gamma.offshore.single.dens, 2))
  # reformatted_table$Effective_R0[idx] <- as.character(round(R0.offshore.single.dens, 2))
  # 
  # # Midchannel Single Frequency
  # idx <- which(reformatted_table$Simulation == "Midchannel" & reformatted_table$Host == "Single")
  # reformatted_table$Betas[idx] <- as.character(round(beta.midchannel.single.dens, 2))
  # reformatted_table$Gammas[idx] <- as.character(round(gamma.midchannel.single.dens, 2))
  # reformatted_table$Effective_R0[idx] <- as.character(round(R0.midchannel.single.dens, 2))
  # 
  # # Nearshore Single Frequency
  # idx <- which(reformatted_table$Simulation == "Nearshore" & reformatted_table$Host == "Single")
  # reformatted_table$Betas[idx] <- as.character(round(beta.nearshore.single.dens, 2))
  # reformatted_table$Gammas[idx] <- as.character(round(gamma.nearshore.single.dens, 2))
  # reformatted_table$Effective_R0[idx] <- as.character(round(R0.nearshore.single.dens, 2))
  # 
  # # Near -> Off Single (uses Nearshore parameters)
  # idx <- which(reformatted_table$Simulation == "Near -> Off" & reformatted_table$Host == "Single")
  # reformatted_table$Betas[idx] <- as.character(round(beta.near.to.off.single.dens, 2))
  # reformatted_table$Gammas[idx] <- as.character(round(gamma.nearshore.single.dens, 2))
  # reformatted_table$Effective_R0[idx] <- as.character(round(R0.near.to.off.single.dens, 2))
  # 
  # # Off -> Near Single (uses Offshore parameters)
  # idx <- which(reformatted_table$Simulation == "Off -> Near" & reformatted_table$Host == "Single")
  # reformatted_table$Betas[idx] <- as.character(round(beta.off.to.near.single.dens, 2))
  # reformatted_table$Gammas[idx] <- as.character(round(gamma.offshore.single.dens, 2))
  # reformatted_table$Effective_R0[idx] <- as.character(round(R0.off.to.near.single.dens, 2))
  # 
  # # Offshore Multi Mixed
  # idx <- which(reformatted_table$Simulation == "Offshore" & reformatted_table$Host == "Multi")
  # reformatted_table$Betas[idx] <- paste(round(beta.LS.offshore.multi, 2), round(beta.MS.offshore.multi, 2), round(beta.HS.offshore.multi, 2), sep = ", ")
  # reformatted_table$Gammas[idx] <- paste(round(gamma.LS.offshore.multi, 2), round(gamma.MS.offshore.multi, 2), round(gamma.HS.offshore.multi, 2), sep = ", ")
  # reformatted_table$Effective_R0[idx] <- paste(round(R0s.offshore.multi, 2), collapse = ", ")
  # 
  # # Midchannel Multi Mixed
  # idx <- which(reformatted_table$Simulation == "Midchannel" & reformatted_table$Host == "Multi")
  # reformatted_table$Betas[idx] <- paste(round(beta.LS.midchannel.multi, 2), round(beta.MS.midchannel.multi, 2), round(beta.HS.midchannel.multi, 2), sep = ", ")
  # reformatted_table$Gammas[idx] <- paste(round(gamma.LS.midchannel.multi, 2), round(gamma.MS.midchannel.multi, 2), round(gamma.HS.midchannel.multi, 2), sep = ", ")
  # reformatted_table$Effective_R0[idx] <- paste(round(R0s.midchannel.multi, 2), collapse = ", ")
  # 
  # # Nearshore Multi Mixed
  # idx <- which(reformatted_table$Simulation == "Nearshore" & reformatted_table$Host == "Multi")
  # reformatted_table$Betas[idx] <- paste(round(beta.LS.nearshore.multi, 2), round(beta.MS.nearshore.multi, 2), round(beta.HS.nearshore.multi, 2), sep = ", ")
  # reformatted_table$Gammas[idx] <- paste(round(gamma.LS.nearshore.multi, 2), round(gamma.MS.nearshore.multi, 2), round(gamma.HS.nearshore.multi, 2), sep = ", ")
  # reformatted_table$Effective_R0[idx] <- paste(round(R0s.nearshore.multi, 2), collapse = ", ")
  # 
  # # Near -> Off Multi (uses Nearshore parameters)
  # idx <- which(reformatted_table$Simulation == "Near -> Off" & reformatted_table$Host == "Multi")
  # reformatted_table$Betas[idx] <- paste(round(beta.LS.nearshore.multi, 2), round(beta.MS.nearshore.multi, 2), round(beta.HS.nearshore.multi, 2), sep = ", ")
  # reformatted_table$Gammas[idx] <- paste(round(gamma.LS.nearshore.multi, 2), round(gamma.MS.nearshore.multi, 2), round(gamma.HS.nearshore.multi, 2), sep = ", ")
  # reformatted_table$Effective_R0[idx] <- paste(round(R0s.nearshore.multi, 2), collapse = ", ")
  # 
  # # Off -> Near Multi (uses Offshore parameters)
  # idx <- which(reformatted_table$Simulation == "Off -> Near" & reformatted_table$Host == "Multi")
  # reformatted_table$Betas[idx] <- paste(round(beta.LS.offshore.multi, 2), round(beta.MS.offshore.multi, 2), round(beta.HS.offshore.multi, 2), sep = ", ")
  # reformatted_table$Gammas[idx] <- paste(round(gamma.LS.offshore.multi, 2), round(gamma.MS.offshore.multi, 2), round(gamma.HS.offshore.multi, 2), sep = ", ")
  # reformatted_table$Effective_R0[idx] <- paste(round(R0s.offshore.multi, 2), collapse = ", ")
  # 
  # # Print the final table
  # print(reformatted_table)  
  # 
  # #write to CSV for easy dumping of data into Word table in manuscript
  # output_path <- here("output", "tablex2_density.csv")
  # 
  # # Write the table to a CSV file without row names
  # write.csv(reformatted_table, output_path, row.names = FALSE)
  # 
  # # Print the file path for confirmation
  # cat("Table saved to: ", output_path, "\n")
  # # Print the final table
  # print(reformatted_table)
  # 
  # 
  ################################## Save output ##################################
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "tables_figures_workspace.RData"))
  