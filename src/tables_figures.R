  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)

  # NOTE - the coral cover conversions done here are not ratios based on particular taxa. this could be improved by re-analyzing CPCe output,
  #         but may not be worth it for this study. the main thing it would affect is presentation of this table, and some interpretation
  #         of the community model
  
  #import workspace from upstream script
  load(here("output/data_processing_workspace.RData"))
  
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
    left_join(susceptible_ref_unique, by = c("site" = "Site")) %>%
    mutate(
      # cover = tottiss / SA.cover.ratio / 100,  # Calculate cover using the site-specific ratio
      cover = tottiss / SA.cover.ratio,  #no need to divide by 100 since this is for cover % presentation purposes
      species_name = recode(spp, !!!species_names)  # Map species codes to full names
    ) %>%
    # Fill missing 'susc' values by taking the first non-NA 'susc' value within the same 'site' and 'species_name' group
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
      total_cover = sum(cover, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step 2: Format the cover in (X, X, X) format for each species and site.
  final_table_formatted <- final_table_summary %>%
    group_by(susc, species_name) %>%
    reframe(
      cover = paste(
        # round(total_cover[site == "off"], 2), 
        # round(total_cover[site == "mid"], 2), 
        # round(total_cover[site == "near"], 2), 
        
        # Apply conditional check for each site to ensure print-out does not reduce very small cover values to 0
        ifelse(round(total_cover[site == "off"], 2) == 0, "<0.01", round(total_cover[site == "off"], 2)), 
        ifelse(round(total_cover[site == "mid"], 2) == 0, "<0.01", round(total_cover[site == "mid"], 2)), 
        ifelse(round(total_cover[site == "near"], 2) == 0, "<0.01", round(total_cover[site == "near"], 2)), 
        sep = ", "
      )
    ) %>%
    arrange(susc, species_name)
  
  # Step 3: Add the gross pathology information based on susceptibility group
  final_table_output <- final_table_formatted %>%
    mutate(
    ) %>%
    # Add the "Total" row
    bind_rows(
      final_table_summary %>%
        group_by(susc) %>%
        summarise(
          species_name = "Total",
          cover = paste(
            # sum(total_cover[site == "off"], na.rm = TRUE), 
            # sum(total_cover[site == "mid"], na.rm = TRUE), 
            # sum(total_cover[site == "near"], na.rm = TRUE),
            
            # Apply rounding to "Total" values
            ifelse(round(sum(total_cover[site == "off"], na.rm = TRUE), 2) == 0, "<0.01", round(sum(total_cover[site == "off"], na.rm = TRUE), 2)),
            ifelse(round(sum(total_cover[site == "mid"], na.rm = TRUE), 2) == 0, "<0.01", round(sum(total_cover[site == "mid"], na.rm = TRUE), 2)),
            ifelse(round(sum(total_cover[site == "near"], na.rm = TRUE), 2) == 0, "<0.01", round(sum(total_cover[site == "near"], na.rm = TRUE), 2)),
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
            ifelse(round(sum(total_cover[site == "off"], na.rm = TRUE), 2) == 0, "<0.01", round(sum(total_cover[site == "off"], na.rm = TRUE), 2)),
            ifelse(round(sum(total_cover[site == "mid"], na.rm = TRUE), 2) == 0, "<0.01", round(sum(total_cover[site == "mid"], na.rm = TRUE), 2)),
            ifelse(round(sum(total_cover[site == "near"], na.rm = TRUE), 2) == 0, "<0.01", round(sum(total_cover[site == "near"], na.rm = TRUE), 2)),
            sep = ", "
          ),
          .groups = "drop"
        )
    ) %>%
    # Sort so that "Total" and "Grand Total" are at the bottom of each 'susc' group
    arrange(susc, species_name == "Total", species_name == "Grand Total", species_name)
  
  #option to remove comma separation
  final_table_output2 <- final_table_output %>%
    mutate(
      cover = gsub(", ", "\t", cover)  # Replace commas with tabs
    )
  
  final_table_output3 = final_table_output %>%
    select(cover)
  
  # Use the "here" package to specify the file path
  output_path <- here("output", "formatted_table.csv")
  
  # Write the table to a CSV file without row names
  write.csv(final_table_output, output_path, row.names = FALSE)
  
  # Print the file path for confirmation
  cat("Table saved to: ", output_path, "\n")  
  # Print the final table
  print(final_table_output)
  
  
  
  
  
  
  
