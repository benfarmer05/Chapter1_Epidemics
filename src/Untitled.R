


# NOTE - sort of failed attempt at editing the starting values of infectious tissue on day 1 at each site.
#   was doing this because of weirdness in multi-host model
#     - 19 feb 2025

# # Create a modified version of only the affected rows
# obs.model_updated <- obs.model %>%
#   filter(Compartment == "Infected", Category == "Total") %>%
#   group_by(Site) %>%
#   mutate(min_tissue = min(tissue[tissue > 0], na.rm = TRUE)) %>%
#   arrange(date) %>%
#   mutate(tissue = case_when(
#     row_number() == min(which(tissue > 0)) ~ min_tissue / case_when(
#       Site == "near" ~ polyp_SA.minimizer.nearshore,
#       Site == "mid" ~ polyp_SA.minimizer.midchannel,
#       Site == "off" ~ polyp_SA.minimizer.offshore
#     ),
#     TRUE ~ tissue  
#   )) %>%
#   ungroup() %>%
#   select(-min_tissue)  
# 
# # Update obs.model in-place without changing row count
# obs.model2 <- obs.model %>%
#   rows_update(obs.model_updated, by = c("Site", "date", "Compartment", "Category"))
# 
# 
# differences <- obs.model2 %>%
#   inner_join(obs.model, by = c("Site", "date", "Compartment", "Category")) %>%
#   filter(tissue.x != tissue.y) %>%
#   select(Site, date, Compartment, Category, tissue_before = tissue.x, tissue_after = tissue.y)
# 
# print(differences)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Create a modified version of only the affected rows
# obs.model_updated <- obs.model %>%
#   filter(Compartment == "Infected") %>%  # Consider all categories for min_tissue
#   group_by(Site) %>%
#   mutate(min_tissue = min(tissue[tissue > 0], na.rm = TRUE)) %>%
#   ungroup() %>%
#   arrange(date) %>%
#   mutate(tissue = case_when(
#     (Category %in% c("Total", "High")) & row_number() == min(which(tissue > 0)) ~ min_tissue / case_when(
#       Site == "near" ~ polyp_SA.minimizer.nearshore,
#       Site == "mid" ~ polyp_SA.minimizer.midchannel,
#       Site == "off" ~ polyp_SA.minimizer.offshore
#     ),
#     TRUE ~ tissue  
#   )) %>%
#   select(-min_tissue)  
# 
# # Update obs.model in-place without changing row count
# obs.model2 <- obs.model %>%
#   rows_update(obs.model_updated, by = c("Site", "date", "Compartment", "Category"))
# 
# # Check differences
# differences <- obs.model2 %>%
#   inner_join(obs.model, by = c("Site", "date", "Compartment", "Category")) %>%
#   filter(tissue.x != tissue.y) %>%
#   select(Site, date, Compartment, Category, tissue_before = tissue.x, tissue_after = tissue.y)
# 
# print(differences)
# 
# 
# 
# 











# Find the smallest nonzero tissue value across all categories within each Site
obs.model_updated <- obs.model %>%
  filter(Compartment == "Infected") %>%
  group_by(Site) %>%
  mutate(min_tissue = min(tissue[tissue > 0], na.rm = TRUE)) %>%
  ungroup()

# Apply changes only to "Total" and "High"
obs.model_updated <- obs.model_updated %>%
  arrange(Site, date) %>%
  group_by(Site, Category) %>%  # Ensure min(which(tissue > 0)) is calculated within each Category
  mutate(tissue = case_when(
    Category %in% c("Total", "High") & row_number() == min(which(tissue > 0)) ~ min_tissue / case_when(
      Site == "near" ~ polyp_SA.minimizer.nearshore,
      Site == "mid" ~ polyp_SA.minimizer.midchannel,
      Site == "off" ~ polyp_SA.minimizer.offshore
    ),
    TRUE ~ tissue
  )) %>%
  ungroup() %>%
  select(-min_tissue)

# Update obs.model in-place without changing row count
obs.model2 <- obs.model %>%
  rows_update(obs.model_updated, by = c("Site", "date", "Compartment", "Category"))

# Check differences
differences <- obs.model2 %>%
  inner_join(obs.model, by = c("Site", "date", "Compartment", "Category")) %>%
  filter(tissue.x != tissue.y) %>%
  select(Site, date, Compartment, Category, tissue_before = tissue.x, tissue_after = tissue.y)

print(differences)






# 
# #establish a scalar with which to minimize the initial infectious tissue on the first infection day
# polyp_SA.minimizer = 5
# 
# 
# # Find the smallest nonzero tissue value across all categories within each Site
# obs_updated <- obs %>%
#   filter(Compartment == "Infected") %>%
#   group_by(Site) %>%
#   mutate(min_tissue = min(tissue[tissue > 0], na.rm = TRUE)) %>%
#   ungroup()
# 
# 
# 
# # Apply changes only to "Total" and "High"
# obs_updated <- obs_updated %>%
#   arrange(Site, date) %>%
#   group_by(Site, Category) %>%  # Ensure min(which(tissue > 0)) is calculated within each Category
#   mutate(tissue = case_when(
#     Category %in% c("Total", "High") & row_number() == min(which(tissue > 0)) ~ min_tissue / case_when(
#       Site == "Nearshore" ~ polyp_SA.minimizer, #polyp_SA.minimizer.nearshore
#       Site == "Midchannel" ~ polyp_SA.minimizer, #polyp_SA.minimizer.midchannel
#       Site == "Offshore" ~ polyp_SA.minimizer #polyp_SA.minimizer.offshore
#     ),
#     TRUE ~ tissue  
#   )) %>%
#   ungroup() %>%
#   select(-min_tissue)
# 
# # Update obs in-place without changing row count
# obs2 <- obs %>%
#   rows_update(obs_updated, by = c("Site", "date", "Compartment", "Category"))
# 
# # Check differences
# differences <- obs2 %>%
#   inner_join(obs, by = c("Site", "date", "Compartment", "Category")) %>%
#   filter(tissue.x != tissue.y) %>%
#   select(Site, date, Compartment, Category, tissue_before = tissue.x, tissue_after = tissue.y)
# 
# print(differences)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #establish a scalar with which to minimize the initial infectious tissue on the first infection day
# polyp_SA.minimizer = 5
# 
# 
# # Find the smallest nonzero tissue value within "Total" category for each Site
# min_tissue_total <- obs %>%
#   filter(Compartment == "Infected", Category == "Total") %>%
#   group_by(Site) %>%
#   summarise(min_tissue = min(tissue[tissue > 0], na.rm = TRUE), .groups = "drop")
# 
# # Apply changes only to "Total" and "High"
# obs_updated <- obs %>%
#   left_join(min_tissue_total, by = "Site") %>%
#   arrange(Site, date) %>%
#   group_by(Site, Category) %>%
#   mutate(tissue = case_when(
#     Compartment == "Infected" & Category %in% c("Total") & row_number() == min(which(tissue > 0)) ~ min_tissue / case_when(
#       Site == "Nearshore" ~ polyp_SA.minimizer, 
#       Site == "Midchannel" ~ polyp_SA.minimizer, 
#       Site == "Offshore" ~ polyp_SA.minimizer 
#     ),
#     TRUE ~ tissue  
#   )) %>%
#   ungroup() %>%
#   select(-min_tissue)
# 
# # Update obs in-place without changing row count
# obs2 <- obs %>%
#   rows_update(obs_updated, by = c("Site", "date", "Compartment", "Category"))
# 
# # Check differences
# differences <- obs2 %>%
#   inner_join(obs, by = c("Site", "date", "Compartment", "Category")) %>%
#   filter(tissue.x != tissue.y) %>%
#   select(Site, date, Compartment, Category, tissue_before = tissue.x, tissue_after = tissue.y)
# 
# print(differences)
# 
