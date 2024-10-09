  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  
  #this data comes from: https://coralreefwatch.noaa.gov/product/index.php
  # - specifically, 'florida_keys.txt' is Regional Virtual Station data, specific to the Florida region
  # - from link above, navigate to point 3 (Regional Virtual Stations and Bleaching Alert Emails), and then click 'home page' for
  #     'Regional Virtual Stations (5km)'. this will open a page with a world map.
  # - navigate on the map to Florida Keys region (not Southeast Florida region). click 'Time Series Graphs & Data' on the prompt, which will open a new page
  # - click 'Time Series Data' for Florida Keys. this will open the 'florida_keys.txt' ASCII data directly for viewing
  # - you can go back and right click the link itself to save the link as a text file to your local directory
  file_content <- readLines(here("data", "florida_keys.txt"))
  metadata <- file_content[1:21]
  
  # NOTE - there is also Single-Pixel Virtual Station data that is closer to the Summerland Key area of the study site
  #         - https://coralreefwatch.noaa.gov/product/vs_single_pixel_exp/florida_keys.php
  #         - these data actually indicate less DHWs. that said, going with the regional source that was already used by Williams 2021
  #         - can always discuss this in manuscript revisions as needed
  #      - monthly composites from CRW may also be worth a look
  
  # Extract data
  data_lines <- file_content[22:length(file_content)]
  
  # Read the data into a data frame
  DHW.CRW <- read.table(text = data_lines, header = TRUE)
  
  # Convert metadata to a named list
  metadata_list <- list(
    Name = sub("Name: ", "", metadata[1]),
    Polygon_Middle_Longitude = as.numeric(sub("Polygon Middle Longitude: ", "", metadata[3])),
    Polygon_Middle_Latitude = as.numeric(sub("Polygon Middle Latitude: ", "", metadata[5])),
    Averaged_Maximum_Monthly_Mean = as.numeric(sub("Averaged Maximum Monthly Mean: ", "", metadata[7])),
    Averaged_Monthly_Mean = as.numeric(unlist(strsplit(sub("Averaged Monthly Mean (Jan-Dec): ", "", metadata[9]), " "))),
    First_Valid_DHW_Date = as.Date(paste(unlist(strsplit(sub("First Valid DHW Date: ", "", metadata[11]), " ")), collapse = "-"), format = "%Y-%d-%m"),
    First_Valid_BAA_Date = as.Date(paste(unlist(strsplit(sub("First Valid BAA Date: ", "", metadata[13]), " ")), collapse = "-"), format = "%Y-%d-%m")
  )
  
  #pull Williams 2021 environmental data to compare
  # - these data are identical to Coral Reef Watch data above, but binned to the study timepoints
  # - 'sst' column is 'SST.90th_HS'
  # - 'bleachalert' is 'BAA_7day_max'
  # - 'dhw' is 'DHW_from_90th_HS.1'
  # - 't.anom' is 'SSTA.90th_HS' (SST anomaly)
  # - not sure what 'tl.means', 'tl.se' are
  # - 'all.incidence.newInc' is from their disease data I believe
  DHW.Williams = read.csv(here("data/extended_envdisFig3.csv"))

  #filter datasets
  DHW.Williams = DHW.Williams %>%
    mutate(X = as.Date(gsub("^X", "", X), format = "%m.%d.%y")) %>%
    rename(date = X) %>%
    filter(date >= as.Date("2018-05-01") & date <= as.Date("2019-12-06")) #include first day of surveying in the field
    # filter(date >= as.Date("2018-08-17") & date <= as.Date("2019-12-06")) #skip to first SCTLD day
  
  DHW.CRW <- DHW.CRW %>%
    mutate(date = as.Date(paste(YYYY, MM, DD, sep = "-"))) %>%
    filter(date >= as.Date("2018-05-01") & date <= as.Date("2019-12-06")) #include first day of surveying in the field
    # filter(date >= as.Date("2018-08-17") & date <= as.Date("2019-12-06")) #skip to first SCTLD day
    select(-YYYY, -MM, -DD) %>%
    select(date, everything())
  