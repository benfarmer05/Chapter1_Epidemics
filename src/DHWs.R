# a goal is to use this script as a starting point to better version control and open-source my R scripting
#   (while also accomplishing the goal of bringing in DHW data!)
#
#   1.) Incorporate guidance for 'here' library
#   2.) Deal with 'rm(list=ls())'
#   3.) Get this working in VSCode
#   4.) Get git-pushing functional


library(here)

file_content <- readLines(here("data", "florida_keys.txt"))
metadata <- file_content[1:21]

# Extract data
data_lines <- file_content[22:length(file_content)]

# Read the data into a data frame
data <- read.table(text = data_lines, header = TRUE)


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


reefwatch = read.table(here("data", "florida_keys.txt"))
# survey = read.table(here("data", "SCTLD_END_Vpub_ts.csv"))
# cover = read.csv(here("data", "cover_long.csv"))
# prograte = read.csv(here("data", "SWG_SCTLDprogrates.csv"))


# would like to:
#   - pull the .txt files for the regional and single-pixel station data straight from the internet to here in the script
#   - also could look into the full data product suite (especially the monthly composites!): https://coralreefwatch.noaa.gov/product/index.php
#       - kind of curious if those composites would be the most helpful, or if the virtual stations are "good enough"
