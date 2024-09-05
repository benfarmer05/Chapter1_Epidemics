
  # .rs.restartR(clean = TRUE)

  rm(list=ls())

  
  library(tidyverse)
  library(data.table)
  library(mgcv)
  library(gratia)
  library(here)

  survey = read.csv(here("data", "SCTLD_END_Vpub_ts.csv"))
  cover = read.csv(here("data", "cover_long.csv"))
  prograte = read.csv(here("data", "SWG_SCTLDprogrates.csv"))
  
  #   LOW SUSCEPTIBILITY (LS)
  #         - Slow onset/slow tissue loss [SSID] -> high cover
  #         - Moderate onset/slow tissue loss [SINT, PCLI] -> high cover
  #   MODERATE SUSCEPTIBILITY (MS)
  #         - Moderate onset/moderate tissue loss [MCAV] -> moderate cover
  #         - Fast onset/slow tissue loss [OANN, OFAV] -> low cover
  #         - Fast onset/moderate tissue loss [SBOU] -> low cover
  #   HIGH SUSCEPTIBILITY (HS)
  #         - Fast onset/fast tissue loss [DSTO, CNAT, PSTR, DLAB] -> moderate cover
  
  #pull wide-format temporal data from lower Florida Keys (Williams et al. 2021), which recorded individual coral colonies; categorize species
  #   by susceptibility group
  survey = survey %>%
    mutate(
      Site_type = case_when(
        Site == 1 ~ "Midchannel", #Site 1, transects 23 and 25 is Wonderland site
        Site == 2 ~ "Offshore", #Site 2, transects 27 and 28 is ACER site
        Site == 3 ~ "Nearshore" #Site 3, transects 45 and 47 is N. Birthday site
      ),
      Sps = as.factor(Sps),
      Sus_Cat = case_when(
        Sps %in% c('PCLI', 'SINT', 'SSID') ~ 'LS', #lowest rates of progression - however, OANN/OFAB have short disease onset time
        Sps %in% c('OANN', 'OFAV', 'MCAV', 'SBOU') ~ 'MS', #moderate rates of progression - but MCAV slow disease onset! SBOU quick onset
        Sps %in% c('CNAT', 'DLAB', 'DSTO', 'MMEA', 'MYCE', 'PSTR') ~ 'HS', #quickest rates of  progression and shortest disease onset times
        Sps %in% c('AAGA', 'ACER', 'OCUL', 'ODIF', 'PAST', 'PDIV', 'PPOR', 'SRAD') ~ 'Unaffected'
      )
    )
  
  #pull middle Florida Keys (Sharp et al. 2020) dataset, which had similar data as 'survey' but measured colonies in 3 dimensions
  #   rather than one. using this dataset to create a statistical association between maximum diameter & colony surface area (SA) by assuming
  #   a hemi-ellipsoid (Holstein et al. 2015), and then predicting SA in 'survey' dataset using maximum diameter alone
  Sharp2020 = read.csv(here("data", "Sharp_Maxwell_2020.csv"))
  
  Sharp2020 = select(Sharp2020, Plot, Transect, Coral, Species, Width, Width_2, Height)
  
  #filter meaningless data
  Sharp2020[Sharp2020==-99] = NA 
  Sharp2020 = Sharp2020[complete.cases(Sharp2020),]
  
  #setting '0' cm heights as an arbitrarily small amount (cm) to allow the hemi-ellipsoid estimation to function correctly. coral recruits 
  #   have very little tangible height, and were recorded underwater as 0 cm height
  Sharp2020$Height[which(Sharp2020$Height == 0)] = 0.01
  
  #filter out repeats (not interested in a time series here)
  Sharp2020 = distinct(Sharp2020)
  
  #set up the variables for the hemi-ellipsoid estimation. p is a dimensionless constant; all else in square cm
  a = Sharp2020$Height; b = Sharp2020$Width_2/2; c = Sharp2020$Width/2; p = 1.6075
  Sharp2020$SA = 2*pi*((((a*b)^p + (a*c)^p + (b*c)^p)/3)^(1/p))
  Sharp2020$SA = Sharp2020$SA/10000 #in sq m
  
  #construct dataframe of estimated surfacea area, as well as measured maximum diameter ('width' in this dataset)
  replen = length(Sharp2020$Width)
  SA = data.frame(x = Sharp2020$Width, y = Sharp2020$SA)
  
  # #a measure to prevent over-predicting surface area of the largest corals from small sample size [decided as unnecessary]
  # SA = SA %>% filter(x < 145)

  plot(SA$x, SA$y)
  
  #Generalized additive model (GAM) to predict colony SA in lower Florida Keys
  SA.GAM = gam(data = SA, y ~ s(x), family = gaussian, method = "REML")
  summary(SA.GAM)
  AIC(SA.GAM)
  plot(SA.GAM, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  draw(SA.GAM, select = 1)
  gam.check(SA.GAM)
  appraise(SA.GAM)
  
  x = survey$Max_width
  x.GAM = as.data.frame(x)
  x.GAM$x = as.numeric(x.GAM$x)
  y = predict(SA.GAM, newdata = x.GAM, type = "response")
  y[y<0] = 0.00005 #ensure GAM-predicted surface area does not dip below zero for small recruits. arbitrary value
  survey$Tissue = y
  plot(x,y, xlim = c(0, 460), xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  sum(survey$Tissue) #quick summary to assess overall tissue SA across all sites
  
  #incorporate pre-SCTLD old mortality from lower Florida Keys plots (same as the ones in 'survey'). last observation of old mortality
  #   before SCTLD arrived was 17 August 2018
  #   SCTLD reached the three sites:
  #   - Offshore: 30 October 2018 [74 days after last observation of old mortality]
  #   - Midchannel: 29 November 2018 [30 days after Offshore]
  #   - Nearshore: 8 February 2019 [71 days after Midchannel; 101 days after Offshore]
  old_mortality = read.csv(here("data", "Updated_June19_SWG.csv"))
  
  #replace 'FWRI database' instances with 'NA' string for clarity in relevant 'New_death' columns
  old_mortality$New_death_110918 = gsub('FWRI Database', 'NA', old_mortality$New_death_110918)
  old_mortality$New_death_103018 = gsub('FWRI Database', 'NA', old_mortality$New_death_103018)
  
  #filter to only most recent recordings of mortality pre-SCTLD
  #   - normal health states / 'Found' statuses for 10 May 2018. but lots of notes (corals turned over especially)
  #   - no new mortality recorded 1 June or 21 June 2018; those are excluded. normal health states / notes / 'Found' statuses as well
  #   - normal health states / 'Found' statuses for 16 July 2018
  old_mortality = old_mortality %>% select(Plot, Sps, Max_width, Coral_ID, New_death_051018, Old_death_051018, New_death_060118,
                                           New_death_062118, New_death_071618, Notes_071618, New_death__081718, Old_death_081718,
                                           Health_state_081718, Notes_081718, Found_081718, SCTLD_081718)

  #create new column that indicates total mortality pre-SCTLD. a tally of new & old mortality columns
  #   must first replace NA values in 'Old_death_081718' with 0 so that summation can be performed correctly
  #     NOTE - assuming that NA actually means there wasn't old mortality. only affects 3 corals, so should not be a big deal but should
  #            be considered. this makes calculations of tissue loss later on more simple and avoids NA problems
  #               - corals are:
  #     - 1_p24_t8_s5_c21_MCAV
  #     - 2_p27_t8_s5_c30_DSTO
  #     - 3_p46_t4_s0_c30_SINT
  old_mortality$Old_death_081718[is.na(old_mortality$Old_death_081718)] = 0
  old_mortality$tot_mortality = old_mortality$New_death__081718 + old_mortality$Old_death_081718
  old_mortality$tot_mortality = as.numeric(old_mortality$tot_mortality)
  
  #simplify dataframe
  old_mortality = old_mortality %>% select(Plot, Sps, Max_width, Coral_ID, tot_mortality)
  
  #label the rows which match plots later monitored for SCTLD
  old_mortality = old_mortality %>%
    mutate(
      Site_type = case_when(
        Plot == 23 | Plot == 25 ~ "Midchannel",
        Plot == 27 | Plot == 28 ~ "Offshore",
        Plot == 45 | Plot == 47 ~ "Nearshore"
      )
    )
  
  #filter down to only plots monitored for SCTLD. there were extra plots initially part of monitoring!
  old_mortality = old_mortality %>%
    filter(
      Site_type == "Midchannel" | Site_type == "Offshore" | Site_type == "Nearshore"
    )
  
  #mismatched data summary (just for peace of mind):
  #   - 1_p25_t9_s5_c4_PAST. unclear what the mismatch is
  #   - 2_p27_t10_s0_c1_OFAV. 'Sps' column labeled as SINT in 'survey'. Brought this up in Mote call
  #   - 2_p27_t10_s0_c3_OFAV. 'Sps' column labeled as SINT in 'survey'. Brought this up in Mote call
  #   *-* 2_p28_t4_s5_c8_SINT. in 'survey' only; no coordinates info
  #   - 3_p47_t1_s0_c1_MCAV. 'Max_width' misentered somewhere. It is 50 in 'survey', but 25 in 'old_mortality'. Brought this up in Mote call
  #   *-* 3_p47_t5_s0_c24_PAST. in 'survey' only; no coordinates info
  # mismatches = anti_join(survey, old_mortality, by = c("Coral_ID")) # *-*
  mismatches = anti_join(survey, old_mortality) #ignore all the ones due to '(1)' or similar in the name. they're fine
  
  # Find and edit the column/row indices where 'Sps' needs to be updated for mismatched data
  row_index = which(survey$Coral_ID == '2_p27_t10_s0_c1_OFAV')
  survey$Sps[row_index] = 'OFAV'
  
  row_index = which(survey$Coral_ID == '2_p27_t10_s0_c3_OFAV')
  survey$Sps[row_index] = 'OFAV'
  
  row_index = which(old_mortality$Coral_ID == '3_p47_t1_s0_c1_MCAV')
  old_mortality$Max_width[row_index] = 50

  #function to replace ' (1)', ' (2)', or ' (3)' with ''
  replace_string <- function(data) {
    for (i in 1:ncol(data)) {
      if (is.character(data[, i])) {  # Check if the column is character
        data[, i] <- gsub(' \\((1|2|3)\\)', '', data[, i])
      }
    }
    return(data)
  }
  
  #replace ' (1)', ' (2)', or ' (3)' with '' in the entire dataframe
  survey = replace_string(survey)
  old_mortality = replace_string(old_mortality)
  
  #merge in the pre-SCTLD old mortality data
  survey = left_join(
    survey,
    old_mortality,
    by = c("Coral_ID", "Plot", "Site_type", 'Sps', 'Max_width'),
    copy = FALSE,
    suffix = c(".x", ".y"),
    keep = NULL
  )
  
  #these two corals are in 'survey' only, so will also assume that they started the SCTLD outbreak with zero old mortality
  row_index = which(survey$Coral_ID == '2_p28_t4_s5_c8_SINT')
  survey$tot_mortality[row_index] = 0
  row_index = which(survey$Coral_ID == '3_p47_t5_s0_c24_PAST')
  survey$tot_mortality[row_index] = 0
  
  #this coral for some reason did not join correctly (possibly a subtle issue with the 'char' string), but does indeed have 10% old mortality
  row_index = which(grepl('p25_t9_s5_c4_P', survey$Coral_ID))
  survey$tot_mortality[row_index] = 10
  
  ###pivot survey data to long format
  survey_nearshore = survey %>%
    subset(
      Site_type == "Nearshore"
    ) %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "date"
    )
  
  survey_midchannel = survey %>%
    subset(
      Site_type == "Midchannel"
    ) %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "date"
    )
  
  survey_offshore = survey %>%
    subset(
      Site_type == "Offshore"
    ) %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "date"
    )
  
  survey_long = rbind(survey_nearshore, survey_midchannel, survey_offshore)
  survey_long$tot_diseased = as.factor(survey_long$tot_diseased)
  levels(survey_long$tot_diseased)
  
  #convert character-based coral index to numeric-based index. there are 2012 corals, and now they have a unique numeric value across all
  # timepoints
  hold = survey_long %>%
    group_by(Coral_ID) %>%
    mutate(coral_numID = cur_group_id())
  survey_long$coral_numID = hold$coral_numID
  rm(hold)
  survey_long = survey_long %>% select(coral_numID, everything())
  
  #append the same coral numerical ID to 'prograte', a dataframe containing time-series lesion progression rates of each diseased coral colony
  prograte = prograte %>%
    mutate(coral_numID = survey_long$coral_numID[match(Coral_ID, survey_long$Coral_ID)])
  
  #add index to each dataframe row to make breaking it apart and then merging & joining later easier
  survey_long$ID = seq.int(nrow(survey_long))
  survey_long = survey_long %>% select(ID, everything())
  
  ###extract instantaneous lesion areal extents (m2) for specific coral host/date combinations from 'progrates' dataframe
  # reformat progrates
  names(prograte)[names(prograte) == 'X'] <- '_X'#have to rename this one since it starts with an X
  for(i in 1:ncol(prograte)){
    if(substr(colnames(prograte[i]), start = 1, stop = 1) == "X"){
      newdatename = as.character(as.POSIXct(str_sub(colnames(prograte[i]), 2), format = "%m.%d.%y"))
      colnames(prograte)[i] = newdatename
    }
  }
  
  survey_long$percloss = NA # percentage loss of live tissue between timepoints
  survey_long$percloss = as.numeric(survey_long$percloss)
  survey_long$progdays = NA # percentage loss of live tissue between timepoints
  survey_long$progdays = as.numeric(survey_long$progdays)
  survey_long$percinf = NA #estimated instantaneous (daily) %loss of live tissue - proxy of instantaneous %infectious tissue
  survey_long$percinf = as.numeric(survey_long$percinf)
  survey_long$starttiss = NA
  survey_long$starttiss = as.numeric(survey_long$starttiss)
  survey_long$inftiss = NA
  survey_long$inftiss = as.numeric(survey_long$inftiss)
  survey_long$remaintiss = NA
  survey_long$remaintiss = as.numeric(survey_long$remaintiss)
  survey_long$removedtiss = NA
  survey_long$removedtiss = as.numeric(survey_long$removedtiss)
  
  #initialize amount of tissue on each coral colony
  survey_long$starttiss = 1-(survey_long$tot_mortality/100) #tissue percentage [scalar] present on the colony at beginning of SCTLD outbreak, after accounting for old mortality

  #convert date format to POSIXct
  survey_long$date = str_sub(survey_long$date, 2)
  survey_long$date = as.factor(survey_long$date)
  survey_long$date = as.POSIXct(survey_long$date, format = "%m.%d.%y")
  
  #sort by coral numerical ID and date to ensure no unintended behavior when calculating tissue loss from lesions through time
  survey_long = survey_long %>%
    arrange(coral_numID, date)
  
  #filter down to post-SCTLD introduction (October 30th 2018) infections
  surveydiseased = survey_long %>%
    filter(
      tot_diseased == "Dis",
      date >= as.POSIXct("2018-10-30")
    ) %>%
    arrange(coral_numID, date) #failsafe to ensure preservation of ID-sorting
  
  #remove exact rows from 'survey_long' as are being extracted for 'surveydiseased', to easily rbind back together after data wrangling
  survey_trimmed = anti_join(survey_long, surveydiseased) %>%
    arrange(coral_numID, date)
  
  # 'SURPRISE DEAD CORALS'
  #
  #there were 57 corals that died suddenly between timepoints but were otherwise observed as healthy. these were all quite small and likely died from
  # SCTLD during the first and second infection waves. to account for this cryptic tissue loss, convert those final timepoints into 100% mortality from
  # disease
  # note - to find these, filter by 2019-12-06 and tot_diseased = 'Health' & value = 'Dead'
  #   - almost all of these corals were small, HS corals which makes sense. important to keep track of their loss between timepoints
  #   - also, applying this modification now means that a small 11 cm-diameter CNAT  becomes our patient zero - it was already dead before the DSTO on
  #       October 30th, 2018.
  
  #   - Surprise-dead corals which had 100% old mortality in pre-SCTLD (mortality) dataset (so, are being *removed* entirely before simulation)
  #       - 2_p27_t9_s5_c2_DSTO - Maybe a discrepancy in mortality dataset; it suddenly died 09-17-2019 in new dataset. said 100% recent mortality on 08-17-2018 in old dataset
  #             - SW: Probably dead-dead on 8-17-2018. Based on field notes
  #       - 2_p27_t8_s5_c21_DSTO - Maybe a discrepancy in mortality dataset; it suddenly died 05-02-2019. 100% old mortality suddenly on 08-17-2018 in old dataset. VERY small coral (2 cm)
  #             - SW: dead-dead on 8-17-2018.
  #       - 1_p23_t1_s0_c6_CNAT - Maybe a discrepancy in mortality dataset; it suddenly died 08-30-2018. 75% old mortality suddenly on 05-10-2018, then suddenly 100% old mortality or not found on 08-17-2018.
  #             _ SW: dead-dead on 8-17-2018. was 75% dead by 5-10-2018, if that is useful
  
  #filter surprise-dead corals (read details above), by their health condition on the last date of surveying (2019-12-06)
  surprise.dead = survey_trimmed %>%
    filter(date == as.POSIXct("2019-12-06"), tot_diseased == 'Health', value == 'Dead') %>%
    arrange(coral_numID, date)

  #pull the full T1 - T26 rows for each surprise-dead coral
  first.dead.full = survey_trimmed %>%
    filter(coral_numID %in% surprise.dead$coral_numID) %>%
    arrange(coral_numID, date)

  #filter to the date that the surprise-dead coral was first documented as dead
  first.dead = first.dead.full %>%
    filter(grepl("\\Dead", value)) %>%
    group_by(coral_numID) %>%
    top_n(n=1, wt=desc(date)) %>%
    ungroup() %>%
    mutate(dead_date = date, #create duplicate column to be used later, and indicate the first documented date of complete mortality
           percloss = 100) %>% #set mortality from SCTLD to '100%'
    arrange(coral_numID, date)
  
  first.dead.full = first.dead.full %>% left_join(
    first.dead %>% select(coral_numID, dead_date),
    by = c("coral_numID")
  ) %>%
    arrange(coral_numID, date)
  
  # Update percloss in first.dead.full based on values from first.dead
  # Perform the left join to bring in percloss and date from first.dead
  first.dead.full <- first.dead.full %>%
    left_join(
      first.dead %>% select(coral_numID, date, percloss),
      by = c("coral_numID", "date")) %>%
    mutate(
      # Update percloss where available from first.dead
      percloss = coalesce(percloss.y, percloss.x),
    ) %>%
    select(-percloss.x, -percloss.y) %>%
    arrange(coral_numID, date) # Sort by coral_numID and date
  
  # Initialize prior_date column & assign any straggling NA's
  first.dead.full$prior_date <- NA
  first.dead.full$percloss <- ifelse(is.na(first.dead.full$percloss), NA, first.dead.full$percloss)
  
  # Loop through the data frame
  for (i in 2:nrow(first.dead.full)) {
    if (first.dead.full$date[i] == first.dead.full$dead_date[i]) {
      # Set percloss to 100% for the current dead_date
      first.dead.full$percloss[i] <- 100
      
      # Calculate percloss for the previous row (prior_date)
      first.dead.full$prior_date[i] <- first.dead.full$date[i-1]
      first.dead.full$percloss[i-1] <- 100 / as.numeric(difftime(first.dead.full$dead_date[i], first.dead.full$date[i-1], units = "days"))
    }
  }
  
  # Drop the prior_date column if you don't need it anymore
  first.dead.full$prior_date <- NULL
  
  #update original dataframe with surprise-dead coral mortality
  survey_trimmed = survey_trimmed %>%
    left_join(
      first.dead.full %>% select(coral_numID, date, dead_date, percloss),
      by = c("coral_numID", "date")) %>%
    mutate(
      percloss.x = coalesce(percloss.y, percloss.x)
    ) %>%
    rename(percloss = percloss.x) %>%
    select(-percloss.y) %>%
    arrange(coral_numID, date)

  #variable name changes for clarity
  names(survey_trimmed)[names(survey_trimmed) == 'tot_mortality'] = 'old_mortality'
  names(survey_trimmed)[names(survey_trimmed) == 'Sps.x'] = 'Sps'
  names(survey_trimmed)[names(survey_trimmed) == 'Max_width.x'] = 'Max_width'
  names(surveydiseased)[names(surveydiseased) == 'tot_mortality'] = 'old_mortality'
  names(surveydiseased)[names(surveydiseased) == 'Sps.x'] = 'Sps'
  names(surveydiseased)[names(surveydiseased) == 'Max_width.x'] = 'Max_width'
  
  # Prepare a dataframe for calculating infected tissue in surprise-dead corals (T1 - T26)
  # only calculating instanteous infected tissue (tissue loss / sloughing within 24 hours) for corals that haven't completed died already
  # NOTE - the way 'percloss' is handled here is different than the loop for confirmed diseased corals below. here, it was assigned
  #         manually be me, to backtrack daily percentage loss back to the last timepoint (and this translates directly to 
  #         instantaneous infected tissue in a day). but below, percloss was a measure of how much tissue was lost between timepoints
  #         in confirmed SCTLD-infected corals. that's why it needs be backtracked and then divided by the number of days between
  #         timepoints. I might have benefited from using a different term than 'percloss' for the surprise-dead corals - could maybe
  #         go back and call it 'instantloss' or 'dailyloss'. similar nomenclature may be useful below as well
  surprise.dead.infections <- survey_trimmed %>%
    filter(coral_numID %in% surprise.dead$coral_numID) %>%
    group_by(coral_numID) %>%
    filter(!all(old_mortality == 100)) %>%  # Exclude coral_numIDs with 100% old mortality before SCTLD surveying
    ungroup() %>%
    mutate(
      percinf = if_else(!is.na(percloss) & percloss > 0 & percloss != 100, percloss, NA_real_), #exclude corals that are 100% dead (from SCTLD)
      inftiss = if_else(!is.na(percinf),
                        (percinf / 100) * Tissue * starttiss,
                        NA_real_)
    ) %>%
    arrange(coral_numID, date)
    

  #update original dataframe with surprise-dead coral mortality
  survey_trimmed = survey_trimmed %>%
    left_join(
      surprise.dead.infections %>% select(coral_numID, date, percinf, inftiss),
      by = c("coral_numID", "date")) %>%
    mutate(
      percinf.x = coalesce(percinf.y, percinf.x),
      inftiss.x = coalesce(inftiss.y, inftiss.x)
    ) %>%
    rename(percinf = percinf.x, inftiss = inftiss.x) %>%
    select(-percinf.y, -inftiss.y) %>%
    arrange(coral_numID, date)
  #
  # 'SURPRISE DEAD CORALS'
  
  
  
  # STOPPING POINT 
  #   - 4 September 2024: okay this is driving me nuts and I need to set it aside for a second, but basically what I need to do is
  #                     identify the patient zero corals for each site, and then make sure there is a timepoint 30 days preceding which
  #                     I can backtrack the infection to (and if not, add one - becoming the new patient zero date). this approach will
  #                     then be recapitulated downstreamm in the script to apply to all corals that get infected. the "Dead" day for a coral
  #                     is final, though - there is backtracked infection, but not instantaneous infection like with any ongoing affliction
  # - Offshore patient zero: 2_p27_t2_s0_c1_DSTO, 10-30-2018, 90% loss. backtrack a week? this is only a 5 inch coral
  # - Midchannel patient zero: 1_p25_t2_s0_c22_DSTO, 11-29-2018, 10% loss
  # - Nearshore patient zero: 3_p47_t3_s0_c8_PSTR (10% loss) and 3_p47_t4_s5_c15_PSTR (5% loss), 2019-02-08
  #
  # This ensures that I am properly backtracking the initial loss of tissue. Currently, I am assuming that infection was progressing
  #   on the coral all the way up until the moment the surveyor dove up down logged an amount of loss, but I'd rather assume that the
  #   coral became infected right as the surveyor was ascending back to the surface (because it is convenient)
  
  
  
  #steps:
  # 1.) backtrack percloss of current timepoint to last timepoint (account for patient zeros somehow): percloss/progdays
  # 2.) apply that backtracked value (inftiss) to the last timepoint, not the current one
  # 3.) the last inftiss should be on the date preceding the final day of new percloss. whether that is lesion cessation or host death
  
  # 'CORALS DOCUMENTED AS DISEASED'
  #
  #prepare a dataframe for calculating infected tissue in corals documented as diseased in the field
  observed.infected = surveydiseased %>%
    mutate(percloss = NA_real_, progdays = NA_real_, percinf = NA_real_, inftiss = NA_real_) %>%
    arrange(coral_numID, date) 
  
  # Initialize error log
  error_log = data.frame(
    prevdate = character(),
    curr_coral_ID = numeric(),
    curr_ID = numeric(),
    currdate = character(),
    stringsAsFactors = FALSE
  )
  
  #loop to calculate instantaneous infected tissue for confirmed SCTLD-infected corals
  for(i in 1:nrow(observed.infected)){
    
    # #test
    # i = 589 #single out the problem coral (2_p27_t2_s0_c1_DSTO). this was a patient zero
    # # i = 5 #first infected coral in dataframe
    
    # Extract values from currIDsdates
    curr_values <- observed.infected %>%
      slice(i) %>%
      mutate(
        coral_numID = as.numeric(coral_numID),
        ID = as.numeric(ID),
        date = as.character(date)
      )
    
    curr_coral_ID = curr_values$coral_numID
    curr_ID = curr_values$ID
    currdate = curr_values$date
    
    dateind = match(curr_values$date, colnames(prograte))
    percloss = as.numeric(as.matrix(setDT(prograte)[coral_numID == curr_coral_ID, ])[dateind])
    IDslice = surveydiseased %>%
      filter(ID == curr_ID)
    availtiss = IDslice$starttiss

    # Only calculate if percloss is valid
    if(percloss > 0 & !is.na(percloss)){
      
      # Update percloss in observed.infected
      observed.infected[i, "percloss"] <- percloss

      # NOTE - the below assumes 'prograte' date columns are in order, which they are. but another approach may pivot
      #         'prograte' to long format earlier in the script for easier sorting
      prevdate = prograte %>%
        filter(
          coral_numID == curr_coral_ID
        ) %>%
        select(0,grep(currdate, colnames(prograte))-1) %>% # NOTE - line of interest, regarding above comment
      colnames()
      
      # Check if prevdate is an invalid date format
      if (length(prevdate) == 0 || is.na(tryCatch(as.Date(prevdate, format="%Y-%m-%d"), error = function(e) NA))){ # Log the error
        error_log = rbind(error_log, data.frame(
          prevdate = paste(prevdate, collapse = ", "),
          curr_coral_ID = curr_coral_ID,
          curr_ID = curr_ID,
          currdate = currdate,
          stringsAsFactors = FALSE
        ))
        
        # Skip to the next iteration. this is an error checker for any problem corals (i.e., they don't have a time point before 10-30-2018)
        next
      }
      
      # Calculate progdays and update in observed.infected
      progdays = as.numeric(difftime(currdate, prevdate, units = "days"))
      observed.infected[i, "progdays"] <- progdays
      
      # Calculate percinf and update in observed.infected
      percinf = percloss / progdays
      observed.infected[i-1, "percinf"] <- percinf
      
      # Calculate inftiss and update in observed.infected
      tissue = IDslice$Tissue
      inftiss = (percinf / 100) * tissue * availtiss
      observed.infected[i-1, "inftiss"] <- inftiss
    }
  }
  
  # Print the error log if any errors were logged
  #   - this is particularly for '2_p27_t2_s0_c1_DSTO', since it is the one patient zero coral with a special date
  if (nrow(error_log) > 0) {
    print("Error log:")
    print(error_log)
  }
  
  #update original dataframe with confirmed diseased coral infections
  survey_trimmed = survey_trimmed %>%
    left_join(
      observed.infected %>% select(coral_numID, date, percloss, progdays, percinf, inftiss),
      by = c("coral_numID", "date")) %>%
    mutate(
      percloss.x = coalesce(percloss.y, percloss.x),
      progdays.x = coalesce(progdays.y, progdays.x),
      percinf.x = coalesce(percinf.y, percinf.x),
      inftiss.x = coalesce(inftiss.y, inftiss.x)
    ) %>%
    rename(percloss = percloss.x, progdays = progdays.x, percinf = percinf.x, inftiss = inftiss.x) %>%
    select(-percloss.y, -progdays.y, -percinf.y, -inftiss.y) %>%
    arrange(coral_numID, date)
  
  
  
  
  #update main dataframe with diseased tissue information
  survey_tissue = rbind(survey_trimmed, surveydiseased) %>%
    arrange(coral_numID, date)
  
  
  
  #
  # 'CORALS DOCUMENTED AS DISEASED'
  
  # CALCULATING REMOVED / REMAINING TISSUE THROUGH TIME
  #
  #   - removed tissue is not instantaneous in the same way that 'inftiss' is. removed tissue represents an accumulated lost tissue
  #
  num_IDs = max(survey_tissue$coral_numID)

  #loop to calculate removed tissue per coral
  for(i in 1:num_IDs){

    # #test
    # # i = 1822 #coral is 3_p47_t5_s0_c12_MCAV (good for testing calculations on corals w/ intermittent infections)
    # i = 1534 #coral is 3_p45_t9_s5_c27_SBOU (good for testing both intermittent infection AND total mortality of a coral)

    curr_coral = survey_tissue %>%
      filter(coral_numID %in% i) %>%
      arrange(date)

    availtiss = curr_coral[1,]$starttiss
    tissue = curr_coral[1,]$Tissue

#     #test
#     # j = 15 #date is 2019-04-11, the first day of infection for 3_p47_t5_s0_c12_MCAV
#     j = 15 #date is 2019-04-11, the first day of infection for 3_p45_t9_s5_c27_SBOU
    
    #loop through timepoints of monitoring period
    for(j in 1:nrow(curr_coral)){

      percloss = curr_coral[j,]$percloss #perc. of coral lost to SCTLD (active infected or dead); no distinction made in field
      remaintiss = curr_coral[j,]$remaintiss
      curr_remaintiss = NA
      curr_removedtiss = NA
      
      # recent_inftiss = if (any(!is.na(curr_coral[1:(j-1),]$inftiss))) {
      #   max(na.omit(curr_coral[1:(j-1),]$inftiss))
      # } else {
      #   0
      # }
      
      #set the amount of infected tissue in the prior timepoint
      if (j > 1) { # Skip the first timepoint (j == 1) since it's not relevant
        # If the previous timepoint's 'inftiss' is not NA, use that value
        # Otherwise, set 'recent_inftiss' to 0
        recent_inftiss <- if (!is.na(curr_coral[j-1,]$inftiss)) {
          curr_coral[j-1,]$inftiss
        } else {
          0
        }
      } else {
        # Set 'recent_inftiss' to 0 if j is 1, as no previous timepoint exists
        recent_inftiss <- 0
      }
      
      
    
      # logical tree:
      # is there an active infection [is 'percloss' nonzero / not NA]?
      #   if NO --> move on to next iteration and do not change the 'remaintiss' value. there is nothing to update.
      #   if YES --> is there any value in 'remaintiss'?
      #     if NO --> this is the first observation of infection! set 'remaintiss' of the *CURRENT* timepoint as a subtraction of 'inftiss'
      #       from 'availtiss*tissue'. remaintiss' must be zero if 'percloss' was 100%. otherwise, 'remaintiss' is nonzero.
      #     if YES --> this is an active infection that *continued* after the first-infected timepoint! set 'remaintiss' of the *CURRENT*
      #       timepoint as a subtraction of 'inftiss' from percloss*availtiss*tissue. this accounts for all accumulated recent tissue loss
      
      if(percloss > 0 & !is.na(percloss)){
        #there is an active infection

        if(remaintiss > 0 & !is.na(remaintiss)){
          #this is a *continued* active infection

          if(!is.na(curr_remaintiss) && curr_remaintiss < 0){
            
            #the coral has died. this is the last day of infection
            #   NOTE - since we cannot know exactly which day the coral officially died, we are making the assumption that SCTLD was actively
            #          eating away at tissue until the day a diver arrived at the colony and documented its death. a more nuanced approach
            #          might interpolate a distribution of times of "true" mortality for all colonies in the study, but this would introduce
            #          more assumptions and possibly overcomplicate the study
            curr_remaintiss = 0
            curr_removedtiss = availtiss*tissue - curr_remaintiss
            curr_coral[j:nrow(curr_coral),]$remaintiss = curr_remaintiss
            curr_coral[j:nrow(curr_coral),]$removedtiss = curr_removedtiss
            
          } else{
            #the coral has not yet died
    
            
            # STOPPING POINT
            #   - 15 august 2024
            #   - this part is broken, need to properly account for the final-mortality day. also, second-guessing
            #     how I am treating "death" here - do I actually need to do some backtracking between timepoints
            #     to properly separate compartments? I think so. ugh
            curr_remaintiss = remaintiss - availtiss*tissue*(percloss/100) + curr_coral[j,]$inftiss - recent_inftiss
            curr_removedtiss = availtiss*tissue - curr_remaintiss
            curr_coral[j:nrow(curr_coral),]$remaintiss = curr_remaintiss
            curr_coral[j:nrow(curr_coral),]$removedtiss = curr_removedtiss
            
          }

        } else{
          
          #this is the first observation of infection. set remaining tissue at the current timepoint based on percentage lost.
          #   the current instantaneous infected tissue is added back in, so that there isn't a "doubling up" between the infected
          #   and removed compartments
          curr_remaintiss = availtiss*tissue - availtiss*tissue*(percloss/100) + curr_coral[j,]$inftiss
          curr_removedtiss = availtiss*tissue - curr_remaintiss
          curr_coral[j:nrow(curr_coral),]$remaintiss = curr_remaintiss
          curr_coral[j:nrow(curr_coral),]$removedtiss = curr_removedtiss
          
        }
      } else{
        
        #no active infection. simply accumulate removal of remaining instantaneous infected tissue from prior timepoint (if any)
        curr_remaintiss = remaintiss - recent_inftiss
        curr_removedtiss = availtiss*tissue - curr_remaintiss
        curr_coral[j:nrow(curr_coral),]$remaintiss = curr_remaintiss
        curr_coral[j:nrow(curr_coral),]$removedtiss = curr_removedtiss
        
      }
      
      
      
      
      
      
      
    }
    
    
    
    survey_tissue[which(survey_tissue$coral_numID%in%i), which(colnames(survey_tissue) %in% 'remaintiss')] = curr_coral$remaintiss
    survey_tissue[which(survey_tissue$coral_numID%in%i), which(colnames(survey_tissue) %in% 'removedtiss')] = curr_coral$removedtiss
  }
  
  
  # STOPPING POINT
  # - 14 August 2024. see below:
  #
  # NOTE - need to make sure that I actually correctly coded this for a coral that indeed dies
  # also NOTE - the below important note has now been coded, but I think there is a little bug still nestled in there. in a case
  #               where there is an intermittent infection, the "dangling" 'inftiss' from the last infected timepoint does not#
  #               correctly propagate to the intermittent timepoint, but instead gets passed to the next infected timepoint.
  #               this doesn't make sense for a few reasons: 1.) cases where the coral never gets infected again and survives,
  #               because now there is some tissue that is never accounted for in the 'removed' compartment, and 2.) cases where
  #               a coral indeed does get infected again, but the 'removed' compartment is updated weeks to months too late.
  #
  #               I think the solution may be 'knowing' when there is an infection gap either preceding the current timepoint
  #               when currently on an infected timepoint and adding 'recent_inftiss' - and then *not* adding that 'recent_inftiss'.
  #               The requirement there, is that when there is 'no active infection' in the if/else tree, need to check if there is
  #               a dangling 'inftiss' from the *immediately preceding* timepoint, and go ahead and add it to 'removed' (and thus,
  #               'remaining')
  
  # NOTE IMPORTANT - I think I caught something minor, but important. currently, the removed tissue added at a timepoint includes the
  #         current instantaneous 'inftiss'. this is not correct, since a swath of tissue cannot be both actively infected
  #         and also dead at the same time. so, I need to subtract 'inftiss' of the current timepoint from the removed tissue
  #         at that timepoint. this is a little tricky because I then have to go back to the previous timepoint, at least in
  #         the if/else case where there is an ONGOING infection, and ADD the 'inftiss' from the previous timepoint to the
  #         removed tissue of the current timepoint. very doable, though.
  #
  
  
  
  
  
  #make special case for patient zero Dichocoenia stokesi on 10-30-2018 since loops don't handle its exception. assume its 90%
  # tissue loss occurred over 3 weeks since we don't have info except 3 months prior when colony was healthy
  # note - filter for '2_p27_t2_s0_c1_DSTO'
  # NOTE - may be causing the infectious tissue "spark" to be too high at Offshore site. could project the infection out to like 30 days instead?
  # 011: workshop some ideas to get this initial infection started off in a way that's a little more natural
  survey_tissue[survey_tissue$Coral_ID == '2_p27_t2_s0_c1_DSTO' & survey_tissue$date == '2018-10-30', ]$progdays = 21
  survey_tissue[survey_tissue$Coral_ID == '2_p27_t2_s0_c1_DSTO' & survey_tissue$date == '2018-10-30', ]$percinf = 90/21
  survey_tissue[survey_tissue$Coral_ID == '2_p27_t2_s0_c1_DSTO' & survey_tissue$date == '2018-10-30', ]$inftiss = (90/21/100)*(0.02288368)*(.98)
  survey_tissue[survey_tissue$Coral_ID == '2_p27_t2_s0_c1_DSTO' & survey_tissue$date == '2018-10-30', ]$removedtiss = ((90/21/100)*(0.02288368)*(.98))*21
  
  ### NOTE - code above assumes that tissue loss for 'unknown' reason is due to SCTLD (filter by value = 'Unknown' and look at percloss that is nonzero)
  ###           - there are 17 corals which experienced this, mostly LS/MS corals
  ###           - I checked all of them, and each ended up also having SCTLD in near or immediately adjacent timepoints - fair assumption that all loss was SCTLD
  ###           - Solenastrea might be a bit more susceptible than I'd thought. a few cases of total mortality in 3-4 weeks. also seeing evidence that 
  ###           -   PCLI can sustain long infections; might match up well with susceptibility knowledge from VI
  
  # # NOTE - after bringing in the old mortality information, the below-described coral was actually already dead prior to the start of
  #           confirmed SCTLD monitoring. while it is certainly possible it actually died of SCTLD, since we can't say for sure, it is not
  #           adding any infectious or removed tissue to the SIR model.
  # #another special case for "true" patient zero Colpophyllia natans on 10-30-2018 which had suddenly died by that survey date. assume its 100%
  # # tissue loss occurred over 3 weeks as well
  # # NOTE - may be causing the infectious tissue "spark" to be too high, particularly at Midchannel! look back at this
  # survey_tissue[survey_tissue$Coral_ID == '1_p23_t1_s0_c6_CNAT' & survey_tissue$date == '2018-10-30', ]$progdays = 21
  # survey_tissue[survey_tissue$Coral_ID == '1_p23_t1_s0_c6_CNAT' & survey_tissue$date == '2018-10-30', ]$percinf = 100/21
  # survey_tissue[survey_tissue$Coral_ID == '1_p23_t1_s0_c6_CNAT' & survey_tissue$date == '2018-10-30', ]$inftiss = (100/21/100)*(0.01773432)
  
  ### NOTE - there are a few corals with 'SCTLD' for their health status even after 100% percloss:
  ###         - 3_p47_t3_s0_c8_PSTR
  ###         - 3_p47_t3_s0_c15_CNAT
  ###         - 3_p47_t7_s0_c2_PSTR
  ###         - 1_p23_t1_s0_c5_DSTO
  ###      - May have been a result of data cleaning from the original 2021 analysis, where active 'SCTLD' infection was noted in the
  #           field, then in post-analysis, it was determined the coral was actually already dead (can check on this further if needed)
  #
  # CALCULATING REMOVED / REMAINING TISSUE THROUGH TIME

  

  
  
  
  # STOPPING POINT - 7 may 2024
  #   why not just do the below earlier in the script?
  # STOPPING POINT - 7 may 2024
  
  #add column for old mortality-adjusted starting tissue values. this should have been done earlier in script, but it works fine
  survey_tissue$starttiss = survey_tissue$Tissue*(1-((survey_tissue$old_mortality)/100))
  
  ###sum tissue cover values across colonies into bins (Site, Sus_Cat, Date)
  survey_tissue$date_factor = as.factor(survey_tissue$date)
  levels(survey_tissue$date_factor)
  
  survey_tissue = within(survey_tissue, {
    Timepoint = paste("T", factor(date_factor, labels=seq(unique(date_factor))), sep="")
  })
  
  survey_tissue$fill = survey_tissue$Tissue*(1-survey_tissue$old_mortality/100)
  indices = which(is.na(survey_tissue$remaintiss))
  survey_tissue$remaintiss[is.na(survey_tissue$remaintiss)] = survey_tissue$fill[indices]
  survey_tissue = survey_tissue[,-which(colnames(survey_tissue) %in% 'fill')]
  
  #clean up dataframe
  survey_tissue$Site_type = as.factor(survey_tissue$Site_type)
  survey_tissue$Timepoint = as.factor(survey_tissue$Timepoint)
  survey_tissue = subset(survey_tissue, select = -date_factor)
  
  #remove the SCTLD-unaffected corals from the simulation entirely. they should not be included in the infection "pool"
  survey_tissue = survey_tissue[!grepl('Unaffected', survey_tissue$Sus_Cat),]
  
  #summary table for infected tissue carpet at plot level
  cols = c("Site_type", "date", "Timepoint", "Plot", "Sus_Cat")
  tissue.summary.plots = survey_tissue %>%
    group_by(across(all_of(cols))) %>%
    summarize(plot.sustiss = sum(remaintiss, na.rm = TRUE),
              plot.inftiss = sum(inftiss, na.rm = TRUE),
              plot.remtiss = sum(removedtiss, na.rm = TRUE),
              num.inf = across(starts_with("inf"), ~ sum(.x != 0, na.rm = TRUE)),
              census = n()
    )
  
  #silly column naming procedure because of R quirk
  temp = as.vector(tissue.summary.plots$num.inf)
  temp = temp[[1]]
  tissue.summary.plots = subset(tissue.summary.plots, select = -num.inf)
  tissue.summary.plots$num.inf = temp
  
  # NOTE - including surprise-dead unaffected corals (there were 2 OCUL & 1 PAST as infected, but not including unaffected corals in
  #         census count since so few of them ever "died" of SCTLD
  cols = c("Site_type", "date", "Timepoint")
  tissue.summary = tissue.summary.plots %>%
    group_by(across(all_of(cols))) %>%
    summarize(
      
      tot.sustiss = sum(plot.sustiss[Sus_Cat=="LS"]) + sum(plot.sustiss[Sus_Cat=="MS"]) + sum(plot.sustiss[Sus_Cat=="HS"]), #excluding unaffected
      LS.sustiss = sum(plot.sustiss[Sus_Cat=="LS"]),
      MS.sustiss = sum(plot.sustiss[Sus_Cat=="MS"]),
      HS.sustiss = sum(plot.sustiss[Sus_Cat=="HS"]),

      tot.inftiss = sum(plot.inftiss[Sus_Cat=="LS"]) + sum(plot.inftiss[Sus_Cat=="MS"]) + sum(plot.inftiss[Sus_Cat=="HS"]), #excluding (very few) infected "unaffected"
      LS.inftiss = sum(plot.inftiss[Sus_Cat=="LS"]),
      MS.inftiss = sum(plot.inftiss[Sus_Cat=="MS"]),
      HS.inftiss = sum(plot.inftiss[Sus_Cat=="HS"]),

      tot.remtiss = sum(plot.remtiss[Sus_Cat=="LS"]) + sum(plot.remtiss[Sus_Cat=="MS"]) + sum(plot.remtiss[Sus_Cat=="HS"]),
      LS.remtiss = sum(plot.remtiss[Sus_Cat=="LS"]),
      MS.remtiss = sum(plot.remtiss[Sus_Cat=="MS"]),
      HS.remtiss = sum(plot.remtiss[Sus_Cat=="HS"]),

      tot.susnum = sum(census[Sus_Cat=="LS"]) + sum(census[Sus_Cat=="MS"]) + sum(census[Sus_Cat=="HS"]), #excluding unaffected
      LS.susnum = sum(census[Sus_Cat=="LS"]),
      MS.susnum = sum(census[Sus_Cat=="MS"]),
      HS.susnum = sum(census[Sus_Cat=="HS"]),

      tot.infnum = sum(num.inf[Sus_Cat=="LS"]) + sum(num.inf[Sus_Cat=="MS"]) + sum(num.inf[Sus_Cat=="HS"]), #excluding (very few) infected "unaffected"
      LS.infnum = sum(num.inf[Sus_Cat=="LS"]),
      MS.infnum = sum(num.inf[Sus_Cat=="MS"]),
      HS.infnum = sum(num.inf[Sus_Cat=="HS"]),
    )
  
  #incorporate days since first infection for each site
  sites = levels(tissue.summary$Site_type)
  tissue.summary['Day'] = NA
  tissue.summary$Day[is.na(tissue.summary$Day)] = 0
  list(time_list <- vector("list", length = length(sites)))
  for(i in 1:length(sites)){
    site = sites[i]
    # site = 'Midchannel' #for testing
    
    curr.timepoint = ''
    if(site == 'Nearshore'){
      curr.timepoint = 'T12' #timepoint with first documented infection  [site specific]
      days.backtrack = 21 #14 #first infections were a 120-cm PSTR (10% loss), 105-cm PSTR (5%), and 12-cm OANN (5%)
    } else if(site == "Midchannel"){
      curr.timepoint = 'T8'
      days.backtrack = 21 #7 #first infection was a 31-cm DSTO, 10% loss
    } else{ #Offshore
      curr.timepoint = 'T6'
      days.backtrack = 21 #28 #10 #first infection was a 13-cm DSTO, 90% loss
    }
    
    currdate = as.POSIXct(min(tissue.summary$date[tissue.summary$Site_type==site & tissue.summary$Timepoint==curr.timepoint]))
    date.backtrack = currdate - days(days.backtrack)
    
    uniquedates = sort(unique(tissue.summary$date[tissue.summary$Site_type==site & tissue.summary$Timepoint!="T0"]))
    enddate = max(uniquedates)
    numdays = as.numeric(difftime(enddate, date.backtrack))
    # time = seq(0, numdays, by = 1) #days of simulated epidemic, to match empirical observations
    
    for (j in 1:length(uniquedates)){
      timediff = as.numeric(uniquedates[j]-date.backtrack)
      tissue.summary[which(tissue.summary$date%in%uniquedates & tissue.summary$Site_type==site),][j,]$Day = timediff
    }
    
    tissue.summary$Day = round(tissue.summary$Day)
    tissue.summary$Day[tissue.summary$Day < 2] = NA
    
    time_list[[i]] = seq(0, numdays, by = 1) #days of simulated epidemic, to match empirical observations
  }
  
  #most readable table
  totals.only.summary = tissue.summary %>% select(Site_type, date, Timepoint, tot.sustiss, tot.inftiss, tot.remtiss, tot.susnum, tot.infnum)
  totals.only.summary$prevtiss = (totals.only.summary$tot.inftiss/totals.only.summary$tot.sustiss)*100
  totals.only.summary$prevnum = (totals.only.summary$tot.infnum/totals.only.summary$tot.susnum)*100
  totals.only.summary$sustiss.host.ratio = totals.only.summary$tot.sustiss/totals.only.summary$tot.susnum
  totals.only.summary$inftiss.host.ratio = totals.only.summary$tot.inftiss/totals.only.summary$tot.infnum
  
  ### NOTE - also consider the Healthy / Unknown (at last time point) corals - did any of these actually die? (64 of these)
  ###           - I think safe to assume none of them died or suffered tissue loss. just weird lookin' corals
  ### NOTE - and the Dis / Unknown - did any of these actually die? (13 of these)
  ###           - they did not! these are all corals which were infected with SCTLD and experienced lesion cessation after prolonged (month+) infection. very cool that the data shows this - they are all fairly large, LS and MS corals too
  ### NOTE - also check out Dis / Healthy (11 of these)
  ###           - similar story here! in fact, we may want to consider the SCTLD-infected corals with "0" tissue loss to at least have minimum 1% loss somehow - to mimic long-lasting but dormant infections
  ### NOTE - and Dis / SCTLD (37 of these)
  ###           - these are all active infections at the end of sampling (sometimes initiated near the very end of sampling, sometimes simply extremely prolonged for 6+ months)
  ###           - almost exclusively large LS/MS corals. probably the last survivors & a second wave. a few of them had long sustained SCLTD but with 0% tissue loss - could force these to be "active". a couple examples of HS corals with sustained infections - hardy genotype etc.? go back and sample them??
  ### NOTE - finally, Dis / Dead (115 of these)
  ###           - vast majority of these are DSTO, PSTR, & CNAT which die in just a matter of weeks - interesting how consistent this is.
  
  ### plots ###
  
  #Offshore
  site = 'Offshore'
  prev.timepoint = 'T5'
  sample.day.list = pull(subset(tissue.summary, Site_type==site)[, "tot.infnum"])
  indices.offshore = which(sample.day.list > 0) #record indices for timepoints including and after the first infection day
  indices.offshore = full_seq(indices.offshore, 1) #ensure that if #/infected dips briefly to 0, that  timepoint is still considered for fitting
  days = tissue.summary %>% na.omit() %>% filter(Site_type == site) %>% pull(Day) 
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  
  #could edit this to include dates
  obs.sus.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.sustiss"])[indices.offshore], Category = "LS", Compartment = "Susceptible", Site = site)
  obs.sus.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.sustiss"])[indices.offshore], Category = "MS", Compartment = "Susceptible", Site = site)
  obs.sus.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.sustiss"])[indices.offshore], Category = "HS", Compartment = "Susceptible", Site = site)
  obs.sus.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.sustiss"])[indices.offshore], Compartment = "Susceptible", Site = site)
  
  obs.inf.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.inftiss"])[indices.offshore], Category = "LS", Compartment = "Infected", Site = site)
  obs.inf.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.inftiss"])[indices.offshore], Category = "MS", Compartment = "Infected", Site = site)
  obs.inf.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.inftiss"])[indices.offshore], Category = "HS", Compartment = "Infected", Site = site)
  obs.inf.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.inftiss"])[indices.offshore], Compartment = "Infected", Site = site)
  
  obs.rem.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.remtiss"])[indices.offshore], Category = "LS", Compartment = "Dead", Site = site)
  obs.rem.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.remtiss"])[indices.offshore], Category = "MS", Compartment = "Dead", Site = site)
  obs.rem.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.remtiss"])[indices.offshore], Category = "HS", Compartment = "Dead", Site = site)
  obs.rem.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.remtiss"])[indices.offshore], Compartment = "Dead", Site = site)
  
  obs.offshore = rbind(obs.sus.LS.tiss, obs.sus.MS.tiss, obs.sus.HS.tiss,
               obs.inf.LS.tiss, obs.inf.MS.tiss, obs.inf.HS.tiss,
               obs.rem.LS.tiss, obs.rem.MS.tiss, obs.rem.HS.tiss)
  
  obs.offshore.basic = rbind(obs.sus.tiss, obs.inf.tiss, obs.rem.tiss)
  
  obs.sus.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.sustiss"]/S0.snapshot$LS.sustiss)[indices.offshore], Category = "LS", Compartment = "Susceptible", Site = site)
  obs.sus.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.sustiss"]/S0.snapshot$MS.sustiss)[indices.offshore], Category = "MS", Compartment = "Susceptible", Site = site)
  obs.sus.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.sustiss"]/S0.snapshot$HS.sustiss)[indices.offshore], Category = "HS", Compartment = "Susceptible", Site = site)
  obs.sus.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.sustiss"]/S0.snapshot$tot.sustiss)[indices.offshore], Compartment = "Susceptible", Site = site)
  
  obs.inf.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.inftiss"]/S0.snapshot$LS.sustiss)[indices.offshore], Category = "LS", Compartment = "Infected", Site = site)
  obs.inf.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.inftiss"]/S0.snapshot$MS.sustiss)[indices.offshore], Category = "MS", Compartment = "Infected", Site = site)
  obs.inf.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.inftiss"]/S0.snapshot$HS.sustiss)[indices.offshore], Category = "HS", Compartment = "Infected", Site = site)
  obs.inf.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.inftiss"]/S0.snapshot$tot.sustiss)[indices.offshore], Compartment = "Infected", Site = site)
  
  obs.rem.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.remtiss"]/S0.snapshot$LS.sustiss)[indices.offshore], Category = "LS", Compartment = "Dead", Site = site)
  obs.rem.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.remtiss"]/S0.snapshot$MS.sustiss)[indices.offshore], Category = "MS", Compartment = "Dead", Site = site)
  obs.rem.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.remtiss"]/S0.snapshot$HS.sustiss)[indices.offshore], Category = "HS", Compartment = "Dead", Site = site)
  obs.rem.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.remtiss"]/S0.snapshot$tot.sustiss)[indices.offshore], Compartment = "Dead", Site = site)
  
  obs.offshore.scale = rbind(obs.sus.LS.tiss.scale, obs.sus.MS.tiss.scale, obs.sus.HS.tiss.scale,
                             obs.inf.LS.tiss.scale, obs.inf.MS.tiss.scale, obs.inf.HS.tiss.scale,
                             obs.rem.LS.tiss.scale, obs.rem.MS.tiss.scale, obs.rem.HS.tiss.scale)
  
  obs.offshore.scale.basic = rbind(obs.sus.tiss.scale, obs.inf.tiss.scale, obs.rem.tiss.scale)
  
  #Midchannel
  site = 'Midchannel'
  prev.timepoint = 'T7'
  sample.day.list = pull(subset(tissue.summary, Site_type==site)[, "tot.infnum"])
  indices.midchannel = which(sample.day.list > 0) #record indices.midchannel for timepoints including and after the first infection day
  indices.midchannel = full_seq(indices.midchannel, 1) #ensure that if #/infected dips briefly to 0, that  timepoint is still considered for fitting
  days = tissue.summary %>% na.omit() %>% filter(Site_type == site) %>% pull(Day) 
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)

  obs.sus.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.sustiss"])[indices.midchannel], Category = "LS", Compartment = "Susceptible", Site = site)
  obs.sus.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.sustiss"])[indices.midchannel], Category = "MS", Compartment = "Susceptible", Site = site)
  obs.sus.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.sustiss"])[indices.midchannel], Category = "HS", Compartment = "Susceptible", Site = site)
  obs.sus.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.sustiss"])[indices.midchannel], Compartment = "Susceptible", Site = site)
  
  obs.inf.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.inftiss"])[indices.midchannel], Category = "LS", Compartment = "Infected", Site = site)
  obs.inf.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.inftiss"])[indices.midchannel], Category = "MS", Compartment = "Infected", Site = site)
  obs.inf.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.inftiss"])[indices.midchannel], Category = "HS", Compartment = "Infected", Site = site)
  obs.inf.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.inftiss"])[indices.midchannel], Compartment = "Infected", Site = site)
  
  obs.rem.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.remtiss"])[indices.midchannel], Category = "LS", Compartment = "Dead", Site = site)
  obs.rem.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.remtiss"])[indices.midchannel], Category = "MS", Compartment = "Dead", Site = site)
  obs.rem.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.remtiss"])[indices.midchannel], Category = "HS", Compartment = "Dead", Site = site)
  obs.rem.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.remtiss"])[indices.midchannel], Compartment = "Dead", Site = site)
  
  obs.midchannel = rbind(obs.sus.LS.tiss, obs.sus.MS.tiss, obs.sus.HS.tiss,
                       obs.inf.LS.tiss, obs.inf.MS.tiss, obs.inf.HS.tiss,
                       obs.rem.LS.tiss, obs.rem.MS.tiss, obs.rem.HS.tiss)
  
  obs.midchannel.basic = rbind(obs.sus.tiss, obs.inf.tiss, obs.rem.tiss)
  
  obs.sus.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.sustiss"]/S0.snapshot$LS.sustiss)[indices.midchannel], Category = "LS", Compartment = "Susceptible", Site = site)
  obs.sus.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.sustiss"]/S0.snapshot$MS.sustiss)[indices.midchannel], Category = "MS", Compartment = "Susceptible", Site = site)
  obs.sus.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.sustiss"]/S0.snapshot$HS.sustiss)[indices.midchannel], Category = "HS", Compartment = "Susceptible", Site = site)
  obs.sus.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.sustiss"]/S0.snapshot$tot.sustiss)[indices.midchannel], Compartment = "Susceptible", Site = site)
  
  obs.inf.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.inftiss"]/S0.snapshot$LS.sustiss)[indices.midchannel], Category = "LS", Compartment = "Infected", Site = site)
  obs.inf.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.inftiss"]/S0.snapshot$MS.sustiss)[indices.midchannel], Category = "MS", Compartment = "Infected", Site = site)
  obs.inf.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.inftiss"]/S0.snapshot$HS.sustiss)[indices.midchannel], Category = "HS", Compartment = "Infected", Site = site)
  obs.inf.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.inftiss"]/S0.snapshot$tot.sustiss)[indices.midchannel], Compartment = "Infected", Site = site)
  
  obs.rem.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.remtiss"]/S0.snapshot$LS.sustiss)[indices.midchannel], Category = "LS", Compartment = "Dead", Site = site)
  obs.rem.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.remtiss"]/S0.snapshot$MS.sustiss)[indices.midchannel], Category = "MS", Compartment = "Dead", Site = site)
  obs.rem.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.remtiss"]/S0.snapshot$HS.sustiss)[indices.midchannel], Category = "HS", Compartment = "Dead", Site = site)
  obs.rem.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.remtiss"]/S0.snapshot$tot.sustiss)[indices.midchannel], Compartment = "Dead", Site = site)
  
  obs.midchannel.scale = rbind(obs.sus.LS.tiss.scale, obs.sus.MS.tiss.scale, obs.sus.HS.tiss.scale,
                             obs.inf.LS.tiss.scale, obs.inf.MS.tiss.scale, obs.inf.HS.tiss.scale,
                             obs.rem.LS.tiss.scale, obs.rem.MS.tiss.scale, obs.rem.HS.tiss.scale)
  
  obs.midchannel.scale.basic = rbind(obs.sus.tiss.scale, obs.inf.tiss.scale, obs.rem.tiss.scale)
  
  #Nearshore
  site = 'Nearshore'
  prev.timepoint = 'T11'
  sample.day.list = pull(subset(tissue.summary, Site_type==site)[, "tot.infnum"])
  indices.nearshore = which(sample.day.list > 0) #record indices.nearshore for timepoints including and after the first infection day
  indices.nearshore = full_seq(indices.nearshore, 1) #ensure that if #/infected dips briefly to 0, that  timepoint is still considered for fitting
  days = tissue.summary %>% na.omit() %>% filter(Site_type == site) %>% pull(Day) 
  S0.snapshot = subset(tissue.summary, Site_type==site & Timepoint == prev.timepoint)
  
  obs.sus.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.sustiss"])[indices.nearshore], Category = "LS", Compartment = "Susceptible", Site = site)
  obs.sus.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.sustiss"])[indices.nearshore], Category = "MS", Compartment = "Susceptible", Site = site)
  obs.sus.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.sustiss"])[indices.nearshore], Category = "HS", Compartment = "Susceptible", Site = site)
  obs.sus.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.sustiss"])[indices.nearshore], Compartment = "Susceptible", Site = site)
  
  obs.inf.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.inftiss"])[indices.nearshore], Category = "LS", Compartment = "Infected", Site = site)
  obs.inf.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.inftiss"])[indices.nearshore], Category = "MS", Compartment = "Infected", Site = site)
  obs.inf.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.inftiss"])[indices.nearshore], Category = "HS", Compartment = "Infected", Site = site)
  obs.inf.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.inftiss"])[indices.nearshore], Compartment = "Infected", Site = site)
  
  obs.rem.LS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.remtiss"])[indices.nearshore], Category = "LS", Compartment = "Dead", Site = site)
  obs.rem.MS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.remtiss"])[indices.nearshore], Category = "MS", Compartment = "Dead", Site = site)
  obs.rem.HS.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.remtiss"])[indices.nearshore], Category = "HS", Compartment = "Dead", Site = site)
  obs.rem.tiss = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.remtiss"])[indices.nearshore], Compartment = "Dead", Site = site)
  
  obs.nearshore = rbind(obs.sus.LS.tiss, obs.sus.MS.tiss, obs.sus.HS.tiss,
                         obs.inf.LS.tiss, obs.inf.MS.tiss, obs.inf.HS.tiss,
                         obs.rem.LS.tiss, obs.rem.MS.tiss, obs.rem.HS.tiss)
  
  obs.nearshore.basic = rbind(obs.sus.tiss, obs.inf.tiss, obs.rem.tiss)
  
  obs.sus.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.sustiss"]/S0.snapshot$LS.sustiss)[indices.nearshore], Category = "LS", Compartment = "Susceptible", Site = site)
  obs.sus.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.sustiss"]/S0.snapshot$MS.sustiss)[indices.nearshore], Category = "MS", Compartment = "Susceptible", Site = site)
  obs.sus.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.sustiss"]/S0.snapshot$HS.sustiss)[indices.nearshore], Category = "HS", Compartment = "Susceptible", Site = site)
  obs.sus.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.sustiss"]/S0.snapshot$tot.sustiss)[indices.nearshore], Compartment = "Susceptible", Site = site)
  
  obs.inf.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.inftiss"]/S0.snapshot$LS.sustiss)[indices.nearshore], Category = "LS", Compartment = "Infected", Site = site)
  obs.inf.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.inftiss"]/S0.snapshot$MS.sustiss)[indices.nearshore], Category = "MS", Compartment = "Infected", Site = site)
  obs.inf.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.inftiss"]/S0.snapshot$HS.sustiss)[indices.nearshore], Category = "HS", Compartment = "Infected", Site = site)
  obs.inf.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.inftiss"]/S0.snapshot$tot.sustiss)[indices.nearshore], Compartment = "Infected", Site = site)
  
  obs.rem.LS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "LS.remtiss"]/S0.snapshot$LS.sustiss)[indices.nearshore], Category = "LS", Compartment = "Dead", Site = site)
  obs.rem.MS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "MS.remtiss"]/S0.snapshot$MS.sustiss)[indices.nearshore], Category = "MS", Compartment = "Dead", Site = site)
  obs.rem.HS.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "HS.remtiss"]/S0.snapshot$HS.sustiss)[indices.nearshore], Category = "HS", Compartment = "Dead", Site = site)
  obs.rem.tiss.scale = data.frame(days, prop = pull(subset(tissue.summary, Site_type==site)[, "tot.remtiss"]/S0.snapshot$tot.sustiss)[indices.nearshore], Compartment = "Dead", Site = site)
  
  obs.nearshore.scale = rbind(obs.sus.LS.tiss.scale, obs.sus.MS.tiss.scale, obs.sus.HS.tiss.scale,
                               obs.inf.LS.tiss.scale, obs.inf.MS.tiss.scale, obs.inf.HS.tiss.scale,
                               obs.rem.LS.tiss.scale, obs.rem.MS.tiss.scale, obs.rem.HS.tiss.scale)
  
  obs.nearshore.scale.basic = rbind(obs.sus.tiss.scale, obs.inf.tiss.scale, obs.rem.tiss.scale)
  
  obs = rbind(obs.offshore, obs.midchannel, obs.nearshore)
  obs$Category = as.factor(obs$Category)
  obs$Compartment = as.factor(obs$Compartment)
  obs$Site = as.factor(obs$Site)
  
  obs.basic = rbind(obs.offshore.basic, obs.midchannel.basic, obs.nearshore.basic)
  obs.basic$Compartment = as.factor(obs.basic$Compartment)
  obs.basic$Site = as.factor(obs.basic$Site)
  
  obs.scale = rbind(obs.offshore.scale, obs.midchannel.scale, obs.nearshore.scale)
  obs.scale$Category = as.factor(obs.scale$Category)
  obs.scale$Compartment = as.factor(obs.scale$Compartment)
  obs.scale$Site = as.factor(obs.scale$Site)
  
  obs.scale.basic = rbind(obs.offshore.scale.basic, obs.midchannel.scale.basic, obs.nearshore.scale.basic)
  obs.scale.basic$Compartment = as.factor(obs.scale.basic$Compartment)
  obs.scale.basic$Site = as.factor(obs.scale.basic$Site)
  
  #arbitrary value - took lowest proportional 'firstday' infection, a 31-cm Midchannel DSTO w/ 10% loss, divided by 3. seemed to work well
  polyp_SA.midchannel = 1.67e-05 #this one consistently looks a little high if I leave it the same as the other two
  polyp_SA.nearshore = 3e-05  # 5e-05 
  polyp_SA.offshore = 3e-05 # 5e-05

  # MULTI-GROUP SIR
  # scaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.scale = ggplot(data = obs.scale %>% filter(Site == "Offshore"), aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.offshore.scale = ggplot(data = obs.scale %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.scale = ggplot(data = obs.scale %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()

  #Midchannel
  p.SIR.midchannel.scale = ggplot(data = obs.scale %>% filter(Site == "Midchannel"), aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.scale = ggplot(data = obs.scale %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel.scale = ggplot(data = obs.scale %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.scale = ggplot(data = obs.scale %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.scale = ggplot(data = obs.scale %>% filter(Site == "Nearshore"), aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.nearshore.scale = ggplot(data = obs.scale %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.scale = ggplot(data = obs.scale %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # unscaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore = ggplot(data = obs %>% filter(Site == "Offshore"), aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.offshore = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.offshore = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore = ggplot(data = obs %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel"), aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel = ggplot(data = obs %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore"), aes(days, prop, colour = Compartment, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore = ggplot(data = obs %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, prop, linetype = Category)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # BASIC SIR
  # scaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Offshore"), aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.offshore.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Midchannel"), aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.midchannel.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Nearshore"), aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.I.nearshore.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.scale.basic = ggplot(data = obs.scale.basic %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Proportion of tissue") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # unscaled plots
  #Offshore
  # display.brewer.all(colorblindFriendly = TRUE)
  p.SIR.offshore.basic = ggplot(data = obs.basic %>% filter(Site == "Offshore"), aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.offshore.basic = ggplot(data = obs.basic %>% filter(Site == "Offshore", Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.offshore.basic = ggplot(data = obs.basic %>% filter(Site == "Offshore", Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.offshore.basic = ggplot(data = obs.basic %>% filter(Site == "Offshore", Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Offshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Midchannel
  p.SIR.midchannel.basic = ggplot(data = obs.basic %>% filter(Site == "Midchannel"), aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.midchannel.basic = ggplot(data = obs.basic %>% filter(Site == "Midchannel", Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.midchannel.basic = ggplot(data = obs.basic %>% filter(Site == "Midchannel", Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.midchannel.basic = ggplot(data = obs.basic %>% filter(Site == "Midchannel", Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Midchannel'), collapse="")) +
    geom_line() +
    theme_classic()
  
  #Nearshore
  p.SIR.nearshore.basic = ggplot(data = obs.basic %>% filter(Site == "Nearshore"), aes(days, prop, colour = Compartment)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
    theme_classic()
  
  p.S.nearshore.basic = ggplot(data = obs.basic %>% filter(Site == "Nearshore", Compartment == "Susceptible"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.I.nearshore.basic = ggplot(data = obs.basic %>% filter(Site == "Nearshore", Compartment == "Infected"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  p.R.nearshore.basic = ggplot(data = obs.basic %>% filter(Site == "Nearshore", Compartment == "Dead"), aes(days, prop)) +
    xlab("Day of observation period") +
    ylab("Tissue Surface Area (m2)") +
    ggtitle(paste(c("", 'Nearshore'), collapse="")) +
    geom_line() +
    theme_classic()
  
  # # Save workspace
  # save.image(file = "FLKEYS_workspace.RData")
  