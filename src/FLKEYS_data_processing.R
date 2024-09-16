
  # .rs.restartR(clean = TRUE)
  rm(list=ls())

  library(tidyverse)
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
  
  #   PATIENT ZERO CORALS
  #   - Offshore: 2_p27_t2_s0_c1_DSTO, 10-30-2018, 90% loss. backtrack a week? this is only a 5 inch coral
  #   - Midchannel: 1_p25_t2_s0_c22_DSTO, 11-29-2018, 10% loss
  #   - Nearshore: 3_p47_t3_s0_c8_PSTR (10% loss) and 3_p47_t4_s5_c15_PSTR (5% loss), 2019-02-08
  
  #pull wide-format temporal data from lower Florida Keys (Williams et al. 2021), which recorded individual coral colonies; categorize species
  #   by susceptibility group
  survey = survey %>%
    rename(site = Site, spp = Sps) %>%
    mutate(
      site = case_when(
        site == 1 ~ "mid", #Site 1, transects 23 and 25 is Wonderland site (midchannel)
        site == 2 ~ "off", #Site 2, transects 27 and 28 is ACER site (offshore)
        site == 3 ~ "near" #Site 3, transects 45 and 47 is N. Birthday site (nearshore)
      ),
      spp = as.factor(spp),
      susc = case_when(
        spp %in% c('PCLI', 'SINT', 'SSID') ~ 'LS', #lowest rates of progression - however, OANN/OFAB have short disease onset time
        spp %in% c('OANN', 'OFAV', 'MCAV', 'SBOU') ~ 'MS', #moderate rates of progression - but MCAV slow disease onset! SBOU quick onset
        spp %in% c('CNAT', 'DLAB', 'DSTO', 'MMEA', 'MYCE', 'PSTR') ~ 'HS', #quickest rates of  progression and shortest disease onset times
        spp %in% c('AAGA', 'ACER', 'OCUL', 'ODIF', 'PAST', 'PDIV', 'PPOR', 'SRAD') ~ 'Unaffected'
      )
    ) %>%
    select(-coords_x, -coords_y, -total, -total_bin, -days_dis, -weeks_dis) # NOTE - could further analyze days/weeks diseased for stats

  #pull middle Florida Keys (Sharp et al. 2020) dataset, which had similar data as 'survey' but measured colonies in 3 dimensions
  #   rather than one. use this dataset to create a statistical association between maximum diameter & colony surface area (SA) by assuming
  #   a hemi-ellipsoid (Holstein et al. 2015), and then predicting SA in 'survey' dataset using maximum diameter alone
  # Sharp2020 = read.csv(here("data", "Sharp_Maxwell_2020.csv")) %>%
  #   select(Plot, Transect, Coral, Species, Width, Width_2, Height)
  # Sharp2020[Sharp2020==-99] = NA
  # Sharp2020 = Sharp2020[complete.cases(Sharp2020),]
  
  Sharp2020 <- read.csv(here("data", "Sharp_Maxwell_2020.csv")) %>%
    select(Plot, Transect, Coral, Species, Width, Width_2, Height) %>%
    mutate(across(where(is.numeric), ~na_if(., -99))) %>% #convert -99 fill values to NA
    drop_na() #remove NAs from dataset
  
  #set '0' cm heights as an arbitrarily small amount (cm) to allow the hemi-ellipsoid estimation to function correctly. coral recruits 
  #   have very little tangible height, and were recorded underwater as 0 cm height
  Sharp2020 = Sharp2020 %>%
    mutate(Height = ifelse(Height == 0, 0.01, Height)) %>%
    distinct() #also filter out repeats (not interested in a time series here)
  
  #hemi-ellipsoid estimation. p is a dimensionless constant; all else in square cm
  Sharp2020 <- Sharp2020 %>%
    # Calculate the surface area
    mutate(
      a = Height,
      b = Width_2 / 2,
      c = Width / 2,
      p = 1.6075,
      SA = 2 * pi * (((a * b)^p + (a * c)^p + (b * c)^p) / 3)^(1 / p),
      SA = SA / 10000  # Convert from square cm to square meters
    )
  
  # Construct dataframe of estimated surface area and measured maximum diameter ('Width' in this dataset)
  #   NOTE - ignoring the mortality data from Sharp study, as the attempt is to predict dimensions of the original coral,
  #          including any areas of bare skeleton
  SA <- Sharp2020 %>%
    select(x = Width, y = SA) 
  # plot(SA$x, SA$y)
  
  #a measure to prevent over-predicting surface area of the largest corals from small sample size [decided as unnecessary]
  SA = SA %>% filter(x < 122) #145
  
  #STOPPING POINT - 16 Sep 2024
  #   - attempts at fancy model comparisons below. one fitting issue may be differently scaled x and y, and also repeated x values
  #   - consider how well the GAMs are doing at predicting "correct" SAs - cross-reference between Sharp and Williams datasets
  #   - values dipping below zero for smallest corals is a problem for GAMs here - see below as well
  SA.linear <- lm(y ~ x + 0, data = SA) #trying to make sure predictions don't go below zero - struggling with this
  SA.GAM <- gam(y ~ s(x), data = SA, family = gaussian, method = "REML")
  SA.GAM_identity = gam(y ~ s(x, bs = "cr"), data = SA, family = gaussian(link = "identity")) #GAM model with a zero intercept constraint
  SA.GAM_linear <- gam(y ~ s(x) + x, data = SA, family = gaussian, method = "REML")
  SA.GAM_cr <- gam(y ~ s(x, bs = "cr", k = 20), data = SA, family = gaussian, method = "REML", gamma = 1.4)
  SA.GAM_gamma <- gam(y ~ s(x), data = SA, family = Gamma(link = "log"), method = "REML")
  SA$x_scaled <- scale(SA$x)  # rescale x
  SA.GAM_scaled <- gam(y ~ s(x_scaled, bs = "cr", k = 20), data = SA, family = gaussian, method = "REML")
  
  # #diagnostics
  # summary(SA.GAM)
  # draw(SA.GAM, select = 1)
  # gam.check(SA.GAM)
  # appraise(SA.GAM)
  
  # Generate predictions from each model
  SA_predictions <- SA %>%
    mutate(
      pred_linear = predict(SA.linear),
      pred_gam = predict(SA.GAM),
      pred_gam_identity = predict(SA.GAM_identity, newdata = SA),
      pred_gam_linear = predict(SA.GAM_linear),
      pred_gam_cr = predict(SA.GAM_cr),
      pred_gam_gamma = predict(SA.GAM_gamma, type = "response"),
      pred_gam_scaled = predict(SA.GAM_scaled, newdata = SA)
    )
  
  ggplot(SA_predictions, aes(x = x, y = y)) +
    geom_point(color = "black", size = 2, alpha = 0.5, show.legend = FALSE) + # Points are black
    # geom_line(aes(y = pred_linear, color = "Linear"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam, color = "GAM"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_identity, color = "GAM with Zero Intercept"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_linear, color = "GAM with Linear Term"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_cr, color = "GAM with CR Basis"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_gamma, color = "GAM with Gamma"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_scaled, color = "GAM with Scaled x"), linewidth = 1, alpha = 0.7) +
    labs(
      x = "Max diameter (cm)",
      y = "Surface area (m2)",
      title = "Model Fits for Surface Area vs. Max Diameter",
      color = "Model"
    ) +
    scale_color_manual(
      values = c(
        "Linear" = "blue",
        "GAM" = "red",
        "GAM with Zero Intercept" = "pink",
        "GAM with Linear Term" = "green",
        "GAM with CR Basis" = "purple",
        "GAM with Gamma" = "orange",
        "GAM with Scaled x" = "brown"
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  #Akaike comparison
  AIC(SA.linear, SA.GAM, SA.GAM_identity, SA.GAM_linear, SA.GAM_cr, SA.GAM_gamma, SA.GAM_scaled)
  
  # # just linear prediction-fits
  # ggplot(SA_predictions, aes(x = x, y = y)) +
  #   geom_point() +
  #   geom_line(aes(y = pred_linear), color = "red") +
  #   labs(x = "Max diameter (cm)", y = "Surface area (m2)", title = "Linear Model Fit") +
  #   theme_minimal()
  # 
  # # just identity GAM prediction-fits
  # ggplot(SA_predictions, aes(x = x, y = y)) +
  #   geom_point() +
  #   geom_line(aes(y = pred_gam_identity), color = "red") +
  #   labs(x = "Max diameter (cm)", y = "Surface area (m2)", title = "GAM Model Fit") +
  #   theme_minimal()
  # 
  # # just GAM prediction-fits
  # ggplot(SA_predictions, aes(x = x, y = y)) +
  #   geom_point() +
  #   geom_line(aes(y = pred_gam), color = "red") +
  #   labs(x = "Max diameter (cm)", y = "Surface area (m2)", title = "GAM Model Fit") +
  #   theme_minimal()
  # 
  # #all prediction-fits
  # plot(SA.GAM, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  # plot(SA.GAM_identity, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  # plot(SA.GAM_linear, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  # plot(SA.GAM_cr, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  # plot(SA.GAM_gamma, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  # plot(SA.GAM_scaled, scale = 0, all.terms = TRUE, shade = T, shade.col="lightpink", xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  
  # #an attempt to correct for GAM dipping below zero. not perfect
  # x = survey$Max_width
  # x.GAM = as.data.frame(x)
  # x.GAM$x = as.numeric(x.GAM$x)
  # y = predict(SA.GAM, newdata = x.GAM, type = "response")
  # y[y<0] = 0.00005 #ensure GAM-predicted surface area does not dip below zero for small recruits. arbitrary value
  # survey$tissue = y
  # plot(x,y, xlim = c(0, 460), xlab = "Max diameter (cm)", ylab = "Surface area (m2)")
  # sum(survey$tissue) #quick summary to assess overall tissue SA across all sites

  # Generate predictions for the survey data
  survey_predictions <- survey %>%
    mutate(
      pred_gam = predict(SA.GAM, newdata = data.frame(x = Max_width)),
      pred_gam_identity = predict(SA.GAM_identity, newdata = data.frame(x = Max_width)),
      pred_gam_linear = predict(SA.GAM_linear, newdata = data.frame(x = Max_width)),
      pred_gam_cr = predict(SA.GAM_cr, newdata = data.frame(x = Max_width)),
      pred_gam_gamma = predict(SA.GAM_gamma, newdata = data.frame(x = Max_width)),
      pred_gam_scaled = predict(SA.GAM_scaled, newdata = data.frame(x_scaled = scale(Max_width)))
    )
  
  # Plot predictions with ggplot
  ggplot(survey_predictions, aes(x = Max_width)) +
    
    # Lines for each model
    # geom_line(aes(y = pred_gam, color = "GAM"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_identity, color = "GAM with Zero Intercept"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_linear, color = "GAM with Linear Term"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_cr, color = "GAM with CR Basis"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_gamma, color = "GAM with Gamma"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_scaled, color = "GAM with Scaled x"), linewidth = 1, alpha = 0.7) +
    
    # Points for each model's predictions
    # geom_point(aes(y = pred_gam, color = "GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(y = pred_gam_identity, color = "GAM with Zero Intercept"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_linear, color = "GAM with Linear Term"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_cr, color = "GAM with CR Basis"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_gamma, color = "GAM with Gamma"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_scaled, color = "GAM with Scaled x"), size = 2, alpha = 0.7, show.legend = FALSE) +
    
    #points from Sharp dataset for comparison
    geom_point(data = SA_predictions, aes(x = x, y = y), color = "black", size = 2, alpha = 0.5, show.legend = FALSE) + # Points are black
    
    # Labels and theme
    labs(
      x = "Max diameter (cm)",
      y = "Surface area (m2)",
      title = "Model Predictions for Surface Area vs. Max Diameter (Survey Data)",
      color = "Model"
    ) +
    scale_color_manual(
      values = c(
        "GAM" = "red",
        "GAM with Zero Intercept" = "pink",
        "GAM with Linear Term" = "green",
        "GAM with CR Basis" = "purple",
        "GAM with Gamma" = "orange",
        "GAM with Scaled x" = "brown"
      )
    ) +
    # xlim(0, max(SA_predictions$x)) +
    # ylim(0, max(SA_predictions$y)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  
  # stopping point - 16 sep 2024
  # - I think the GAM comparisons look great, and if anything are further support for using a basic GAM
  # - 11 sq meters is the largest predicted coral surface area - and that's about the size of a parking lot space. seems reasonable to me
  # - add tick marks to show observations in 'survey' - helps to understand level of extrapolation for outlier large corals
  # - now just move forward with finalizing the script cleaning below, re-running analysis with bleaching info
  
  
  #incorporate pre-SCTLD old mortality from lower Florida Keys plots (same as the ones in 'survey'). last observation of old mortality
  # also replace 'FWRI database' instances with 'NA's, for clarity, in relevant 'New_death' columns
  #   before SCTLD arrived was 17 August 2018
  #   SCTLD reached the three sites:
  #   - Offshore: 30 October 2018 [74 days after last observation of old mortality]
  #   - Midchannel: 29 November 2018 [30 days after Offshore]
  #   - Nearshore: 8 February 2019 [71 days after Midchannel; 101 days after Offshore]
  old_mortality = read.csv(here("data", "Updated_June19_SWG.csv")) %>%
    mutate(New_death_110918 = str_replace(New_death_110918, 'FWRI Database', 'NA'),
           New_death_103018 = str_replace(New_death_103018, 'FWRI Database', 'NA'))
  
  #filter to only most recent recordings of mortality pre-SCTLD
  #   - there were normal health states / 'Found' statuses for 10 May 2018. but lots of notes (corals turned over especially)
  #   - there were no new mortality recorded 1 June or 21 June 2018; those are excluded. normal health states / notes / 'Found' statuses as well
  #   - there were normal health states / 'Found' statuses for 16 July 2018
  old_mortality = old_mortality %>% select(Plot, Sps, Max_width, Coral_ID, New_death_051018, Old_death_051018, New_death_060118,
                                           New_death_062118, New_death_071618, Notes_071618, New_death__081718, Old_death_081718,
                                           Health_state_081718, Notes_081718, Found_081718, SCTLD_081718)
  
  # Create a new column indicating total mortality pre-SCTLD
  old_mortality <- old_mortality %>%
    mutate(
      # Replace NA values in 'Old_death_081718' with 0 (NOTE - assume NA actually means no old mortality - only the case for 3 corals)
      #   Those three corals:
      #     - 1_p24_t8_s5_c21_MCAV
      #     - 2_p27_t8_s5_c30_DSTO
      #     - 3_p46_t4_s0_c30_SINT
      Old_death_081718 = replace_na(Old_death_081718, 0),
      OM = New_death__081718 + Old_death_081718, #calculate total mortality by summing new and old mortality columns
    ) %>%
    select(Plot, Sps, Max_width, Coral_ID, OM) %>%
    rename(spp = Sps)
  
  # STOPPING POINT - 13 SEP 2024
  #   - go down below and consider whether this summing of mortality from 08-17-2018 should be reconsidered!! it might 
  #     conflict with how I am backtracking infections!
  
  #label the rows which match plots later monitored for SCTLD
  old_mortality = old_mortality %>%
    mutate(
      site = case_when(
        Plot == 23 | Plot == 25 ~ "mid",
        Plot == 27 | Plot == 28 ~ "off",
        Plot == 45 | Plot == 47 ~ "near"
      )
    ) %>%
    #keep only the plots continuously monitored for SCTLD. there were extra plots initially part of monitoring!
    filter(site == "mid" | site == "off" | site == "near")
  
  #replace ' (1)', ' (2)', or ' (3)' with '' in dataframes. there were some strange coral_ID entries
  replace_string <- function(data) {
    for (i in 1:ncol(data)) {
      if (is.character(data[, i])) {  # Check if the column is character
        data[, i] <- gsub(' \\((1|2|3)\\)', '', data[, i])
      }
    }
    return(data)
  }
  survey = replace_string(survey)
  old_mortality = replace_string(old_mortality)
  
  # NOTE - mismatched data summary (just for peace of mind):
  #   - 1_p25_t9_s5_c4_PAST. unclear what the mismatch is
  #   - 2_p27_t10_s0_c1_OFAV. 'Sps' column labeled as SINT in 'survey'. Brought this up in Mote call
  #   - 2_p27_t10_s0_c3_OFAV. 'Sps' column labeled as SINT in 'survey'. Brought this up in Mote call
  #   *-* 2_p28_t4_s5_c8_SINT. in 'survey' only; no coordinates info
  #   - 3_p47_t1_s0_c1_MCAV. 'Max_width' misentered somewhere. It is 50 in 'survey', but 25 in 'old_mortality'. Brought this up in Mote call
  #   *-* 3_p47_t5_s0_c24_PAST. in 'survey' only; no coordinates info
  mismatches = anti_join(survey, old_mortality)
  
  # Update species and max width for mismatched data
  survey <- survey %>%
    mutate(spp = case_when(
      Coral_ID == '2_p27_t10_s0_c1_OFAV' ~ 'OFAV',
      Coral_ID == '2_p27_t10_s0_c3_OFAV' ~ 'OFAV',
      TRUE ~ spp  # Keep existing values unchanged
    ))
  
  old_mortality <- old_mortality %>%
    mutate(Max_width = case_when(
      Coral_ID == '3_p47_t1_s0_c1_MCAV' ~ 50,
      TRUE ~ Max_width  # Keep existing values unchanged
    ))
  
  #merge in the pre-SCTLD old mortality data
  survey = left_join(
      survey,
      old_mortality,
      by = c("Coral_ID", "Plot", "site", 'spp', 'Max_width'),
      copy = FALSE,
      suffix = c(".x", ".y"),
      keep = NULL
    )
  
  #clean up the metadata after joining
  survey = survey %>%
    rename(diam = Max_width) %>%
    select(-Plot)
  
  # Update tot_mortality for specific Coral_ID values
  survey <- survey %>%
    mutate(OM = case_when(
      Coral_ID == '2_p28_t4_s5_c8_SINT' ~ 0, #was in 'survey' only, assume started SCTLD outbreak with zero old mortality
      Coral_ID == '3_p47_t5_s0_c24_PAST' ~ 0, #was in 'survey' only, assume started SCTLD outbreak with zero old mortality
      grepl('p25_t9_s5_c4_PAST', Coral_ID) ~ 10, #possible issue with 'char' string, ensure it has 10% old mortality
      TRUE ~ OM  # Keep existing values unchanged
    ))
  
  #pivot survey to long format
  survey_long <- survey %>%
    # Filter based on Site_type
    filter(site %in% c("near", "mid", "off")) %>%
    # Pivot to long format
    pivot_longer(
      cols = starts_with("X"),
      names_to = "date"
    ) %>%
    # Convert 'tot_diseased' to a factor
    mutate(tot_diseased = as.factor(tot_diseased)) %>%
    rename(D = tot_diseased, status_field = value) %>%
    mutate(
      D = case_when(
        D == "Dis" ~ "Y",
        D == "Health" ~ "N",
      ),
      status_model = status_field
    )
  
  #create unique numeric value for each coral
  survey_long <- survey_long %>%
    group_by(Coral_ID) %>%
    mutate(coral_numID = cur_group_id()) %>%
    ungroup() %>%
    select(coral_numID, everything())
  
  #append the same coral numerical ID to 'prograte', a dataframe containing time-series lesion progression rates of each diseased coral colony
  # NOTE - better approach would have been to join prograte with survey dataframe at beginning of script. it's fine though
  prograte = prograte %>%
    mutate(coral_numID = survey_long$coral_numID[match(Coral_ID, survey_long$Coral_ID)])

  #add index to each dataframe row (standard convention for IDs, but also makes breaking it apart and then merging & joining later easier)
  survey_long <- survey_long %>%
    mutate(ID = row_number()) %>%
    select(ID, everything())
  
  #format progrates date columns
  names(prograte)[names(prograte) == 'X'] <- '_X'#this old column happens to start with an X, like the date columns - so, rename it
  for(i in 1:ncol(prograte)){
    if(substr(colnames(prograte[i]), start = 1, stop = 1) == "X"){
      newdatename = as.character(as.POSIXct(str_sub(colnames(prograte[i]), 2), format = "%m.%d.%y"))
      colnames(prograte)[i] = newdatename
    }
  }

  #initialize the dataframe with all required columns for calculating tissue loss through time
  survey_long <- survey_long %>%
    # Initialize columns with NA values and convert to numeric
    mutate(
      start_perctiss = NA_real_,
      percloss = NA_real_, #percentage loss of live tissue between timepoints
      cum_percloss = NA_real_,
      percinf = NA_real_, #estimated instantaneous (daily) %loss of live tissue - proxy of instantaneous %infectious tissue
      progdays = NA_real_, #amount of days between timepoints
      remaintiss = NA_real_,
      inftiss = NA_real_,
      cum_tissloss = NA_real_,
      start_perctiss = 1 - (OM / 100), #initialize amount of tissue on each coral colony
      date = str_sub(date, 2) %>%
        as.factor() %>%
        as.POSIXct(format = "%m.%d.%y"), #cnvert date format to POSIXct
      surpdead = 'N', #prep dataset for new column indicating surprise-dead (Y/N) and "ever died" ('died') statuses
      died = 'N'
    ) %>%
    rename(coral_long_ID = Coral_ID) %>%
    arrange(coral_numID, date) #sort by coral ID/date to ensure proper temporal analysis downstream
  
  #filter down to post-SCTLD introduction (October 30th 2018) infections
  # NOTE / STOPPING POINT - make sure this doesn't mess up the backtracking
  surveydiseased = survey_long %>%
    filter(
      D == "Y",
      date >= as.POSIXct("2018-10-30")
    ) %>%
    arrange(coral_numID, date) #failsafe to ensure preservation of ID-sorting
  
  #remove exact rows from 'survey_long' as are being extracted for 'surveydiseased', to easily rbind back together after data wrangling
  survey_trimmed = anti_join(survey_long, surveydiseased) %>%
    arrange(coral_numID, date)
  
  #prep original dataset for new column indicating surprise-dead (Y/N) and "ever died" ('died') statuses
  survey_trimmed = survey_trimmed %>%
    mutate(surpdead = 'N', died = 'N') %>%
    arrange(coral_numID, date)
  
  # 'SURPRISE DEAD CORALS'
  #
  #there were 57 corals that died suddenly between timepoints but were otherwise observed as healthy. these were all quite small and likely died from
  # SCTLD during the first and second infection waves. to account for this cryptic tissue loss, convert those final timepoints into 100% mortality from
  # disease
  # note - to find these, filter by 2019-12-06 and D = 'N' & status = 'Dead'
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
    filter(date == as.POSIXct("2019-12-06"), D == 'N', status_field == 'Dead') %>%
    mutate(surpdead = 'Y', died = 'Y', D = 'Y') %>% #update presumed disease history to SCTLD-positive
    arrange(coral_numID, date)

  #pull the full T1 - T26 rows for each surprise-dead coral
  first.dead.full = survey_trimmed %>%
    filter(coral_numID %in% surprise.dead$coral_numID) %>%
    mutate(surpdead = 'Y', died = 'Y') %>%
    arrange(coral_numID, date)

  #filter to the date that the surprise-dead coral was first documented as dead
  first.dead = first.dead.full %>%
    filter(grepl("\\Dead", status_field)) %>%
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
  
  #calculate accumulated percloss
  first.dead.full <- first.dead.full %>%
    group_by(coral_numID) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(cum_percloss = cumsum(replace_na(percloss, 0))) %>%
    ungroup()
  
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
  
  #update backtracked SCTLD-infected status for the timepoint prior to the surprise-dead date
  first.dead.full <- first.dead.full %>%
    mutate(
      status_model = ifelse(!is.na(percloss) & percloss < 100, 'backtracked_SCTLD', status_model)
    )
  
  
  
  # STOPPING POINT - final touch-up: insert correct status for backtracked infections. maybe implement as new column to distinguish
  #   from status as defined by the field findings? or just supersede it.
  

  
  # Drop the prior_date column
  first.dead.full$prior_date <- NULL
  
  #update original dataframe with surprise-dead coral mortality
  survey_trimmed = survey_trimmed %>%
    left_join(
      first.dead.full %>% select(coral_numID, date, status_model, surpdead, died, dead_date, percloss, cum_percloss),
      by = c("coral_numID", "date")) %>%
    mutate(
      status_model.x = coalesce(status_model.y, status_model.x),
      percloss.x = coalesce(percloss.y, percloss.x),
      surpdead.x = coalesce(surpdead.y, surpdead.x),
      died.x = coalesce(died.y, died.x),
      cum_percloss.x = coalesce(cum_percloss.y, cum_percloss.x)
    ) %>%
    rename(status_model = status_model.x, 
           percloss = percloss.x,
           surpdead = surpdead.x,
           died = died.x,
           cum_percloss = cum_percloss.x) %>%
    select(-status_model.y, -percloss.y, -surpdead.y, -died.y, -cum_percloss.y) %>%
    relocate(dead_date, .after = last_col()) %>%
    arrange(coral_numID, date)

  # Prepare a dataframe for calculating infected tissue in surprise-dead corals (T1 - T26)
  # only calculating instantaneous infected tissue (tissue loss / sloughing within 24 hours) for corals that haven't completed died already
  # NOTE - the way 'percloss' is handled here is different than the loop for confirmed diseased corals below. here, it was assigned
  #         manually by me, to backtrack daily percentage loss back to the last timepoint (and this translates directly to 
  #         instantaneous infected tissue in a day). but below, percloss was a measure of how much tissue was lost between timepoints
  #         in confirmed SCTLD-infected corals. that's why it needs be backtracked and then divided by the number of days between
  #         timepoints. I might have benefited from using a different term than 'percloss' for the surprise-dead corals - could maybe
  #         go back and call it 'instantloss' or 'dailyloss'. similar nomenclature may be useful below as well
  surprise.dead.infections <- survey_trimmed %>%
    filter(coral_numID %in% surprise.dead$coral_numID) %>%
    group_by(coral_numID) %>%
    filter(!all(OM == 100)) %>%  # Exclude coral_numIDs with 100% old mortality before SCTLD surveying
    ungroup() %>%
    mutate(
      percinf = if_else(!is.na(percloss) & percloss > 0 & percloss != 100, percloss, NA_real_), #exclude corals that are 100% dead (from SCTLD)
      inftiss = if_else(!is.na(percinf),
                        (percinf / 100) * tissue * start_perctiss,
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
  
  # 'CORALS DOCUMENTED AS DISEASED'
  #
  #steps:
  # 1.) backtrack percloss of current timepoint to last timepoint (account for patient zeros somehow): percloss/progdays
  #      - assumes that the coral became infected right as the surveyor was ascending back to the surface (convenient)
  # 2.) apply that backtracked value (inftiss) to the last timepoint, not the current one
  # 3.) the last inftiss is on the date preceding the final day of new percloss. whether that is lesion cessation or host death
  
  #prepare a dataframe for calculating infected tissue in corals documented as diseased in the field
  observed.infected = surveydiseased %>%
    mutate(percloss = NA_real_,
           progdays = NA_real_,
           percinf = NA_real_,
           inftiss = NA_real_) %>%
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
    
    # dateind = match(curr_values$date, colnames(prograte))
    # percloss = as.numeric(as.matrix(setDT(prograte)[coral_numID == curr_coral_ID, ])[dateind])
    
    percloss <- prograte %>%
      filter(coral_numID == curr_coral_ID) %>%
      select(all_of(curr_values$date)) %>%
      pull() %>%
      { if (length(.) == 0) NA_real_ else as.numeric(.) } #set as NA if no match in prograte dataframe
    
    IDslice = surveydiseased %>%
      filter(ID == curr_ID)
    availtiss = IDslice$start_perctiss

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
        
        # Skip to the next iteration. this is an error checker for any problem corals (i.e., they don't have a timepoint before 10-30-2018)
        next
      }
      
      # Calculate progdays and update in observed.infected
      progdays = as.numeric(difftime(currdate, prevdate, units = "days"))
      observed.infected[i, "progdays"] <- progdays
      
      # Calculate percinf and update in observed.infected
      percinf = percloss / progdays
      observed.infected[i-1, "percinf"] <- percinf
      
      # Calculate inftiss and update in observed.infected
      tissue = IDslice$tissue
      inftiss = (percinf / 100) * tissue * availtiss
      observed.infected[i-1, "inftiss"] <- inftiss
      
      # Update backtracked infection status
      if (observed.infected[i-1, "status_model"] == "Healthy") {
        observed.infected[i-1, "status_model"] <- "backtracked_SCTLD"
      }
    }
  }
  
  # Print the error log if any errors were logged
  #   - this is particularly for '2_p27_t2_s0_c1_DSTO', since it is the one patient zero coral with a special date
  if (nrow(error_log) > 0) {
    print("Error log:")
    print(error_log)
  }
  
  #update original dataframe with confirmed diseased coral infections
  surveydiseased = surveydiseased %>%
    left_join(
      observed.infected %>% select(coral_numID, date, status_model, percloss, progdays, percinf, inftiss),
      by = c("coral_numID", "date")) %>%
    mutate(
      status_model.x = coalesce(status_model.y, status_model.x),
      percloss.x = coalesce(percloss.y, percloss.x),
      progdays.x = coalesce(progdays.y, progdays.x),
      percinf.x = coalesce(percinf.y, percinf.x),
      inftiss.x = coalesce(inftiss.y, inftiss.x)
    ) %>%
    rename(status_model = status_model.x,
           percloss = percloss.x,
           progdays = progdays.x,
           percinf = percinf.x,
           inftiss = inftiss.x) %>%
    select(-status_model.y, -percloss.y, -progdays.y, -percinf.y, -inftiss.y) %>%
    arrange(coral_numID, date)

  #calculate accumulated percloss
  surveydiseased <- surveydiseased %>%
    group_by(coral_numID) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(cum_percloss = cumsum(replace_na(percloss, 0))) %>%
    ungroup()
  
  #assign dead dates (first date when accumulated percloss is 100)
  first_death_dates <- surveydiseased %>%
    filter(cum_percloss >= 100) %>%
    group_by(coral_numID) %>%
    summarize(dead_date = min(date)) %>%
    ungroup()
  
  surveydiseased <- surveydiseased %>%
    left_join(first_death_dates, by = "coral_numID") %>%
    arrange(coral_numID, date)
  
  #assemble main dataframe with diseased tissue information
  survey_tissue = rbind(survey_trimmed, surveydiseased) %>%
    arrange(coral_numID, date)
  
  #remove, rearrange, and rename columns for clarity
  # NOTE - there may be timepoint-specific bleaching (thermal stress) data for each colony somewhere. can come up later in results revisions
  survey_tissue <- survey_tissue %>%
    rename(BL = tot_stressed, BL_D = tot_both) %>%
    mutate(
      BL = case_when(
        BL == "S" ~ "Y",
        BL == "NS" ~ "N"
      ),
      BL_D = case_when(
        BL_D == "both" ~ "Y",
        BL_D == "meh" ~ "N"
      )
    ) %>%
    filter(date >= as.Date("2018-08-17") & OM != 100) #remove pre-SCTLD temporal data and any corals that were 100% dead before SCTLD
  #
  # 'CORALS DOCUMENTED AS DISEASED'
  
  # REMOVED / REMAINING TISSUE CALCULATION
  #
  # STOPPING POINT - 13 SEP 2024
  # NOTE - need to figure out getting the dead date (and cum_percloss etc.) to properly propagate before 10-30-2018 - or also I can just delete that timepoint
  
  # 2_p28_t1_s0_c9_MCAV for test coral
  # STOPPING POINT / NOTE - 13 Sep 2024
  #   - seeing here that 'backtracked_SCTLD' does not always fill correctly. check out '2_p28_t4_s0_c5_SSID'
  #   - could also clean up columns: S, I, and R for tissue and percentage should be named and arranged accordingly
  survey_tissue <- survey_tissue %>%
    mutate(
      cum_tissloss = (cum_percloss/100) * (tissue * start_perctiss),
      remaintiss = (tissue * start_perctiss) - cum_tissloss,
    )
  
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
  # REMOVED / REMAINING TISSUE CALCULATION
  
  #make special case for patient zero 2_p27_t2_s0_c1_DSTO on 10-30-2018 since loops don't handle its exception. assume its 90%
  # tissue loss occurred over 14 days since we don't have info except 3 months prior when colony was healthy
  progdays.new <- 14
  reference_date <- ymd('2018-10-30')
  new_date <- reference_date - days(progdays.new) #backtrack 14 days from surprise-dead day of 10-30-2018

  # Calculate percinf (90% percloss over the time difference)
  percinf.new <- 90 / progdays.new
  
  # Extract the existing row for the coral with coral_longID '2_p27_t2_s0_c1_DSTO'
  coral_row <- survey_tissue %>%
    filter(coral_long_ID == '2_p27_t2_s0_c1_DSTO' & date == reference_date) #%>%
    # select(-date) # Remove the date column to replace it with the new date
  
  # Calculate inftiss
  inftiss.new <- (percinf.new / 100) * coral_row$tissue * coral_row$start_perctiss
  
  # Create the new row with the calculated values
  new_row <- coral_row %>%
    mutate(
      date = new_date,
      progdays = progdays.new,
      percinf = percinf.new,
      inftiss = inftiss.new,
      ID_SIR = coral_row$ID,
      ID = NA  # Set ID to NA for the new row
    )
  
  # Insert the new row into survey_tissue and update `ID_SIR` column to be identical to `ID` and increment the values after the new row
  survey_tissue <- survey_tissue %>%
    mutate(ID_SIR = NA) %>%
    add_row(new_row) %>%
    arrange(coral_numID, date) %>%
    mutate(
      ID_SIR = ifelse(is.na(ID_SIR), ID, ID_SIR),  # Keep the correct `ID_SIR`
      ID_SIR = ifelse(
        ID >= new_row$ID_SIR & !is.na(ID), ID_SIR + 1, ID_SIR  # Increment after new row
      )
    )
  
  # STOPPING POINT - Sep 13 2024
  #   - almost there!!!!
  #   - I don't think I was handling the patient zero coral quite right. something is awry there with the two '90' percloss values
  #   - in addition to above stopping points, need to add in something that double-checks 'BL', 'D', 'BL_D', and 'status' actually make
  #       sense with the surprise dead corals included in the dataset as dead from SCTLD. I think just 'D' and 'status' need to be
  #       looked at [update - yes, 'BL' is fully from field data that I cannot corroborate, so take it as truth, and update 'D'
  #       perhaps as 'D_model', and rename 'D' to 'D_field'. 'BL_D' should be updated accordingly. 'status' is taken care of]
  #

  # AGGREGATE DATA INTO SUMMARY TABLES
  #
  ###sum tissue cover values across colonies into bins (Site, Sus_Cat, Date)
  survey_tissue$date_factor = as.factor(survey_tissue$date)
  levels(survey_tissue$date_factor)
  
  # STOPPING POINT - should decide if 2018-08-17 timepoint (T1 now) is required at all - I don't think so. patient zero 2018-10-16 should
  # be all that is needed
  survey_tissue = within(survey_tissue, {
    Timepoint = paste("T", factor(date_factor, labels=seq(unique(date_factor))), sep="")
  })
  
  survey_tissue$fill = survey_tissue$tissue*(1-survey_tissue$OM/100)
  indices = which(is.na(survey_tissue$remaintiss))
  survey_tissue$remaintiss[is.na(survey_tissue$remaintiss)] = survey_tissue$fill[indices]
  survey_tissue = survey_tissue[,-which(colnames(survey_tissue) %in% 'fill')]
  
  #clean up dataframe
  survey_tissue$Site_type = as.factor(survey_tissue$site)
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
  