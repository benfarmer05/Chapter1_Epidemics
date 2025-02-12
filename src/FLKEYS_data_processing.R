  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(mgcv)
  library(gratia)
  
  survey = read.csv(here("data", "SCTLD_END_Vpub_ts.csv"))
  cover = read.csv(here("data", "cover_long.csv"))
  prograte = read.csv(here("data", "SWG_SCTLDprogrates.csv"))
  cover = read.csv(here("data/cover_long.csv"))
  
  #   LOW SUSCEPTIBILITY (LS)
  #         - Slow onset/slow tissue loss [SSID] -> high cover
  #         - moderate onset/slow tissue loss [SINT, PCLI] -> high cover
  #   MODERATE SUSCEPTIBILITY (MS)
  #         - moderate onset/moderate tissue loss [MCAV] -> moderate cover
  #         - Fast onset/slow tissue loss [OANN, OFAV] -> low cover
  #         - Fast onset/moderate tissue loss [SBOU] -> low cover
  #   HIGH SUSCEPTIBILITY (HS)
  #         - Fast onset/fast tissue loss [DSTO, CNAT, PSTR, DLAB] -> moderate cover
  
  #   PATIENT ZERO CORALS
  #   - Offshore: 2_p27_t2_s0_c1_DSTO, 10-30-2018, 90% loss. backtrack a week? this is only a 5 inch coral
  #   - Midchannel: 1_p25_t2_s0_c22_DSTO, 11-29-2018, 10% loss
  #   - Nearshore: 3_p47_t3_s0_c8_PSTR (10% loss) and 3_p47_t4_s5_c15_PSTR (5% loss), 2019-02-08
  
  #define scalars to minimize starting infection amounts (for patient zeros) downstream
  # NOTE - a better solution may be seeing if the models can use a standard amount of starting tissue (e.g., 1e-4) so that R0 is
  #         more standardized and can be better projected onto different populations
  polyp_SA.minimizer.nearshore = 1
  # polyp_SA.minimizer.nearshore = 0.2 #a test to actually amplify polyp_SA since nearshore is having issues
  # polyp_SA.minimizer.nearshore = 10 #10X less than the smallest inftiss value across patient zero's time series (NOTE - should be coded in accordance with above)
  polyp_SA.minimizer.midchannel = 3
  # polyp_SA.minimizer.midchannel = 25 #25X - did this because the first loss of tissue at midchannel is quite high without a strong dampening applied
  polyp_SA.minimizer.offshore = 5
  # polyp_SA.minimizer.offshore = 66 #this was tweaked to be about 10X less than the smallest inftiss value across patient zero's time series. could (should) be coded in
  
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
        spp %in% c('PCLI', 'SINT', 'SSID') ~ 'low', #lowest rates of progression - however, OANN/OFAB have short disease onset time
        spp %in% c('OANN', 'OFAV', 'MCAV', 'SBOU') ~ 'moderate', #moderate rates of progression - but MCAV slow disease onset! SBOU quick onset
        spp %in% c('CNAT', 'DLAB', 'DSTO', 'MMEA', 'MYCE', 'PSTR') ~ 'high', #quickest rates of  progression and shortest disease onset times
        spp %in% c('AAGA', 'ACER', 'OCUL', 'ODIF', 'PAST', 'PDIV', 'PPOR', 'SRAD') ~ 'Unaffected'
      )
    ) %>%
    # filter(!grepl('Unaffected', susc)) %>% #remove presumed SCTLD-unaffected corals, mainly PAST
    select(-coords_x, -coords_y, -total, -total_bin, -days_dis, -weeks_dis) # NOTE - could further analyze days/weeks diseased for stats

  #my understanding from Williams 2021 is that there are 100 photos per site, with CPCe-based % cover being the value associated with
  #   each photo. so, below I have taken the mean, per site, of all 100 photos to get site-level % cover
  cover = cover %>%
    mutate(Site = case_when(
      Site == "Midchannel" ~ "mid",
      Site == "Nearshore"  ~ "near",
      Site == "Offshore"   ~ "off"
    )) %>%
    filter(timept == 't1') %>% #filter to pre-SCTLD assessment of coral cover
    filter(plotnum %in% c('p23', 'p25', 'p27', 'p28', 'p45', 'p47')) %>%  #filter to only plots included in SCTLD monitoring
    group_by(Site) %>%
    summarize(mean.percent.cover = mean(percent.cover, na.rm = TRUE), .groups = 'drop')
  
  #pull middle Florida Keys (Sharp et al. 2020) dataset, which had similar data as 'survey' but measured colonies in 3 dimensions
  #   rather than one. use this dataset to create a statistical association between maximum diameter & colony surface area (SA) by assuming
  #   a hemi-ellipsoid (Holstein et al. 2015), and then predicting SA in 'survey' dataset using maximum diameter alone
  Sharp2020 = read.csv(here("data", "Sharp_Maxwell_2020.csv")) %>%
    select(Plot, Transect, Coral, Species, Width, Width_2, Height) %>%
    mutate(across(where(is.numeric), ~na_if(., -99))) %>% #convert -99 fill values to NA
    drop_na() #remove NAs from dataset
  
  #set '0' cm heights as an arbitrarily small amount (cm) to allow the hemi-ellipsoid estimation to function correctly. coral recruits 
  #   have very little tangible height, and were recorded underwater as 0 cm height
  Sharp2020 = Sharp2020 %>%
    mutate(Height = ifelse(Height == 0, 0.01, Height)) %>%
    distinct() #also filter out repeats (not interested in a time series here)
  
  #hemi-ellipsoid estimation. p is a dimensionless constant; all else in square cm
  Sharp2020 = Sharp2020 %>%
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
  SA = Sharp2020 %>%
    select(x = Width, y = SA) 
  # plot(SA$x, SA$y)
  
  # #a measure to prevent over-predicting surface area of the largest corals from small sample size [decided as unnecessary]
  # SA = SA %>% filter(x < 122) #145
  
  # NOTE - GAM comparisons below. one fitting issue may be differently scaled x and y, and also repeated x values
  #   - values dipping below zero for smallest corals is a problem for GAMs here. Gamma distribution addresses this fairly well
  SA.linear = lm(y ~ x + 0, data = SA) #trying to make sure predictions don't go below zero - struggling with this
  SA.GAM = gam(y ~ s(x), data = SA, family = gaussian, method = "REML")
  SA.GAM_identity = gam(y ~ s(x, bs = "cr"), data = SA, family = gaussian(link = "identity")) #GAM model with a zero intercept constraint
  SA.GAM_linear = gam(y ~ s(x) + x, data = SA, family = gaussian, method = "REML")
  SA.GAM_cr = gam(y ~ s(x, bs = "cr", k = 20), data = SA, family = gaussian, method = "REML", gamma = 2.0)
  SA.GAM_gamma = gam(y ~ s(x), data = SA, family = Gamma(link = "log"), method = "REML")
  SA$x_scaled = scale(SA$x)  # rescale x
  SA.GAM_scaled = gam(y ~ s(x_scaled, bs = "cr", k = 20), data = SA, family = gaussian, method = "REML")
  
  # More GAM variations with additional smooth terms, penalties, and log adjustments
  #
  # Combination of smooth terms and varying degrees of freedom
  SA.GAM_complex = gam(y ~ s(x, bs = "tp", k = 25) + te(x), data = SA, family = gaussian, method = "REML")
  # SA.GAM_complex = gam(y ~ s(x, bs = "cr", k = 10) + s(x, bs = "tp", k = 10), data = SA)
  
  # Penalize large extrapolation (use gamma > 1)
  SA.GAM_penalized = gam(y ~ s(x, bs = "tp", k = 25), data = SA, family = gaussian, method = "REML", gamma = 1.4)
  
  # Applying log transformation to the dependent variable and fitting a model
  SA.GAM_log_y = gam(log(y) ~ s(x, bs = "cr", k = 20), data = SA, family = gaussian(link = "identity"), method = "REML")
  
  # Introducing additional model with adaptive smooth
  SA.GAM_adaptive = gam(y ~ s(x, bs = "ad", k = 25), data = SA, family = gaussian, method = "REML")
  
  # Modifying model to use a different number of knots and different basis functions
  SA.GAM_varying_knots = gam(y ~ s(x, bs = "cr", k = 30), data = SA, family = gaussian, method = "REML")
  
  #non-wiggly gamma GAMs
  SA.GAM_gamma = gam(y ~ s(x), data = SA, family = Gamma(link = "log"), method = "REML")
  SA.GAM_gamma_1.5 = gam(y ~ s(x), data = SA, family = Gamma(link = "log"), method = "REML", gamma = 3)
  SA.GAM_gamma_tp = gam(y ~ s(x, bs = "tp", k = 20), data = SA, family = Gamma(link = "log"), method = 'REML', gamma = 1.5)
  SA.GAM_gamma_ps = gam(y ~ s(x, bs = "ps", k = 20), data = SA, family = Gamma(link = "log"), method = 'REML', gamma = 1.5)
  
  # #diagnostics
  # summary(SA.GAM)
  # plot(SA.GAM, residuals = TRUE)
  # draw(SA.GAM, select = 1)
  # gam.check(SA.GAM)
  # appraise(SA.GAM)
  # 
  # summary(SA.GAM_gamma_tp)
  # plot(SA.GAM_gamma_tp, residuals = TRUE)
  # draw(SA.GAM_gamma_tp, select = 1)
  # gam.check(SA.GAM_gamma_tp)
  # appraise(SA.GAM_gamma_tp)
  
  # Generate predictions from each model
  SA_predictions = SA %>%
    mutate(
      pred_linear = predict(SA.linear),
      pred_gam = predict(SA.GAM),
      pred_gam_identity = predict(SA.GAM_identity), #, newdata = SA
      pred_gam_linear = predict(SA.GAM_linear),
      pred_gam_cr = predict(SA.GAM_cr),
      pred_gam_gamma = predict(SA.GAM_gamma, type = "response"),
      pred_gam_scaled = predict(SA.GAM_scaled), #, newdata = SA
      pred_gam_complex = predict(SA.GAM_complex),
      pred_gam_penalized = predict(SA.GAM_penalized),
      pred_gam_log_y = exp(predict(SA.GAM_log_y)),  # Undo log-transformation for predictions
      pred_gam_adaptive = predict(SA.GAM_adaptive),
      pred_gam_varying_knots = predict(SA.GAM_varying_knots),
      pred_gam_gamma_1.5 = predict(SA.GAM_gamma_1.5, type = 'response'),
      pred_gam_gamma_tp = predict(SA.GAM_gamma_tp, type = 'response'),
      pred_gam_gamma_ps = predict(SA.GAM_gamma_ps, type = 'response')
    )
  
  ggplot(SA_predictions, aes(x = x, y = y)) +
    geom_point(color = "black", size = 2, alpha = 0.5, show.legend = FALSE) + # Points are black
    # geom_line(aes(y = pred_linear, color = "Linear"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam, color = "GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_identity, color = "GAM with Zero Intercept"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_linear, color = "GAM with Linear Term"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_cr, color = "GAM with CR Basis"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_gamma, color = "GAM with Gamma"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_scaled, color = "GAM with Scaled x"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_complex, color = "Complex GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_penalized, color = "Penalized GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_log_y, color = "Log-transformed GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_adaptive, color = "Adaptive GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_varying_knots, color = "Varying Knots GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_gamma_1.5, color = "GAM with Gamma 1.5"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_gamma_tp, color = "GAM with Tensor Product"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_gamma_ps, color = "GAM with P-splines"), linewidth = 1, alpha = 0.7) +
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
        "GAM with Scaled x" = "brown",
        "Complex GAM" = "cyan",
        "Penalized GAM" = "magenta",
        "Log-transformed GAM" = "yellow",
        "Adaptive GAM" = "darkgreen",
        "Varying Knots GAM" = "darkblue",
        "GAM with Gamma 1.5" = 'darkorange',
        "GAM with Tensor Product" = 'darkmagenta',
        "GAM with P-splines" = 'darkred'
      )
    ) +
    # xlim(0, 20) +
    # ylim(0, 0.10) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  #Akaike comparison
  aic_values = AIC(SA.linear, SA.GAM, SA.GAM_identity, SA.GAM_linear, SA.GAM_cr, SA.GAM_gamma, SA.GAM_scaled, 
      SA.GAM_complex, SA.GAM_penalized, SA.GAM_log_y, SA.GAM_adaptive, SA.GAM_varying_knots,
      SA.GAM_gamma_1.5, SA.GAM_gamma_tp, SA.GAM_gamma_ps)
  
  # Create a tibble with model names, degrees of freedom, and AIC values
  model_names = c("SA.linear", "SA.GAM", "SA.GAM_identity", "SA.GAM_linear", "SA.GAM_cr", 
                   "SA.GAM_gamma", "SA.GAM_scaled", "SA.GAM_complex", "SA.GAM_penalized", 
                   "SA.GAM_log_y", "SA.GAM_adaptive", "SA.GAM_varying_knots", 
                   "SA.GAM_gamma_1.5", "SA.GAM_gamma_tp", "SA.GAM_gamma_ps")
  
  aic_df = tibble(
    Model = model_names,
    df = aic_values$df,
    AIC = aic_values$AIC
  )
  
  # Sort by AIC values in ascending order; view "best" models
  sorted_aic = aic_df %>%
    arrange(AIC)  
  sorted_aic
  
  # Generate predictions for the survey data
  survey_predictions = survey %>%
    mutate(
      pred_gam = predict(SA.GAM, newdata = data.frame(x = Max_width)),
      pred_gam_identity = predict(SA.GAM_identity, newdata = data.frame(x = Max_width)),
      pred_gam_linear = predict(SA.GAM_linear, newdata = data.frame(x = Max_width)),
      pred_gam_cr = predict(SA.GAM_cr, newdata = data.frame(x = Max_width)),
      pred_gam_gamma = predict(SA.GAM_gamma, newdata = data.frame(x = Max_width), type = 'response'),
      pred_gam_scaled = predict(SA.GAM_scaled, newdata = data.frame(x_scaled = scale(Max_width))),
      pred_gam_complex = predict(SA.GAM_complex, newdata = data.frame(x = Max_width)),
      pred_gam_penalized = predict(SA.GAM_penalized, newdata = data.frame(x = Max_width)),
      pred_gam_log_y = exp(predict(SA.GAM_log_y, newdata = data.frame(x = Max_width))),  # Undo log-transformation for predictions
      pred_gam_adaptive = predict(SA.GAM_adaptive, newdata = data.frame(x = Max_width)),
      pred_gam_varying_knots = predict(SA.GAM_varying_knots, newdata = data.frame(x = Max_width)),
      pred_gam_gamma_1.5 = predict(SA.GAM_gamma_1.5, newdata = data.frame(x = Max_width), type = 'response'),
      pred_gam_gamma_tp = predict(SA.GAM_gamma_tp, newdata = data.frame(x = Max_width), type = 'response'),
      pred_gam_gamma_ps = predict(SA.GAM_gamma_ps, newdata = data.frame(x = Max_width), type = 'response')
    )
  
  # Plot predictions with ggplot
  ggplot(survey_predictions, aes(x = Max_width)) +
    
    # Lines for each model
    geom_line(aes(y = pred_gam, color = "GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_identity, color = "GAM with Zero Intercept"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_linear, color = "GAM with Linear Term"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_cr, color = "GAM with CR Basis"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_gamma, color = "GAM with Gamma"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_scaled, color = "GAM with Scaled x"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_complex, color = "Complex GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_penalized, color = "Penalized GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_log_y, color = "Log-transformed GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_adaptive, color = "Adaptive GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_varying_knots, color = "Varying Knots GAM"), linewidth = 1, alpha = 0.7) +
    # geom_line(aes(y = pred_gam_gamma_1.5, color = "GAM Gamma 1.5"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_gamma_tp, color = "GAM Tensor Product"), linewidth = 1, alpha = 0.7) +
    geom_line(aes(y = pred_gam_gamma_ps, color = "GAM P-splines"), linewidth = 1, alpha = 0.7) +
    
    # Points for each model's predictions
    geom_point(aes(y = pred_gam, color = "GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_identity, color = "GAM with Zero Intercept"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_linear, color = "GAM with Linear Term"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_cr, color = "GAM with CR Basis"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_gamma, color = "GAM with Gamma"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_scaled, color = "GAM with Scaled x"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_complex, color = "Complex GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_penalized, color = "Penalized GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_log_y, color = "Log-transformed GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_adaptive, color = "Adaptive GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_varying_knots, color = "Varying Knots GAM"), size = 2, alpha = 0.7, show.legend = FALSE) +
    # geom_point(aes(y = pred_gam_gamma_1.5, color = "GAM Gamma 1.5"), size = 2, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(y = pred_gam_gamma_tp, color = "GAM Tensor Product"), size = 2, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(y = pred_gam_gamma_ps, color = "GAM P-splines"), size = 2, alpha = 0.7, show.legend = FALSE) +
    
    #points from Sharp dataset for comparison
    geom_point(data = SA_predictions, aes(x = x, y = y), color = "black", size = 2, alpha = 0.5, show.legend = FALSE) + # Points are black

    # # Add density plot or histogram for observation density
    # geom_density(data = survey, aes(x = Max_width), color = "blue", size = 1, alpha = 0.3, linetype = "dashed") +  # Density plot
    # # geom_histogram(data = survey, aes(x = Max_width, y = ..density..), bins = 30, fill = "blue", alpha = 0.3) +  # Histogram
    
    # Add tick marks on the bottom
    geom_rug(data = survey, aes(x = Max_width), sides = "b", color = "grey", alpha = 0.5) +
    
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
        "GAM with Scaled x" = "brown",
        "Complex GAM" = "cyan",
        "Penalized GAM" = "magenta",
        "Log-transformed GAM" = "yellow",
        "Adaptive GAM" = "darkgreen",
        "Varying Knots GAM" = "darkblue",
        "GAM Gamma 1.5" = "darkred",
        "GAM Tensor Product" = "darkorange",
        "GAM P-splines" = "darkcyan"        )
    ) +
    # xlim(0, max(SA_predictions$x)) +
    # ylim(0, max(SA_predictions$y)) +
    # xlim(0, 20) +
    # ylim(0, 0.10) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Generate predictions and confidence intervals
  plot_data <- smooth_estimates(SA.GAM_gamma_tp, smooth = "s(x)", overall = TRUE)
  
  # Add confidence intervals
  plot_data <- plot_data %>%
    mutate(
      lower_ci = .estimate - 1.96 * .se,  # 95% CI lower bound
      upper_ci = .estimate + 1.96 * .se   # 95% CI upper bound
    )
  
  # Plot with confidence intervals
  ggplot(plot_data, aes(x = x, y = .estimate)) +
    geom_line(color = "blue", size = 1) +  # GAM prediction line
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3) +  # Confidence intervals
    labs(
      title = "GAM Prediction with Uncertainty",
      x = "Maximum Diameter (x)",
      y = "Surface Area (SA)"
    ) +
    theme_minimal()  
  
  # - NOTE - The GAM gamma p-spline and tensor product predictions both look pretty good. They account well for the right-skewed,
  #           positive, low-variance nature of the dataset - or seem to - because the small corals are predicted well. but, there is 
  #           still some over-prediction at the extremes, something a regular GAM performs better with and avoids. so, the below:
  #       - but, still manually editing the one or two largest extrapolated corals to a reasonable size
  #
  # Retrieve row indices in 'survey' for the two largest corals
  largest_width_rows = survey %>%
    slice_max(Max_width, n = 2, with_ties = FALSE)
  largest_width_indices = which(survey$Max_width %in% largest_width_rows$Max_width)
  
  # Update survey with predicted SA from the "best" GAM (gamma with log-link, tensor product). extrapolated values from the two largest
  #   corals as predicted by the simple GAM are pasted in, since a simple GAM does not over-predict as much at the extremes
  survey = survey %>%
    left_join(survey_predictions %>% select(Coral_ID, pred_gam_gamma_tp, pred_gam), by = "Coral_ID") %>%
    mutate(colony_SA = pred_gam_gamma_tp)
  survey$colony_SA[largest_width_indices] = survey_predictions$pred_gam[largest_width_indices]
  survey = survey %>% select(-pred_gam_gamma_tp, -pred_gam)
  
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
  old_mortality = old_mortality %>%
    mutate(
      # Replace NA values in 'Old_death_081718' with 0 (NOTE - assume NA actually means no old mortality - only the case for 3 corals)
      #   Those three corals:
      #     - 1_p24_t8_s5_c21_MCAV
      #     - 2_p27_t8_s5_c30_DSTO
      #     - 3_p46_t4_s0_c30_SINT
      Old_death_081718 = replace_na(Old_death_081718, 0),
      percOM = New_death__081718 + Old_death_081718, #calculate total mortality by summing new and old mortality columns
    ) %>%
    select(Plot, Sps, Max_width, Coral_ID, percOM) %>%
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
  replace_string = function(data) {
    for (i in 1:ncol(data)) {
      if (is.character(data[, i])) {  # Check if the column is character
        data[, i] = gsub(' \\((1|2|3)\\)', '', data[, i])
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
  survey = survey %>%
    mutate(
      spp = case_when(
        Coral_ID == '2_p27_t10_s0_c1_OFAV' ~ 'OFAV',
        Coral_ID == '2_p27_t10_s0_c3_OFAV' ~ 'OFAV',
        TRUE ~ spp  # Keep existing values unchanged
      ),
      susc = case_when(
        Coral_ID %in% c('2_p27_t10_s0_c1_OFAV', '2_p27_t10_s0_c3_OFAV') ~ 'moderate',
        TRUE ~ susc  # Keep existing values unchanged
      )
    )
  
  old_mortality = old_mortality %>%
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
  ) %>%
  group_by(Coral_ID) %>%
  filter(!all(percOM == 100)) %>% #exclude corals with 100% old mortality before SCTLD surveying
  ungroup() %>%
  rename(diam = Max_width) %>%
  select(-Plot)
  
  # Update tot_mortality for specific errant corals
  survey = survey %>%
    mutate(percOM = case_when(
      Coral_ID == '2_p28_t4_s5_c8_SINT' ~ 0, #was in 'survey' only, assume started SCTLD outbreak with zero old mortality
      Coral_ID == '3_p47_t5_s0_c24_PAST' ~ 0, #was in 'survey' only, assume started SCTLD outbreak with zero old mortality
      grepl('p25_t9_s5_c4_PAST', Coral_ID) ~ 10, #possible issue with 'char' string, ensure it has 10% old mortality
      TRUE ~ percOM  # Keep existing values unchanged
    ))
  
  #pivot survey to long format
  survey_long = survey %>%
    # Filter based on Site_type
    filter(site %in% c("near", "mid", "off")) %>%
    # Pivot to long format
    pivot_longer(
      cols = starts_with("X"),
      names_to = "date"
    ) %>%
    # Convert 'tot_diseased' to a factor
    mutate(tot_diseased = as.factor(tot_diseased)) %>%
    rename(D_field = tot_diseased, status_field = value) %>%
    mutate(
      D_field = case_when(
        D_field == "Dis" ~ "Y",
        D_field == "Health" ~ "N",
      ),
      D_model = D_field,
      status_model = status_field
    ) %>%
    select(1:which(names(.) == "D_field"), D_model, everything())  #reorder to place D_model after D_field
  
  #create unique numeric value for each coral
  survey_long = survey_long %>%
    group_by(Coral_ID) %>%
    mutate(coral_numID = cur_group_id()) %>%
    ungroup() %>%
    select(coral_numID, everything())
  
  #append the same coral numerical ID to 'prograte', a dataframe containing time-series lesion progression rates of each diseased coral colony
  # NOTE - better approach would have been to join prograte with survey dataframe at beginning of script. it's fine though
  prograte = prograte %>%
    mutate(coral_numID = survey_long$coral_numID[match(Coral_ID, survey_long$Coral_ID)])

  #add index to each dataframe row (standard convention for IDs, but also makes breaking it apart and then merging & joining later easier)
  survey_long = survey_long %>%
    mutate(ID = row_number()) %>%
    select(ID, everything())
  
  #format progrates date columns
  names(prograte)[names(prograte) == 'X'] = '_X'#this old column happens to start with an X, like the date columns - so, rename it
  for(i in 1:ncol(prograte)){
    if(substr(colnames(prograte[i]), start = 1, stop = 1) == "X"){
      newdatename = as.character(as.POSIXct(str_sub(colnames(prograte[i]), 2), format = "%m.%d.%y"))
      colnames(prograte)[i] = newdatename
    }
  }
  
  #initialize the dataframe with all required columns for calculating tissue loss through time
  survey_long = survey_long %>%
    # Initialize columns with NA values and convert to numeric
    mutate(
      start_perctiss = NA_real_, #starting percentage of whole-colony SA that is live tissue
      percloss = NA_real_, #percentage loss of live tissue between timepoints
      cum_percloss = NA_real_,
      percinf = NA_real_, #estimated instantaneous (daily) %loss of live tissue - proxy of instantaneous %infectious tissue
      progdays = NA_real_, #amount of days between timepoints
      starttiss = NA_real_, #starting amount of tissue (SA in m2)
      remaintiss = NA_real_,
      inftiss = NA_real_,
      cum_tissloss = NA_real_,
      start_perctiss = 100 - percOM, #initialize amount of tissue on each coral colony
      starttiss = colony_SA * (start_perctiss / 100),
      date = str_sub(date, 2) %>%
        as.factor() %>%
        as.POSIXct(format = "%m.%d.%y"), #convert date format to POSIXct
      surpdead = 'N', #prep dataset for new column indicating surprise-dead (Y/N) and "ever died" ('died') statuses
      died = 'N'
    ) %>%
    rename(coral_long_ID = Coral_ID) %>%
    filter(date >= as.Date("2018-08-17")) %>% #remove pre-SCTLD timepoints except for one, to easily establish pre-SCTLD susceptible tissue
    # filter(date >= as.Date("2018-10-30")) %>% #remove pre-SCTLD timepoints
    arrange(coral_numID, date) #sort by coral ID/date to ensure proper temporal analysis downstream
  
  #feed DHW data into survey data
  #this data comes from: https://coralreefwatch.noaa.gov/product/index.php
  # - specifically, 'florida_keys.txt' is Regional Virtual Station data, specific to the Florida region
  # - from link above, navigate to point 3 (Regional Virtual Stations and Bleaching Alert Emails), and then click 'home page' for
  #     'Regional Virtual Stations (5km)'. this will open a page with a world map.
  # - navigate on the map to Florida Keys region (not Southeast Florida region). click 'Time Series Graphs & Data' on the prompt, which will open a new page
  # - click 'Time Series Data' for Florida Keys. this will open the 'florida_keys.txt' ASCII data directly for viewing
  # - you can go back and right click the link itself to save the link as a text file to your local directory
  file_content = readLines(here("data", "florida_keys.txt"))
  data_lines = file_content[22:length(file_content)]
  DHW.CRW = read.table(text = data_lines, header = TRUE)
  DHW.CRW.full = DHW.CRW %>%
    mutate(date = as.Date(paste(YYYY, MM, DD, sep = "-"))) %>%
    select(-YYYY, -MM, -DD) %>%
    select(date, everything())
  DHW.CRW = DHW.CRW.full %>%
    # filter(date >= as.Date("2018-05-01") & date <= as.Date("2019-12-06")) %>% #include first day of surveying in the field
    filter(date >= as.Date("2018-08-17") & date <= as.Date("2019-12-06")) #skip to first SCTLD day
  survey_long <- survey_long %>%
    left_join(DHW.CRW, by = "date")
  
  #prepare dataframe for calculating infections of confirmed-diseased corals
  # NOTE / STOPPING POINT - make sure this doesn't mess up the backtracking
  surveydiseased = survey_long %>%
    filter(
      D_field == "Y",
      date >= as.POSIXct("2018-10-30") # 'prograte' only begins at 10-30-2018
    ) %>%
    arrange(coral_numID, date) #failsafe to ensure preservation of ID-sorting. recapitulated throughout script
  
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
  # note - to find these, filter by 2019-12-06 and D_field = 'N' & status = 'Dead'
  #   - almost all of these corals were small, HS corals which makes sense. important to keep track of their loss between timepoints
  
  #   - Surprise-dead corals which had 100% old mortality in pre-SCTLD (mortality) dataset (so, are being *removed* entirely before simulation)
  #       - 2_p27_t9_s5_c2_DSTO - Maybe a discrepancy in mortality dataset; it suddenly died 09-17-2019 in new dataset. said 100% recent mortality on 08-17-2018 in old dataset
  #             - SW: Probably dead-dead on 8-17-2018. Based on field notes
  #       - 2_p27_t8_s5_c21_DSTO - Maybe a discrepancy in mortality dataset; it suddenly died 05-02-2019. 100% old mortality suddenly on 08-17-2018 in old dataset. VERY small coral (2 cm)
  #             - SW: dead-dead on 8-17-2018.
  #       - 1_p23_t1_s0_c6_CNAT - Maybe a discrepancy in mortality dataset; it suddenly died 08-30-2018. 75% old mortality suddenly on 05-10-2018, then suddenly 100% old mortality or not found on 08-17-2018.
  #             _ SW: dead-dead on 8-17-2018. was 75% dead by 5-10-2018, if that is useful
  
  #filter surprise-dead corals (read details above), by their health condition on the last date of surveying (2019-12-06)
  surprise.dead = survey_trimmed %>%
    group_by(coral_numID) %>%
    filter(D_field == 'N', status_field == 'Dead') %>% #filter the corals that were never marked as being diseased but died unexpectedly
    ungroup() %>%
    filter(date == as.POSIXct("2019-12-06")) %>% #filter to final day of surveying
    mutate(surpdead = 'Y', died = 'Y', D_model = 'Y') %>%  # update presumed disease history to SCTLD-positive
    arrange(coral_numID, date)

  #pull the full T1 - T26 rows for each surprise-dead coral
  first.dead.full = survey_trimmed %>%
    filter(coral_numID %in% surprise.dead$coral_numID) %>%
    mutate(D_model = 'Y', surpdead = 'Y', died = 'Y') %>%
    arrange(coral_numID, date)
  
  #filter to the date that the surprise-dead coral was first documented as dead and document it
  first.dead = first.dead.full %>%
    filter(grepl("\\Dead", status_field)) %>%
    group_by(coral_numID) %>%
    top_n(n=1, wt=desc(date)) %>%
    ungroup() %>%
    mutate(dead_date = date) %>% #create duplicate column to be used later, and indicate the first documented date of complete mortality
    arrange(coral_numID, date)
  
  first.dead.full = first.dead.full %>% left_join(
    first.dead %>% select(coral_numID, dead_date),
    by = c("coral_numID")
  ) %>%
    arrange(coral_numID, date)
  
  #backtrack infections
  for (i in 2:nrow(first.dead.full)) { #start at 2 to accommodate reverse indexing
    if (first.dead.full$date[i] == first.dead.full$dead_date[i]) { #when the date of death is the current date
      
      # Set backtracked %infected to the prior date
      first.dead.full$percinf[i-1] = 100 / as.numeric(difftime(first.dead.full$dead_date[i], first.dead.full$date[i-1], units = "days"))
      
      # Set 100% loss to SCTLD for current (dead) date
      first.dead.full$percloss[i] = 100
    }
  }
  
  #calculate accumulated tissue loss percentage
  first.dead.full = first.dead.full %>%
    group_by(coral_numID) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(cum_percloss = cumsum(replace_na(percloss, 0))) %>%
    mutate(cum_percloss = replace(cum_percloss, cum_percloss == 0, NA)) %>% #replace 0's pre-SCTLD with NA
    ungroup()
  
  #update backtracked SCTLD-infected status for the timepoint prior to the surprise-dead date
  first.dead.full = first.dead.full %>%
    mutate(
      status_model = ifelse(!is.na(percinf) & percinf > 0, 'backtracked_SCTLD', status_model)
    )
  
  #update original dataframe with surprise-dead coral mortality
  survey_trimmed = survey_trimmed %>%
    left_join(
      first.dead.full %>% select(coral_numID, date, D_model, status_model, surpdead, died, dead_date, percloss, cum_percloss, percinf),
      by = c("coral_numID", "date")) %>%
    mutate(
      D_model.x = coalesce(D_model.y, D_model.x),
      status_model.x = coalesce(status_model.y, status_model.x),
      percloss.x = coalesce(percloss.y, percloss.x),
      surpdead.x = coalesce(surpdead.y, surpdead.x),
      died.x = coalesce(died.y, died.x),
      cum_percloss.x = coalesce(cum_percloss.y, cum_percloss.x),
      percinf.x = coalesce(percinf.y, percinf.x)
    ) %>%
    rename(D_model = D_model.x,
           status_model = status_model.x, 
           percloss = percloss.x,
           surpdead = surpdead.x,
           died = died.x,
           cum_percloss = cum_percloss.x,
           percinf = percinf.x) %>%
    select(-D_model.y, -status_model.y, -percloss.y, -surpdead.y, -died.y, -cum_percloss.y, -percinf.y) %>%
    relocate(dead_date, .after = last_col()) %>%
    arrange(coral_numID, date)

  # calculate backtracked instantaneous (tissue loss / sloughing within 24 hours) infected tissue in surprise-dead corals 
  surprise.dead.infections = survey_trimmed %>%
    filter(coral_numID %in% surprise.dead$coral_numID) %>%
    mutate(
      inftiss = if_else(!is.na(percinf),
                        (percinf / 100) * starttiss,
                        NA_real_)
    ) %>%
    arrange(coral_numID, date)
  
  #update original dataframe with surprise-dead coral mortality
  survey_trimmed = survey_trimmed %>%
    left_join(
      surprise.dead.infections %>% select(coral_numID, date, inftiss),
      by = c("coral_numID", "date")) %>%
    mutate(
      inftiss.x = coalesce(inftiss.y, inftiss.x)
    ) %>%
    rename(inftiss = inftiss.x) %>%
    select(-inftiss.y) %>%
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
    # # i = 589 #single out the problem coral (2_p27_t2_s0_c1_DSTO). this was a patient zero
    # i = 5 #first infected coral in dataframe
    
    # Extract values from currIDsdates
    curr_values = observed.infected %>%
      slice(i) %>%
      mutate(
        coral_numID = as.numeric(coral_numID),
        ID = as.numeric(ID),
        date = as.character(date)
      )
    
    curr_coral_ID = curr_values$coral_numID
    curr_ID = curr_values$ID
    currdate = curr_values$date
    
    percloss = prograte %>%
      filter(coral_numID == curr_coral_ID) %>%
      select(all_of(curr_values$date)) %>%
      pull() %>%
      { if (length(.) == 0) NA_real_ else as.numeric(.) } #set as NA if no match in prograte dataframe
    
    IDslice = surveydiseased %>%
      filter(ID == curr_ID)
    starttiss = IDslice$starttiss

    # Only calculate if percloss is valid
    if(percloss > 0 & !is.na(percloss)){
      
      # Update percloss in observed.infected
      observed.infected[i, "percloss"] = percloss

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
      
      # Calculate #/days it took for the amount of loss observed to accumulate
      progdays = as.numeric(difftime(currdate, prevdate, units = "days"))
      observed.infected[i, "progdays"] = progdays
      
      # Calculate backtracked % of coral tissue that was infected, each day, from prior timepoint until current observation
      percinf = percloss / progdays
      observed.infected[i-1, "percinf"] = percinf
      
      # Populate backtracked amount of tissue (SA) infected, from percentage value (percinf)
      inftiss = (percinf / 100) * starttiss
      observed.infected[i-1, "inftiss"] = inftiss
      
      # Update backtracked infection status
      if (observed.infected[i-1, "status_model"] == "Healthy") {
        observed.infected[i-1, "status_model"] = "backtracked_SCTLD"
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
  surveydiseased = surveydiseased %>%
    group_by(coral_numID) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      cum_percloss = cumsum(replace_na(percloss, 0)),
      cum_percloss = if_else(cum_percloss == 0, NA_real_, cum_percloss)  # Replace 0s with NA
    ) %>%
    ungroup()
  
  #assign dead dates (first date when accumulated percloss is 100)
  first_death_dates = surveydiseased %>%
    filter(cum_percloss >= 100) %>%
    group_by(coral_numID) %>%
    summarize(dead_date = min(date)) %>%
    ungroup()
  
  surveydiseased = surveydiseased %>%
    left_join(first_death_dates, by = "coral_numID") %>%
    arrange(coral_numID, date)
  
  #assemble main dataframe with diseased tissue information
  survey_tissue = rbind(survey_trimmed, surveydiseased) %>%
    arrange(coral_numID, date)
  
  #remove, rearrange, and rename columns for clarity; update disease statuses for model
  # NOTE - 3_p47_t1_s0_c4_MCAV, 3_p47_t6_s0_c17_CNAT, and 2_p28_t4_s0_c5_SSID are useful corals to check to make sure
  #         intermittent/lapsed/recurrent infection statuses are populated correctly
  # NOTE - there may be timepoint-specific bleaching (thermal stress) data for each colony somewhere. can come up later in results revisions
  survey_tissue = survey_tissue %>%
    rename(BL_field = tot_stressed, BL_D_model = tot_both) %>%
    group_by(coral_numID) %>%
    mutate(
      dead_date = if(any(!is.na(dead_date))) first(na.omit(dead_date)) else dead_date, #populate 'dead_date' before 10-30-2018
      died = if(any(!is.na(dead_date))) "Y" else died, #mark coral as infected
      status_model = if_else(!is.na(cum_percloss) & cum_percloss == 100, "Dead", status_model), # set status_model to 'Dead' if cum_percloss is 100 (redundant, for peace of mind)
      status_model = if_else( #set 'status_model' to 'SCTLD' if there is active infection (redundant, for peace of mind)
        status_model != "backtracked_SCTLD" & 
          status_model != "SCTLD" & 
          !is.na(cum_percloss) & 
          cum_percloss > 0 & 
          cum_percloss < 100, 
        "SCTLD", 
        status_model
      ), 
      status_model = if_else( #corals not actively infected (but marked 'Unknown' in field) are marked healthy
        status_model == "Unknown", 
        "Healthy", 
        status_model
      ),
      # status_model = if_else( #corals marked with % loss in field, but no *active* infection as defined by model, marked as 'halted_SCTLD'
      #   !is.na(percloss) &
      #     is.na(percinf) &
      #     cum_percloss != 100,
      #   "Halted_SCTLD", 
      #   status_model
      # )
      status_model = if_else( #corals marked with % loss in field, but no *active* infection as defined by model, marked as 'halted_SCTLD'
        !is.na(cum_percloss) &
          cum_percloss != 100 &
          is.na(percinf),
        "Halted_SCTLD", 
        status_model
      )
    ) %>%
    ungroup() %>%
    mutate(
      BL_field = case_when(
        BL_field == "S" ~ "Y",
        BL_field == "NS" ~ "N"
      ),
      BL_D_field = case_when( #mark coral as either 'Y' (cormorbid for disease & bleaching across timepoints) or 'N' (no comorbidity)
        BL_field == "Y" & D_field == "Y" ~ "Y", 
        BL_field == "N" | D_field == "N" ~ "N"
      ),
      BL_D_model = case_when( #mark coral as either 'Y' (cormorbid for disease & bleaching across timepoints) or 'N' (no comorbidity)
        BL_field == "Y" & D_model == "Y" ~ "Y", 
        BL_field == "N" | D_model == "N" ~ "N"
      )
    ) %>%
    select(1:which(names(.) == "D_model"), BL_D_field, everything())  #reorder to place BL_D_field after D_model
  #
  ### NOTE - code above assumes that tissue loss for 'unknown' reason is due to SCTLD (filter by value = 'Unknown' and look at percloss that is nonzero)
  ###           - there are 17 corals which experienced this, mostly LS/MS corals
  ###           - I checked all of them, and each ended up also having SCTLD in near or immediately adjacent timepoints - fair assumption that all loss was SCTLD
  ###           - Solenastrea might be a bit more susceptible than I'd thought. a few cases of total mortality in 3-4 weeks. also seeing evidence that 
  ###           -   PCLI can sustain long infections; might match up well with susceptibility knowledge from VI
  #
  # 'CORALS DOCUMENTED AS DISEASED'

  # REMOVED / REMAINING TISSUE CALCULATION
  #
  survey_tissue = survey_tissue %>%
    mutate(
      cum_tissloss = (cum_percloss/100) * starttiss,
      remaintiss = starttiss - ifelse(is.na(inftiss), 0, inftiss) - ifelse(is.na(cum_tissloss), 0, cum_tissloss), #subtract non-NA infected and accumulated dead tissue from starting tissue to get remaining tissue
      remaintiss = coalesce(remaintiss, starttiss)  #anywhere 'remaintiss' is NA (no tissue loss has begun), revert to starting tissue value
    )
  #
  # REMOVED / REMAINING TISSUE CALCULATION
  
  # PATIENT ZERO SPECIAL CASES
  #
  #true patient zero is 2_p27_t2_s0_c1_DSTO on 10-30-2018. assume its 90% tissue loss occurred over 14 days since we don't have info
  #   except 3 months prior when colony was healthy
  progdays.new = 14
  reference_date = ymd('2018-10-30')
  new_date = reference_date - days(progdays.new) #backtrack 14 days from surprise-dead day of 10-30-2018

  # Calculate percinf (90% percloss over the time difference)
  percinf.new = 90 / progdays.new
  
  # Extract the existing row for the offshore coral with coral_longID '2_p27_t2_s0_c1_DSTO'
  patientzero.offshore.row = survey_tissue %>%
    filter(coral_long_ID == '2_p27_t2_s0_c1_DSTO' & date == reference_date) #%>%
    # select(-date) # Remove the date column to replace it with the new date
  
  # Calculate inftiss
  inftiss.new = (percinf.new / 100) * patientzero.offshore.row$starttiss
  
  #shrink the patient zero coral's infection surface area considerably, to approximate only a few polyps being infected in the very
  #   earliest stages of gross lesion presentation. desirable for priming epidemic model
  #     NOTE - the minimizing scalar is arbitrary, and can be tweaked to whatever works well for initializing the epidemic model runs
  inftiss.new = inftiss.new / polyp_SA.minimizer.offshore
  
  # Create the new row with the calculated values
  new_row = patientzero.offshore.row %>%
    mutate(
      date = new_date,
      status_field = NA_character_,
      status_model = 'patientzero_SCTLD',
      percloss = NA_real_,
      cum_percloss = NA_real_,
      progdays = NA_real_,
      # progdays = progdays.new,
      percinf = percinf.new,
      inftiss = inftiss.new,
      remaintiss = starttiss - inftiss,
      cum_tissloss = NA_real_,
      ID_SIR = patientzero.offshore.row$ID,
      ID = NA  # Set ID to NA for the new row
    )
  
  # Insert the new row into survey_tissue and update `ID_SIR` column to be identical to `ID` and increment the values after the new row
  survey_tissue = survey_tissue %>%
    mutate(ID_SIR = NA) %>%
    add_row(new_row) %>%
    arrange(coral_numID, date) %>%
    mutate(
      ID_SIR = ifelse(is.na(ID_SIR), ID, ID_SIR),  # Keep the correct `ID_SIR`
      ID_SIR = ifelse(
        ID >= new_row$ID_SIR & !is.na(ID), ID_SIR + 1, ID_SIR  # Increment after new row
      ),
      progdays = ifelse(coral_long_ID == '2_p27_t2_s0_c1_DSTO' & #update forward-tracked days of progression since backtracked prior date
                          date == '2018-10-30',
                        progdays.new,
                        progdays)
    ) %>%
    select(1:which(names(.) == "ID"), ID_SIR, everything())
  
  #similarly update site-specific patient zero corals' infected surface area for the two sites (midchannel and nearshore) which contracted
  #   SCTLD later on than offshore
  min_inftiss_midchannel <- survey_tissue %>%
    filter(coral_long_ID == '1_p25_t2_s0_c22_DSTO') %>%
    filter(!is.na(inftiss) & inftiss > 0) %>%
    arrange(date) %>%  # Sort by date to find the earliest entry
    # # optionally, just take the first infection day instead of a minimum of all days
    # slice(1) %>%  # Select the first (earliest) row #this was deprecated in favor of pulling the most minimal inftiss value in that patient zero's time series (which sometimes is not the initial value)
    # pull(ID_SIR)
    filter(inftiss == min(inftiss)) %>%  # Find the row with the minimum inftiss
    pull(inftiss)  # Get the minimum inftiss value from the earliest row
  
  # Get the earliest date that has valid inftiss values
  earliest_date_with_valid_inftiss <- survey_tissue %>%
    filter(coral_long_ID == '1_p25_t2_s0_c22_DSTO' & !is.na(inftiss) & inftiss > 0) %>%
    summarize(earliest_date = min(date)) %>%
    pull(earliest_date)
  
  # Update only the earliest date row for the specific coral_long_ID
  survey_tissue <- survey_tissue %>%
    mutate(inftiss = if_else(coral_long_ID == '1_p25_t2_s0_c22_DSTO' & date == earliest_date_with_valid_inftiss,
                             min_inftiss_midchannel / polyp_SA.minimizer.midchannel,
                             inftiss))  
  
  #the same as above, for nearshore site
  min_inftiss_nearshore <- survey_tissue %>%
    filter(coral_long_ID %in% c('3_p47_t3_s0_c8_PSTR', '3_p47_t4_s5_c15_PSTR')) %>%  # two patient zero corals
    filter(!is.na(inftiss) & inftiss > 0) %>%  # Keep rows with non-NA, non-zero 'inftiss'
    group_by(coral_long_ID) %>%  # Group by coral_long_ID
    summarise(min_inftiss = min(inftiss), .groups = 'drop')  # Get the minimum inftiss for each group
  
  # Get the earliest date that has valid inftiss values
  earliest_date_with_valid_inftiss <- survey_tissue %>%
    filter(coral_long_ID %in% c('3_p47_t3_s0_c8_PSTR', '3_p47_t4_s5_c15_PSTR') & !is.na(inftiss) & inftiss > 0) %>%
    group_by(coral_long_ID) %>%
    summarise(earliest_date = min(date), .groups = 'drop')
  
  # Update only the earliest date row for the specific coral_long_ID
  survey_tissue <- survey_tissue %>%
    left_join(min_inftiss_nearshore, by = "coral_long_ID") %>%  # Join to get the min_inftiss for each coral
    left_join(earliest_date_with_valid_inftiss, by = "coral_long_ID") %>%  # Join to get the earliest_date for each coral
    mutate(inftiss = case_when(
      coral_long_ID == '3_p47_t3_s0_c8_PSTR' & date == earliest_date ~ min_inftiss / polyp_SA.minimizer.nearshore,
      coral_long_ID == '3_p47_t4_s5_c15_PSTR' & date == earliest_date ~ min_inftiss / polyp_SA.minimizer.nearshore,
      TRUE ~ inftiss  # Keep the original inftiss if conditions don't match
    )) %>%
    select(-min_inftiss, -earliest_date)  # Remove the joined columns after use
  #
  # PATIENT ZERO SPECIAL CASES
  
  # AGGREGATE DATA INTO SUMMARY TABLES
  #
  summary_grouped = survey_tissue %>%
    # filter(!grepl('Unaffected', susc)) %>% #remove presumed SCTLD-unaffected corals, mainly PAST
    mutate(TP = sprintf("%02d", dense_rank(date))) %>%
    group_by(site, date, TP, susc) %>%
    mutate(
      is_susceptible = ifelse(is.na(inftiss) & (is.na(cum_percloss) | cum_percloss != 100), 1, 0),  # Temporary column for susceptible corals
      is_infected = ifelse(inftiss > 0, 1, 0),  # Temporary column for infected corals
      is_dead = ifelse(cum_percloss == 100, 1, 0)  # Temporary column for dead corals
    ) %>%
    summarise(
      tottiss = sum(remaintiss, inftiss, cum_tissloss, na.rm = TRUE),
      sustiss = sum(remaintiss, na.rm = TRUE),
      inftiss = sum(inftiss, na.rm = TRUE),
      deadtiss = sum(cum_tissloss, na.rm = TRUE),
      totnum = n(),  # Count of all entries in the group
      susnum = sum(is_susceptible, na.rm = TRUE),  # Count susceptible corals
      infnum = sum(is_infected, na.rm = TRUE),  # Count infected corals
      deadnum = sum(is_dead, na.rm = TRUE),  # Count dead corals
      .groups = 'drop'  # Ungroup after summarizing
    ) %>%
    complete(site, TP, susc, fill = list(tottiss = NA, sustiss = NA, inftiss = NA, deadtiss = NA, totnum = NA, susnum = NA,
                                         infnum = NA, deadnum = NA)) %>%
    group_by(site, susc) %>% 
    mutate( #fill NA values for backtracked patient zero timepoint using subsequent timepoint. requires further cleaning, as below
      tottiss = ifelse(TP == "02" & is.na(tottiss), lead(tottiss, default = NA), tottiss),
      sustiss = ifelse(TP == "02" & is.na(sustiss), lead(sustiss, default = NA), sustiss),
      inftiss = ifelse(TP == "02" & is.na(inftiss), lead(inftiss, default = NA), inftiss),
      deadtiss = ifelse(TP == "02" & is.na(deadtiss), lead(deadtiss, default = NA), deadtiss),
      totnum = ifelse(TP == "02" & is.na(totnum), lead(totnum, default = NA), totnum),
      susnum = ifelse(TP == "02" & is.na(susnum), lead(susnum, default = NA), susnum),
      infnum = ifelse(TP == "02" & is.na(infnum), lead(infnum, default = NA), infnum),
      deadnum = ifelse(TP == "02" & is.na(deadnum), lead(deadnum, default = NA), deadnum)
    ) %>%
    ungroup() %>%
    mutate( #properly update 'tottiss' and 'sustiss' values for patient zero timepoint
      #handle high-susceptibles (which contain patient zero)
      tottiss = ifelse(site == "off" & TP == "02" & susc == "high",
                       tottiss[site == "off" & TP == "03" & susc == "high"], 
                       tottiss),
      sustiss = ifelse(site == "off" & TP == "02" & susc == "high",
                       tottiss[site == "off" & TP == "03" & susc == "high"] - inftiss[site == "off" & TP == "02" & susc == "high"], 
                       sustiss),
      totnum = ifelse(site == "off" & TP == "02" & susc == "high",
                      totnum[site == "off" & TP == "03" & susc == "high"],
                      totnum),
      susnum = ifelse(site == "off" & TP == "02" & susc == "high",
                      totnum[site == "off" & TP == "03" & susc == "high"] - infnum[site == "off" & TP == "02" & susc == "high"],
                      susnum),
      #handle moderate-susceptibles (which do not contain patient zero but do get infected in the subsequent timepoint)
      sustiss = ifelse(site == "off" & TP == "02" & susc == "moderate", 
                       tottiss, 
                       sustiss),
      inftiss = ifelse(site == "off" & TP == "02" & susc == "moderate", 
                       0, 
                       inftiss),
      susnum = ifelse(site == "off" & TP == "02" & susc == "moderate", 
                      totnum, 
                      susnum),
      infnum = ifelse(site == "off" & TP == "02" & susc == "moderate", 
                      0, 
                      infnum)
    )
  # NOTE - the above '02' and '03 hard-coding is not ideal, and this is an area to look to if there are bugs in the future

  #populate dates for patient zero timepoint
  summary_grouped$date[summary_grouped$TP == "02"] = as.POSIXct(toString(new_date))

  #include DHW data, with days that match the surveying dates (this is important for meshing with the SIR model downstream)
  days.DHW <- survey_tissue %>%
    filter(!is.na(inftiss) & inftiss > 0) %>%
    arrange(date) %>%
    slice(1) %>%
    pull(date) %>%
    { tibble(date = seq(., max(survey_tissue$date, na.rm = TRUE), by = "day"),
             day = seq(0, as.numeric(max(survey_tissue$date, na.rm = TRUE) - .))) } %>%
    mutate(
      date = as.Date(date)
    )
  
  DHW.CRW = DHW.CRW %>%
    left_join(days.DHW, by = 'date') %>%
    select(date, day, everything())

  summary_grouped = summary_grouped %>%
    left_join(DHW.CRW, by = "date")
  
  #further summary
  # NOTE - infections are always zero for the final timepoint of the study. this is because there is no further timepoint with which to
  #         backtrack infections to. it is a limitation of the approach, but I think the best choice for the available data
  summary = summary_grouped %>%
    group_by(site, date, TP) %>%
    summarize(
      
      #low-susceptibility
      low.sustiss = sum(sustiss[susc=="low"]),
      low.inftiss = sum(inftiss[susc=="low"]),
      low.deadtiss = sum(deadtiss[susc=="low"]),
      
      low.susnum = sum(susnum[susc=="low"]),
      low.infnum = sum(infnum[susc=="low"]),
      low.deadnum = sum(deadnum[susc=="low"]),
      
      #moderate-susceptibility
      moderate.sustiss = sum(sustiss[susc=="moderate"]),
      moderate.inftiss = sum(inftiss[susc=="moderate"]),
      moderate.deadtiss = sum(deadtiss[susc=="moderate"]),
      
      moderate.susnum = sum(susnum[susc=="moderate"]),
      moderate.infnum = sum(infnum[susc=="moderate"]),
      moderate.deadnum = sum(deadnum[susc=="moderate"]),
      
      #high-susceptibility
      high.sustiss = sum(sustiss[susc=="high"]),
      high.inftiss = sum(inftiss[susc=="high"]),
      high.deadtiss = sum(deadtiss[susc=="high"]),
      
      high.susnum = sum(susnum[susc=="high"]),
      high.infnum = sum(infnum[susc=="high"]),
      high.deadnum = sum(deadnum[susc=="high"]),
      
      #unaffected
      unaffected.tiss = sum(sustiss[susc=="Unaffected"]),
      unaffected.num = sum(susnum[susc=="Unaffected"]),
      
      #total
      tot.sustiss = low.sustiss + moderate.sustiss + high.sustiss,
      tot.inftiss = low.inftiss + moderate.inftiss + high.inftiss,
      tot.deadtiss = low.deadtiss + moderate.deadtiss + high.deadtiss,
      
      tot.susnum = low.susnum + moderate.susnum + high.susnum,
      tot.infnum = low.infnum + moderate.infnum + high.infnum,
      tot.deadnum = low.deadnum + moderate.deadnum + high.deadnum,
      
      SA.with.unaffected = low.sustiss + moderate.sustiss + high.sustiss + unaffected.tiss,
      count.with.unaffected = low.susnum + moderate.susnum + high.susnum + unaffected.num,
      
      #DHW data
      SST_MIN = first(SST_MIN),
      SST_MAX = first(SST_MAX),
      SST.90th_HS = first(SST.90th_HS),
      SSTA.90th_HS = first(SSTA.90th_HS),
      X90th_HS.0 = first(X90th_HS.0),
      DHW_from_90th_HS.1 = first(DHW_from_90th_HS.1),
      BAA_7day_max = first(BAA_7day_max)
    )
  
  #easily legible table for SCTLD-affected corals
  tot.summary = summary %>%
    select(site, date, TP, tot.sustiss, tot.inftiss, tot.deadtiss, tot.susnum, tot.infnum, tot.deadnum, SA.with.unaffected, count.with.unaffected) %>%
    mutate(
      perc.prevalence.tiss = (tot.inftiss / tot.sustiss) * 100,
      perc.prevalence.num = (tot.infnum / tot.susnum) * 100,
      sustiss.host.ratio = tot.sustiss / tot.susnum,
      inftiss.host.ratio = tot.inftiss / tot.infnum
    )
  #
  # AGGREGATE DATA INTO SUMMARY TABLES
  
  # TRACK DAYS SINCE FIRST INFECTION(S)
  #
  # Identify the first timepoint and first infection date for each site
  summary = summary %>%
    group_by(site) %>%
    mutate( 
      first.infTP = TP[which(tot.infnum > 0)[1]], # Find the first 'TP' (timepoint) where 'tot.infnum' is non-zero
      first.infdate.site = min(date[tot.infnum > 0], na.rm = TRUE), # Get the date of the first infection occurrence for this site
      days.inf.site = as.numeric(difftime(date, first.infdate.site, units = "days")) # Calculate days since the first infection for this specific site (days.inf.site)
    ) %>%
    ungroup()
  
  # Calculate the earliest infection date across all sites
  global_first_infection_date = summary %>%
    filter(tot.infnum > 0) %>%
    summarise(first_infection = min(date, na.rm = TRUE)) %>%
    pull(first_infection)
  
  # Calculate the 'days.inf' column, which tracks days since the first infection across all sites
  summary = summary %>%
    mutate(
      # Calculate the number of days since the very first infection across any site
      days.survey = as.numeric(difftime(date, global_first_infection_date, units = "days"))
    )
  
  # Clean up: round both 'days.inf' and 'days.inf.site', and set any negative values to NA
  summary = summary %>%
    mutate(
      days.survey = round(days.survey),
      days.inf.site = round(days.inf.site),
      # Optional: Replace negative values with NA if necessary
      days.survey = ifelse(days.survey < 0, NA, days.survey),
      days.inf.site = ifelse(days.inf.site < 0, NA, days.inf.site)
    )
  
  # NOTE - in the above code, corals marked as SCTLD-affected in the field but with 0% tissue loss / active infection are considered as
  #         Halted_SCTLD in the context of our model. it may be worth considering that these are actually still "infected" but dormant,
  #         or otherwise not entirely in either a susceptible or infected compartment. also, due to the nature of coral diseases, 
  #         especially SCTLD, a coral is probably always still "susceptible" even when it is infected...perhaps another argument for
  #         using a tissue model. but it will be interesting to plot the rise and fall of susceptible hosts, even as susceptible
  #         tissue SA does nothing but stagnate or fall (by necessity and design)
  #
  # TRACK DAYS SINCE FIRST INFECTION(S)

  ### PLOTTING PREPARATIONS
  #
  # Get the list of infected counts for all sites
  sample_day_list = summary %>%
    group_by(site) %>%
    summarize(inf_count = list(tot.infnum), .groups = 'drop')
  
  # Create a list to store the indices for each site
  indices__list = sample_day_list %>%
    mutate(indices = map(inf_count, ~ which(.x > 0))) %>%
    ungroup()
  
  #ensure the final timepoint is included, even though the model assumes no infections there (to account for newly removed tissue)
  indices__list = indices__list %>%
    mutate(indices = map2(inf_count, indices, ~ {
      first_nonzero = min(which(.x > 0))  # Find the index of the first non-zero
      subsequent_zeros = which(.x == 0 & seq_along(.x) >= first_nonzero)  # Keep zeros after the first non-zero
      full_seq(c(.y, subsequent_zeros), 1)  # Extend the sequence
    }))
  
  # Merge back into summary to prepare for observation creation
  obs_summary = summary %>%
    left_join(indices__list %>% select(site, indices), by = "site") %>%
    group_by(site) %>%
    mutate(valid_rows = list(slice(., indices[[1]]))) %>%
    ungroup()
  
  # Extract the 'Susceptible' values at TP == 1 for each Category
  susceptible_ref = obs_summary %>%
    filter(TP == '01') %>%
    pivot_longer(cols = c(low.sustiss, moderate.sustiss, high.sustiss, tot.sustiss, SA.with.unaffected),
                 names_to = "Category",
                 values_to = "tissue_ref") %>%
    pivot_longer(cols = c(low.susnum, moderate.susnum, high.susnum, tot.susnum, count.with.unaffected),
                 names_to = "Category_count",
                 values_to = "count_ref") %>%
    mutate(
      Category = case_when(
        str_detect(Category, "low") ~ "Low",
        str_detect(Category, "moderate") ~ "Moderate",
        str_detect(Category, "high") ~ "High",
        str_detect(Category, "tot") ~ "Total",
        str_detect(Category, "unaffected") ~ "Total_with_Unaffected"
      ),
      Category_count = case_when(
        str_detect(Category_count, "low") ~ "Low",
        str_detect(Category_count, "moderate") ~ "Moderate",
        str_detect(Category_count, "high") ~ "High",
        str_detect(Category_count, "tot") ~ "Total",
        str_detect(Category_count, "unaffected") ~ "Total_with_Unaffected"
      )
    ) %>%
    filter(Category == Category_count) %>%
    select(Site = site, Category, tissue_ref, count_ref)

  ### ESTIMATION OF SITE AND GROUP-LEVEL PERCENT COVER
  #
  # here, the CPCe-derived observations of coral cover from Williams et al. (2021) are used to estimate percent cover
  #   - Williams' cover values included all corals at the site, including those unaffected by SCTLD. it also is not provided by species,
  #     but at a site-level
  #   - so, we create a relationship (simple ratio), for each site, between site-level 2-dimensional coral cover and site-level 
  #       3-dimensional surface area. this ratio is assumed to be applicable to, and identical for, any given species of coral
  #   - this assumption is not taking into account the differing morphologies and sizes inherent to each coral among susceptibility groups,
  #       which in reality may modify that SA:cover ratio, but we are okay with that since the differences may not be that large
  #   - with that assumption in mind, the ratio (which is pretty similar but varies a bit for each site) is applied to
  #       susceptibility-group-level 3-D surface area to produce an estimated group-level cover for SCTLD-susceptible corals
  #   - this same process is applied at the site level, excluding only the SCTLD-unaffected corals
  site_ref = susceptible_ref %>%
    filter(Category == "Total") %>%
    select(Site, tissue_ref) %>%
    rename(N.site = tissue_ref)
  
  site_ref_unaffected = susceptible_ref %>%
    filter(Category == "Total_with_Unaffected") %>%
    select(Site, tissue_ref) %>%
    rename(N.site.with.unaffected = tissue_ref) 
  
  # Join the cover data into susceptible_ref
  susceptible_ref <- susceptible_ref %>%
    left_join(cover, by = "Site") %>%
    rename(site.cover.Williams = mean.percent.cover) #%>%

  susceptible_ref = susceptible_ref %>%
    left_join(site_ref, by = "Site") %>%
    left_join(site_ref_unaffected, by = "Site") %>%
    mutate(rectangle.area.site = 200) %>% #each site had two 10 x 10 m quadrats (Williams et al. 2021). so, 200 m2 2-dimensional area
    group_by(Site) %>%
    mutate(SA.cover.ratio = tissue_ref[Category == "Total_with_Unaffected"] / site.cover.Williams) %>% #SA is the 3-dimensional GAM-estimated surface area of all corals in each site
    # mutate(SA.cover.ratio.prior = 3.0) %>% #before, I was simply using a crude 3:1 ratio which happened to approximate cover well for all three sites compared to Williams 2021
    ungroup() %>%
    mutate(
      cover.site = N.site / SA.cover.ratio / 100,
      # cover.site.prior = N.site / rectangle.area.site / SA.cover.ratio.prior,
      cover_ref = tissue_ref / SA.cover.ratio / 100,
      # cover_ref.prior = tissue_ref / rectangle.area.site / SA.cover.ratio.prior
    )
  #
  ### ESTIMATION OF SITE AND GROUP-LEVEL PERCENT COVER
  
  ### LAST STEPS TO PREP FOR PLOTTING AND DOWNSTREAM MODELING SCRIPTS
  #
  # Create a function to generate the observations
  create_observations = function(data, comp, tissue_col, count_col, cat) {
    data %>%
      transmute(
        TP = TP,
        days.inf.site = days.inf.site,
        days.survey = days.survey,
        tissue = !!sym(tissue_col), #bang bang operator to convert wide format to long
        count = !!sym(count_col),
        Category = cat,
        Compartment = comp,
        Site = site,
        SST_MIN = SST_MIN,
        SST_MAX = SST_MAX,
        SST_90th_HS = SST.90th_HS,
        SSTA_90th_HS = SSTA.90th_HS,
        X90th_HS_0 = X90th_HS.0,
        DHW_from_90th_HS_1 = DHW_from_90th_HS.1,
        BAA_7day_max = BAA_7day_max
      )
  }
  
  # Create data frames for all categories, including count columns
  obs_sus = bind_rows(
    create_observations(obs_summary, "Susceptible", "low.sustiss", "low.susnum", "Low"),
    create_observations(obs_summary, "Susceptible", "moderate.sustiss", "moderate.susnum", "Moderate"),
    create_observations(obs_summary, "Susceptible", "high.sustiss", "high.susnum", "High"),
    create_observations(obs_summary, "Susceptible", "tot.sustiss", "tot.susnum", "Total")
  )
  
  obs_inf = bind_rows(
    create_observations(obs_summary, "Infected", "low.inftiss", "low.infnum", "Low"),
    create_observations(obs_summary, "Infected", "moderate.inftiss", "moderate.infnum", "Moderate"),
    create_observations(obs_summary, "Infected", "high.inftiss", "high.infnum", "High"),
    create_observations(obs_summary, "Infected", "tot.inftiss", "tot.infnum", "Total")
  )
  
  obs_dead = bind_rows(
    create_observations(obs_summary, "Dead", "low.deadtiss", "low.deadnum", "Low"),
    create_observations(obs_summary, "Dead", "moderate.deadtiss", "moderate.deadnum", "Moderate"),
    create_observations(obs_summary, "Dead", "high.deadtiss", "high.deadnum", "High"),
    create_observations(obs_summary, "Dead", "tot.deadtiss", "tot.deadnum", "Total")
  )
  
  # Combine all observations into a single data frame
  obs = bind_rows(obs_sus, obs_inf, obs_dead)
  
  # Add scaled tissue and count columns to the existing 'obs' dataframe
  obs = obs %>%
    left_join(susceptible_ref, by = c("Site", "Category")) %>%
    mutate(
      tissue.scaled = tissue / tissue_ref,  # Scale by the reference tissue
      count.scaled = count / count_ref,     # Scale by the reference count
    ) %>%
    select(-tissue_ref, -count_ref)  # Remove the reference columns after scaling
  
  # Factorize and add dates
  summary_unique = summary %>%
    group_by(TP) %>%
    summarize(date = first(date), .groups = "drop")
  obs = obs %>%
    mutate(across(c(TP, Category, Compartment, Site), as.factor)) %>%
    left_join(summary_unique, by = "TP") %>%
    select(TP, date, everything())  # This will put 'date' as the second column
  
  # Identify patient zero tissue and count in 'obs' (Compartment == 'Infected')
  patient_zero = obs %>%
    filter(Compartment == 'Infected' & !is.na(tissue) & tissue > 0 & !is.na(count) & count > 0) %>%
    group_by(Site, Category) %>%
    arrange(date) %>%
    summarize(
      patientzero_tissue = first(tissue),
      patientzero_count = first(count),
      .groups = 'drop'
    ) %>%
    mutate(
      patientzero_tissue = ifelse(Category %in% c('Low', 'Moderate'), NA, patientzero_tissue),
      patientzero_count = ifelse(Category %in% c('Low', 'Moderate'), NA, patientzero_count)
    )
  
  # Join patient zero information to 'susceptible_ref'
  susceptible_ref = susceptible_ref %>%
    left_join(patient_zero, by = c("Site", "Category"))
  
  obs = obs %>%
    mutate(
      Site = case_when( # Map site names to new values
        Site == "off" ~ "Offshore",
        Site == "mid" ~ "Midchannel",
        Site == "near" ~ "Nearshore",
        TRUE ~ Site
      )
    )
  #
  ### LAST STEPS TO PREP FOR PLOTTING AND DOWNSTREAM MODELING SCRIPTS
  
  # ### CREATE ERROR EVALUATION DATAFRAME
  # # Define combinations of sites, hosts, and model types
  # sites <- c("off", "mid", "near")
  # hosts <- c("Single-host", "Multi-host")
  # types <- c("Fitted", "DHW", "Projected")
  # 
  # # Create an empty Error_eval dataframe with all combinations
  # error_eval <- expand.grid(site = sites, 
  #                           host = hosts, 
  #                           type = types, 
  #                           SSR = NA, 
  #                           TSS = NA, 
  #                           R_squared = NA, 
  #                           stringsAsFactors = FALSE)
  # ### CREATE ERROR EVALUATION DATAFRAME
  
  library(tibble)
  
  # Define combinations of sites, hosts, and model types
  sites <- c("off", "mid", "near")
  hosts <- c("Single-host", "Multi-host")
  types <- c("Fitted", "DHW", "Projected")
  waves = c("Pre-heat", "Post-heat", "Both", "Full")
  
  # Create an empty tibble with all combinations and list columns for vectors
  error_eval <- expand.grid(site = sites, 
                            host = hosts, 
                            type = types, 
                            wave = waves,
                            stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(
      SST_threshold = NA,
      DHW_threshold = NA,
      date_thresh = as.Date(NA),
      # SSR_optimizer = NA_real_, #residual sum of squares used for fitting optimization
      SSR = NA_real_, #residual sum of squared residuals
      TSS = NA_real_, #total sum of squares for entire outbreak
      R_squared = NA_real_, #R-squared (1 - SSR / TSS) for entire outbreak
      sim_inf = vector("list", n()),  # List-column for simulated infected
      sim_rem = vector("list", n()),  # List-column for simulated removed
      obs_inf = vector("list", n()),  # List-column for observed infected
      obs_rem = vector("list", n())   # List-column for observed removed
    )
  
  # #save workspace for use in subsequent scripts
  # save.image(file = here("output", "data_processing_workspace.RData"))
  