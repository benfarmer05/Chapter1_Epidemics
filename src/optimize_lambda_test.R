 

 # NOTE - run this while working within plots_basic.R

  
  library(dplyr)
  
  # Define the search space
  lambda_values <- seq(0, 10, length.out = 3000)
  
  best_r_squared <- -Inf
  best_lambda <- NA
  
  for (lambda in lambda_values) {
    
    # Compute offset
    offset = 1 - 1 / (1 + exp(-lambda * 1.0))
    
    # Run the model
    output.basic.offshore.transfer = data.frame(ode(c(S = S.offshore, I = I.offshore, R = R.offshore),
                                                    days.model.offshore, SIR, c(b = beta.nearshore, g = gamma.nearshore,
                                                                                N = N.offshore,
                                                                                C = cover.offshore,
                                                                                l = lambda)))
    
    
    # Compute R-squared
    sim.rem.total = output.basic.offshore.transfer[which(output.basic.offshore.transfer$time %in% days.obs),
                                                     which(colnames(output.basic.offshore.transfer) %in% c('R'))]
    
    if (length(obs.rem.total) > length(sim.rem.total)) {
      obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
    }
    
    diff.rem.total = (sim.rem.total - obs.rem.total)
    sum_diff.total = sum(diff.rem.total^2)
    mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
    tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
    r_squared.near.to.off.basic = 1 - (sum_diff.total / tss_rem.total)
    
    # Check if this is the best so far
    if (r_squared.near.to.off.basic > best_r_squared) {
      best_r_squared <- r_squared.near.to.off.basic
      best_lambda <- lambda
    }
  }
  
  # Print the best parameters
  print(best_lambda)
  print(best_r_squared)




  # NOTE - run this while working within plots_multi.R
  
  library(dplyr)
  
  # Define the search space
  # lambda_values <- seq(0, 50, length.out = 10)  # Small grid for speed
  # lambda_values <- seq(0, 20, length.out = 10)  # Small grid for speed
  lambda_values <- (seq(0, 1, length.out = 20))^4 * 20 #power-based sequence with sharper resolution near zero
  
  best_r_squared <- -Inf
  best_lambda <- c(lambda.LS = NA, lambda.MS = NA, lambda.HS = NA)
  
  for (lambda.LS in lambda_values) {
    for (lambda.MS in lambda_values) {
      for (lambda.HS in lambda_values) {
  
        # Compute offsets
        offset.LS = 1 - 1 / (1 + exp(-lambda.LS * 1.0))
        offset.MS = 1 - 1 / (1 + exp(-lambda.MS * 1.0))
        offset.HS = 1 - 1 / (1 + exp(-lambda.HS * 1.0))
  
        # Run the model
        output.near.to.off.multi = data.frame(ode(c(S.LS = S.LS.offshore, I.LS = I.LS.offshore, R.LS = R.LS.offshore,
                                                    S.MS = S.MS.offshore, I.MS = I.MS.offshore, R.MS = R.MS.offshore,
                                                    S.HS = S.HS.offshore, I.HS = I.HS.offshore, R.HS = R.HS.offshore),
                                                  days.model.offshore, SIR.multi, c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
                                                                                    b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
                                                                                    b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
                                                                                    N.LS = N.LS.offshore, N.MS = N.MS.offshore, N.HS = N.HS.offshore,
                                                                                    C = cover.offshore,
                                                                                    C.LS = cover.offshore.LS, C.MS = cover.offshore.MS, C.HS = cover.offshore.HS,
                                                                                    l = lambda)))
  
        # Compute R-squared
        sim.rem.total = rowSums(output.near.to.off.multi[which(output.near.to.off.multi$time %in% days.obs),
                                                         which(colnames(output.near.to.off.multi) %in% c('R.LS', 'R.MS', 'R.HS'))],
                                na.rm = TRUE)
  
        if (length(obs.rem.total) > length(sim.rem.total)) {
          obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
        }
  
        diff.rem.total = (sim.rem.total - obs.rem.total)
        sum_diff.total = sum(diff.rem.total^2)
        mean_obs_rem.total = mean(obs.rem.total, na.rm = TRUE)
        tss_rem.total = sum((obs.rem.total - mean_obs_rem.total)^2)
        r_squared.near.to.off.multi = 1 - (sum_diff.total / tss_rem.total)
  
        # Check if this is the best so far
        if (r_squared.near.to.off.multi > best_r_squared) {
          best_r_squared <- r_squared.near.to.off.multi
          best_lambda <- c(lambda.LS, lambda.MS, lambda.HS)
        }
      }
    }
  }
  
  # Print the best parameters
  print(best_lambda)
  print(best_r_squared)
  
  # 
  # # version that is meant to show progress but doesn't seem to work as well
  # library(dplyr)
  # library(pbapply)  # Progress bar for tracking
  # 
  # # Define parameter grid
  # lambda_values <- seq(0, 20, length.out = 10)  # Adjust resolution for speed
  # # lambda_values <- seq(0, 50, length.out = 10)  # Adjust resolution for speed
  # param_grid <- expand.grid(lambda.LS = lambda_values,
  #                           lambda.MS = lambda_values,
  #                           lambda.HS = lambda_values)
  # 
  # # Function to evaluate a given set of lambda values
  # evaluate_lambda <- function(params) {
  #   lambda.LS <- params[1]
  #   lambda.MS <- params[2]
  #   lambda.HS <- params[3]
  # 
  #   offset.LS <- 1 - 1 / (1 + exp(-lambda.LS * 1.0))
  #   offset.MS <- 1 - 1 / (1 + exp(-lambda.MS * 1.0))
  #   offset.HS <- 1 - 1 / (1 + exp(-lambda.HS * 1.0))
  # 
  #   output <- data.frame(ode(c(S.LS = S.LS.offshore, I.LS = I.LS.offshore, R.LS = R.LS.offshore,
  #                              S.MS = S.MS.offshore, I.MS = I.MS.offshore, R.MS = R.MS.offshore,
  #                              S.HS = S.HS.offshore, I.HS = I.HS.offshore, R.HS = R.HS.offshore),
  #                            days.model.offshore, SIR.multi,
  #                            c(b.LS = beta.nearshore.LS, g.LS = gamma.nearshore.LS,
  #                              b.MS = beta.nearshore.MS, g.MS = gamma.nearshore.MS,
  #                              b.HS = beta.nearshore.HS, g.HS = gamma.nearshore.HS,
  #                              N.LS = N.LS.offshore, N.MS = N.MS.offshore, N.HS = N.HS.offshore,
  #                              C = cover.offshore,
  #                              C.LS = cover.offshore.LS, C.MS = cover.offshore.MS, C.HS = cover.offshore.HS,
  #                              l = lambda)))
  # 
  #   sim.rem.total <- rowSums(output[which(output$time %in% days.obs),
  #                                   which(colnames(output) %in% c('R.LS', 'R.MS', 'R.HS'))], na.rm = TRUE)
  # 
  #   obs.rem.total <- obs.model %>%
  #     filter(Site == curr.site, Category == "Total", Compartment == "Dead") %>%
  #     slice(head(row_number(), n()-DHW.modifier)) %>%
  #     pull(tissue)
  # 
  #   if (length(obs.rem.total) > length(sim.rem.total)) {
  #     obs.rem.total <- obs.rem.total[(length(obs.rem.total) - length(sim.rem.total) + 1):length(obs.rem.total)]
  #   }
  # 
  #   diff.rem.total <- (sim.rem.total - obs.rem.total)
  #   sum_diff.total <- sum(diff.rem.total^2)
  #   mean_obs_rem.total <- mean(obs.rem.total, na.rm = TRUE)
  #   tss_rem.total <- sum((obs.rem.total - mean_obs_rem.total)^2)
  # 
  #   r_squared <- 1 - (sum_diff.total / tss_rem.total)
  # 
  #   return(list(lambda.LS = lambda.LS, lambda.MS = lambda.MS, lambda.HS = lambda.HS, r_squared = r_squared))
  # }
  # 
  # # Run optimization with progress tracking
  # results <- pblapply(1:nrow(param_grid), function(i) evaluate_lambda(param_grid[i, ]))
  # 
  # # Extract best parameters
  # best_result <- results[[which.max(sapply(results, function(x) x$r_squared))]]
  # best_lambda.LS <- best_result$lambda.LS
  # best_lambda.MS <- best_result$lambda.MS
  # best_lambda.HS <- best_result$lambda.HS
  # 
  # cat("Best parameters found:\n")
  # cat(sprintf("lambda.LS = %.2f, lambda.MS = %.2f, lambda.HS = %.2f\n",
  #             best_lambda.LS, best_lambda.MS, best_lambda.HS))
  # 
