library(tseries)
library(dplyr)
library(urca)
library(forecast)
library(ggplot2)
library(dplyr)
library(parallel)
rm(list = ls())

NUM_OBSERVATIONS <- 300
# Fxn to geneate random ARIMA data and model
generate_random_arima <- function(max_p, max_d, max_q, num_observations, seed = NULL) {

  # Setting a seed for reproducibility (optional, can be removed or modified)
  if(!is.null(seed)){
    print(seed)
    set.seed(seed)
  }
  # Randomly select the values for p, d, and q
  repeat
  {
    p <- sample(0:max_p, 1)
    d <- sample(0:max_d, 1)
    q <- sample(0:max_q, 1)
    # Generate random AR and MA coefficients
    ar_coefs <- if(p > 0)
    {
      ar_coefs <- runif(p, min = -0.4, max = 0.5)
    }
    else
    {
      numeric(0)  # Handles case when p=0
    }
    ma_coefs <- if(q > 0)
    {
      runif(q, min = -1, max = 1)
    }
    else
    {
      numeric(0)  # Handles case when q=0
    }
    # Define the ARIMA model
    rand_model <- list(ar = ar_coefs, ma = ma_coefs, order = c(p, d, q))

    # Generate and return the ARIMA data, the error we are concerned about is "Error: AR is not stationary, so we repeat the process until a given interation has a stationary AR model
    rand_arima <- tryCatch({
      arima.sim(n = num_observations, model = rand_model)
    }, error = function(e) {
      # If an error occurs, print a message and return NULL
      cat("Error in generating ARIMA data: ", e$message, "\n")
      return(NULL)
    })

    # If arima.sim was successful and did not return NULL, break the repeat loop
    if (!is.null(rand_arima)) {
      break
    }
  }
  model_and_data <<- list(rand_arima, rand_model)
  return(model_and_data)
}

# Fxn to fit ARIMA models
fit_arima_models <- function(random_arima_data, max_p, max_q) {
  model_list <- list()
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      model_list[[paste0("ARIMA_", p, "_1_", q)]] <- arima(random_arima_data, order = c(p, 1, q))
    }
  }
  return(model_list)
}

# Fxn to calculate signficiancen p_values for parameters.
get_significant_p_values <- function(model) {
  # Extract coefficients and their standard errors
  coefficients <- coef(model)
  std_errors <- sqrt(diag(vcov(model)))
  # Calculate z-scores
  z_scores <- coefficients / std_errors
  # Calculate two-tailed p-values
  p_values <- 2 * (1 - pnorm(abs(z_scores)))

  model_p_values <- list(
    model = model,
    p_value = p_values
  )
  return(model_p_values)
}

# Find optimal model with least sum of p_values.
find_best_model_least_p <- function(model_list) {
  p_sw_arr <- array(NaN, dim = length(model_list))
  p_lbq_arr <- array(NaN, dim = length(model_list))
  p_ttest_arr <- array(NaN, dim = length(model_list))

  for (index in seq_along((model_list))) {
    model <- model_list[[index]]
    residuals <- model$residuals
    p_sw_arr[index] <- shapiro.test(residuals)$p.value
    p_lbq_arr[index] <- Box.test(residuals, lag = 20, type = "Ljung-Box")$p.value
    p_ttest_arr[index] <- t.test(residuals)$p.value
  }

  arr_sum_p_vals <- p_sw_arr + p_lbq_arr + p_ttest_arr
  min_p_val <- min(arr_sum_p_vals)
  index_of_min_p_val <- which.min(arr_sum_p_vals)

  best_model_min_p <- model_list[[index_of_min_p_val]]

  return(list("min_p_val" = min_p_val, "index_of_min_p_val" = index_of_min_p_val, "best_model_min_p" = best_model_min_p))
}

# Function to find the model with the minimum sum of AIC and BIC
find_best_model_least_aic_bic <- function(model_list) {
  best_aic_bic <- Inf  # Initialize with infinity, to ensure any real sum will be less
  best_model <- NULL
  # Iterate over the model list and compute the sum of AIC and BIC
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    model_name_here <- model_name
    aic_bic_sum <- AIC(model) + BIC(model)

    # If this model has a lower sum than the current best, update the best values
    if (aic_bic_sum < best_aic_bic) {
      best_aic_bic <- aic_bic_sum
      best_model <- model
      best_model_name <- model_name  # Optionally keep track of the model name
    }
  }

  # Return the best model and its AIC and BIC sum
  list("best_model" = best_model, "aic_bic_sum" = best_aic_bic, "model_name" = best_model_name)
}

# Chose model with least sum of signifiane p_values and least aic and bic:
# Function to find the best model based on minimum sum of p-values and least sum of AIC and BIC
find_best_model_least_ic_high_prop_sig_p_vals <- function(model_list) {
  model_list_test <- model_list
  best_combined_score <- Inf
  best_model <- NULL
  best_model_name <- NULL

  for (model_name in names(model_list)) {
    model_name <- model_name
    model <- model_list[[model_name]]

    # Calculate the sum of AIC and BIC
    aic_bic_sum <- AIC(model) + BIC(model)
    # Get p-values and calculate the proportion of significant p-values
    coefficients <- coef(model)
    std_errors <- sqrt(diag(vcov(model)))
    # Calculate z-scores
    z_scores <- coefficients / std_errors
    # Calculate two-tailed p-values
    p_values <- 2 * (1 - pnorm(abs(z_scores)))
    significant_p_values <- p_values[p_values < 0.05]  # threshold for significance
    if(length(p_values) != 0)
    {
      prop_significant <- length(significant_p_values) / (length(p_values))
    } else {
      prop_significant <- 0
    }
    # Combine AIC/BIC sum and proportion of significant p-values
    # Since higher proportion is better, subtract it from the combined score
    combined_score <- aic_bic_sum - prop_significant

    # Check if this model has a better (lower) score than the current best
    if (combined_score < best_combined_score) {
      best_combined_score <- combined_score
      best_model <- model
      best_model_name <- model_name
    }
  }

  list("best_model" = best_model, "model_name" = best_model_name)
}

# Function that selects model based on least aic and bic and least sum of 3 tests
# Function to find the model with the least sum of AIC, BIC, and p-values from sw_test, t_test, and lbq_test
find_best_model_least_ic_3_p_vals <- function(model_list) {
  best_score <- Inf
  best_model <- NULL
  best_model_name <- NULL

  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]

    # Calculate AIC and BIC sum
    aic_bic_sum <- AIC(model) + BIC(model)

    # Calculate p-values for the Shapiro-Wilk test, Ljung-Box test, and t-test
    residuals <- residuals(model)
    p_sw <- shapiro.test(residuals)$p.value
    p_lbq <- Box.test(residuals, lag = 20, type = "Ljung-Box")$p.value
    p_t <- t.test(residuals)$p.value

    # Sum of p-values
    sum_p_values <- p_sw + p_lbq + p_t

    # Combined score
    combined_score <- aic_bic_sum + sum_p_values

    # Selecting the model with the lowest combined score
    if (combined_score < best_score) {
      best_score <- combined_score
      best_model <- model
      best_model_name <- model_name
    }
  }

  list("best_model" = best_model, "model_name" = best_model_name)
}










analyse_significant_parameters <- function(model, significance_level = 0.05) {
  coefficients <- coef(model)
  std_errors <- sqrt(diag(vcov(model)))
  z_scores <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_scores)))

  # Check if any p-value exceeds the significance level
  if (any(p_values > significance_level)) {
    # Flag this model as having less significant coefficients
    less_significant_parameters <- list(
      model = model,
      p_values = p_values,
      less_significant_parameters = names(p_values)[p_values > significance_level],
      count_less_significant_parameters = sum(p_values > significance_level)
    )
  } else {
    less_significant_parameters <- list(
      message = "All parameters are significant"
    )
  }

  return(less_significant_parameters)
}


# Fxn to split the data into training and testing sets:
split_data <- function(data, train_frac = 0.8) {
  train_size <- round(length(data) * train_frac)
  list(
    train = data[1:train_size],
    test = data[(train_size + 1):length(data)]
  )
}

forecast_model <- function(model, h) {
  # print(paste("h is: ", h))
  # print(paste("model is: ", model$coef))
  forecast_result <- forecast::forecast(model, h = h)
  print(forecast_result)
}

# Function to calculate performance metrics
evaluate_performance <- function(actual, forecasted) {
  mae <- mean(abs(forecasted - actual))
  rmse <- sqrt(mean((forecasted - actual)^2))
  mape <- mean(abs((forecasted - actual) / actual)) * 100
  list(MAE = mae, RMSE = rmse, MAPE = mape)
}



# Cross validation:
# Function for time series cross-validation
perform_ts_cross_validation <- function(random_arima_data, fitted_model, initial_train_frac = 0.5) {
  # Initial training set size
  initial_train_size <- round(length(random_arima_data) * initial_train_frac)

  # Total length of the time series
  total_length <- length(random_arima_data)

  # Store forecast errors
  errors <- data.frame(MAE = numeric(), RMSE = numeric())

  # Cross-validation loop
  for (i in initial_train_size:(total_length - 1)) {
    # Define the test set
    test_set <- random_arima_data[(i + 1):(i + 1)]  # One step ahead forecast

    # Forecast using the provided model
    forecast <- forecast::forecast(fitted_model, h = 1)

    # Calculate errors
    mae <- mean(abs(test_set - forecast$mean))
    rmse <- sqrt(mean((test_set - forecast$mean)^2))

    # Store errors
    errors <- rbind(errors, data.frame(MAE = mae, RMSE = rmse))
  }

  # Average errors
  average_errors <- colMeans(errors)
  return(average_errors)
}
# Example usage
# Assuming `random_arima_data` is your time series data and `best_model_auto.arima` is a function that fits your selected model




run_arima_analysis <- function(num_observations, random_arima_data) {
  num_observations <- 400
  # Generate random ARIMA data and model
  if (is.null(random_arima_data[[1]])) {
    stop("ARIMA data generation failed.")
  }
  # Fit ARIMA models
  model_list <- fit_arima_models(random_arima_data, max_p = 3, max_q = 3)

  # Find best models based on different criteria
  best_model_least_p <- find_best_model_least_p(model_list)$best_model
  best_model_auto_arima <- auto.arima(random_arima_data, allowdrift = TRUE, approximation = FALSE, stepwise = FALSE)
  best_model_least_aic_bic <- find_best_model_least_aic_bic(model_list)$best_model
  best_model_least_ic_high_prop_sig_p_vals <- find_best_model_least_ic_high_prop_sig_p_vals(model_list)$best_model
  best_model_least_ic_3_p_vals <- find_best_model_least_ic_3_p_vals(model_list)$best_model

  # Analyze models
  models_analysis_least_p <- analyse_significant_parameters(best_model_least_p, significance_level = 0.05)
  models_analysis_auto_arima <- analyse_significant_parameters(best_model_auto_arima, significance_level = 0.05)
  models_analysis_least_aic_bic <- analyse_significant_parameters(best_model_least_aic_bic, significance_level = 0.05)
  models_analysis_least_ic_high_prop_sig_p_vals <- analyse_significant_parameters(best_model_least_ic_high_prop_sig_p_vals, significance_level = 0.05)
  models_analysis_least_ic_3_p_vals <- analyse_significant_parameters(best_model_least_ic_3_p_vals, significance_level = 0.05)


  # Split data and forecast
  split <- split_data(random_arima_data)
  forecast_least_p <- forecast::forecast(best_model_least_p, length(split$test))
  print("forecast_tes working")
  forecast_auto_arima <- forecast::forecast(best_model_auto_arima, length(split$test))
  forecast_least_aic_bic <- forecast::forecast(best_model_least_aic_bic, length(split$test))
  forecast_least_ic_high_prop_sig_p_vals <- forecast::forecast(best_model_least_ic_high_prop_sig_p_vals, length(split$test))
  forecast_least_ic_3_p_vals <- forecast::forecast(best_model_least_ic_3_p_vals, length(split$test))

  # Perform cross-validation
  average_errors_least_p <- perform_ts_cross_validation(random_arima_data, best_model_least_p)
  average_errors_auto_arima <- perform_ts_cross_validation(random_arima_data, best_model_auto_arima)
  average_errors_least_aic_bic <- perform_ts_cross_validation(random_arima_data, best_model_least_aic_bic)
  average_errors_least_ic_high_prop_sig_p_vals <- perform_ts_cross_validation(random_arima_data, best_model_least_ic_high_prop_sig_p_vals)
  average_errors_least_ic_3_p_vals <- perform_ts_cross_validation(random_arima_data, best_model_least_ic_3_p_vals)





  # Evaluate performance
  performance_least_p <- evaluate_performance(split$test, forecast_least_p$mean)
  performance_auto_arima <- evaluate_performance(split$test, forecast_auto_arima$mean)
  performance_least_aic_bic <- evaluate_performance(split$test, forecast_least_aic_bic$mean)
  performance_least_ic_high_prop_p_vals <- evaluate_performance(split$test, forecast_least_ic_high_prop_sig_p_vals$mean)
  performance_least_ic_3_p_vals <- evaluate_performance(split$test, forecast_least_ic_3_p_vals$mean)


  model_list_for_plotting <- list(
    "Least p-value Model" = list(
      model = best_model_least_p,
      overfit_analysis = models_analysis_least_p,
      cross_validation_errors = average_errors_least_p,
      forecast = forecast_least_p,
      performance = performance_least_p
    ),
    "Auto.arima Model" = list(
      model = best_model_auto_arima,
      overfit_analysis = models_analysis_auto_arima,
      cross_validation_errors = average_errors_auto_arima,
      forecast = forecast_auto_arima,
      performance = performance_auto_arima
    ),
    "Least AIC/BIC Model" = list(
      model = best_model_least_aic_bic,
      overfit_analysis = models_analysis_least_aic_bic,
      cross_validation_errors = average_errors_least_aic_bic,
      forecast = forecast_least_aic_bic,
      performance = performance_least_aic_bic
    ),
    "High Prop Sig p-vals Model" = list(
      model = best_model_least_ic_high_prop_sig_p_vals,
      overfit_analysis = models_analysis_least_ic_high_prop_sig_p_vals,
      cross_validation_errors = average_errors_least_ic_high_prop_sig_p_vals,
      forecast = forecast_least_ic_high_prop_sig_p_vals,
      performance = performance_least_ic_high_prop_p_vals
    ),
    "Least IC 3 p-vals Model" = list(
      model = best_model_least_ic_3_p_vals,
      overfit_analysis = models_analysis_least_ic_3_p_vals,
      cross_validation_errors = average_errors_least_ic_3_p_vals,
      forecast = forecast_least_ic_3_p_vals,
      performance = performance_least_ic_3_p_vals
    )
  )

}

# rank_models_rmse <- function (analysis){
# model_names <- names(analysis)
# rmse_values <- sapply(analysis, function(x) x$performance$RMSE)
# # Creating a data frame with model names and their RMSE values
# rmse_df <- data.frame(Model = model_names, RMSE = rmse_values)
# # Using dplyr to rank the models with dense_rank()
# rmse_df <- rmse_df %>%
# arrange(RMSE) %>%
# mutate(Rank = dense_rank(RMSE))
# return(rmse_df)
# }
# rank_models_rmse_combined <- function(analysis) {
#   model_names <- names(analysis)
#   # Extract RMSE values from cross-validation
#   rmse_values_cv <- sapply(analysis, function(x) x$cross_validation_errors[[2]])
#   # Extract RMSE values from evaluate_performance
#   rmse_values_ep <- sapply(analysis, function(x) x$performance$RMSE)
#   # Creating data frames for each RMSE type
#   rmse_df_cv <- data.frame(Model = model_names, RMSE_CV = rmse_values_cv)
#   rmse_df_ep <- data.frame(Model = model_names, RMSE_EP = rmse_values_ep)
#   # Ranking the models based on cross-validation RMSE
#   if(is.null(rmse_df_cv))
#   {
#     print("rmse_df_cv is null")
#   }
#   rmse_df_cv <- rmse_df_cv %>%
#     arrange(RMSE_CV) %>%
#     mutate(Rank_CV = dense_rank(RMSE_CV))
#   # Ranking the models based on evaluate_performance RMSE
#   rmse_df_ep <- rmse_df_ep %>%
#     arrange(RMSE_EP) %>%
#     mutate(Rank_EP = dense_rank(RMSE_EP))
#   # Merging the two data frames
#   rmse_df_combined <- merge(rmse_df_cv, rmse_df_ep, by = "Model")
#   return(rmse_df_combined)
# }

# Function to rank models based on chosen RMSE type or both
rank_models_rmse <- function(analysis, use_rmse_type = "cv") {
  model_names <- names(analysis)

  # Initialize an empty data frame
  rmse_df <- data.frame(Model = model_names)

  # Choose which RMSE to use based on the use_rmse_type
  if (use_rmse_type == "cv" || use_rmse_type == "both") {
    # Use RMSE from cross-validation
    rmse_values_cv <- sapply(analysis, function(x) x$cross_validation_errors[[2]])
    rmse_df$RMSE_CV <- rmse_values_cv
  }
  if (use_rmse_type == "ep" || use_rmse_type == "both") {
    # Use RMSE from evaluate_performance
    rmse_values_ep <- sapply(analysis, function(x) x$performance$RMSE)
    rmse_df$RMSE_EP <- rmse_values_ep
  }

  # Ranking the models based on chosen RMSE type(s)
  if (use_rmse_type == "cv") {
    rmse_df <- rmse_df %>%
      arrange(RMSE_CV) %>%
      mutate(Rank_CV = dense_rank(RMSE_CV))
  } else if (use_rmse_type == "ep") {
    rmse_df <- rmse_df %>%
      arrange(RMSE_EP) %>%
      mutate(Rank_EP = dense_rank(RMSE_EP))
  } else if (use_rmse_type == "both") {
    rmse_df <- rmse_df %>%
      arrange(RMSE_CV, RMSE_EP) %>%
      mutate(Rank_CV = dense_rank(RMSE_CV),
             Rank_EP = dense_rank(RMSE_EP))
  }
  return(rmse_df)
}

# Example usage:
# To rank using RMSE from cross-validation



# Run the analysis
print("Running ARIMA analysis...")
num_observations <- NUM_OBSERVATIONS
# random_arima_data <- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations, seed = NULL)[[1]]
# analysis <- run_arima_analysis(NUM_OBSERVATIONS, random_arima_data)
# # Assuming 'analysis' is your list
# ranked_models_cv <- rank_models_rmse(analysis, use_rmse_type = "cv")
# Function to run a single iteration of ARIMA analysis and ranking
run_single_iteration <- function(num_observations) {
  random_arima_data <- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations, seed = NULL)[[1]]
  analysis <- run_arima_analysis(NUM_OBSERVATIONS, random_arima_data)
  # Rank models based on cross-validation RMSE
  ranked_models_cv <- rank_models_rmse(analysis, use_rmse_type = "cv")
  return(ranked_models_cv)
}
# Number of iterations to perform

get_all_rankings_paralell <- function(num_iterations, num_observations) {
  # Function definitions (e.g., run_single_iteration) should be included here or sourced if they are defined in separate files
  # Detect the number of cores and set up a parallel cluster
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)

  # Export necessary functions and objects to the cluster
  clusterExport(cl, c("run_single_iteration", "NUM_OBSERVATIONS", "generate_random_arima", "run_arima_analysis",
                      "rank_models_rmse", "fit_arima_models", "analyse_significant_parameters", "split_data",
                      "forecast_model", "evaluate_performance", "perform_ts_cross_validation",
                      "find_best_model_least_ic_3_p_vals", "find_best_model_least_ic_high_prop_sig_p_vals",
                      "find_best_model_least_aic_bic", "find_best_model_least_p"))

  # Load required libraries in each node of the cluster
  clusterEvalQ(cl, {
    library(forecast)
    library(dplyr)
    # Load other necessary libraries here
  })
  # Run the analysis in parallel using parLapply
  all_rankings <- parLapply(cl, 1:num_iterations, function(x) run_single_iteration(num_observations))
  # Shut down the cluster
  stopCluster(cl)
  return(all_rankings)
}
# get_all_rankings <- function (num_iterations, num_observations, all_rankings) {
# # # Run the analysis and save rankings for each iteration
# for (i in 1:num_iterations) {
#   all_rankings[[i]] <- run_single_iteration(num_observations)
# }
# }
# Function to calculate average rankings
calculate_average_rankings <- function(all_rankings) {
  # Combine all data frames into one
  combined_rankings <- do.call(rbind, all_rankings)
  # Calculate average RMSE for each model
  average_rmse <- aggregate(RMSE_CV ~ Model, data = combined_rankings, mean)
  # Order models based on average RMSE in ascending order
  average_rmse <- average_rmse[order(average_rmse$RMSE_CV), ]
  # Add a new rank column based on the sorted average RMSE
  average_rmse$Rank <- seq_along(average_rmse$RMSE_CV)
  return(average_rmse)
}
#Calculate average rankings


num_iterations <- 10 # Adjust this as needed
all_rankings <- list()
all_rankings <- get_all_rankings_paralell(num_iterations, NUM_OBSERVATIONS)
# all_rankings <- get_all_rankings(num_iterations, NUM_OBSERVATIONS, all_rankings)
# for (i in 1:num_iterations) {
#   all_rankings[[i]] <- run_single_iteration(num_observations)
# }

average_rankings <- calculate_average_rankings(all_rankings)

# To rank using RMSE from evaluate_performance
# ranked_models_ep <- rank_models_rmse(analysis, use_rmse_type = "ep")
# To rank using both RMSE CV and RMSE EP


# ranked_models_both <- rank_models_rmse(analysis, use_rmse_type = "both")

# Defining and usage of error terms:
# MAE: Use when we care about average of errors and not much about the outliers.
# RMSE: Use when we care about outliers, since it is squared, large errors are given more weight. Good if we care about occasional large errors
#  MAPE: Perentage of error, easy of interpret and scale independent; helpful to understand the relative size of errors.




# plot_actual_vs_fitted <- function(actual, fitted, title) {
#   # Convert 'ts' objects to regular vectors
#   actual_vector <- as.vector(actual)
#   fitted_vector <- as.vector(fitted)
#
#   # Create a data frame for the actual and fitted values
#   data_combined <- data.frame(
#     Time = seq_along(actual_vector),
#     Value = c(actual_vector, fitted_vector),
#     Type = rep(c('Actual', 'Fitted'), each = length(actual_vector))
#   )
#
#   # Create the plot with improved aesthetics
#   gg <- ggplot(data_combined, aes(x = Time, y = Value, color = Type, linetype = Type)) +
#     geom_line(size = 1) +
#     geom_point(data = subset(data_combined, Type == 'Actual'), size = 2) +
#     scale_color_manual(values = c("Actual" = "blue", "Fitted" = "red")) +
#     scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "dashed")) +
#     labs(title = title, x = "Time", y = "Value") +
#     theme_minimal() +
#     theme(legend.title = element_blank(), legend.position = "bottom")
#
#   print(gg)
# }

# Example usage:
# Extract fitted values from the models and convert them to regular vectors if necessary
# plot_graphs <- function (){
#
# fitted_values_least_p <- as.vector(fitted(best_model_least_p))
# fitted_values_auto_arima <- as.vector(fitted(best_model_auto.arima))
#
# # Plot for least p-value model
# plot_actual_vs_fitted(actual = random_arima_data,
# fitted = fitted_values_least_p,
# title = "Actual vs Fitted Values: Least p-value Model")
#
# # Plot for auto.arima model
# plot_actual_vs_fitted(actual = random_arima_data,
# fitted = fitted_values_auto_arima,
# title = "Actual vs Fitted Values: Auto.arima Model")
# }
# plot_graphs()

# models_analysis_auto_arima$less_significant_parameters
# model_analysis_least_aic_bic$less_significant_parameters
# rm(list = ls())