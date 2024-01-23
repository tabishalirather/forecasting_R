library(tseries)
library(dplyr)
library(urca)
library(forecast)
library(ggplot2)
library(parallel)
rm(list = ls())

NUM_OBSERVATIONS <- 300
MAX_P <- 5
MAX_Q <- 5
MAX_D <- 0
# Fxn to geneate random ARIMA data and model
generate_random_arima <- function(max_p, max_d, max_q, num_observations, seed = NULL)
{
  # Setting a seed for reproducibility (optional, can be removed or modified)
  if(!is.null(seed)){
    print(seed)
    set.seed(seed)
  }

  # Randomly select the values for p, d, and q
  # p_for_now <- 4
  # q_for_now <- 3
  repeat
  {
    p <- sample(0:max_p, 1)
    # d <- sample(0:max_d, 1)
    d <- 0
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
  if(MAX_D == 0){
    do_include_constant <- TRUE
  }else{
    do_include_constant <- FALSE
  }
  model_list <- list()
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      model_list[[paste0("ARIMA_", p, "_0_", q)]] <- Arima(random_arima_data, order = c(p, 0, q), include.constant = do_include_constant)
    }
  }
  return(model_list)
}

# Fxn to calculate signficiancen p_values for parameters.
get_significant_p_values <- function(model) {
  # Extract coefficients and their standard errors
  coefficients <- coef(model)
  vcov_matrix <- vcov(model)
  vcov_matrix[vcov_matrix < 0] <- 0
  std_errors <- sqrt(diag(vcov_matrix))
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
  p_sw_arr <<- array(NaN, dim = length(model_list))
  p_lbq_arr <<- array(NaN, dim = length(model_list))
  p_ttest_arr <<- array(NaN, dim = length(model_list))
  for (index in seq_along((model_list))) {
    model <- model_list[[index]]
    residuals <- model$residuals
    p_sw_arr[index] <<- shapiro.test(residuals)$p.value
    p_lbq_arr[index] <<- Box.test(residuals, lag = 20, type = "Ljung-Box")$p.value
    p_ttest_arr[index] <<- t.test(residuals)$p.value
  }
  arr_sum_p_vals <- p_sw_arr + p_lbq_arr + p_ttest_arr
  # min_p_val <- min(arr_sum_p_vals)
  min_p_val <- max(arr_sum_p_vals)
  index_of_min_p_val <- which.max(arr_sum_p_vals)
  best_model_min_p <- model_list[[index_of_min_p_val]]
  order <- c(best_model_min_p$arma[[1]], best_model_min_p$arma[[6]], best_model_min_p$arma[[2]])
  return(list("min_p_val" = min_p_val, "index_of_min_p_val" = index_of_min_p_val, "best_model_min_p" = best_model_min_p, "order" = order))
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
  order <- c(best_model$arma[[1]], best_model$arma[[6]], best_model$arma[[2]])
  # Return the best model and its AIC and BIC sum
  list("best_model" = best_model, "aic_bic_sum" = best_aic_bic, "model_name" = best_model_name, "order" = order)
}

# Chose model with least sum of signifiane p_values and least aic and bic:
# Function to find the best model based on minimum sum of p-values and least sum of AIC and BIC
find_best_model_least_ic_high_prop_sig_p_vals <- function(model_list) {
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
    # std_errors <- sqrt(diag(vcov(model)))
    vcov_matrix <- vcov(model)
    vcov_matrix[vcov_matrix < 0] <- 0
    std_errors <- sqrt(diag(vcov_matrix))
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
  order <- c(best_model$arma[[1]], best_model$arma[[6]], best_model$arma[[2]])
  list("best_model" = best_model, "model_name" = best_model_name, "order" = order)
}


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
  order <- c(best_model$arma[[1]], best_model$arma[[6]], best_model$arma[[2]])
  list("best_model" = best_model, "model_name" = best_model_name, "order" = order)
}


analyse_significant_parameters <- function(model, significance_level = 0.05) {
  coefficients <- coef(model)
  vcov_matrix <- vcov(model)
  vcov_matrix[vcov_matrix < 0] <- 0
  std_errors <- sqrt(diag(vcov_matrix))
  z_scores <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_scores)))
  # Flag this model as having less significant coefficients
  if(length(p_values) != 0)
  {
    prop_significant <- sum(p_values < significance_level) / (length(p_values))
  } else {
    prop_significant <- 0
  }
  less_significant_parameters <- list(
    model = model,
    p_values = p_values,
    less_significant_parameters = names(p_values)[p_values > significance_level],
    count_less_significant_parameters = sum(p_values > significance_level),
    prop_sig_params = prop_significant
  )
  return(less_significant_parameters)
}


calculate_proportion_significant <- function(model_analysis, significance_level = 0.05) {
  # print("THis funciton is being called")
  model_analysis_in_calc <<- model_analysis
  if (is.null(model_analysis$p_values)) {
    return(NA)  # Return NA if p_values are not present
  }
  p_values <- model_analysis$p_values
  # Calculate the proportion of significant p_values
  num_sig_params <- sum(p_values < significance_level)
  proportion_significant <- num_sig_params / length(p_values)
  return(list(
    proportion_significant =  proportion_significant,
    num_sig_params = num_sig_params
  ))
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

  signifance_prop_least_p <- calculate_proportion_significant(models_analysis_least_p)
  signifance_prop_auto_arima <- calculate_proportion_significant(models_analysis_auto_arima)
  signifance_prop_least_aic_bic <- calculate_proportion_significant(models_analysis_least_aic_bic)
  signifance_prop_least_ic_high_prop_sig_p_vals <- calculate_proportion_significant(models_analysis_least_ic_high_prop_sig_p_vals)
  signifance_prop_least_ic_3_p_vals <- calculate_proportion_significant(models_analysis_least_ic_3_p_vals)


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


  model_list_for_plotting <<- list(
    "Least p-value Model" = list(
      model = best_model_least_p,
      overfit_analysis = models_analysis_least_p,
      cross_validation_errors = average_errors_least_p,
      forecast = forecast_least_p,
      performance = performance_least_p,
      significance_prop = signifance_prop_least_p
    ),
    "Auto.arima Model" = list(
      model = best_model_auto_arima,
      overfit_analysis = models_analysis_auto_arima,
      cross_validation_errors = average_errors_auto_arima,
      forecast = forecast_auto_arima,
      performance = performance_auto_arima,
      significance_prop = signifance_prop_auto_arima
    ),
    "Least AIC/BIC Model" = list(
      model = best_model_least_aic_bic,
      overfit_analysis = models_analysis_least_aic_bic,
      cross_validation_errors = average_errors_least_aic_bic,
      forecast = forecast_least_aic_bic,
      performance = performance_least_aic_bic,
      significance_prop = signifance_prop_least_aic_bic
    ),
    "High Prop Sig p-vals Model" = list(
      model = best_model_least_ic_high_prop_sig_p_vals,
      overfit_analysis = models_analysis_least_ic_high_prop_sig_p_vals,
      cross_validation_errors = average_errors_least_ic_high_prop_sig_p_vals,
      forecast = forecast_least_ic_high_prop_sig_p_vals,
      performance = performance_least_ic_high_prop_p_vals,
      significance_prop = signifance_prop_least_ic_high_prop_sig_p_vals
    ),
    "Least IC 3 p-vals Model" = list(
      model = best_model_least_ic_3_p_vals,
      overfit_analysis = models_analysis_least_ic_3_p_vals,
      cross_validation_errors = average_errors_least_ic_3_p_vals,
      forecast = forecast_least_ic_3_p_vals,
      performance = performance_least_ic_3_p_vals,
      significance_prop = signifance_prop_least_ic_3_p_vals

    )
  )
}


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


refit_least_p_model <- function (order_comparison) {
  new_rmse_vals <- list()
  refitted_model_analysis_least_p <- list()
  for(i in seq_along(order_comparison))
  {
    overfit_coeffs_least_p <- order_comparison[[i]]$models_analysis_least_p$less_significant_parameters
    print(paste("overfit_coeffs_least_p is: ", overfit_coeffs_least_p))
    original_p <- order_comparison[[i]]$least_p_order[[1]]
    print(paste("original_p in least p is: ", original_p))
    original_q <- order_comparison[[i]]$least_p_order[[3]]
    print(paste("original_q in least p is: ", original_q))

    adjusted_p <- original_p - sum(overfit_coeffs_least_p %in% paste0("ar", 1:original_p))
    print(paste("adjusted_p in least p is: ", adjusted_p))
    adjusted_q <- original_q - sum(overfit_coeffs_least_p %in% paste0("ma", 1:original_q))
    print(paste("adjusted_q in least p is: ", adjusted_q))

    exclude_intercept <- "intercept" %in% overfit_coeffs_least_p

    current_data <- order_comparison[[i]]$data

    refitted_model <- Arima(current_data, order=c(adjusted_p, 0, adjusted_q), include.mean=!exclude_intercept)
    # print(paste("refitted_model is: ", refitted_model))
    new_rmse_vals[[i]] <- perform_ts_cross_validation(current_data, refitted_model)[[2]]
    refitted_model_analysis_least_p[[i]] <- analyse_significant_parameters(refitted_model)
  }
  avg_new_rms_vals_p_model <- mean(unlist(new_rmse_vals))
  return(list(avg_new_rms_vals_p_model, refitted_model_analysis_least_p))
}


refit_auto_arima_model <- function (order_comparison) {
  new_rmse_vals <- list()
  refitted_model_analysis_auto_arima <- list()
  for(i in seq_along(order_comparison))
  {
    overfit_coeffs_auto_arima <- order_comparison[[i]]$models_analysis_auto_arima$less_significant_parameters
    p <- order_comparison[[i]]$auto_arima_order[[1]]
    q <- order_comparison[[i]]$auto_arima_order[[3]]
  print(paste("overfit_coeffs_auto_arima ", overfit_coeffs_auto_arima))
    print(paste("original_p ", p))
    print(paste("original_q ", q))

    p <- p - sum(overfit_coeffs_auto_arima %in% paste0("ar", 1:p))
    q <- q - sum(overfit_coeffs_auto_arima %in% paste0("ma", 1:q))

    exclude_intercept <- "intercept" %in% overfit_coeffs_auto_arima

    current_data <- order_comparison[[i]]$data
  # print(paste("c(adjusted_p, 0, adjusted_q) ", c(adjusted_p, 0, adjusted_q)))
    print(paste("adjusted_p ", p))
    print(paste("adjusted_q ", q))
    refitted_model_auto_arima <- Arima(current_data, order=c(p, 0, q), include.mean=!exclude_intercept)
    new_rmse_vals[[i]] <- perform_ts_cross_validation(current_data, refitted_model_auto_arima)[[2]]
    refitted_model_analysis_auto_arima[[i]] <- analyse_significant_parameters(refitted_model_auto_arima)

  }
  avg_new_rms_vals_auto_arima <- mean(unlist(new_rmse_vals))
  return(list(avg_new_rms_vals_auto_arima, refitted_model_analysis_auto_arima))
}


compare_arima_order <- function (random_arima_data_for_testing)
{
  model_list <- fit_arima_models(random_arima_data_for_testing, max_p = 3, max_q = 3)

  least_p_model <- find_best_model_least_p(model_list)
  auto_arima_model <- auto.arima(random_arima_data_for_testing, allowdrift = TRUE, approximation = FALSE, stepwise = FALSE )
  auto_arima_model$order <- c(auto_arima_model$arma[[1]], auto_arima_model$arma[[6]], auto_arima_model$arma[[2]])
  least_aic_bic_model <- find_best_model_least_aic_bic(model_list)
  least_ic_high_prop_sig_p_vals_model <- find_best_model_least_ic_high_prop_sig_p_vals(model_list)
  least_ic_3_p_vals_model <- find_best_model_least_ic_3_p_vals(model_list)
  # least_p_model$order <- c(least_p_model$best_model_min_p$arma[[1]], least_p_model$best_model_min_p$arma[[2]], least_p_model$best_model_min_p$arma[[6]])
  actual_order <- model_and_data[[2]]$order
  auto_arima_order <- auto_arima_model$order
  least_p_order <- least_p_model$order
  # return(list("least_p_model" = least_p_model, "auto_arima_model" = auto_arima_model, "least_aic_bic_model" = least_aic_bic_model, "least_ic_high_prop_sig_p_vals_model" = least_ic_high_prop_sig_p_vals_model, "least_ic_3_p_vals_model" = least_ic_3_p_vals_model, "actual_order" = actual_order))
  # model_list <- fit_arima_models(random_arima_data_for_testing, max_p = 3, max_q = 3)
  # least_p_model <- find_best_model_least_p(model_list)
  avg_rmse_least_p_model <- perform_ts_cross_validation(random_arima_data_for_testing,least_p_model$best_model_min_p)
  avg_rmse_auto_arima_model <- perform_ts_cross_validation(random_arima_data_for_testing,auto_arima_model)

  models_analysis_least_p <- analyse_significant_parameters(least_p_model$best_model_min_p)
  models_analysis_auto_arima <- analyse_significant_parameters(auto_arima_model)

  # Run the analysis once and store the result
  arima_analysis_result <- run_arima_analysis(NUM_OBSERVATIONS, random_arima_data_for_testing)

  # Extract the count of insignificant parameters for both models
  count_insig_params_p_model <- arima_analysis_result$`Least p-value Model`$overfit_analysis$count_less_significant_parameters
  count_insig_params_auto_arima_model <- arima_analysis_result$`Auto.arima Model`$overfit_analysis$count_less_significant_parameters


  comparison_results <- list(
    auto_arima_order_match = all(auto_arima_order == actual_order),
    auto_arima_order = auto_arima_order,
    least_p_order_match = all(least_p_order == actual_order),
    least_p_order = least_p_order,
    avg_rmse_least_p_model = avg_rmse_least_p_model,
    avg_rmse_auto_arima_model = avg_rmse_auto_arima_model,
    count_insig_params_p_model = count_insig_params_p_model,
    count_insig_params_auto_arima_model = count_insig_params_auto_arima_model,
    least_p_model = least_p_model$best_model_min_p,
    auto_arima_model = auto_arima_model,
    models_analysis_least_p = models_analysis_least_p,
    models_analysis_auto_arima = models_analysis_auto_arima
  )
  # }
  return(comparison_results)

}


count_order_matches <- function (order_comparison){ # Counting TRUE values for auto_arima_order_match
  # calc avg statistics for the order comparison.
  auto_arima_true_count <- sum(sapply(order_comparison, function(x) x$auto_arima_order_match))
  # TODO: get the indicies of true matches
  # TODO: get the indicies of true matchesth
  # Counting FALSE values for auto_arima_order_match
  auto_arima_false_count <- sum(sapply(order_comparison, function(x) !x$auto_arima_order_match))
  # Counting TRUE values for least_p_order_match
  least_p_true_count <- sum(sapply(order_comparison, function(x) x$least_p_order_match))
  # Counting FALSE values for least_p_order_match
  least_p_false_count <- sum(sapply(order_comparison, function(x) !x$least_p_order_match))

  # calc-avg_rmse
  avg_rmse_least_p_model_values <- mean(sapply(order_comparison, function(x) x$avg_rmse_least_p_model))
  avg_rmse_auto_arima_model <- mean(sapply(order_comparison, function(x) x$avg_rmse_auto_arima_model))

  # avg non_sig parameters
  avg_insig_params_least_p_model <- mean(sapply(order_comparison, function(x) x$count_insig_params_p_model))
  avg_insig_params_auto_arima_model <- mean(sapply(order_comparison, function(x) x$count_insig_params_auto_arima_model))


  from_refit_least_p_model <- refit_least_p_model(order_comparison)
  adjusted_least_p_avg_rmse <- from_refit_least_p_model[[1]]
  model_analysis_refit_p <<- from_refit_least_p_model[[2]]

  from_refit_auto_arima <- refit_auto_arima_model(order_comparison)
  adjusted_auto_arima_avg_rmse <-from_refit_auto_arima[[1]]
  model_analysis_refit_auto_arima <<-from_refit_auto_arima[[2]]

  avg_insig_params_least_p_model_refit <- mean(sapply(model_analysis_refit_p, function(x) x$count_less_significant_parameters))
  print(paste("avg_insig_params_least_p_model_refit is: ", avg_insig_params_least_p_model_refit))
  avg_insig_params_auto_arima_model_refit <- mean(sapply(model_analysis_refit_auto_arima, function(x) x$count_less_significant_parameters))
  print(paste("avg_insig_params_auto_arima_model_refit is: ", avg_insig_params_auto_arima_model_refit))

  # Create a DF:
  test_data <- data.frame(
    actual_order = paste(model_and_data[[2]]$order, collapse = " "),
    num_obser = NUM_OBSERVATIONS,
    auto_arima_true_count = auto_arima_true_count,
    auto_arima_false_count = auto_arima_false_count,
    least_p_true_count = least_p_true_count,
    least_p_false_count = least_p_false_count,
    avg_rmse_updated_least_p_model = adjusted_least_p_avg_rmse,
    avg_rmse_least_p_model_values = avg_rmse_least_p_model_values,
    avg_rmse_updated_auto_arima = adjusted_auto_arima_avg_rmse,
    avg_rmse_auto_arima_model = avg_rmse_auto_arima_model,
    avg_insig_params_least_p_model_refit = avg_insig_params_least_p_model_refit,
    avg_insig_params_least_p_model = avg_insig_params_least_p_model,
    avg_insig_params_auto_arima_model_refit = avg_insig_params_auto_arima_model_refit,
    avg_insig_params_auto_arima_model = avg_insig_params_auto_arima_model
  )
  # Output the counts
  output_file_path <- "C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/TestingR/Week 1+/output_2.csv"
  write.table(test_data, file = output_file_path, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(output_file_path), quote = FALSE)
}
# Run the analysis
print("Running ARIMA analysis...")
num_observations <- NUM_OBSERVATIONS
all_analysis <- list()
run_single_iteration <- function(num_observations) {
  # all_analysis[[iteration]] <- iteration
  random_arima_data <- generate_random_arima(max_p = 5, max_d = 0, max_q = 5, num_observations, seed = NULL)[[1]]
  analysis <<- run_arima_analysis(NUM_OBSERVATIONS, random_arima_data)
  # Rank models based on cross-validation RMSE
  ranked_models_cv <- rank_models_rmse(analysis, use_rmse_type = "both")
  return(list(ranked_models_cv, analysis))

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
                      "find_best_model_least_aic_bic", "find_best_model_least_p", "calculate_proportion_significant", "all_analysis"))

  # Load required libraries in each node of the cluster
  clusterEvalQ(cl, {
    library(forecast)
    library(dplyr)
    # Load other necessary libraries here
  })
  stopCluster(cl)
  return(all_rankings)
}
# Function to calculate average rankings
calculate_average_rankings <- function(all_rankings) {
  # Combine all data frames into one
  combined_rankings <- do.call(rbind, all_rankings)
  # Calculate average RMSE for each model
  average_rmse <- aggregate(cbind(RMSE_CV) ~ Model, data = combined_rankings, mean)
  # Order models based on average RMSE in ascending order
  average_rmse <- average_rmse[order(average_rmse$RMSE_CV), ]
  # Add a new rank column based on the sorted average RMSE
  average_rmse$Rank <- seq_along(average_rmse$RMSE_CV)
  return(average_rmse)
}
num_iterations <- 1 # Adjust this as needed
all_rankings <- list()
for (i in 1:num_iterations) {
  # all_rankings[[i]] <- run_single_iteration(num_observations)[[1]]
}

# cmmt these out for now:
add_parameter_info <- function(analysis, ranked_models_df) {
  # Initialize new columns with NAs
  ranked_models_df$Count_Less_Significant_Parameters <- NA
  ranked_models_df$Prop_Significant_Params <- NA
  # Loop through each row in ranked_models_df
  for (i in seq_len(nrow(ranked_models_df))) {
    model_name <- ranked_models_df$Model[i]
    # Check if the model_name exists in analysis
    if (model_name %in% names(analysis)) {
      # Extract the count and proportion for the matching model
      ranked_models_df$Count_Less_Significant_Parameters[i] <- analysis[[model_name]]$overfit_analysis$count_less_significant_parameters
      ranked_models_df$Prop_Significant_Params[i] <- analysis[[model_name]]$overfit_analysis$prop_sig_params
    }
  }

  return(ranked_models_df)
}

order_comparison <- list()
# comparison_results <- list()
num_iters_csv <- 20
for(i in seq_along(1:num_iters_csv))
{
  print(i)
  random_arima_data_for_testing <- generate_random_arima(max_p = MAX_P, max_d = MAX_D, max_q = MAX_Q, num_observations = NUM_OBSERVATIONS, seed = NULL)[[1]]
  order_comparison[[i]] <- compare_arima_order(random_arima_data_for_testing)
  order_comparison[[i]]$data <- random_arima_data_for_testing
}

count_order_matches(order_comparison)
