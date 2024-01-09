library(tseries)
library(dplyr)
library(urca)
library(forecast)
library(ggplot2)

rm(list = ls())
NUM_OBSERVATIONS <- 400
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
      ar_coefs <- runif(p, min = -0.8, max = 0.8)
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
  model_and_data <- list(rand_arima, rand_model)
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
  iq <<- 0
  # Iterate over the model list and compute the sum of AIC and BIC
  for (model_name in names(model_list)) {
    iq <<- iq + 1
    model <- model_list[[model_name]]
    model_name_here <<- model_name
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
  model_list_test <<- model_list
  best_combined_score <- Inf
  best_model <- NULL
  best_model_name <- NULL
  i <<- 0
  for (model_name in names(model_list)) {
    model_name <<- model_name
    model <<- model_list[[model_name]]
    i <<- i + 1

    # Calculate the sum of AIC and BIC
    aic_bic_sum <- AIC(model) + BIC(model)
    # Get p-values and calculate the proportion of significant p-values
    coefficients <- coef(model)
    std_errors <- sqrt(diag(vcov(model)))
  # Calculate z-scores
    z_scores <- coefficients / std_errors
  # Calculate two-tailed p-values
    p_values <<- 2 * (1 - pnorm(abs(z_scores)))
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










analyze_p_values <- function(model, significance_level = 0.05) {
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
  forecast::forecast(model, h = h)
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

















random_arima_data <- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = NUM_OBSERVATIONS, seed = NULL)[[1]]
random_arima_model <- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = NUM_OBSERVATIONS, seed = NULL)[[2]]
model_list <- fit_arima_models(random_arima_data, max_p = 3, max_q = 3)

best_model_least_p <- find_best_model_least_p(model_list)$best_model_min_p
models_analysis_least_p <- analyze_p_values(best_model_least_p, significance_level = 0.05)
print("Best model with least p_values")
best_model_least_p$coef
models_analysis_least_p$less_significant_parameters

best_model_auto.arima <- auto.arima(random_arima_data, allowdrift = TRUE, approximation = FALSE, stepwise = FALSE )
models_analysis_auto_arima <- analyze_p_values(best_model_auto.arima, significance_level = 0.05)
print("Best model with auto.arima")
best_model_auto.arima$coef
models_analysis_auto_arima$less_significant_parameters

best_model_least_aic_bic <- find_best_model_least_aic_bic(model_list)
model_analysis_least_aic_bic <- analyze_p_values(best_model_least_aic_bic$best_model, significance_level = 0.05)
print(paste("Best model based on AIC and BIC"))
print(best_model_least_aic_bic$best_model$coef)
model_analysis_least_aic_bic$less_significant_parameters

best_model_least_ic_3_p_vals <- find_best_model_least_ic_3_p_vals(model_list)
model_analysis_least_ic_3_p_vals <- analyze_p_values(best_model_least_ic_3_p_vals$best_model, significance_level = 0.05)
print(paste("Best combined model based on least ic and 3 p_vals"))
print(best_model_least_ic_3_p_vals$best_model$coef)
model_analysis_least_ic_3_p_vals$less_significant_parameters



find_best_model_least_ic_high_prop_sig_p_vals <- find_best_model_least_ic_high_prop_sig_p_vals(model_list)
print(paste("Best model based on low aic, high prop sig coeff is:"))
find_best_model_least_ic_high_prop_sig_p_vals$best_model$coef

# Calculate significant p_values:

# models_analysis_least_p <- analyze_p_values(best_model_least_p, significance_level = 0.05)
# models_analysis_auto_arima <- analyze_p_values(best_model_auto.arima, significance_level = 0.05)
# model_analysis_least_aic_bic <- analyze_p_values(best_model_least_aic_bic$best_model, significance_level = 0.05)
signifiance_p_values_best_model_least_p <- get_significant_p_values(best_model_least_p)
significance_p_values_auto_arima <- get_significant_p_values(best_model_auto.arima)
signifiance_p_values_best_model_least_aic_bic <- get_significant_p_values(best_model_least_aic_bic$best_model)
signifiance_p_values_least_ic_high_sig_prop <- get_significant_p_values(find_best_model_least_ic_high_prop_sig_p_vals$best_model)
signifiance_p_values_best_model_least_ic_3_p_vals <- get_significant_p_values(best_model_least_ic_3_p_vals$best_model)

# Cross validation testing:
average_errors_least_p <- perform_ts_cross_validation(random_arima_data, best_model_least_p)
print("Average errors for least p-value model")
print(average_errors_least_p)
average_errors_auto_arima <- perform_ts_cross_validation(random_arima_data, best_model_auto.arima)
print("Average errors for auto.arima model")
print(average_errors_auto_arima)
average_errors_least_aic_bic <- perform_ts_cross_validation(random_arima_data, best_model_least_aic_bic$best_model)
print("Average errors for least aic and bic model")
print(average_errors_least_aic_bic)
average_errors_least_ic_significant_p_vals <- perform_ts_cross_validation(random_arima_data, find_best_model_least_ic_high_prop_sig_p_vals$best_model)
print("Average errors for least ic, high prop sig parameters")
print(average_errors_least_ic_significant_p_vals)
average_errors_least_ic_3_p_vals <- perform_ts_cross_validation(random_arima_data, best_model_least_ic_3_p_vals$best_model)
print("Average errors for least ic and 3 p_vals model")
print(average_errors_least_ic_3_p_vals)








# Split the data into training and testing sets and forecasting and evaluating performance
data <- random_arima_data
split <- split_data(data)
model_least_p <- best_model_least_p
model_auto_arima <- best_model_auto.arima

forecast_least_p <- forecast_model(model_least_p, length(split$test))
forecast_auto_arima <- forecast_model(model_auto_arima, length(split$test))
forecast_least_aic_bic <- forecast_model(best_model_least_aic_bic$best_model, length(split$test))
forecast_least_ic_high_prop_sig_p_vals <- forecast_model(find_best_model_least_ic_high_prop_sig_p_vals$best_model, length(split$test))
forecast_least_ic_3_p_vals <- forecast_model(best_model_least_ic_3_p_vals$best_model, length(split$test))

performance_least_p <- evaluate_performance(split$test, forecast_least_p$mean)
performance_auto_arima <- evaluate_performance(split$test, forecast_auto_arima$mean)
performance_least_aic_bic <- evaluate_performance(split$test, forecast_least_aic_bic$mean)
performance_least_ic_high_prop_p_vals <- evaluate_performance(split$test, forecast_least_ic_high_prop_sig_p_vals$mean)
performance_least_ic_3_p_vals <- evaluate_performance(split$test, forecast_least_ic_3_p_vals$mean)


# Defining and usage of error terms:
# MAE: Use when we care about average of errors and not much about the outliers.
# RMSE: Use when we care about outliers, since it is squared, large errors are given more weight. Good if we care about occasional large errors
#  MAPE: Perentage of error, easy of interpret and scale independent; helpful to understand the relative size of errors.




plot_actual_vs_fitted <- function(actual, fitted, title) {
  # Convert 'ts' objects to regular vectors
  actual_vector <- as.vector(actual)
  fitted_vector <- as.vector(fitted)

  # Create a data frame for the actual and fitted values
  data_combined <- data.frame(
    Time = seq_along(actual_vector),
    Value = c(actual_vector, fitted_vector),
    Type = rep(c('Actual', 'Fitted'), each = length(actual_vector))
  )

  # Create the plot with improved aesthetics
  gg <- ggplot(data_combined, aes(x = Time, y = Value, color = Type, linetype = Type)) +
    geom_line(size = 1) +
    geom_point(data = subset(data_combined, Type == 'Actual'), size = 2) +
    scale_color_manual(values = c("Actual" = "blue", "Fitted" = "red")) +
    scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "dashed")) +
    labs(title = title, x = "Time", y = "Value") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "bottom")

  print(gg)
}

# Example usage:
# Extract fitted values from the models and convert them to regular vectors if necessary
plot_graphs <- function (){

fitted_values_least_p <- as.vector(fitted(best_model_least_p))
fitted_values_auto_arima <- as.vector(fitted(best_model_auto.arima))

# Plot for least p-value model
plot_actual_vs_fitted(actual = random_arima_data,
fitted = fitted_values_least_p,
title = "Actual vs Fitted Values: Least p-value Model")

# Plot for auto.arima model
plot_actual_vs_fitted(actual = random_arima_data,
fitted = fitted_values_auto_arima,
title = "Actual vs Fitted Values: Auto.arima Model")
}
# plot_graphs()

# models_analysis_auto_arima$less_significant_parameters
# model_analysis_least_aic_bic$less_significant_parameters
