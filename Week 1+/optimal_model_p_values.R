library(forecast)
library(foreach)
library(doParallel)
rm(list = ls())

NUM_OBSERVATIONS <- 300
MAX_P <- 5
MAX_Q <- 5
MAX_D <- 1
OUTPUT_PATH <- "C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/TestingR/Week 1+/output_5_non_zero_d.csv"
NUM_TEST_MODELS <- 10
# Fxn to geneate random ARIMA data and model
generate_random_arima <- function(max_p, max_d, max_q, num_observations, seed = NULL)
{
  # Setting a seed for reproducibility (optional, can be removed or modified)
  if(!is.null(seed)){
    print(seed)
    set.seed(seed)
  }
  # Randomly select the values for p, d, and q
  repeat
  {
    p <- sample(0:max_p, 1)
    # d <- sample(0:max_d, 1)
    d <- max_d
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
  d <- MAX_D
  if(MAX_D == 0){
    do_include_constant <- TRUE
  }else{
    do_include_constant <- FALSE
  }
  model_list <<- list()
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      model_list[[paste0("ARIMA_", p, "_0_", q)]] <<- Arima(random_arima_data, order = c(p, d, q), include.constant = do_include_constant)
    }
  }
  return(model_list)
}

# Find optimal model with least sum of p_values.
find_best_model_max_p <- function(model_list) {
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
  # min_p_val <- min(arr_sum_p_vals)
  min_p_val <- max(arr_sum_p_vals)
  index_of_min_p_val <- which.max(arr_sum_p_vals)
  best_model_min_p <- model_list[[index_of_min_p_val]]
  order <- c(best_model_min_p$arma[[1]], best_model_min_p$arma[[6]], best_model_min_p$arma[[2]])
  return(list("min_p_val" = min_p_val, "index_of_min_p_val" = index_of_min_p_val, "best_model_min_p" = best_model_min_p, "order" = order))
}

# Function to find the model with the least sum of AIC, BIC, and p-values from sw_test, t_test, and lbq_test

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

refit_least_p_model <- function (order_comparison) {
  new_rmse_vals <- list()
  refitted_model_analysis_least_p <- list()
  for(i in seq_along(order_comparison))
  {
    overfit_coeffs_least_p <- order_comparison[[i]]$models_analysis_least_p$less_significant_parameters
    # print(paste("overfit_coeffs_least_p is: ", overfit_coeffs_least_p))
    original_p <- order_comparison[[i]]$least_p_order[[1]]
    # print(paste("original_p in least p is: ", original_p))
    original_q <- order_comparison[[i]]$least_p_order[[3]]
    # print(paste("original_q in least p is: ", original_q))

    adjusted_p <- original_p - sum(overfit_coeffs_least_p %in% paste0("ar", 1:original_p))
    print(paste("adjusted_p in least p is: ", adjusted_p))
    adjusted_q <- original_q - sum(overfit_coeffs_least_p %in% paste0("ma", 1:original_q))
    print(paste("adjusted_q in least p is: ", adjusted_q))

    exclude_intercept <- "intercept" %in% overfit_coeffs_least_p
    print(paste("exclude_intercept in least p is: ", exclude_intercept))
    current_data <- order_comparison[[i]]$data

    refitted_model <- Arima(current_data, order=c(adjusted_p, MAX_D, adjusted_q), include.mean=!exclude_intercept)
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
    # print(paste("overfit_coeffs_auto_arima ", overfit_coeffs_auto_arima))
    # print(paste("original_p ", p))
    # print(paste("original_q ", q))

    p <- p - sum(overfit_coeffs_auto_arima %in% paste0("ar", 1:p))
    q <- q - sum(overfit_coeffs_auto_arima %in% paste0("ma", 1:q))

    exclude_intercept <- "intercept" %in% overfit_coeffs_auto_arima

    current_data <- order_comparison[[i]]$data
    # print(paste("c(adjusted_p, 0, adjusted_q) ", c(adjusted_p, 0, adjusted_q)))
    print(paste("adjusted_p ", p))
    print(paste("adjusted_q ", q))
    refitted_model_auto_arima <- Arima(current_data, order=c(p, MAX_D, q), include.mean=!exclude_intercept)
    new_rmse_vals[[i]] <- perform_ts_cross_validation(current_data, refitted_model_auto_arima)[[2]]
    refitted_model_analysis_auto_arima[[i]] <- analyse_significant_parameters(refitted_model_auto_arima)

  }
  avg_new_rms_vals_auto_arima <- mean(unlist(new_rmse_vals))
  return(list(avg_new_rms_vals_auto_arima, refitted_model_analysis_auto_arima))
}

calculate_p_values <- function(auto_arima_model) {
  # Calculate p-values
  sw_p_value <- shapiro.test(auto_arima_model$residuals)$p.value
  lbq_test_p_value <- Box.test(auto_arima_model$residuals, lag = 20, type = "Ljung-Box")$p.value
  t_test_p_value <- t.test(auto_arima_model$residuals)$p.value

  # Return p-values as a list
  list(sw_p_value = sw_p_value, lbq_test_p_value = lbq_test_p_value, t_test_p_value = t_test_p_value)
}

average_p_values <- function(order_comparison) {
  # Calculate the average of each p-value separately
  sw_p_values <<- lapply(order_comparison, function(x) x$p_values_list$sw_p_value)
  lbq_test_p_values <<- lapply(order_comparison, function(x) x$p_values_list$lbq_test_p_value)
  t_test_p_values <<- lapply(order_comparison, function(x) x$p_values_list$t_test_p_value)
  avg_sw_p_value <- mean(unlist(sw_p_values))
  avg_lbq_test_p_value <- mean(unlist(lbq_test_p_values))
  avg_t_test_p_value <- mean(unlist(t_test_p_values))
  # Return the average p-values as a list
  list(avg_sw_p_value = avg_sw_p_value, avg_lbq_test_p_value = avg_lbq_test_p_value, avg_t_test_p_value = avg_t_test_p_value)
}

calc_avg_validation_auto_arima <- function(order_comparison)
{
  avg_validation_score <- mean(unlist(lapply(order_comparison, function(x) x$validation_score)))
  # count_validated_models <- sum(unlist(lapply(order_comparison, function(x) x$validation_score)))
  return(avg_validation_score)

}

assign_validatin_auto_arima <- function(p_values_list) {

  # Check if all p-values are greater than 0.05
  if(all(unlist(p_values_list) > 0.05)) {
    validation_score <- 1
  } else {
    validation_score <- 0
  }

  # Store the binary variable in a list
  validation_score_list <- list(validation_score = validation_score)

  return(validation_score_list)
}

get_info_about_models <- function (random_arima_data_for_testing)
{
  model_list <- fit_arima_models(random_arima_data_for_testing, max_p = MAX_P, max_q = MAX_D)

  least_p_model <- find_best_model_max_p(model_list)
  auto_arima_model <- auto.arima(random_arima_data_for_testing, allowdrift = TRUE, approximation = FALSE, stepwise = FALSE )
  auto_arima_model$order <- c(auto_arima_model$arma[[1]], auto_arima_model$arma[[6]], auto_arima_model$arma[[2]])

  p_values_list <- calculate_p_values(auto_arima_model)
  validation_score <- assign_validatin_auto_arima(p_values_list)




  actual_order <- model_and_data[[2]]$order
  auto_arima_order <- auto_arima_model$order
  least_p_order <- least_p_model$order

  avg_rmse_least_p_model <- perform_ts_cross_validation(random_arima_data_for_testing,least_p_model$best_model_min_p)
  avg_rmse_auto_arima_model <- perform_ts_cross_validation(random_arima_data_for_testing,auto_arima_model)

  models_analysis_least_p <- analyse_significant_parameters(least_p_model$best_model_min_p)
  models_analysis_auto_arima <- analyse_significant_parameters(auto_arima_model)

  # Run the analysis once and store the result
  # arima_analysis_result <- run_arima_analysis(random_arima_data_for_testing)
  split <- split_data(random_arima_data_for_testing)
  forecast_least_p <- forecast::forecast(least_p_model$best_model_min_p, length(split$test))
  print("forecast_tes working")
  forecast_auto_arima <- forecast::forecast(auto_arima_model, length(split$test))
  # Extract the count of insignificant parameters for both models
  # count_insig_params_p_model <- arima_analysis_result$`Least p-value Model`$overfit_analysis$count_less_significant_parameters
  performance_least_p <- evaluate_performance(split$test, forecast_least_p$mean)
  performance_auto_arima <- evaluate_performance(split$test,forecast_auto_arima$mean)

  count_insig_params_p_model <- models_analysis_least_p$count_less_significant_parameters

  # count_insig_params_auto_arima_model <- arima_analysis_result$`Auto.arima Model`$overfit_analysis$count_less_significant_parameters
  count_insig_params_auto_arima_model <- models_analysis_auto_arima$count_less_significant_parameters


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
    models_analysis_auto_arima = models_analysis_auto_arima,

    performance_least_p = performance_least_p,
    performance_auto_arima = performance_auto_arima,

    validation_score = validation_score,
    p_values_list = p_values_list
  )
  # }
  return(comparison_results)

}



# Let's breakdown count_order_matches



gen_avg_values_for_orgl_and_refit <- function (order_comparison, model_and_data){ # Counting TRUE values for auto_arima_order_match
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
  actual_order <<- model_and_data[[2]]$order
  # Create a DF:
  test_data <- data.frame(
    actual_order = paste(model_and_data[[2]]$order, collapse = " "),
    # actual_order = 204,
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
  output_file_path <- OUTPUT_PATH
  write.table(test_data, file = output_file_path, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

  invisible("Process complete")
}
append_to_file <- function(avg_p_vals, avg_validation_score, output_file_path) {
  # Read the file as a vector of lines
  lines <- readLines(output_file_path)

  # Extract the values from the list
  avg_sw_p_value <- avg_p_vals$avg_sw_p_value
  avg_lbq_test_p_value <- avg_p_vals$avg_lbq_test_p_value
  avg_t_test_p_value <- avg_p_vals$avg_t_test_p_value

  # Add the new values to the last line
  last_line <- lines[length(lines)]
  new_line <- paste(last_line, avg_sw_p_value, avg_lbq_test_p_value, avg_t_test_p_value, avg_validation_score, sep = ",")

  # Replace the last line with the new line
  lines[length(lines)] <- new_line

  # Write the lines back to the file
  writeLines(lines, output_file_path)
}
# Run the analysis

print("Running ARIMA analysis...")

# Function for parallel processing
run_parallel_processing <- function(num_iters_csv) {
  print("Start pll processing:")
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  all_objects <- ls(globalenv())
  order_comparison <- list()
  model_and_data <- list()
  foreach(i = 1:num_iters_csv, .combine = 'c', .packages = c("forecast"), .multicombine = TRUE, .export = all_objects) %dopar% {
    model_and_data[[i]] <- generate_random_arima(max_p = MAX_P, max_d = MAX_D, max_q = MAX_Q, num_observations = NUM_OBSERVATIONS, seed = NULL)
    random_arima_data_for_testing <- model_and_data[[i]][[1]]
    result <- get_info_about_models(random_arima_data_for_testing)
    result$data <- random_arima_data_for_testing
    result$model_and_data <- model_and_data[[i]]
    order_comparison[[i]] <- list(result)
  }
  stopCluster(cl)
  gen_avg_values_for_orgl_and_refit(order_comparison, model_and_data)
  return(list(order_comparison = order_comparison, model_and_data = model_and_data))
}


# Function for sequential processing
run_sequential_processing <- function(num_iters_csv, MAX_P, MAX_D, MAX_Q, NUM_OBSERVATIONS) {
  order_comparison <- list()
  for(i in seq_along(1:num_iters_csv)) {
    model_and_data <- generate_random_arima(max_p = MAX_P, max_d = MAX_D, max_q = MAX_Q, num_observations = NUM_OBSERVATIONS, seed = NULL)
    print(paste("model_and_data is: ", model_and_data))
    random_arima_data_for_testing <- model_and_data[[1]]
    order_comparison[[i]] <- get_info_about_models(random_arima_data_for_testing)
    order_comparison[[i]]$data <- random_arima_data_for_testing
    order_comparison[[i]]$model <- model_and_data[[2]]$order
    # order_comparison[[i]]$model_and_data <- model_and_data
  }
  gen_avg_values_for_orgl_and_refit(order_comparison, model_and_data)
  # process <- lapply(order_comparison, function(x) gen_avg_values_for_orgl_and_refit(order_comparison, x$model_and_data))
  return(order_comparison)
}

# Function to choose between parallel and sequential processing
run_processing <- function(use_parallel_processing, num_iters_csv, MAX_P, MAX_D, MAX_Q, NUM_OBSERVATIONS) {
  if (use_parallel_processing) {
    run_parallel_processing(num_iters_csv)
  } else {
    run_sequential_processing(num_iters_csv, MAX_P, MAX_D, MAX_Q, NUM_OBSERVATIONS)
  }
}

num_iters_csv <- NUM_TEST_MODELS
order_comparison <- run_processing(use_parallel_processing = FALSE, num_iters_csv, MAX_P, MAX_D, MAX_Q, NUM_OBSERVATIONS)

avg_p_vals <- average_p_values(order_comparison)
avg_validation_score <- calc_avg_validation_auto_arima(order_comparison)

# Read the file as a vector of lines
output_file_path <- OUTPUT_PATH
append_to_file(avg_p_vals, avg_validation_score, output_file_path)