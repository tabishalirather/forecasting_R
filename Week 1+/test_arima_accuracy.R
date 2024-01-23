# clear the console.
rm(list = ls())
# load required libraries
library(forecast)
library(astsa)
library(doParallel)
library(foreach)
library(progress)
# define global variables
NUM_TESTS <- 10
NUM_OBSERVATIONS <- 1000
IC <- "aic"  # Can be "aic", "aicc", "bic"
# Define the base seed for reproducibility as a random number.
BASE_SEED <- as.integer(Sys.time())
print(paste("BASE_SEED:", BASE_SEED))



# Fxn to generate random arima data
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

# DoneTODO: automate this testing. Conduct 100 tests for different number of observations and record how well auto.arima performs

# test_arima_accuracy <- function(num_tests, max_p, max_d, max_q, num_observations) {
#   true_postive_count <- 0
#   actual_orders <- list()
#   chosen_orders <- list()
#   comparison_results <- list()
#
#   for (i in 1:num_tests)
#   {
#     iter_seed <- BASE_SEED + i
#     model_and_data <- generate_random_arima(max_p, max_d, max_q, num_observations, iter_seed)
#     actual_arima_order <- model_and_data[[2]]$order
#     random_arima_data <- model_and_data[[1]]
#     best_auto_arima_model <<- auto.arima(random_arima_data, max.p = max_p, max.d = max_d, max.q = max_q, approximation = FALSE, stepwise = FALSE, trace = TRUE, ic = "bic")
#     chosen_arima_order <<- c(best_auto_arima_model$arma[1], best_auto_arima_model$arma[6], best_auto_arima_model$arma[2])
#   # print(i)
#     actual_orders[[i]] <- actual_arima_order
#     chosen_orders[[i]] <- chosen_arima_order
#     comparison_results[[i]] <- list(actual_order = actual_arima_order, chosen_order = chosen_arima_order)
#   }
#  true_positive_count <- sum(mapply(identical, actual_orders, chosen_orders))
#   # print(true_positive_count)
#   accuracy_percentage <-  (true_positive_count / num_tests) * 100
#
#   return(list("true_positive_count" = true_positive_count,
#               "accuracy_percentage" = accuracy_percentage,
#               "actual_orders" = actual_orders,
#               "chosen_orders" = chosen_orders,
#               "comparison_results" = comparison_results))
# }

# comp_test_arima_accuracy <- function ()
# {
#   print("This is being called")
# start_time <- Sys.time()
# result <- test_arima_accuracy(num_tests = NUM_TESTS, max_p = 3, max_d = 1, max_q = 3, num_observations = NUM_OBSERVATIONS)
# end_time <- Sys.time()
# time_taken <- as.numeric(end_time - start_time, units = "secs")
#   print(paste0(result$accuracy_percentage, "% of the time, auto.arima chose the correct order"))
# print(paste0("Time taken to test in parallel: ", time_taken))
# comparison_results <- result$comparison_results
# }

initiate_parallel_processing <- function ()
{
  # parallel accuracy testing for auto.arima fxn
num_cores <- detectCores() - 1  # Leave one core free for system processes
# Create a cluster and register it
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Export the function and required variables to the cluster
clusterExport(cl, varlist = c("generate_random_arima", "BASE_SEED"))
  return(list("cl" = cl))
}

from_initiate_parallel_processing <- initiate_parallel_processing()

# Fxn to run the test using parallel processing.
test_arima_accuracy_parallel <- function(num_tests, max_p, max_d, max_q, num_observations, ic)
{
  print("I am being called")
  results<- foreach(i = 1:num_tests, .packages = c("forecast", "astsa")) %dopar%
  {
    iter_seed <- BASE_SEED + i
    model_and_data <- generate_random_arima(max_p, max_d, max_q, num_observations, iter_seed)
    actual_arima_order <- model_and_data[[2]]$order
    random_arima_data <- model_and_data[[1]]
    best_auto_arima_model <- auto.arima(random_arima_data, max.p = max_p, max.d = max_d, max.q = max_q, approximation = FALSE, stepwise = FALSE, trace = TRUE, ic=ic)
    chosen_arima_order <- c(best_auto_arima_model$arma[1], best_auto_arima_model$arma[6], best_auto_arima_model$arma[2])
    return(list(actual_order = actual_arima_order, chosen_order = chosen_arima_order, model = best_auto_arima_model))
  }

  true_postive_count <- sum(mapply(function(actual, chosen) identical(actual, chosen),
                                   actual_orders <- lapply(results, `[[`, "actual_order"),
                                   chosen_orders <- lapply(results, `[[`, "chosen_order")))

  error_percentage <- (1 - (true_postive_count / num_tests)) * 100

  return(list("true_positive_count" = true_postive_count,
              "error_percentage" = error_percentage,
              "actual_orders" = actual_orders,
              "chosen_orders" = chosen_orders,
              "comparison_results_pll" = results))
}


log_message <- function(message, filepath) {
  cat(message, file = filepath, append = TRUE, sep = "\n")
}


run_parallel_accuracy_test <- function()
{
  log_file <- "C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/TestingR/Week 1+/accuracy_test_results.txt"  # Define the log file name
  print("Performing accuracy test in parallel.....")
  start_time <- Sys.time()
  result_pll <- test_arima_accuracy_parallel(num_tests = NUM_TESTS, max_p = 3, max_d = 1, max_q = 3, num_observations = NUM_OBSERVATIONS, ic = IC)
  end_time <- Sys.time()
  time_taken_parallel <- as.numeric(end_time - start_time, units = "secs")
  # result_pll
  stopCluster(from_initiate_parallel_processing$cl)
  # print(result_pll)
  accuracy_decimal <- result_pll$error_percentage
  num_observations <- NUM_OBSERVATIONS
  num_tests <- NUM_TESTS
  ic_used <- IC
  # comparison_results_pll <<- result_pll$comparison_results_pll
  # Create a data frame to store the results

  results_df <- data.frame(
      accuracy = accuracy_decimal,
      num_observations = num_observations,
      num_tests = num_tests,
      ic_used = ic_used,
      stringsAsFactors = FALSE
    )
  # Write the results to the log file
  write.table(results_df, file = log_file, append = TRUE, row.names = FALSE, sep = " | ")

  print(paste0(result_pll$error_percentage, "% of the time, auto.arima chose the wrong order"))
  print(paste0(result_pll$true_positive_count, " out of ", NUM_TESTS, " tests were successful"))
  print(paste0("Time taken to test in parallel: ", time_taken_parallel))

  results_from_paralllel_test <- result_pll

}
results_pll <- run_parallel_accuracy_test()
# TODO: Look into which part of ARIMA does autoarima get wrong most often. Is it the AR or MA part?

compare_orders <- function(comparison_results_pll) {
  diff_count <- 0
  ar_diff_count <- 0
  diff_diff_count <- 0
  ma_diff_count <- 0
  # Iterate over the comparison_results_pll list
  for(i in seq_along(comparison_results_pll)) {
    actual_order <- comparison_results_pll[[i]]$actual_order
    chosen_order <- comparison_results_pll[[i]]$chosen_order

    # Check if the actual_order and chosen_order are not the same
    if(!identical(actual_order, chosen_order)) {
      # Compare each element of the orders
      print("-------------------------------------=")
      for(j in seq_along(actual_order)) {
        if(actual_order[j] != chosen_order[j]) {
          # Print out the differences
          # print("New line")
          if(j == 1) {
            print(paste("Testing", i, "AR Part: Actual =", actual_order[j], ", Chosen =", chosen_order[j]))
            ar_diff_count <- ar_diff_count + 1
          } else if(j == 2) {
            print(paste("Testing", i, "Differencing Part: Actual =", actual_order[j], ", Chosen =", chosen_order[j]))
            diff_diff_count <- diff_diff_count + 1
          } else if(j == 3) {
            print(paste("Testing", i, "MA Part: Actual =", actual_order[j], ", Chosen =", chosen_order[j]))
            ma_diff_count <- ma_diff_count + 1
          }
        }
      }
    }
  }
  diff_count <- list(ar_diff_count, diff_diff_count, ma_diff_count)
  return(diff_count)
}

# full_list_comparison <- list()
# for(i in 1:5) {
#   full_list_comparison[[i]] <-compare_orders(comparison_results_pll)
#
# }


# Done_TODO: Compare the the results of this functions for different ics
# Now, let's write the manual hyndman, khandakar algo to test how well that works.
# TODO: now use the test to compare the accuracy of manual auto.arima with different number of observations
# Now, let's write the manual hyndman, khandakar algo to test how well that works.

# random_arima_data <<- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = 1000, NULL)[[1]]
# random_arima_model <<- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = 1000, NULL)[[2]]

fit_arima_models <- function(random_arima_data, max_p, max_q) {
  model_list <- list()
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      model_list[[paste0("ARIMA_", p, "_1_", q)]] <- arima(random_arima_data, order = c(p, 1, q))
    }
  }
  return(model_list)
}

  AICC <- function(model) {
  n <- length(model_list$residuals)
  k <- length(coef(model))
  AIC(model) + 2*k*(k+1) / (n-k-1)
}
find_optimal_model_ic <- function(model_list) {
  aic_array <- sapply(model_list, function(model) model$aic)
  # model$bic <- sapply(model_list, function(model) BIC(model))
  bic_array <- sapply(model_list, function(model) model$bic <- BIC(model))
  aicc <- sapply(model_list, function (model) AICC(model))
  # model$aicc <- aicc
  aicc_array <- sapply(model_list, function(model) model$aicc <- AICC(model))

  optimal_model_aic <- model_list[[which.min(aic_array)]]
  optimal_model_bic <- model_list[[which.min(bic_array)]]
  optimal_model_bic$bic <- min(bic_array)
  optimal_model_aicc <- model_list[[which.min(aicc_array)]]
  optimal_model_aicc$aicc <- min(aicc_array)
  overall_optimal_model <- model_list[[which.min(aic_array + bic_array)]]

  return(list(optimal_model_aic = optimal_model_aic,
              optimal_model_bic = optimal_model_bic,
              optimal_model_aicc = optimal_model_aicc,
              overall_optimal_model = overall_optimal_model))
}
# This fxn find separate optimal optimal models for each ic. Done_TODO: modify this to find the optimal model for which the sum of ics is the least
#TODO:Consider using parameter count to check for overfitting.
# TODO: test the accuracy of overall optimal arima model and other aic, and bic models, and comapre with auto.arima
chosen_order_aic_list <- list()
chosen_order_overall_list <- list()
actual_order_list <- list()
match_count <- 0
comp_ord <- list()
# for(i in 1:30){
#   # Generate new ARIMA data for each iteration
#   random_arima_data <<- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = 1000, NULL)[[1]]
#   random_arima_model <<- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = 1000, NULL)[[2]]
#
#   model_list <- fit_arima_models(random_arima_data, max_p = 3, max_q = 3)
#   actual_arima_order <- random_arima_model$order
#   optimal_models <- find_optimal_model(model_list)
#   chosen_arima_order_overall <- c(optimal_models$overall_optimal_model$arma[1], optimal_models$overall_optimal_model$arma[6], optimal_models$overall_optimal_model$arma[2])
#   chosen_arima_order_aic <- c(optimal_models$optimal_model_aic$arma[1], optimal_models$optimal_model_aic$arma[6], optimal_models$optimal_model_aic$arma[2])
#   chosen_order_aic_list[[i]] <- chosen_arima_order_aic
#   chosen_order_overall_list[[i]] <- chosen_arima_order_overall
#   actual_order_list[[i]] <- actual_arima_order
#   # Check if chosen order is equal to actual order and keep count of that.
#   print(paste("Overall optimal model:", paste(chosen_arima_order_overall, collapse = " ")))
#   print(paste("AIC optimal model:", paste(chosen_arima_order_aic, collapse = " ")))
#   print(paste("Actual optimal model:", paste(actual_arima_order, collapse = " ")))
#   sublist <- list(
#    chosen_order_aic_list[[i]],
#     chosen_order_overall_list[[i]],
#     actual_order_list[[i]]
#   )
#   comp_ord[i] <- sublist
#   print("   ")
#   if(identical(chosen_arima_order_overall, actual_arima_order )) {
#     match_count <- match_count + 1
#   }
# }
#   print(match_count)
# TODO: Look for how much the difference is in the arimas that are chosen by auto.arima and the actual arima model.
comparison_results_pll <-  results_pll$comparison_results_pll
compare_ar_ma_d_accuracy <- function(comparison_results_pll) {
  ar_error_count <- 0
  ma_error_count <- 0
  d_error_count <- 0
  total_models <- length(comparison_results_pll)

  for (i in seq_along(comparison_results_pll)) {
    actual_order <- comparison_results_pll[[i]]$actual_order
    chosen_order <- comparison_results_pll[[i]]$chosen_order

    # Compare AR part
    if (actual_order[1] != chosen_order[1]) {
      ar_error_count <- ar_error_count + 1
    }

    # Compare MA part
    if (actual_order[3] != chosen_order[3]) {
      ma_error_count <- ma_error_count + 1
    }

    #Compare differencing part
    if (actual_order[2] != chosen_order[2]) {
      d_error_count <- d_error_count + 1
    }
  }

  ar_error_percentage <- (ar_error_count / total_models) * 100
  ma_error_percentage <- (ma_error_count / total_models) * 100
  d_error_percentage <- (d_error_count / total_models) * 100

  return(list(ar_error_percentage = ar_error_percentage,
              ma_error_percentage = ma_error_percentage,
              d_error_percentage = d_error_percentage))
}

# Apply this function to your comparison results
# # ar_ma_comparison <- compare_ar_ma_d_accuracy(comparison_results_pll)
# print(paste0("AR error percentage: ", ar_ma_comparison$ar_error_percentage))
# print(paste0("MA error percentage: ", ar_ma_comparison$ma_error_percentage))
# print(paste0("Differencing error percentage: ", ar_ma_comparison$d_error_percentage))
# AR_count_error = 1 + 1 + 0 + 1 + 1 + 0 + 1 +
# MA_count_error =

# Fxn to get p_values for each coefficient`
get_arima_p_values <- function(model) {
    # Extract coefficients and their standard errors
    coefficients <- coef(model)
    std_errors <- sqrt(diag(vcov(model)))
    # TODO: How to estimate standard errors for a model.
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

# Now anaylze p_values
analyze_p_values <- function(p_values_list, significance_level = 0.05) {
    less_significant_models <- list()

    for (i in seq_along(p_values_list)) {
        p_values <- p_values_list[[i]]$p_value
        model <- p_values_list[[i]]$model

        # Check if any p-value exceeds the significance level
        if (any(p_values > significance_level)) {
            # Flag this model as having less significant coefficients
            less_significant_models[[length(less_significant_models) + 1]] <- list(
                model_index = i,
                model = model,
                p_values = p_values
            )
        }
    }

    return(less_significant_models)
}




check_overfitting <- function(comparison_results_pll) {
  overfit_models_order <- list()
  overfit_models_residuals <- list()
  overfit_count <- 0
  total_models <- length(comparison_results_pll)
  not_overfit_models <- list()
  p_values_list <- list()

  for (i in seq_along(comparison_results_pll)) {
    # print(class(comparison_results_pll[[i]]))
    # print(names(comparison_results_pll[[i]]))

    actual_order <- comparison_results_pll[[i]]$actual_order
    chosen_order <- comparison_results_pll[[i]]$chosen_order
    model <- comparison_results_pll[[i]]$model # Assuming the ARIMA model object is stored here

    # Order comparison for overfitting
    order_overfit <- sum(chosen_order > actual_order) > 0
    print(paste("Order overfit: ", order_overfit))

    # Residual analysis (Ljung-Box test)
    residuals_analysis <<- Box.test(residuals(model), type = "Ljung-Box")
    residual_overfit <- residuals_analysis$p.value < 0.05 #Testing the asssumption that if the model is overfit, the residuals will be correlated, they'll be capturing the noise in the data.
    print(paste("Residual overfit: ", residual_overfit))
    # P-value check for coefficients
    p_values_list[[i]] <- get_arima_p_values(model)
    # print(paste("P-values: ", model_p_values$p_value, "for model" , model_p_values$model))


    # Determine if the model is overfitting
    #
    # print(paste("Order overfit: ", order_overfit))
    # print(paste("Residual overfit: ", residual_overfit))
    #
    # print(paste0("Model ", i, " is overfit: ", is_overfit))

    if (order_overfit) {
      overfit_models_order <- c(overfit_models_order, list(comparison_results_pll[[i]]))
    }
    if(residual_overfit) {
      overfit_models_residuals <- c(overfit_models_residuals, list(comparison_results_pll[[i]]))
    }
    is_overfit <- order_overfit || residual_overfit #|| coef_overfit
    if(is_overfit){
      overfit_count <- overfit_count + 1
    }
  }

  overfit_percentage <- (overfit_count / total_models) * 100

  return(list(
    overfit_count = overfit_count,
    overfit_percentage = overfit_percentage,
    total_models = total_models,
    overfit_models_order = overfit_models_order,
    overfit_models_residuals = overfit_models_residuals,
    residual_overfit = residual_overfit,
    order_overfit =  order_overfit,
    p_values_list = p_values_list
  ))
}


overfitting_test_results <- check_overfitting(comparison_results_pll)
p_values_list <- overfitting_test_results$p_values_list

# Analyze p-values
less_significant_models <- analyze_p_values(p_values_list)
# TODO Standard errors for the model using our model.
#  TODO If we cannot calculate standard errors, then how to check for overiftting and calcualtion of p_values
# TODO:

# get_arima_p_values(model)

random_arima_data <- generate_random_arima(max_p = 3, max_d = 1, max_q = 3, num_observations = 1000, NULL)[[1]]
model_list <- fit_arima_models(random_arima_data, max_p = 3, max_q = 3)
find_best_model <- function(model_list) {
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

best_model <- find_best_model(model_list)

auto_arima_model <- auto.arima(random_arima_data, allowdrift = TRUE, approximation = FALSE, stepwise = FALSE )