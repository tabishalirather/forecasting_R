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
IC <- "aicc"  # Can be "aic", "aicc", "bic" or "hqic"
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

test_arima_accuracy <- function(num_tests, max_p, max_d, max_q, num_observations) {
  true_postive_count <- 0
  actual_orders <- list()
  chosen_orders <- list()
  comparison_results <- list()

  for (i in 1:num_tests)
  {
    iter_seed <- BASE_SEED + i
    model_and_data <- generate_random_arima(max_p, max_d, max_q, num_observations, iter_seed)
    actual_arima_order <- model_and_data[[2]]$order
    random_arima_data <- model_and_data[[1]]
    best_auto_arima_model <<- auto.arima(random_arima_data, max.p = max_p, max.d = max_d, max.q = max_q, approximation = FALSE, stepwise = FALSE, trace = TRUE, ic = "bic")
    chosen_arima_order <<- c(best_auto_arima_model$arma[1], best_auto_arima_model$arma[6], best_auto_arima_model$arma[2])
  # print(i)
    actual_orders[[i]] <- actual_arima_order
    chosen_orders[[i]] <- chosen_arima_order
    comparison_results[[i]] <- list(actual_order = actual_arima_order, chosen_order = chosen_arima_order)
  }
 true_positive_count <- sum(mapply(identical, actual_orders, chosen_orders))
  # print(true_positive_count)
  accuracy_percentage <- (true_positive_count / num_tests) * 100

  return(list("true_positive_count" = true_positive_count,
              "accuracy_percentage" = accuracy_percentage,
              "actual_orders" = actual_orders,
              "chosen_orders" = chosen_orders,
              "comparison_results" = comparison_results))
}

comp_test_arima_accuracy <- function ()
{
  print("This is being called")
start_time <- Sys.time()
result <- test_arima_accuracy(num_tests = NUM_TESTS, max_p = 3, max_d = 1, max_q = 3, num_observations = NUM_OBSERVATIONS)
end_time <- Sys.time()
time_taken <- as.numeric(end_time - start_time, units = "secs")
  print(paste0(result$accuracy_percentage, "% of the time, auto.arima chose the correct order"))
print(paste0("Time taken to test in parallel: ", time_taken))
comparison_results <- result$comparison_results
}

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
  results<- foreach(i = 1:num_tests, .packages = c("forecast", "astsa")) %dopar%
  {
    iter_seed <- BASE_SEED + i
    model_and_data <- generate_random_arima(max_p, max_d, max_q, num_observations, iter_seed)
    actual_arima_order <- model_and_data[[2]]$order
    random_arima_data <- model_and_data[[1]]
    best_auto_arima_model <- auto.arima(random_arima_data, max.p = max_p, max.d = max_d, max.q = max_q, approximation = FALSE, stepwise = FALSE, trace = TRUE, ic=ic)
    chosen_arima_order <- c(best_auto_arima_model$arma[1], best_auto_arima_model$arma[6], best_auto_arima_model$arma[2])
    return(list(actual_order = actual_arima_order, chosen_order = chosen_arima_order))
  }

  true_postive_count <- sum(mapply(function(actual, chosen) identical(actual, chosen),
                                   actual_orders <- lapply(results, `[[`, "actual_order"),
                                   chosen_orders <- lapply(results, `[[`, "chosen_order")))

  accuracy_percentage <- (true_postive_count / num_tests) * 100

  return(list("true_positive_count" = true_postive_count,
              "accuracy_percentage" = accuracy_percentage,
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
  stopCluster(from_initiate_parallel_processing$cl)

  accuracy_decimal <- result_pll$accuracy_percentage / 100
  num_observations <- NUM_OBSERVATIONS
  num_tests <- NUM_TESTS
  ic_used <- IC

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

  print(paste0(result_pll$accuracy_percentage, "% of the time, auto.arima chose the correct order"))
  print(paste0(result_pll$true_positive_count, " out of ", NUM_TESTS, " tests were successful"))
  print(paste0("Time taken to test in parallel: ", time_taken_parallel))
  comparison_results_pll <<- result_pll$comparison_results_pll
}
run_parallel_accuracy_test()
# TODO: Compare the the results of this functions for different ics
# TODO: now use the test to compare the accuracy of manual auto.arima with different number of observations
