rm(list = ls())
library(forecast)
library(astsa)
library(doParallel)
library(foreach)
print("testing")
# define a global variable
NUM_TESTS <- 10
NUM_OBSERVATIONS <- 1000
BASE_SEED <- as.integer(Sys.time())
print(paste("BASE_SEED:", BASE_SEED))

num_cores <- detectCores() - 1  # Leave one core free for system processes
# Create a cluster and register it
cl <- makeCluster(num_cores)
registerDoParallel(cl)


generate_random_arima <- function(max_p, max_d, max_q, num_observations, seed = NULL) {

  # Setting a seed for reproducibility (optional, can be removed or modified)
  if(!is.null(seed)){
	print(seed)
  	set.seed(seed)
  }
  arima_model <- NULL
  # set.seed(234)
  # Randomly select the values for p, d, and q
  repeat
  {
  p <- sample(0:max_p, 1)
  d <- sample(0:max_d, 1)
  q <- sample(0:max_q, 1)
  #  print(p)
  # print(d)
  # print(q)
  # Generate random AR and MA coefficients
  ar_coefs <- if(p > 0)
  {
  ar_coefs <- runif(p, min = -1, max = 1)
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

  # Generate and return the ARIMA data
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
  # return(value)
}

# make fxns and variables avaialble to the parallell cluster
clusterExport(cl, varlist = c("generate_random_arima", "BASE_SEED"))
# model_and_data <- generate_random_arima(3,1,3,1000)
# auto arima seems to fail for 1000 obersvations
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


    # if (identical(actual_arima_order, chosen_arima_order)) {
    #   true_postive_count <- true_postive_count + 1
    # }

  # return(list(actual_order = actual_orders, chosen_order = chosen_orders))

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
# comp_test_arima_accuracy()
# print(result$true_positive_count)


# 6th is d

# best_auto_arima_model
test_arima_accuracy_parallel <- function(num_tests, max_p, max_d, max_q, num_observations) {
  results<- foreach(i = 1:num_tests, .packages = c("forecast", "astsa")) %dopar%
  {
    iter_seed <- BASE_SEED + i
    model_and_data <- generate_random_arima(max_p, max_d, max_q, num_observations, iter_seed)
    actual_arima_order <- model_and_data[[2]]$order
    random_arima_data <- model_and_data[[1]]
    best_auto_arima_model <- auto.arima(random_arima_data, max.p = max_p, max.d = max_d, max.q = max_q, approximation = FALSE, stepwise = FALSE, trace = TRUE)
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

print("Performing accuracy test in parallel.....")

start_time <- Sys.time()
result_pll <- test_arima_accuracy_parallel(num_tests = NUM_TESTS, max_p = 3, max_d = 1, max_q = 3, num_observations = NUM_OBSERVATIONS)
#
end_time <- Sys.time()
time_taken_parallel <- as.numeric(end_time - start_time, units = "secs")

# Print results
# print(result_pll$true_positive_count)
print(paste0(result_pll$accuracy_percentage, "% of the time, auto.arima chose the correct order"))
print(paste0("Time taken to test in parallel: ", time_taken_parallel))
comparison_results_pll <- result_pll$comparison_results_pll

# Stop the cluster
stopCluster(cl)



