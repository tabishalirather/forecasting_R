# This program explores the idea of overfitting of models.
library(tseries)
library(dplyr)
library(urca)
library(forecast)
stock_data <- read.csv("C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/Data/02_02_2019_o_10_06_2021_MSFT data.csv", header = TRUE)
stock_data <- stock_data[order(as.Date(stock_data$Date)),]
#Training set:
# [1] "Number of observations in training set: 505 (82.11%)"
# Assuming 'stock_data' is the variable holding the entire dataset
# stock_data$daily_returns <- 0
stock_data$daily_returns <- c(0.015770969 , diff(stock_data$Adj.Close) / head(stock_data$Adj.Close,-1))
training_set_size <- floor(0.8211 * length(stock_data$daily_returns))
training_set <- stock_data$daily_returns[1:training_set_size]
# //First element of training set maybe incorrect
testing_set <- stock_data$daily_returns[(training_set_size + 1):length(stock_data$daily_returns)]
detect_significant_lags_acf <- function(time_series, lag.max = NULL, confidence_level = 0.95) {
  N <- length(na.omit(time_series)) # Effective sample size
  acf_values <- acf(time_series, lag.max = lag.max, plot = TRUE, main = "ACF for Training Set")

  # The standard error under the assumption of white noise is 1/sqrt(N)
  standard_error <- 1 / sqrt(N)

  # Calculate the critical value for the given confidence level
  critical_value <- qnorm((1 + confidence_level) / 2) * standard_error

  # Find lags where the absolute ACF value is greater than the critical value
  significant_lags <- which(abs(acf_values$acf) > critical_value)

  # Remove the first element (lag 0) since we are interested in lags greater than 0
  significant_lags_acf <- significant_lags[significant_lags > 1] - 1

  # Return the significant lags
  return(significant_lags_acf)
}
detect_significant_lags_pacf <- function(time_series, lag.max = NULL, confidence_level = 0.95) {
  N <- length(na.omit(time_series)) # Effective sample size
  pacf_values <- pacf(time_series, lag.max = lag.max, plot = TRUE, main = "PACF for Training Set")

  # The standard error under the assumption of white noise is 1/sqrt(N)
  standard_error <- 1 / sqrt(N)

  # Calculate the critical value for the given confidence level
  critical_value <- qnorm((1 + confidence_level) / 2) * standard_error

  # Find lags where the absolute ACF value is greater than the critical value
  significant_lags <- which(abs(pacf_values$acf) > critical_value)

  # Remove the first element (lag 0) since we are interested in lags greater than 0
  significant_lags_pacf <- significant_lags[significant_lags > 1] - 1

  # Return the significant lags
  return(significant_lags_pacf)
}
# Example usage:
signif_lags_acf <- detect_significant_lags_acf(training_set, lag.max = 250)
signif_lags_pacf <- detect_significant_lags_pacf(training_set, lag.max = 250)

# print(signif_lags_acf)
# print(signif_lags_pacf)
training_set <- ts(training_set, frequency = 252)
summary_training_set <- summary(training_set)
decomp_training_set <- decompose(training_set)
plot(decomp_training_set)
# stl_decomposition <- stl(training_set, s.window = "periodic")
# plot(stl_decomposition)
# Now, let's do stationarity testing.
# We peroform tests to check if the time series is stationary(unit root process) or not.  ADF test has null hypothesis as unit root exists and alternative is
# 1) None, zero mean statinary 2) Non-zero mean stationarity, 3) trend stationarity.
# As we previously seen, there is no any upwards/downwards trend in returns, so we won’t test against trend stationarity (for our alternative hypothesis Ha). As returns are oscillating around zero, we won’t test it against a non-zero mean stationarity as well.
# So for our alternative hypothesis we’ll chose parameter None - zero-means stationarity (either an intercept nor a trend is included in the test regression).
adf_result <- adf.test(training_set, alternative = "stationary")
# print(adf_result)

# //manual adf test
manual_adf_test <- function(time_series)
{
	manual_adf_data <- data.frame(training_set = time_series)
	manual_adf_data$lagged_series <- stats::lag(manual_adf_data$training_set, -1)
	manual_adf_data$diff_series <- c(NA,diff(manual_adf_data$training_set, differences = 1))
	manual_adf_data$diff_lag1 <- stats::lag(manual_adf_data$diff_series, -1)
	manual_adf_data <- na.omit(manual_adf_data)
	manual_adf_model <- lm(diff_series ~ lagged_series - 1, data = manual_adf_data)
	summary(manual_adf_model)
}
# adf_test_result <- manual_adf_test(training_set)
# print(adf_test_result)
# ers_test_result <- ur.ers(training_set, type = "DF-GLS", model = "constant", lag.max = 10)
# summary(ers_test_result)
optimal_arima <- auto.arima(training_set, stepwise = FALSE, approximation = FALSE, trace = TRUE)
# Find out how to check for autoarima
summary_optimal_arima <- summary(optimal_arima)
# Extract the variance-covariance matrix
var_coef_matrix <- summary_optimal_arima$var.coef
# Compute the standard errors as the square roots of the diagonal elements
#tod00o: find standard errors
optimal_arima$std_errors <- array(sqrt(diag(var_coef_matrix)))
coefficients <- optimal_arima$coef
standard_errors <- optimal_arima$std_errors
# Print the standard errors
# print(paste("standard errors are:", std_errors))
t_statistics <- coefficients / standard_errors
p_values <- 2 * (1 - pt(abs(t_statistics), df = length(testing_set) - length(coefficients)))
p_values

#  Now, let's check for overfitting:
overfit_one <- arima(x = testing_set, order= c(2,0,4))
var_coef_matrix <- overfit_one$var.coef
overfit_one$std_errors <- array(sqrt(diag(var_coef_matrix)))

overfit_coef <- overfit_one$coef
print(paste("coeff", overfit_coef))
# overfit_coef <- 0.2487
standard_errors <- overfit_one$std_errors
print(paste("errors", standard_errors))
# standard_errors <- 0.1067
t_stats_overfit <- overfit_coef / standard_errors
print(paste("t_stats", t_stats_overfit))
p_val_overfit <-  2*(1 - pt(abs(t_stats_overfit), df = length(testing_set)-length(overfit_coef)))
p_val_overfit
# Now, let's minize the sum of p values of residuls. SW-test, LBQtest and t-test
#sw-test


