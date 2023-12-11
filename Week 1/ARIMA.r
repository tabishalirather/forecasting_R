print("this is working")
library(xts)
library(forecast)
library(tseries)
# arima_r <- function (){
stock_data <- read.csv("C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/Data/ASX.AX.csv", header = TRUE)
duplicates <- duplicated(stock_data)
stock_data <- stock_data[!duplicated(stock_data), ]
class(stock_data$Date)
stock_data$Date <- as.Date(stock_data$Date, format = "%Y-%m-%d")
class(stock_data$Date)
#plot the data
# plot(stock_data$Date, stock_data$Open, type='l')
ts_stock_data <- xts(stock_data$Open, order.by = stock_data$Date)
# dev.off()?

# par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(ts_stock_data)
# null hypothesis in kpss is that the series is stationary. p < 0.05, we reject the null hypothesis, i.e time series is not stationary and we need to difference it to make it stationary.
kpss_test <- kpss.test(stock_data$Open, null = "Trend")
print(kpss_test)
# since p value is less than 0.05, we cannot rejec we'll difference the data to make it strationary
# # Arima(ts_stock_data, order = c(1,1,1)
# differenced_stock_data_open <- diff(ts_stock_data)
ts_stock_data <- ts_stock_data
ts_stock_data <- na.omit(ts_stock_data)
plot(ts_stock_data)


kpss_test <- kpss.test(ts_stock_data, null = "Trend")
print(kpss_test)
# After differencing, p_value is 0.1 > 0.01, which means we don't  need to difference the data again

# ACF Plot
Acf(ts_stock_data, lag.max = 20, main = "ACF for Differenced Data",
    col = "blue", lwd = 2)

Pacf(ts_stock_data, lag.max = 20, main = "PACF for Differenced Data",
     col = "green", lwd = 2)

# Now, let's fit various ARIMA models to our data.

model_list <- list()
max_p <- 3
max_q <- 3
for (p in 0:max_p)
{
  for (q in 0:max_q)
  {
  # model_names <- paste0("ARIMA_", p, "_1_", q)
  model_list[[paste0("ARIMA_", p, "_1_", q)]] <- Arima(ts_stock_data, order = c(p, 1, q))
    # model_names <- array(model_names)
}
  }
  aic_array <- array(NaN, dim = length(model_list))
  model_names <- array(NaN, dim = length(model_list))
for (index in seq_along(model_list))
{
  # print(summary(model))
  model <- model_list[[index]]
  model_name <- names(model_list[index])
  model_names[index] <- model_name
  aic_array[index] <- model$aic

  # print(index)
  print(model_name)
  # print(model$aic)
  print(aic_array)
  # aic_array  <- array(model$aic, dim = length(aic_array)+1)
}
  least_aic <- min(aic_array)
  index_of_least_aic <- which.min(aic_array)
  print(paste(least_aic, ",", index_of_least_aic))

# for (model in model_list)
# {
#   if(model$aic == least_aic)
#   {
#     print(model_names[index_of_least_aic])
#     # print(model_names)
#     print(model)
#     break
#   }
# }
# print(paste("Best model is:", model_names[index_of_least_aic]))
print(( model_list[model_names[index_of_least_aic]]))
best_fit <- model_list[model_names[index_of_least_aic]]
best_fit_test <- model_list[model_names[index_of_least_aic]]

best_fit <- best_fit[[1]]
residuals_best <- best_fit$residuals
checkresiduals(residuals_best)
autoplot(best_fit)
autoplot(forecast(best_fit))

auto_arima_model <- auto.arima(ts_stock_data, allowdrift = TRUE, approximation = FALSE )
autoplot(forecast(auto_arima_model))
autoplot(auto_arima_model)
print(paste('auto_arima chosen model',auto_arima_model, auto_arima_model$aicc))

print(paste('my manually chosen arima', best_fit, best_fit$aicc))
# }
# arima_r()
# let's calculate the residuals for each model and find p_values for swtest, lbqtest and t-test and try to minizie their sum
p_sw_arr <- array(NaN, dim = length(model_list))
p_lbq_arr <- array(NaN, dim = length(model_list))
p_ttest_arr <- array(NaN, dim = length(model_list))

for (index in seq_along((model_list)))
{
  model <- model_list[[index]]
  residuals <- model$residuals
  p_sw_arr[index] <- shapiro.test(residuals)$p.value
  p_lbq_arr[index] <- Box.test(residuals, lag = 20, type = "Ljung-Box")$p.value
  p_ttest_arr[index] <- t.test(residuals)$p.value
}
arr_sum_p_vals <- p_sw_arr + p_lbq_arr + p_ttest_arr
# find min of arr_sum_p_vals
min_p_val <- min(arr_sum_p_vals)
# find index of min of arr_sum_p_vals
index_of_min_p_val <- which.min(arr_sum_p_vals)
print(paste("min p val is:", min_p_val, "at index:", index_of_min_p_val))
# using the min index, find the model in model_list
best_model_min_p <- model_list[[index_of_min_p_val]]
# p_test_val <- shapiro.test(best_fit$residuals)$p.value
# cor_matrix <-
print(model_list)
# Looks like miniziming the sum of p_values is not such a good idea.
