print("this is working")
library(xts)
library(forecast)
library(tseries)
stock_data <- read.csv("C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/Data/ASX.AX.csv", header = TRUE)
duplicates <- duplicated(stock_data)
stock_data <- stock_data[!duplicated(stock_data), ]
class(stock_data$Date)
stock_data$Date <- as.Date(stock_data$Date, format = "%Y-%m-%d")
class(stock_data$Date)
#plot the data
# plot(stock_data$Date, stock_data$Open, type='l')
ts_stock_data <- xts(stock_data$Open, order.by = stock_data$Date)
plot(ts_stock_data)
# null hypothesis in kpss is that the series is stationary. p < 0.05, we reject the null hypothesis, i.e time series is not stationary and we need to difference it to make it stationary.
kpss_test <- kpss.test(stock_data$Open, null = "Trend")
print(kpss_test)
# since p value is less than 0.05, we cannot rejec we'll difference the data to make it strationary
# # Arima(ts_stock_data, order = c(1,1,1)
differenced_stock_data_open <- diff(ts_stock_data)
differenced_stock_data_open <- na.omit(differenced_stock_data_open)
plot(differenced_stock_data_open)


kpss_test <- kpss.test(differenced_stock_data_open, null = "Trend")
print(kpss_test)
# After differencing, p_value is 0.1 > 0.01, which means we don't  need to difference the data again

# ACF Plot
Acf(differenced_stock_data_open, lag.max = 20, main = "ACF for Differenced Data",
    col = "blue", lwd = 2)

Pacf(differenced_stock_data_open, lag.max = 20, main = "PACF for Differenced Data",
     col = "green", lwd = 2)

# Now, let's fit various ARIMA models to our data.

model_list <- list()
max_p <- 4
max_q <- 4
for (p in 0:max_p)
{
  for (q in 0:max_q)
  {
  # model_names <- paste0("ARIMA_", p, "_1_", q)
  model_list[[paste0("ARIMA_", p, "_1_", q)]] <- Arima(differenced_stock_data_open, order = c(p,1,q))
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
model_names[index_of_least_aic]
print(model_list[model_names[index_of_least_aic]])
best_fit <- model_list[model_names[index_of_least_aic]]
residuals_best <- best_fit$ARIMA_0_1_1$residuals
checkresiduals(residuals_best)
autoplot(forecast(best_fit$ARIMA_0_1_1))