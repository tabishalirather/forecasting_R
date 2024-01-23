library(dplyr)
library(ggplot2)

data <- read.csv("Week 1+/output_2.csv")
# Convert columns to numeric
data$auto_arima_true_count <- as.numeric(as.character(data$auto_arima_true_count))
data$least_p_true_count <- as.numeric(as.character(data$least_p_true_count))

data$auto_arima_false_count <- as.numeric(as.character(data$auto_arima_false_count))
data$least_p_false_count <- as.numeric(as.character(data$least_p_false_count))

data$refit_avg_rmse_auto_arima <- as.numeric(as.character(data$refit_avg_rmse_auto_arima))
data$avg_rmse_auto_arima_model <- as.numeric(as.character(data$avg_rmse_auto_arima_model))

data$refit_avg_rmse_least_p_model <- as.numeric(as.character(data$refit_avg_rmse_least_p_model))
data$avg_rmse_least_p_model_values <- as.numeric(as.character(data$avg_rmse_least_p_model_values))

data$refit_avg_insig_params_auto_arima_model <- as.numeric(as.character(data$refit_avg_insig_params_auto_arima_model))
print(paste("refit_avg_insig_params_auto_arima_model", data$refit_avg_insig_params_auto_arima_model))
data$avg_insig_params_auto_arima_model <- as.numeric(as.character(data$avg_insig_params_auto_arima_model))

data$refit_avg_insig_params_least_p_model <- as.numeric(as.character(data$refit_avg_insig_params_least_p_model))
print(paste("refit_avg_insig_params_least_p_model", data$refit_avg_insig_params_least_p_model))
avg_refit_insig_params_least_p_model <- mean(data$refit_avg_insig_params_least_p_model, na.rm = TRUE)
print(paste("avg_refit_insig_params_least_p_model", avg_refit_insig_params_least_p_model))
data$avg_insig_params_least_p_model <- as.numeric(as.character(data$avg_insig_params_least_p_model))

# summary(data)
true_false_counts <- data.frame(
  Model = c("auto_arima", "max_p"),
  True_Count = c(mean(data$auto_arima_true_count), mean(data$least_p_true_count)),
  False_Count = c(mean(data$auto_arima_false_count), mean(data$least_p_false_count))
)

print(true_false_counts)
# Average True Counts for each model
avg_true_auto_arima <- mean(data$auto_arima_true_count, na.rm = TRUE)
avg_true_least_p <- mean(data$least_p_true_count, na.rm = TRUE)

# Average False Counts for each model
avg_false_auto_arima <- mean(data$auto_arima_false_count, na.rm = TRUE)
avg_false_least_p <- mean(data$least_p_false_count, na.rm = TRUE)

# Average RMSE for each model
avg_rmse_auto_arima <- mean(data$avg_rmse_auto_arima_model, na.rm = TRUE)
avg_rmse_least_p <- mean(data$avg_rmse_least_p_model_values, na.rm = TRUE)

refit_avg_rmse_auto_arima <- mean(data$refit_avg_rmse_auto_arima, na.rm = TRUE)
refit_avg_rmse_least_p <- mean(data$refit_avg_rmse_least_p_model, na.rm = TRUE)

# Average Number of Insignificant Parameters for each model
avg_insig_params_auto_arima <- mean(data$avg_insig_params_auto_arima_model, na.rm = TRUE)
avg_insig_params_least_p <- mean(data$avg_insig_params_least_p_model, na.rm = TRUE)

refit_avg_insig_params_auto_arima <- mean(data$refit_avg_insig_params_auto_arima_model, na.rm = TRUE)
refit_avg_insig_params_least_p <- mean(data$refit_avg_insig_params_least_p_model , na.rm = TRUE)

counts_df <- data.frame(
  Model = rep(c("auto_arima", "max_p"), each = 2),
  Count_Type = rep(c("True Count", "False Count"), 2),
  Count = c(avg_true_auto_arima, avg_false_auto_arima, avg_true_least_p, avg_false_least_p)
)

ggplot(counts_df, aes(x = Model, y = Count, fill = Count_Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Model", y = "Count", title = "True vs False Counts for Each Model") +
  scale_fill_brewer(palette = "Set1")


# Prepare data for plotting
rmse_df <- data.frame(
  Model = c("auto_arima", "max_p", 'refit_auto_arima', 'refit_max_p'),
  Avg_RMSE = c(avg_rmse_auto_arima, avg_rmse_least_p, refit_avg_rmse_auto_arima, refit_avg_rmse_least_p)
)
print(rmse_df)
# Plot
ggplot(rmse_df, aes(x = Model, y = Avg_RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Model", y = "Average RMSE", title = "Average RMSE Comparison Between Models")


# Prepare data for plotting
insig_params_df <- data.frame(
  Model = c("auto_arima", "max_p", "refit_auto_arima", "refit_max_p"),
  Avg_Insignificant_Params = c(avg_insig_params_auto_arima, avg_insig_params_least_p, refit_avg_insig_params_auto_arima, refit_avg_insig_params_least_p)
)
print(insig_params_df)

# Plot
ggplot(insig_params_df, aes(x = Model, y = Avg_Insignificant_Params, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Model", y = "Average Number of Insignificant Parameters", title = "Insignificant Parameters Comparison")


# histogram plot
# Histogram for avg_rmse_auto_arima_model
	ggplot(data, aes(x = avg_rmse_auto_arima_model)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.4) +
  labs(x = "Average RMSE (auto_arima)", y = "Density", title = "Histogram of Average RMSE for auto_arima Model") +
  xlim(min(data$avg_rmse_auto_arima_model)-0.2, max(data$avg_rmse_auto_arima_model)+0.2) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data, aes(x = refit_avg_rmse_least_p_model)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.4) +
  labs(x = "Average RMSE (refit_least_p)", y = "Density", title = "Histogram of Average RMSE for refit least p model") +
  xlim(min(data$refit_avg_rmse_least_p_model)-0.2, max(data$refit_avg_rmse_least_p_model)+0.2) +
  theme_minimal() +
  theme(legend.position = "none")
# Calculate IQR
# Q1 <- quantile(data$refit_avg_rmse_least_p_model, 0.25)
# Q3 <- quantile(data$refit_avg_rmse_least_p_model, 0.75)
# IQR <- Q3 - Q1
#
# # Define bounds
# lower_bound <- Q1 - 1.5 * IQR
# upper_bound <- Q3 + 1.5 * IQR
#
# # Extract outliers
# outliers <- data$refit_avg_rmse_least_p_model[data$refit_avg_rmse_least_p_model < lower_bound | data$refit_avg_rmse_least_p_model > upper_bound]
#
# # Print outliers
# outliers
#
# # Calculate IQR for avg_rmse_auto_arima_model
# Q1_auto_arima <- quantile(data$avg_rmse_auto_arima_model, 0.25)
# Q3_auto_arima <- quantile(data$avg_rmse_auto_arima_model, 0.75)
# IQR_auto_arima <- Q3_auto_arima - Q1_auto_arima
#
# # Define bounds for outliers
# lower_bound_auto_arima <- Q1_auto_arima - 1.5 * IQR_auto_arima
# upper_bound_auto_arima <- Q3_auto_arima + 1.5 * IQR_auto_arima
#
# # Extract outliers for avg_rmse_auto_arima_model
# outliers_auto_arima <- data$avg_rmse_auto_arima_model[data$avg_rmse_auto_arima_model < lower_bound_auto_arima | data$avg_rmse_auto_arima_model > upper_bound_auto_arima]
#
# # Print outliers for avg_rmse_auto_arima_model
# outliers_auto_arima

#
# ggplot(data) +
#   geom_histogram(aes(x = avg_rmse_auto_arima_model), binwidth = 0.05, fill = "blue", alpha = 0.5) +
#   geom_histogram(aes(x = avg_rmse_least_p_model_values), binwidth = 0.05, fill = "red", alpha = 0.5) +
#   labs(x = "Average Number of Insignificant Parameters", y = "Frequency",
#        title = "Overlay of Histograms for avg_rmse_least_p and avg_rmse_auto_arima_model") +
#   theme_minimal()




# Histogram for the average number of insignificant parameters
# ggplot(data, aes(x = avg_insig_params_auto_arima_model)) +
#   geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
#   labs(x = "Average Number of Insignificant Parameters for auto arima", y = "Frequency",
#        title = "Distribution of Average Number of Insignificant Parameters") +
#   theme_minimal()


# Histogram for the average number of insignificant parameters
# ggplot(data, aes(x = avg_insig_params_least_p_model)) +
#   geom_histogram(binwidth = 0.07, fill = "blue", color = "black") +
#   labs(x = "Average Number of Insignificant Parameters for least p", y = "Frequency",
#        title = "Distribution of Average Number of Insignificant Parameters") +
#   theme_minimal()
#
# ggplot(data) +
#   geom_histogram(aes(x = avg_insig_params_least_p_model), binwidth = 0.095, fill = "blue", alpha = 0.5) +
#   geom_histogram(aes(x = avg_insig_params_auto_arima_model), binwidth = 0.095, fill = "red", alpha = 0.5) +
#   labs(x = "Average Number of Insignificant Parameters ", y = "Frequency",
#        title = "Overlay of Histograms for avg_insig_params_least_p_model and avg_insig_params_auto_arima_model") +
#   theme_minimal()
