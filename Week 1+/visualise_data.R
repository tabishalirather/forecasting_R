library(dplyr)
library(ggplot2)
rm(list = ls())

data <- read.csv("Week 1+/output_5_non_zero_d.csv")
# Convert columns to numeric
data$auto_arima_true_count <- as.numeric(as.character(data$auto_arima_true_count))
data$auto_arima_false_count <- as.numeric(as.character(data$auto_arima_false_count))

data$least_p_true_count <- as.numeric(as.character(data$least_p_true_count))
data$least_p_false_count <- as.numeric(as.character(data$least_p_false_count))


data$refit_avg_rmse_least_p_model <- as.numeric(as.character(data$refit_avg_rmse_least_p_model))
data$avg_rmse_least_p_model_values <- as.numeric(as.character(data$avg_rmse_least_p_model_values))

data$refit_avg_rmse_auto_arima <- as.numeric(as.character(data$refit_avg_rmse_auto_arima))
data$avg_rmse_auto_arima_model <- as.numeric(as.character(data$avg_rmse_auto_arima_model))

data$refit_avg_insig_params_least_p_model <- as.numeric(as.character(data$refit_avg_insig_params_least_p_model))
# print(paste("refit_avg_insig_params_least_p_model", data$refit_avg_insig_params_least_p_model))
avg_refit_insig_params_least_p_model <- mean(data$refit_avg_insig_params_least_p_model, na.rm = TRUE)
# print(paste("avg_refit_insig_params_least_p_model", avg_refit_insig_params_least_p_model))
data$avg_insig_params_least_p_model <- as.numeric(as.character(data$avg_insig_params_least_p_model))

data$refit_avg_insig_params_auto_arima_model <- as.numeric(as.character(data$refit_avg_insig_params_auto_arima_model))
# print(paste("refit_avg_insig_params_auto_arima_model", data$refit_avg_insig_params_auto_arima_model))
data$avg_insig_params_auto_arima_model <- as.numeric(as.character(data$avg_insig_params_auto_arima_model))

data$avg_sw_p_value <- as.numeric(as.character(data$avg_sw_p_value))
data$avg_lbq_p_value <- as.numeric(as.character(data$avg_lbq_p_value))
data$avg_t_test_p_value <- as.numeric(as.character(data$avg_t_test_p_value))
data$avg_validation_score <- as.numeric(as.character(data$avg_validation_score))


# summary(data)
# true_false_counts <- data.frame(
#   Model = c("auto_arima", "max_p"),
#   True_Count = c(mean(data$auto_arima_true_count), mean(data$least_p_true_count)),
#   False_Count = c(mean(data$auto_arima_false_count), mean(data$least_p_false_count))
# )

# print(true_false_counts)
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

# Function for the first plot
plot_true_vs_false_counts <- function(data) {
  counts_df <- data.frame(
    Model = rep(c("auto_arima", "max_p"), each = 2),
    Count_Type = rep(c("True Count", "False Count"), 2),
    Count = c(mean(data$auto_arima_true_count, na.rm = TRUE), mean(data$auto_arima_false_count, na.rm = TRUE), mean(data$least_p_true_count, na.rm = TRUE), mean(data$least_p_false_count, na.rm = TRUE))
  )

  p <- ggplot(counts_df, aes(x = Model, y = Count, fill = Count_Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Model", y = "Count", title = "True vs False Counts for Each Model") +
    scale_fill_brewer(palette = "Set1")
  print(p)
  return(counts_df)
}

# Function for the second plot
plot_avg_rmse_comparison <- function(data) {
  rmse_df <- data.frame(
    Model = c("auto_arima", "max_p", 'refit_auto_arima', 'refit_max_p'),
    Avg_RMSE = c(mean(data$avg_rmse_auto_arima_model, na.rm = TRUE), mean(data$avg_rmse_least_p_model_values, na.rm = TRUE), mean(data$refit_avg_rmse_auto_arima, na.rm = TRUE), mean(data$refit_avg_rmse_least_p_model, na.rm = TRUE))
  )

  p <- ggplot(rmse_df, aes(x = Model, y = Avg_RMSE, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Model", y = "Average RMSE", title = "Average RMSE Comparison Between Models")
  print(p)
  return(rmse_df)
}

# Function for the third plot
plot_insig_params_comparison <- function(data) {
  insig_params_df <- data.frame(
    Model = c("auto_arima", "max_p", "refit_auto_arima", "refit_max_p"),
    Avg_Insignificant_Params = c(mean(data$avg_insig_params_auto_arima_model, na.rm = TRUE), mean(data$avg_insig_params_least_p_model, na.rm = TRUE), mean(data$refit_avg_insig_params_auto_arima_model, na.rm = TRUE), mean(data$refit_avg_insig_params_least_p_model, na.rm = TRUE))
  )

  p <- ggplot(insig_params_df, aes(x = Model, y = Avg_Insignificant_Params, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Model", y = "Average Number of Insignificant Parameters", title = "Insignificant Parameters Comparison")
  print(p)
  return(insig_params_df)
}

# Function for the fourth plot
plot_histogram_avg_rmse <- function(data) {
  p <- ggplot(data, aes(x = avg_rmse_auto_arima_model)) +
    geom_histogram(fill = "blue", color = "black", alpha = 0.4) +
    labs(x = "Average RMSE (auto_arima)", y = "Density", title = "Histogram of Average RMSE for auto_arima Model") +
    xlim(min(data$avg_rmse_auto_arima_model, na.rm = TRUE)-0.2, max(data$avg_rmse_auto_arima_model, na.rm = TRUE)+0.2) +
    theme_minimal() +
    theme(legend.position = "none")
  print(p)
}

# Function for the fifth plot
plot_histogram_avg_p_values <- function(data) {
  p <- ggplot(data, aes(x = avg_sw_p_value)) +
    geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.5) +
    labs(x = "Average SW P-Value", y = "Frequency", title = "Histogram of Average SW P-Values") +
    theme_minimal()
  print(p)
}

# Function for the sixth plot
plot_avg_validation_score <- function(data) {
  data$avg_validation_score_percent <- data$avg_validation_score * 100
  data$mean_avg_validation_score <- mean(data$avg_validation_score, na.rm = TRUE)

  p <- ggplot(data, aes(x = as.factor(row.names(data)), y = avg_validation_score_percent)) +
    geom_bar(stat = "identity", fill = "purple") +
    geom_hline(yintercept = data$mean_avg_validation_score * 100, color = "red", linetype = "dashed") +
    labs(x = "Test number", y = "Average % of validated models (%)", title = "Average number of validated models as a Percentage for Each Test") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(p)
  return(data$mean_avg_validation_score)
}

count_df <- plot_true_vs_false_counts(data)
# print(count_df)
rmse_df <- plot_avg_rmse_comparison(data)
# print(rmse_df)
insig_params_df <- plot_insig_params_comparison(data)
# print(insig_params_df)
# plot_histogram_avg_rmse(data)
# plot_histogram_avg_p_values(data)
mean_avg_val_score <- plot_avg_validation_score(data)
# print(mean_avg_val_score)