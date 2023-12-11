library(fpp3);
library(tsibble)
# print("Hello World");
# x <- 3
# plot(x-10)
# tssibles are used for time series storage, and are considered time aware.
# tibbles are data frames for general purpose data storage. data frames, are data types build for stastical analysis.
# my_tsibble <- tsibble(
#   year  = 2015:2019,
#   y = c(123,34,56,78,90),
#   # value = 1:12,
#   index = y
# )
#
#   my_data <- tibble(
# 	year = 2015:2019,
# 	y = c(212, 234, 345, 456, 567)
#   ) |>
# 	as_tsibble(index = year)
#
#   my_data
prison <- readr::read_csv("https://OTexts.com/fpp3/extrafiles/prison_population.csv") |>
  mutate(Quarter = yearquarter(date)) |>
  select(-date) |>
  as_tsibble(
	index = Quarter,
	key = c(state, gender, legal, indigenous)
  )
#