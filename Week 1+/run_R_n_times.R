# Function to run an R script multiple times
run_r_script <- function(script_path, times) {
  for (i in 1:times) {
    # output <- system(paste("Rscript", script_path))
    output <- source(script_path)
    print(output)
  }
}

# Path to your main R script
script_path <- "C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/TestingR/Week 1+/optimal_model_p_values.R"

# Run the main R script 5 times
run_r_script(script_path, 3)