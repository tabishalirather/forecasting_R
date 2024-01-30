# Function to run an R script multiple times
run_r_script <- function(script_path, times) {
  for (i in 1:times) {
    print(paste("Running R script", i))
    # Use tryCatch to handle errors
    tryCatch({
      output <- source(script_path)
      print(output)
    }, error = function(e) {
      # Print the error message
      print(paste("Error during run", i, ":", e$message))
    })
  }
}

# Path to your main R script
script_path <- "C:/Users/tabis/OneDrive - Swinburne University/Summer Project 2023/TestingR/Week 1+/optimal_model_p_values.R"

# Run the main R script 10 times
run_r_script(script_path, 2)