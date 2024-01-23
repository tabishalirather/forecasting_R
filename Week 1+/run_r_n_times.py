import subprocess


def run_r_script(times):
    for _ in range(times):
        subprocess.run(["Rscript", "optimal_model_p_values.R"])


# Run the R script 5 times
run_r_script(5)
