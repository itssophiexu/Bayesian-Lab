# Run from the project directory, or use Rscript R/lab_student.R.
# First complete the five expressions in R/ridge_gibbs_student.R.
arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) > 0L) {
  this_file <- normalizePath(sub("^--file=", "", file_argument[1]))
  lab_root <- dirname(dirname(this_file))
} else if (file.exists(file.path("R", "ridge_helpers.R"))) {
  lab_root <- normalizePath(".")
} else if (file.exists(file.path("..", "R", "ridge_helpers.R"))) {
  lab_root <- normalizePath("..")
} else {
  stop("Open the lab project, then source R/lab_student.R.")
}
source(file.path(lab_root, "R", "ridge_helpers.R"))
source(file.path(lab_root, "R", "ridge_gibbs_student.R"))

lab_result <- run_lab_experiment(
  fit_function = ridge_gibbs,
  output_dir = file.path(lab_root, "lab_output"),
  data = simulate_lab_data(n_train = 50, n_test = 300, p = 80,
                           rho = 0.7, signal_variance = 4, seed = 704),
  v_values = c(0.01, 0.1, 1), n_iter = 3000, burn = 1000)

print(lab_result$metrics, digits = 3, row.names = FALSE)
print(lab_result$diagnostics, digits = 3, row.names = FALSE)
cat("\nInspect the PNG plots and compare the independent chains.\n")
cat("The mean-CI response-containment column is not prediction coverage.\n")

