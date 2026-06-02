# ================================================================================
# MAIN SIMULATION EXECUTION SCRIPT
# ================================================================================
# Purpose: Orchestrate the entire simulation study workflow
#          Load modules, run simulations, save results
# ================================================================================

# Clear workspace
rm(list = ls())

# --- STEP 1: LOAD ALL MODULES IN CORRECT ORDER ---

# Module 1: Setup and library loading
source("R/01_setup.R")

# Module 2: Covariance utility functions
source("R/02_covariance_utilities.R")

# Module 3: Data generation functions
source("R/03_data_generation.R")

# Module 4: Simulation execution functions
source("R/04_simulation_functions.R")

# Module 5: Scenario definitions
source("R/05_scenario_definitions.R")

# --- STEP 2: CONFIGURE SIMULATION PARAMETERS ---

# Choose scenario to run
scenario_name <- "Scenario4_TVEM"  # Change this to run different scenarios
scenario <- scenarios_list[[scenario_name]]

# Extract truth structure from scenario (no need to specify separately)
truth_structure <- scenario$truth_structure

# Number of simulation iterations (replicates)
num_simulations <- 50              # For testing; use 3000+ for final results

# Method to run: "STC", "Bucher", or "both"
method <- "both"

# Create output directory if it doesn't exist
if (!dir.exists("results")) {
        dir.create("results")
}

# Generate timestamp for output files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- STEP 3: RUN SIMULATIONS ---

cat("\n")
cat("================================================================================\n")
cat("STARTING SIMULATION\n")
cat("================================================================================\n")
cat(sprintf("Scenario: %s\n", scenario_name))
cat(sprintf("Truth structure: %s\n", truth_structure))
cat(sprintf("Number of iterations: %d\n", num_simulations))
cat(sprintf("Method: %s\n", method))
cat("================================================================================\n\n")

# Initialize result data frames
results_stc <- NULL
results_bucher <- NULL

# Set up progress tracking
pb_total <- progress::progress_bar$new(
        format = "[:bar] :current/:total (:percent) [Elapsed: :elapsedsecs ETA: :etas]",
        total = num_simulations,
        clear = FALSE
)

# Run simulations
start_time <- Sys.time()

for (i in 1:num_simulations) {
        # Set seed for reproducibility
        # Each iteration gets a unique seed based on iteration number
        seed <- 1000 + i

        # Run STC simulation if requested
        if (method %in% c("STC", "both")) {
                result_stc <- run_simulation_stc_ml(
                        seed = seed,
                        scenario = scenario
                )
                results_stc <- rbind(results_stc, result_stc)
        }

        # Run Bucher simulation if requested
        if (method %in% c("Bucher", "both")) {
                result_bucher <- run_simulation_bucher_ml(
                        seed = seed,
                        scenario = scenario
                )
                results_bucher <- rbind(results_bucher, result_bucher)
        }

        # Update progress bar
        pb_total$tick()
}

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")

cat("\n")
cat(sprintf("Simulations completed in %.2f seconds\n", as.numeric(elapsed_time)))
cat("\n")

# --- STEP 4: SAVE RESULTS ---

# Save STC results in RDS format (R native, compressed)
if (!is.null(results_stc)) {
        output_file_stc <- sprintf("results/stc_results_%s_%s.rds",
                scenario_name, timestamp)
        saveRDS(results_stc, output_file_stc, compress = "xz")
        cat(sprintf("STC results saved to: %s\n", output_file_stc))
        cat(sprintf("  File size: %.2f MB\n", file.size(output_file_stc) / 1024^2))
}

# Save Bucher results in RDS format
if (!is.null(results_bucher)) {
        output_file_bucher <- sprintf("results/bucher_results_%s_%s.rds",
                scenario_name, timestamp)
        saveRDS(results_bucher, output_file_bucher, compress = "xz")
        cat(sprintf("Bucher results saved to: %s\n", output_file_bucher))
        cat(sprintf("  File size: %.2f MB\n", file.size(output_file_bucher) / 1024^2))
}

# --- STEP 5: COMPUTE AND DISPLAY SUMMARY STATISTICS ---

cat("\n")
cat("================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

if (!is.null(results_stc)) {
        cat("--- STC RESULTS ---\n\n")

        # Mean and SD of treatment effects
        if (!all(is.na(results_stc$delta_stc_long_mmrm_tvem))) {
                cat("Delta STC Long TVEM:\n")
                cat(sprintf("  Mean: %.4f\n", mean(results_stc$delta_stc_long_mmrm_tvem, na.rm = TRUE)))
                cat(sprintf("  SD:   %.4f\n", sd(results_stc$delta_stc_long_mmrm_tvem, na.rm = TRUE)))
                cat(sprintf("  Bias: %.4f\n", mean(results_stc$delta_stc_long_mmrm_tvem, na.rm = TRUE) - 0))
                cat("\n")
        }

        if (!all(is.na(results_stc$delta_stc_long_mmrm_tiem))) {
                cat("Delta STC Long TIEM:\n")
                cat(sprintf("  Mean: %.4f\n", mean(results_stc$delta_stc_long_mmrm_tiem, na.rm = TRUE)))
                cat(sprintf("  SD:   %.4f\n", sd(results_stc$delta_stc_long_mmrm_tiem, na.rm = TRUE)))
                cat(sprintf("  Bias: %.4f\n", mean(results_stc$delta_stc_long_mmrm_tiem, na.rm = TRUE) - 0))
                cat("\n")
        }

        if (!all(is.na(results_stc$delta_stc_conv_mmrm))) {
                cat("Delta STC Conv:\n")
                cat(sprintf("  Mean: %.4f\n", mean(results_stc$delta_stc_conv_mmrm, na.rm = TRUE)))
                cat(sprintf("  SD:   %.4f\n", sd(results_stc$delta_stc_conv_mmrm, na.rm = TRUE)))
                cat(sprintf("  Bias: %.4f\n", mean(results_stc$delta_stc_conv_mmrm, na.rm = TRUE) - 0))
                cat("\n")
        }

        # Coverage probability (95% CI includes true effect = 0)
        cat("Coverage Probability (95% CI includes 0):\n")
        if (!all(is.na(results_stc$ci_lower_stc_long_mmrm_tvem))) {
                coverage_tvem <- mean((results_stc$ci_lower_stc_long_mmrm_tvem <= 0) &
                        (results_stc$ci_upper_stc_long_mmrm_tvem >= 0), na.rm = TRUE)
                cat(sprintf("  STC Long TVEM: %.1f%%\n", coverage_tvem * 100))
        }
        if (!all(is.na(results_stc$ci_lower_stc_long_mmrm_tiem))) {
                coverage_tiem <- mean((results_stc$ci_lower_stc_long_mmrm_tiem <= 0) &
                        (results_stc$ci_upper_stc_long_mmrm_tiem >= 0), na.rm = TRUE)
                cat(sprintf("  STC Long TIEM: %.1f%%\n", coverage_tiem * 100))
        }
        if (!all(is.na(results_stc$ci_lower_stc_conv_mmrm))) {
                coverage_conv <- mean((results_stc$ci_lower_stc_conv_mmrm <= 0) &
                        (results_stc$ci_upper_stc_conv_mmrm >= 0), na.rm = TRUE)
                cat(sprintf("  STC Conv:      %.1f%%\n", coverage_conv * 100))
        }
        cat("\n")
}

if (!is.null(results_bucher)) {
        cat("--- BUCHER RESULTS ---\n\n")

        cat("Delta Bucher Conv:\n")
        cat(sprintf("  Mean: %.4f\n", mean(results_bucher$delta_bucher_conv_mmrm, na.rm = TRUE)))
        cat(sprintf("  SD:   %.4f\n", sd(results_bucher$delta_bucher_conv_mmrm, na.rm = TRUE)))
        cat(sprintf("  Bias: %.4f\n", mean(results_bucher$delta_bucher_conv_mmrm, na.rm = TRUE) - 0))
        cat("\n")

        coverage_bucher <- mean((results_bucher$ci_lower_bucher_conv_mmrm <= 0) &
                (results_bucher$ci_upper_bucher_conv_mmrm >= 0), na.rm = TRUE)
        cat(sprintf("Coverage Probability (95%% CI includes 0): %.1f%%\n", coverage_bucher * 100))
        cat("\n")
}

cat("================================================================================\n")
cat("SIMULATION COMPLETE\n")
cat("================================================================================\n")

# --- STEP 6: STORE RESULTS IN GLOBAL ENVIRONMENT (Optional) ---

# This allows further exploration in interactive mode
assign("results_stc", results_stc, envir = .GlobalEnv)
assign("results_bucher", results_bucher, envir = .GlobalEnv)

cat("\nResults stored in 'results_stc' and 'results_bucher' objects.\n")
cat("You can now explore results, create plots, and perform additional analyses.\n")

# End of main simulation script
