# ================================================================================
# BATCH SIMULATION: RUN ALL SCENARIOS AUTOMATICALLY
# ================================================================================
# Purpose: Execute simulations for all 24 scenarios (12 TVEM + 12 TIEM) in sequence
#          with configurable number of iterations for each scenario
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

# --- STEP 2: CONFIGURE BATCH PARAMETERS ---

# IMPORTANT: Set these parameters before running
num_simulations <- 3000            # Number of iterations per scenario
                                   # For testing: 50 (5-10 min per scenario)
                                   # For final: 3000 (30-60 min per scenario)

method <- "both"                   # "STC", "Bucher", or "both"

# Scenarios to run
# Options:
#   - "all": Run all 24 scenarios (TVEM and TIEM)
#   - "tvem_only": Run only 12 TVEM scenarios
#   - "tiem_only": Run only 12 TIEM scenarios
scenarios_to_run <- "all"

# Create output directory if it doesn't exist
if (!dir.exists("results")) {
        dir.create("results")
}

# Generate timestamp for output files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- STEP 3: DETERMINE WHICH SCENARIOS TO RUN ---

cat("\n")
cat("================================================================================\n")
cat("BATCH SIMULATION: ALL SCENARIOS\n")
cat("================================================================================\n")
cat(sprintf("Configuration:\n"))
cat(sprintf("  Number of simulations per scenario: %d\n", num_simulations))
cat(sprintf("  Method: %s\n", method))
cat(sprintf("  Scenarios to run: %s\n", scenarios_to_run))
cat("================================================================================\n\n")

# Select scenarios based on filter
scenario_names <- names(scenarios_list)

if (scenarios_to_run == "all") {
        selected_scenarios <- scenario_names
} else if (scenarios_to_run == "tvem_only") {
        selected_scenarios <- grep("_TVEM$", scenario_names, value = TRUE)
} else if (scenarios_to_run == "tiem_only") {
        selected_scenarios <- grep("_TIEM$", scenario_names, value = TRUE)
} else {
        stop("Invalid scenarios_to_run. Use 'all', 'tvem_only', or 'tiem_only'.")
}

cat(sprintf("Will run %d scenarios:\n", length(selected_scenarios)))
for (i in seq_along(selected_scenarios)) {
        cat(sprintf("  %2d. %s\n", i, selected_scenarios[i]))
}
cat("\n")

# --- STEP 4: RUN BATCH SIMULATIONS ---

# Initialize results storage
all_results_stc <- NULL
all_results_bucher <- NULL

# Overall progress tracking
scenario_count <- length(selected_scenarios)
overall_start_time <- Sys.time()

for (scenario_idx in seq_along(selected_scenarios)) {

        scenario_name <- selected_scenarios[scenario_idx]
        scenario <- scenarios_list[[scenario_name]]
        truth_structure <- scenario$truth_structure

        cat("\n")
        cat("================================================================================\n")
        cat(sprintf("SCENARIO %d/%d: %s (%s)\n",
                scenario_idx, scenario_count, scenario_name, truth_structure))
        cat("================================================================================\n")

        # Create scenario-specific progress bar
        pb <- progress::progress_bar$new(
                format = "  [:bar] :current/:total (:percent) [Elapsed: :elapsedsecs ETA: :etas]",
                total = num_simulations,
                clear = FALSE
        )

        # Initialize scenario-specific results
        results_stc_scenario <- NULL
        results_bucher_scenario <- NULL

        # Run simulations for this scenario
        start_time_scenario <- Sys.time()

        for (i in 1:num_simulations) {

                # Set seed for reproducibility
                # Format: scenario_index (3 digits) + iteration (5 digits)
                seed <- (scenario_idx * 100000) + i

                # Run STC simulation if requested
                if (method %in% c("STC", "both")) {
                        result_stc <- run_simulation_stc_ml(
                                seed = seed,
                                scenario = scenario
                        )
                        # Add scenario identifier
                        result_stc$scenario <- scenario_name
                        result_stc$truth_structure <- truth_structure
                        results_stc_scenario <- rbind(results_stc_scenario, result_stc)
                }

                # Run Bucher simulation if requested
                if (method %in% c("Bucher", "both")) {
                        result_bucher <- run_simulation_bucher_ml(
                                seed = seed,
                                scenario = scenario
                        )
                        # Add scenario identifier
                        result_bucher$scenario <- scenario_name
                        result_bucher$truth_structure <- truth_structure
                        results_bucher_scenario <- rbind(results_bucher_scenario, result_bucher)
                }

                # Update progress bar
                pb$tick()
        }

        end_time_scenario <- Sys.time()
        elapsed_scenario <- difftime(end_time_scenario, start_time_scenario, units = "secs")

        cat(sprintf("Completed in %.1f seconds\n\n", as.numeric(elapsed_scenario)))

        # Save scenario results
        if (!is.null(results_stc_scenario)) {
                output_file_stc <- sprintf("results/stc_results_%s_%s.rds",
                        scenario_name, timestamp)
                saveRDS(results_stc_scenario, output_file_stc, compress = "xz")
                cat(sprintf("STC results saved: %s\n", basename(output_file_stc)))
                cat(sprintf("  File size: %.2f MB, Rows: %d\n",
                        file.size(output_file_stc) / 1024^2, nrow(results_stc_scenario)))

                # Append to overall results
                all_results_stc <- rbind(all_results_stc, results_stc_scenario)
        }

        if (!is.null(results_bucher_scenario)) {
                output_file_bucher <- sprintf("results/bucher_results_%s_%s.rds",
                        scenario_name, timestamp)
                saveRDS(results_bucher_scenario, output_file_bucher, compress = "xz")
                cat(sprintf("Bucher results saved: %s\n", basename(output_file_bucher)))
                cat(sprintf("  File size: %.2f MB, Rows: %d\n",
                        file.size(output_file_bucher) / 1024^2, nrow(results_bucher_scenario)))

                # Append to overall results
                all_results_bucher <- rbind(all_results_bucher, results_bucher_scenario)
        }

        cat("\n")
}

overall_end_time <- Sys.time()
overall_elapsed <- difftime(overall_end_time, overall_start_time, units = "secs")

# --- STEP 5: SUMMARY ---

cat("\n")
cat("================================================================================\n")
cat("BATCH SIMULATION COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Total scenarios run: %d\n", scenario_count))
cat(sprintf("Total time elapsed: %.1f seconds (%.1f minutes)\n",
        as.numeric(overall_elapsed), as.numeric(overall_elapsed) / 60))
cat(sprintf("Average time per scenario: %.1f seconds\n",
        as.numeric(overall_elapsed) / scenario_count))

if (!is.null(all_results_stc)) {
        cat(sprintf("\nSTC Results Summary:\n"))
        cat(sprintf("  Total rows: %d\n", nrow(all_results_stc)))
        cat(sprintf("  Scenarios: %s\n", paste(unique(all_results_stc$scenario), collapse = ", ")))
}

if (!is.null(all_results_bucher)) {
        cat(sprintf("\nBucher Results Summary:\n"))
        cat(sprintf("  Total rows: %d\n", nrow(all_results_bucher)))
        cat(sprintf("  Scenarios: %s\n", paste(unique(all_results_bucher$scenario), collapse = ", ")))
}

cat("\n")
cat(sprintf("Timestamp: %s\n", timestamp))
cat("================================================================================\n\n")

# End of batch simulation script
