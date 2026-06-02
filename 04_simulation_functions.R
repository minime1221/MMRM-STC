# ================================================================================
# SIMULATION EXECUTION FUNCTIONS
# ================================================================================
# Purpose: Functions to run STC simulations
#          with MMRM models, extracting treatment effects and confidence intervals
# ================================================================================
# Function: run_simulation_stc_ml
# Description: Execute a single STC simulation iteration
#              Generates AC and BC trial data, fit STC-MMRM, STC-ANCOVA
#
# Parameters:
#   seed (numeric): Random seed for this simulation iteration
#   scenario (list): Scenario list with $AC, $BC, and $truth_structure components
#                    truth_structure must be "TVEM" or "TIEM"
#
# Returns:
#   data.frame with columns for:
#     - seed: The random seed used
#     - miss_*: Missing rates for each time point in AC and BC trials
#     - trt1_*: Treatment effect estimates from different models
#     - var_trt1_*: Variance estimates for treatment effects
#     - delta_*: Difference in treatment effects (A - B)
#     - var_delta_*: Variance of delta estimates
#     - ci_*: 95% confidence interval bounds for delta

run_simulation_stc_ml <- function(seed, scenario) {

        # Extract truth structure from scenario (determines TVEM vs TIEM)
        truth_structure <- scenario$truth_structure

        # Select appropriate data generation function based on truth structure
        if (truth_structure == "TVEM") {
                gen_func <- gen_data_tvem
        } else if (truth_structure == "TIEM") {
                gen_func <- gen_data_tiem
        } else {
                stop("Unknown truth_structure. Scenario must have truth_structure='TVEM' or 'TIEM'.")
        }

        # --- STEP 1: DATA GENERATION ---

        # Generate AC trial data
        AC_data <- gen_func(
                seed = seed, trial_name = "AC", N = scenario$AC$N,
                trt_effect = scenario$AC$delta, beta_prog = scenario$AC$b_prog,
                beta_em = scenario$AC$beta_em,
                beta_em_trt = scenario$AC$beta_em_trt,
                rho = scenario$AC$rho, sd_residual = scenario$AC$sdY,
                meanX = scenario$AC$meanX, covX = scenario$AC$covX,
                dropout = scenario$AC$dropout,
                dropout_model = scenario$AC$dropout_model,
                dropout_inct = scenario$AC$dropout_inct,
                dropout_slope = scenario$AC$dropout_slope,
                scale_cov = scenario$AC$scale_cov
        )

        # Generate BC trial data
        BC_data <- gen_func(
                seed = seed + 1, trial_name = "BC", N = scenario$BC$N,
                trt_effect = scenario$BC$delta, beta_prog = scenario$BC$b_prog,
                beta_em = scenario$BC$beta_em,
                beta_em_trt = scenario$BC$beta_em_trt,
                rho = scenario$BC$rho, sd_residual = scenario$BC$sdY,
                meanX = scenario$BC$meanX, covX = scenario$BC$covX,
                dropout = scenario$BC$dropout,
                dropout_model = scenario$BC$dropout_model,
                dropout_inct = scenario$BC$dropout_inct,
                dropout_slope = scenario$BC$dropout_slope,
                scale_cov = scenario$BC$scale_cov
        )

        # Compute missing rates for each trial
        miss_ac <- compute_missing_rates_for_sim(AC_data)
        miss_bc <- compute_missing_rates_for_sim(BC_data)

        # --- STEP 2: DATA PREPARATION FOR ANALYSIS ---

        # Convert from wide to long format for mixed model analysis
        AC_ipd <- pivot_longer(AC_data, cols = starts_with("y"), names_to = "time", values_to = "y") %>%
                mutate(time = as.numeric(gsub("y", "", time)))
        BC_ipd <- pivot_longer(BC_data, cols = starts_with("y"), names_to = "time", values_to = "y") %>%
                mutate(time = as.numeric(gsub("y", "", time)))

        # Summarize baseline covariates from BC trial
        # These will be used for centering in the AC analysis
        summ_BC <- BC_ipd %>%
                summarize(
                        mean_Prog = mean(Prog, na.rm = TRUE),
                        mean_EM = mean(EM, na.rm = TRUE)
                )

        # Prepare AC data for analysis
        # Center effect modifier using BC trial mean
        # This creates the counterfactual adjustment for AC trial
        AC_ipd <- AC_ipd %>%
                mutate(
                        id = factor(id),
                        trt = factor(trt),
                        time = relevel(factor(time), ref = "4"),  # Set time 4 as reference
                        EM_cen_AC = EM - mean(EM),               # Center using AC mean
                        EM_cen_BC = EM - summ_BC$mean_EM,        # Center using BC mean (counterfactual)
                        Prog_cen_BC = Prog - summ_BC$mean_Prog   # Center prognostic factor using BC mean
                )

        # Prepare BC data for analysis
        BC_ipd <- BC_ipd %>%
                mutate(
                        id = factor(id),
                        trt = factor(trt),
                        time = relevel(factor(time), ref = "4"),
                        EM_cen_BC = EM - summ_BC$mean_EM,
                        Prog_cen_BC = Prog - summ_BC$mean_Prog
                )

        # Extract final time point data for conventional STC approach
        AC_ipd_t4 <- AC_ipd %>% filter(time == 4)

        # --- STEP 3: MODEL FITTING ---
        # NOTE: All 4 models are fit regardless of truth_structure
        # This allows assessment of model misspecification effects
        # (e.g., TVEM data fit with TIEM model, and vice versa)

        # BC trial model - TVEM structure (time-varying effect modifier)
        # This is the "external" or "comparator" trial model with EM*trt*time interaction
        mod_bc_mmrm_tvem <- mmrm(
                y ~ Prog_cen_BC + trt * time + EM_cen_BC * trt * time + us(time | id),
                reml = FALSE,
                data = BC_ipd
        )

        # BC trial model - TIEM structure (time-invariant effect modifier)
        # EM enters as main effect and interaction with trt (no time interaction)
        mod_bc_mmrm_tiem <- mmrm(
                y ~ Prog_cen_BC + trt * time + EM_cen_BC + EM_cen_BC:trt + us(time | id),
                reml = FALSE,
                data = BC_ipd
        )
        # Select the BC aggregate estimate according to the data-generating truth
        # This mimics using one reported aggregate BC estimate per scenario.
        if (scenario$truth_structure == "TVEM") {
                mod_bc_mmrm <- mod_bc_mmrm_tvem
        } else if (scenario$truth_structure == "TIEM") {
                mod_bc_mmrm <- mod_bc_mmrm_tiem
        } else {
                stop("Unknown truth_structure: ", scenario$truth_structure)
        }
        
        # AC trial - STC TVEM model: Longitudinal MMRM with time-varying EM effect
        # BC covariate centering applied (counterfactual adjustment)
        mod_stc_long_tvem <- mmrm(
                y ~ Prog_cen_BC + trt * time + EM_cen_BC * trt * time + us(time | id),
                reml = FALSE,
                data = AC_ipd
        )

        # AC trial - STC TIEM model: Longitudinal MMRM with time-invariant EM effect
        # BC covariate centering applied
        mod_stc_long_tiem <- mmrm(
                y ~ Prog_cen_BC + trt * time + EM_cen_BC + EM_cen_BC:trt + us(time | id),
                reml = FALSE,
                data = AC_ipd
        )

        # AC trial - Conventional STC model (final time point only, no time-treatment interaction)
        # Cross-sectional analysis at time 4 with baseline covariate adjustment
        mod_stc_conv <- lm(
                y ~ Prog_cen_BC + trt * EM_cen_BC,
                data = AC_ipd_t4
        )

        # --- STEP 4: EXTRACT TREATMENT EFFECT ESTIMATES ---
        # Extract treatment effect at reference level (time 4, trt=1)
        # All estimates extracted regardless of truth_structure (allows model comparison)

        # AC trial treatment effects (Longitudinal STC models)
        trt1_stc_long_ac_tvem <- coef(mod_stc_long_tvem)["trt1"]
        trt1_stc_long_ac_tiem <- coef(mod_stc_long_tiem)["trt1"]

        # AC trial treatment effect (Conventional STC - final time point)
        trt1_stc_conv_ac <- coef(mod_stc_conv)["trt1"]

        # BC trial treatment effects (both models fit)
        trt1_mmrm_bc <- coef(mod_bc_mmrm)["trt1"]

        # --- STEP 5: EXTRACT VARIANCE ESTIMATES ---
        # Variance of treatment effect estimates (diagonal of covariance matrix)
        # All variances extracted from all models (enables model misspecification assessment)

        # AC trial variances (Longitudinal STC models)
        var_trt1_stc_long_ac_tvem <- vcov(mod_stc_long_tvem)["trt1", "trt1"]
        var_trt1_stc_long_ac_tiem <- vcov(mod_stc_long_tiem)["trt1", "trt1"]

        # AC trial variance (Conventional STC)
        var_trt1_stc_conv_ac <- vcov(mod_stc_conv)["trt1", "trt1"]

        # BC trial variances (both TVEM and TIEM models)
        var_trt1_mmrm_bc <- vcov(mod_bc_mmrm)["trt1", "trt1"]

        # --- STEP 6: COMPUTE COMPARATIVE EFFECTIVENESS (DELTA) ---
        # Delta = treatment effect in AC - treatment effect in BC
        # This represents the adjusted difference in treatment effects
        # All delta values computed for all method combinations (enables model comparison)

        # STC Long TVEM: AC TVEM vs selected BC aggregate estimate
        delta_stc_long_mmrm_tvem <- trt1_stc_long_ac_tvem - trt1_mmrm_bc
        var_stc_long_mmrm_tvem <- var_trt1_stc_long_ac_tvem + var_trt1_mmrm_bc
        ci_stc_long_mmrm_tvem <- c(
                delta_stc_long_mmrm_tvem - qnorm(0.975) * sqrt(var_stc_long_mmrm_tvem),
                delta_stc_long_mmrm_tvem + qnorm(0.975) * sqrt(var_stc_long_mmrm_tvem)
        )

        # STC Long TIEM: AC TIEM vs selected BC aggregate estimate
        delta_stc_long_mmrm_tiem <- trt1_stc_long_ac_tiem - trt1_mmrm_bc
        var_stc_long_mmrm_tiem <- var_trt1_stc_long_ac_tiem + var_trt1_mmrm_bc
        ci_stc_long_mmrm_tiem <- c(
                delta_stc_long_mmrm_tiem - qnorm(0.975) * sqrt(var_stc_long_mmrm_tiem),
                delta_stc_long_mmrm_tiem + qnorm(0.975) * sqrt(var_stc_long_mmrm_tiem)
        )

        # STC Conventional (ANCOVA): AC final timepoint vs selected BC aggregate estimate
        delta_stc_conv_mmrm <- trt1_stc_conv_ac - trt1_mmrm_bc
        var_stc_conv_mmrm <- var_trt1_stc_conv_ac + var_trt1_mmrm_bc
        ci_stc_conv_mmrm <- c(
                delta_stc_conv_mmrm - qnorm(0.975) * sqrt(var_stc_conv_mmrm),
                delta_stc_conv_mmrm + qnorm(0.975) * sqrt(var_stc_conv_mmrm)
        )

        # --- STEP 7: RETURN RESULTS ---
        # Compile all results into a single data frame
        # All 4 methods have complete delta, variance, and CI estimates
        return(data.frame(
                # ---- Iteration identifier ----
                seed = seed,

                # ---- Missing data rates ----
                miss_ac_y1 = unname(miss_ac["miss_y1"]),
                miss_ac_y2 = unname(miss_ac["miss_y2"]),
                miss_ac_y3 = unname(miss_ac["miss_y3"]),
                miss_ac_y4 = unname(miss_ac["miss_y4"]),
                miss_bc_y1 = unname(miss_bc["miss_y1"]),
                miss_bc_y2 = unname(miss_bc["miss_y2"]),
                miss_bc_y3 = unname(miss_bc["miss_y3"]),
                miss_bc_y4 = unname(miss_bc["miss_y4"]),

                # ---- Individual treatment effect estimates (for reference) ----
                trt1_stc_long_ac_tvem = trt1_stc_long_ac_tvem,
                trt1_stc_long_ac_tiem = trt1_stc_long_ac_tiem,
                trt1_stc_conv_ac = trt1_stc_conv_ac,
                trt1_mmrm_bc = trt1_mmrm_bc,

                # ---- Variance of treatment effect estimates ----
                var_trt1_stc_long_ac_tvem = var_trt1_stc_long_ac_tvem,
                var_trt1_stc_long_ac_tiem = var_trt1_stc_long_ac_tiem,
                var_trt1_stc_conv_ac = var_trt1_stc_conv_ac,
                var_trt1_mmrm_bc = var_trt1_mmrm_bc,
                
                bc_model_used = truth_structure,

                # ---- STC Long TVEM: AC TVEM vs BC TVEM ----
                # Delta represents difference in treatment effects
                delta_stc_long_mmrm_tvem = delta_stc_long_mmrm_tvem,
                var_stc_long_mmrm_tvem = var_stc_long_mmrm_tvem,
                ci_lower_stc_long_mmrm_tvem = ci_stc_long_mmrm_tvem[1],
                ci_upper_stc_long_mmrm_tvem = ci_stc_long_mmrm_tvem[2],

                # ---- STC Long TIEM: AC TIEM vs BC TIEM ----
                delta_stc_long_mmrm_tiem = delta_stc_long_mmrm_tiem,
                var_stc_long_mmrm_tiem = var_stc_long_mmrm_tiem,
                ci_lower_stc_long_mmrm_tiem = ci_stc_long_mmrm_tiem[1],
                ci_upper_stc_long_mmrm_tiem = ci_stc_long_mmrm_tiem[2],

                # ---- STC Conventional (ANCOVA): AC final time point vs BC TVEM ----
                delta_stc_conv_mmrm = delta_stc_conv_mmrm,
                var_stc_conv_mmrm = var_stc_conv_mmrm,
                ci_lower_stc_conv_mmrm = ci_stc_conv_mmrm[1],
                ci_upper_stc_conv_mmrm = ci_stc_conv_mmrm[2]
        ))
}

# ================================================================================

# Function: run_simulation_bucher_ml
# Description: Execute a single Bucher method simulation iteration
#              Generates AC and BC trial data, fits MMRM models for BC and simple
#              treatment model for AC (without covariates), then combines using
#              Bucher's method for comparing treatments
#
# Parameters:
#   seed (numeric): Random seed for this simulation iteration
#   scenario (list): Scenario list with $AC, $BC, and $truth_structure components
#
# Returns:
#   data.frame with treatment effect estimates and confidence intervals

run_simulation_bucher_ml <- function(seed, scenario) {

        # Extract truth structure from scenario (determines TVEM vs TIEM)
        truth_structure <- scenario$truth_structure

        # Select appropriate data generation function based on truth structure
        if (truth_structure == "TVEM") {
                gen_func <- gen_data_tvem
        } else if (truth_structure == "TIEM") {
                gen_func <- gen_data_tiem
        } else {
                stop("Unknown truth_structure. Scenario must have truth_structure='TVEM' or 'TIEM'.")
        }

        # --- STEP 1: DATA GENERATION ---

        # Generate AC trial data
        AC_data <- gen_func(
                seed = seed, trial_name = "AC", N = scenario$AC$N,
                trt_effect = scenario$AC$delta, beta_prog = scenario$AC$b_prog,
                beta_em = scenario$AC$beta_em,
                beta_em_trt = scenario$AC$beta_em_trt,
                rho = scenario$AC$rho, sd_residual = scenario$AC$sdY,
                meanX = scenario$AC$meanX, covX = scenario$AC$covX,
                dropout = scenario$AC$dropout,
                dropout_model = scenario$AC$dropout_model,
                dropout_inct = scenario$AC$dropout_inct,
                dropout_slope = scenario$AC$dropout_slope,
                scale_cov = scenario$AC$scale_cov
        )

        # Generate BC trial data
        BC_data <- gen_func(
                seed = seed + 1, trial_name = "BC", N = scenario$BC$N,
                trt_effect = scenario$BC$delta, beta_prog = scenario$BC$b_prog,
                beta_em = scenario$BC$beta_em,
                beta_em_trt = scenario$BC$beta_em_trt,
                rho = scenario$BC$rho, sd_residual = scenario$BC$sdY,
                meanX = scenario$BC$meanX, covX = scenario$BC$covX,
                dropout = scenario$BC$dropout,
                dropout_model = scenario$BC$dropout_model,
                dropout_inct = scenario$BC$dropout_inct,
                dropout_slope = scenario$BC$dropout_slope,
                scale_cov = scenario$BC$scale_cov
        )

        # --- STEP 2: DATA PREPARATION FOR ANALYSIS ---

        # Convert to long format
        AC_ipd <- pivot_longer(AC_data, cols = starts_with("y"), names_to = "time", values_to = "y") %>%
                mutate(time = as.numeric(gsub("y", "", time)))
        BC_ipd <- pivot_longer(BC_data, cols = starts_with("y"), names_to = "time", values_to = "y") %>%
                mutate(time = as.numeric(gsub("y", "", time)))

        # Summarize covariates from BOTH AC and BC trials
        # For Bucher method, each trial is analyzed on its own scale
        summ_AC <- AC_ipd %>%
                summarize(
                        mean_Prog = mean(Prog, na.rm = TRUE),
                        mean_EM = mean(EM, na.rm = TRUE)
                )
        summ_BC <- BC_ipd %>%
                summarize(
                        mean_Prog = mean(Prog, na.rm = TRUE),
                        mean_EM = mean(EM, na.rm = TRUE)
                )

        # Prepare AC data: center by AC trial means
        AC_ipd <- AC_ipd %>%
                mutate(
                        id = factor(id),
                        trt = factor(trt),
                        time = relevel(factor(time), ref = "4"),
                        EM_cen_AC = EM - summ_AC$mean_EM,
                        Prog_cen_AC = Prog - summ_AC$mean_Prog
                )

        # Prepare BC data: center by BC trial means
        BC_ipd <- BC_ipd %>%
                mutate(
                        id = factor(id),
                        trt = factor(trt),
                        time = relevel(factor(time), ref = "4"),
                        EM_cen_BC = EM - summ_BC$mean_EM,
                        Prog_cen_BC = Prog - summ_BC$mean_Prog
                )

        # Extract final time point for simpler analysis
        AC_ipd_t4 <- AC_ipd %>% filter(time == 4)

        # --- STEP 3: MODEL FITTING ---
        # Both TVEM and TIEM BC models are fit.
        # One BC aggregate estimate is selected according to the
        # data-generating truth and used in the Bucher comparison.
        mod_bc_mmrm_tvem <- mmrm(
                y ~ Prog_cen_BC + trt * time + EM_cen_BC * trt * time + us(time | id),
                reml = FALSE,
                data = BC_ipd
        )

        # BC trial model - TIEM structure: MMRM with time-invariant effect modifier
        mod_bc_mmrm_tiem <- mmrm(
                y ~ Prog_cen_BC + trt * time + EM_cen_BC + EM_cen_BC:trt + us(time | id),
                reml = FALSE,
                data = BC_ipd
        )
        
        # Select the BC aggregate estimate according to the data-generating truth
        # This mimics using one reported aggregate BC estimate per scenario.
        if (scenario$truth_structure == "TVEM") {
                mod_bc_mmrm <- mod_bc_mmrm_tvem
        } else if (scenario$truth_structure == "TIEM") {
                mod_bc_mmrm <- mod_bc_mmrm_tiem
        } else {
                stop("Unknown truth_structure: ", scenario$truth_structure)
        }
        
        # AC trial model: Simple treatment-only model (Bucher method)
        # This represents analysis with minimal covariate adjustment
        # (Only considers treatment effect, ignoring baseline covariates)
        mod_bucher_conv <- lm(y ~ trt, data = AC_ipd_t4)

        # --- STEP 4: EXTRACT TREATMENT EFFECTS ---
        # Treatment effects from all models

        trt1_bucher_conv_ac <- coef(mod_bucher_conv)["trt1"]
        trt1_mmrm_bc <- coef(mod_bc_mmrm)["trt1"]

        # --- STEP 5: EXTRACT VARIANCES ---
        # Variance estimates from all models

        var_trt1_bucher_conv_ac <- vcov(mod_bucher_conv)["trt1", "trt1"]
        var_trt1_mmrm_bc <- vcov(mod_bc_mmrm)["trt1", "trt1"]

        # --- STEP 6: COMPUTE COMPARATIVE EFFECTIVENESS USING BUCHER METHOD ---
        # Bucher method using the selected BC aggregate estimate
        delta_bucher_conv_mmrm <- trt1_bucher_conv_ac - trt1_mmrm_bc
        var_bucher_conv_mmrm <- var_trt1_bucher_conv_ac + var_trt1_mmrm_bc
        ci_bucher_conv_mmrm <- c(
                delta_bucher_conv_mmrm - qnorm(0.975) * sqrt(var_bucher_conv_mmrm),
                delta_bucher_conv_mmrm + qnorm(0.975) * sqrt(var_bucher_conv_mmrm)
        )

        # --- STEP 7: RETURN RESULTS ---
        # Return Bucher results using the selected BC aggregate estimate
        return(data.frame(
                # Iteration identifier
                seed = seed,

                # Treatment effect estimates
                trt1_bucher_conv_ac = trt1_bucher_conv_ac,
                trt1_mmrm_bc = trt1_mmrm_bc,
                
                # Variances
                var_trt1_bucher_conv_ac = var_trt1_bucher_conv_ac,
                var_trt1_mmrm_bc = var_trt1_mmrm_bc,
                
                bc_model_used = truth_structure,

                # Bucher vs BC TVEM
                delta_bucher_conv_mmrm = delta_bucher_conv_mmrm,
                var_bucher_conv_mmrm = var_bucher_conv_mmrm,
                ci_lower_bucher_conv_mmrm = ci_bucher_conv_mmrm[1],
                ci_upper_bucher_conv_mmrm = ci_bucher_conv_mmrm[2]
        ))
}

# End of simulation functions module
