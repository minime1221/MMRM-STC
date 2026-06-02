# ================================================================================
# DATA GENERATION FUNCTIONS
# ================================================================================
# Purpose: Functions for generating synthetic trial data with various configurations
#          Supports both time-varying effect modifier (TVEM) and time-invariant (TIEM)
#          effect modifier structures, with configurable dropout mechanisms
# ================================================================================

# Function: gen_data_tvem
# Description: Generate longitudinal trial data with Time-Varying Effect Modifier (TVEM)
#              The effect modifier's interaction with treatment varies across time points
#
# Parameters:
#   seed (numeric): Random seed for reproducibility
#   trial_name (character): Name of the trial (e.g., "AC", "BC")
#   N (numeric): Sample size
#   allocation (numeric): Treatment allocation ratio (default 0.5 for balanced)
#   time_points (numeric): Number of follow-up time points (default 4)
#   trt_effect (numeric): Vector of treatment effects at each time point
#   placebo_effect (numeric): Vector of placebo effects at each time point
#   beta_prog (numeric): Prognostic factor coefficient
#   beta_em (numeric): Vector of effect modifier main effect coefficients (time-varying)
#   beta_em_trt (numeric): Vector of effect modifier × treatment interaction coefficients
#   rho (numeric): Vector of 6 pairwise correlations
#   sd_residual (numeric): Vector of residual standard deviations per time point
#   meanX (numeric): Vector of mean baseline covariates (length 2)
#   covX (matrix): 2x2 covariance matrix for baseline covariates
#   dropout (character): Type of dropout: "MAR", "NOdrop"
#   dropout_model (character): Dropout mechanism: "logistic" or "threshold"
#   dropout_inct (numeric): Vector of dropout intercepts for each time point
#   dropout_slope (numeric): Slope parameter for logistic dropout
#   dropout_threshold (numeric): Threshold value for threshold-based dropout
#   scale_cov (numeric): Scale factor for covariance matrix

gen_data_tvem <- function(seed, trial_name, N, allocation = 0.5, time_points = 4,
                          trt_effect, placebo_effect = c(0, 0, 0, 0), beta_prog,
                          beta_em, beta_em_trt, rho, sd_residual,
                          meanX, covX, dropout, dropout_model = "logistic",
                          dropout_inct = c(-2, -2, -2),  # t1, t2, t3
                          dropout_slope = -1.5,
                          dropout_threshold = -1.5, scale_cov = 1) {

        # Set random seed for reproducibility
        set.seed(seed)

        # --- STEP 1: ASSIGN TREATMENT AND GENERATE COVARIATES ---

        # Create treatment assignment (1 = treatment, 0 = placebo)
        trt <- c(rep(1, N * allocation), rep(0, N * (1 - allocation)))
        trt <- sample(trt)  # Randomize the order

        # Generate covariance matrix for residuals with scaling factor
        Sigma <- gen_covmat_us(sd = sd_residual, rho = rho) * scale_cov

        # Generate correlated baseline covariates (Prog and EM)
        # Prog = Prognostic factor, EM = Effect modifier
        baseline_covariates <- MASS::mvrnorm(N, mu = meanX, Sigma = covX)
        Prog <- baseline_covariates[, 1]
        EM <- baseline_covariates[, 2]

        # Generate residuals with the specified covariance structure
        # Each subject has correlated residuals across time points
        residuals <- MASS::mvrnorm(N, mu = rep(0, time_points), Sigma = Sigma)

        # --- STEP 2: GENERATE OUTCOME VALUES ---

        # Combine treatment effect and placebo effect for net effect
        trt_effect_net <- trt_effect + placebo_effect
        trt_effect <- trt_effect_net

        # Initialize outcome matrix
        y <- matrix(NA, nrow = N, ncol = time_points)

        # Generate outcome for each time point
        # Model: y[i,j] = Prog*beta_prog + trt[i]*trt_effect[j] +
        #                  (1-trt[i])*placebo_effect[j] +
        #                  EM[i]*beta_em[j] + trt[i]*EM[i]*beta_em_trt[j] + residuals[i,j]
        for (j in 1:time_points) {
                y[, j] <- 0 +
                        Prog * beta_prog +                # Prognostic effect (constant over time)
                        trt * trt_effect[j] +             # Treatment effect (varies by time)
                        (1 - trt) * placebo_effect[j] +   # Placebo effect (varies by time)
                        EM * beta_em[j] +                 # Effect modifier main effect (time-varying)
                        trt * EM * beta_em_trt[j] +       # EM × treatment interaction (time-varying)
                        residuals[, j]                    # Residual error
        }

        # --- STEP 3: IMPLEMENT DROPOUT MECHANISM ---

        # Logit (inverse logit) function for converting linear predictor to probability
        expit <- function(x) {
                exp(x) / (1 + exp(x))
        }

        # Create response indicator matrix (1 = observed, 0 = missing due to dropout)
        r <- matrix(1, nrow = N, ncol = time_points)

        if (dropout == "MAR") {
                # Missing At Random mechanism

                if (dropout_model == "logistic") {
                        # Logistic dropout: probability depends on previous outcome value
                        for (j in 2:time_points) {
                                # Calculate dropout probability based on previous outcome
                                prob_dropout <- expit(dropout_inct[j - 1] + dropout_slope * y[, j - 1])
                                # Once a subject drops out, all subsequent values are missing
                                r[, j] <- r[, j - 1] * rbinom(N, size = 1, prob = 1 - prob_dropout)
                        }
                } else if (dropout_model == "threshold") {
                        # Threshold-based dropout: dropout if outcome falls below threshold
                        for (i in 1:N) {
                                for (j in 1:time_points) {
                                        if (y[i, j] < dropout_threshold) {
                                                # Mark all subsequent time points as missing
                                                if (j < time_points) {
                                                        r[i, (j + 1):time_points] <- 0
                                                }
                                                break
                                        }
                                }
                        }
                } else {
                        stop("Unknown dropout_model. Use 'logistic' or 'threshold'.")
                }

                # Set missing outcome values to NA
                y[r == 0] <- NA

        } else if (dropout == "NOdrop") {
                # No dropout: all subjects complete the study
                # r matrix remains all 1's
                r <- matrix(1, nrow = N, ncol = time_points)
                # No values are set to NA

        } else {
                stop("Unknown dropout type. Use 'MAR' or 'NOdrop'.")
        }

        # --- STEP 4: FORMAT OUTPUT DATA ---

        # Create wide format data frame for output
        wideData <- data.frame(
                id = 1:N,           # Subject identifier
                trial = trial_name, # Trial name (AC, BC, etc.)
                trt = trt,          # Treatment assignment
                Prog = Prog,        # Prognostic factor value
                EM = EM             # Effect modifier value
        )

        # Add outcome columns (y1, y2, y3, y4)
        for (j in 1:time_points) {
                wideData[[paste0("y", j)]] <- y[, j]
        }

        return(wideData)
}

# ================================================================================

# Function: gen_data_tiem
# Description: Generate longitudinal trial data with Time-Invariant Effect Modifier (TIEM)
#              The effect modifier's main effect and treatment interaction are constant
#              across time points (unlike TVEM where they vary by time)
#
# Parameters: (identical to gen_data_tvem except)
#   beta_em (numeric): Single value for effect modifier main effect (constant over time)
#   beta_em_trt (numeric): Single value for EM × treatment interaction (constant)

gen_data_tiem <- function(seed, trial_name, N, allocation = 0.5, time_points = 4,
                          trt_effect, placebo_effect = c(0, 0, 0, 0), beta_prog,
                          beta_em, beta_em_trt, rho, sd_residual,
                          meanX, covX, dropout, dropout_model = "logistic",
                          dropout_inct = c(-2, -2, -2),  # t1, t2, t3
                          dropout_slope = -1.5,
                          dropout_threshold = -1.5, scale_cov = 1) {

        # Set random seed for reproducibility
        set.seed(seed)

        # --- STEP 1: ASSIGN TREATMENT AND GENERATE COVARIATES ---

        # Create treatment assignment (1 = treatment, 0 = placebo)
        trt <- c(rep(1, N * allocation), rep(0, N * (1 - allocation)))
        trt <- sample(trt)  # Randomize the order

        # Generate covariance matrix for residuals with scaling factor
        Sigma <- gen_covmat_us(sd = sd_residual, rho = rho) * scale_cov

        # Generate correlated baseline covariates (Prog and EM)
        baseline_covariates <- MASS::mvrnorm(N, mu = meanX, Sigma = covX)
        Prog <- baseline_covariates[, 1]
        EM <- baseline_covariates[, 2]

        # Generate residuals with the specified covariance structure
        residuals <- MASS::mvrnorm(N, mu = rep(0, time_points), Sigma = Sigma)

        # --- STEP 2: GENERATE OUTCOME VALUES ---

        # Combine treatment effect and placebo effect for net effect
        trt_effect_net <- trt_effect + placebo_effect
        trt_effect <- trt_effect_net

        # Initialize outcome matrix
        y <- matrix(NA, nrow = N, ncol = time_points)

        # Generate outcome for each time point
        # Key difference from TVEM: beta_em and beta_em_trt are scalars (not vectors)
        # Model: y[i,j] = Prog*beta_prog + trt[i]*trt_effect[j] +
        #                  (1-trt[i])*placebo_effect[j] +
        #                  EM[i]*beta_em + trt[i]*EM[i]*beta_em_trt + residuals[i,j]
        for (j in 1:time_points) {
                y[, j] <- 0 +
                        Prog * beta_prog +                # Prognostic effect
                        trt * trt_effect[j] +             # Treatment effect (time-varying)
                        (1 - trt) * placebo_effect[j] +   # Placebo effect (time-varying)
                        EM * beta_em +                    # Effect modifier main effect (time-INVARIANT)
                        trt * EM * beta_em_trt +          # EM × treatment interaction (time-INVARIANT)
                        residuals[, j]                    # Residual error
        }

        # --- STEP 3: IMPLEMENT DROPOUT MECHANISM ---

        # Logit function for converting linear predictor to probability
        expit <- function(x) {
                exp(x) / (1 + exp(x))
        }

        # Create response indicator matrix
        r <- matrix(1, nrow = N, ncol = time_points)

        if (dropout == "MAR") {
                # Missing At Random mechanism

                if (dropout_model == "logistic") {
                        # Logistic dropout based on previous outcome
                        for (j in 2:time_points) {
                                prob_dropout <- expit(dropout_inct[j - 1] + dropout_slope * y[, j - 1])
                                r[, j] <- r[, j - 1] * rbinom(N, size = 1, prob = 1 - prob_dropout)
                        }
                } else if (dropout_model == "threshold") {
                        # Threshold-based dropout
                        for (i in 1:N) {
                                for (j in 1:time_points) {
                                        if (y[i, j] < dropout_threshold) {
                                                if (j < time_points) {
                                                        r[i, (j + 1):time_points] <- 0
                                                }
                                                break
                                        }
                                }
                        }
                } else {
                        stop("Unknown dropout_model. Use 'logistic' or 'threshold'.")
                }

                # Set missing outcome values to NA
                y[r == 0] <- NA

        } else if (dropout == "NOdrop") {
                # No dropout
                r <- matrix(1, nrow = N, ncol = time_points)

        } else {
                stop("Unknown dropout type. Use 'MAR' or 'NOdrop'.")
        }

        # --- STEP 4: FORMAT OUTPUT DATA ---

        # Create wide format data frame for output
        wideData <- data.frame(
                id = 1:N,           # Subject identifier
                trial = trial_name, # Trial name
                trt = trt,          # Treatment assignment
                Prog = Prog,        # Prognostic factor
                EM = EM             # Effect modifier
        )

        # Add outcome columns
        for (j in 1:time_points) {
                wideData[[paste0("y", j)]] <- y[, j]
        }

        return(wideData)
}

# ================================================================================

# Function: compute_missing_rates_for_sim
# Description: Calculate the proportion of missing values for each outcome variable
#              Useful for tracking and reporting dropout rates across time points
#
# Parameters:
#   wide_dat (data.frame): Wide-format data with outcome columns (y1, y2, y3, y4)
#
# Returns:
#   Named numeric vector with missing rates for each time point
#   Names: c("miss_y1", "miss_y2", "miss_y3", "miss_y4")

compute_missing_rates_for_sim <- function(wide_dat) {
        # Find all outcome columns (y1, y2, y3, etc.)
        y_cols <- grep("^y\\d+$", names(wide_dat), value = TRUE)

        if (length(y_cols) == 0) {
                stop("No y* columns found in wide_dat. Expected outcome variables named y1, y2, etc.")
        }

        # Calculate missing rate for each outcome column
        # mean(is.na(x)) = number of NA / total observations
        out <- vapply(y_cols, function(v) mean(is.na(wide_dat[[v]])), numeric(1))

        # Sort outcome columns in numeric order (y1, y2, y3, y4)
        # This ensures consistent ordering regardless of column order in data
        y_cols_ord <- y_cols[order(as.numeric(gsub("y", "", y_cols)))]
        out <- out[y_cols_ord]

        # Rename to include "miss_" prefix for clarity
        names(out) <- paste0("miss_", names(out))

        return(out)
}

# End of data generation functions module
