# ================================================================================
# SCENARIO DEFINITIONS FOR COMPREHENSIVE SIMULATION STUDIES
# ================================================================================
# Purpose: Define 12 core scenarios (× TVEM/TIEM versions) testing different
#          combinations of dropout severity, correlation structure, and baseline
#          covariate overlap between AC and BC trials
# Author: Minhee Seo
# Date: 2024
# ================================================================================
#
# FACTORIAL DESIGN: 3 (Dropout) × 2 (Correlation) × 2 (Overlap) × 2 (TVEM/TIEM)
#
# Dropout Severity (3 levels):
#   - NOdrop: No missing data
#   - Moderate: dropout_inct = c(-3, -2.8, -2.5)
#   - Severe: dropout_inct = c(-3.0, -2.5, 0)
#
# Within-Subject Correlation (2 levels):
#   - Strong: rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68)
#   - Weak: rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25)
#
# Baseline Covariate Overlap (2 levels):
#   - Less overlap: AC meanX = c(0.3, 0.3), BC meanX = c(0.6, 0.6)
#   - More overlap: AC meanX = c(0.5, 0.5), BC meanX = c(0.6, 0.6)
#
# Effect Modifier Structure (2 levels):
#   - TVEM: Time-varying (vectors for beta_em and beta_em_trt)
#   - TIEM: Time-invariant (scalars for beta_em and beta_em_trt)
# ================================================================================

# Helper function to create a scenario template
# Reduces code duplication across the 12 scenarios
make_scenario <- function(dropout = "NOdrop",
                         dropout_inct = NULL,
                         rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68),
                         ac_meanX = c(0.3, 0.3)) {
        if (is.null(dropout_inct)) {
                dropout_inct <- c(-2, -2, -2)
        }

        list(
                AC = list(
                        N = 1000,
                        delta = c(0, 0.5, 1.0, 1.5),
                        b_prog = 0.2,
                        beta_em = c(0.0, 0.2, 0.3, 0.4),
                        beta_em_trt = c(0.0, 0.1, 0.2, 0.3),
                        meanX = ac_meanX,
                        covX = matrix(c(0.04, 0.01, 0.01, 0.04), nrow = 2),
                        rho = rho,
                        sdY = c(1.0, 1.3, 1.7, 2.1),
                        dropout = dropout,
                        dropout_model = "logistic",
                        dropout_inct = dropout_inct,
                        dropout_slope = -1.5,
                        scale_cov = 1
                ),
                BC = list(
                        N = 1000,
                        delta = c(0, 0.5, 1.0, 1.5),
                        b_prog = 0.2,
                        beta_em = c(0.0, 0.2, 0.3, 0.4),
                        beta_em_trt = c(0.0, 0.1, 0.2, 0.3),
                        meanX = c(0.6, 0.6),
                        covX = matrix(c(0.04, 0.01, 0.01, 0.04), nrow = 2),
                        rho = rho,
                        sdY = c(1.0, 1.3, 1.7, 2.1),
                        dropout = dropout,
                        dropout_model = "logistic",
                        dropout_inct = dropout_inct,
                        dropout_slope = -1.5,
                        scale_cov = 1
                )
        )
}

# ================================================================================
# GROUP 1: NO DROPOUT + LESS BASELINE COVARIATE OVERLAP
# ================================================================================

# SCENARIO 1: No dropout, moderate correlation, less overlap
# Key characteristics: Baseline case - no missing data, high correlation
Scenario1 <- c(make_scenario(dropout = "NOdrop", rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68), ac_meanX = c(0.3, 0.3)),
               list(truth_structure = "TVEM"))

# SCENARIO 2: No dropout, low correlation, less overlap
# Key characteristics: Tests robustness to weak correlation structure
Scenario2 <- c(make_scenario(dropout = "NOdrop", rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25), ac_meanX = c(0.3, 0.3)),
               list(truth_structure = "TVEM"))

# ================================================================================
# GROUP 2: MODERATE DROPOUT + LESS BASELINE COVARIATE OVERLAP
# ================================================================================

# SCENARIO 3: Moderate dropout, moderate correlation, less overlap
# Key characteristics: Tests MMRM robustness to moderate missing data (MAR)
Scenario3 <- c(make_scenario(dropout = "MAR",
                             dropout_inct = c(-3, -2.8, -2.5),
                             rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68),
                             ac_meanX = c(0.3, 0.3)),
               list(truth_structure = "TVEM"))

# SCENARIO 4: Moderate dropout, low correlation, less overlap
# Key characteristics: Combines moderate dropout with weak correlation
Scenario4 <- c(make_scenario(dropout = "MAR",
                             dropout_inct = c(-3, -2.8, -2.5),
                             rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25),
                             ac_meanX = c(0.3, 0.3)),
               list(truth_structure = "TVEM"))

# ================================================================================
# GROUP 3: SEVERE DROPOUT + LESS BASELINE COVARIATE OVERLAP
# ================================================================================

# SCENARIO 5: Severe dropout, moderate correlation, less overlap
# Key characteristics: Tests MMRM performance under heavy missingness (MAR)
Scenario5 <- c(make_scenario(dropout = "MAR",
                             dropout_inct = c(-3.0, -2.5, 0),
                             rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68),
                             ac_meanX = c(0.3, 0.3)),
               list(truth_structure = "TVEM"))

# SCENARIO 6: Severe dropout, low correlation, less overlap
# Key characteristics: Most challenging scenario - severe dropout + weak correlation
Scenario6 <- c(make_scenario(dropout = "MAR",
                             dropout_inct = c(-3.0, -2.5, 0),
                             rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25),
                             ac_meanX = c(0.3, 0.3)),
               list(truth_structure = "TVEM"))

# ================================================================================
# GROUP 4: NO DROPOUT + MORE BASELINE COVARIATE OVERLAP
# ================================================================================

# SCENARIO 7: No dropout, moderate correlation, more overlap
# Key characteristics: Better covariate balance between AC and BC trials
Scenario7 <- c(make_scenario(dropout = "NOdrop",
                             rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68),
                             ac_meanX = c(0.5, 0.5)),
               list(truth_structure = "TVEM"))

# SCENARIO 8: No dropout, low correlation, more overlap
# Key characteristics: Good overlap but weak correlation structure
Scenario8 <- c(make_scenario(dropout = "NOdrop",
                             rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25),
                             ac_meanX = c(0.5, 0.5)),
               list(truth_structure = "TVEM"))

# ================================================================================
# GROUP 5: MODERATE DROPOUT + MORE BASELINE COVARIATE OVERLAP
# ================================================================================

# SCENARIO 9: Moderate dropout, moderate correlation, more overlap
# Key characteristics: Moderate dropout with better covariate balance
Scenario9 <- c(make_scenario(dropout = "MAR",
                             dropout_inct = c(-3, -2.8, -2.5),
                             rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68),
                             ac_meanX = c(0.5, 0.5)),
               list(truth_structure = "TVEM"))

# SCENARIO 10: Moderate dropout, low correlation, more overlap
# Key characteristics: Moderate dropout with better overlap but weak correlation
Scenario10 <- c(make_scenario(dropout = "MAR",
                              dropout_inct = c(-3, -2.8, -2.5),
                              rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25),
                              ac_meanX = c(0.5, 0.5)),
                list(truth_structure = "TVEM"))

# ================================================================================
# GROUP 6: SEVERE DROPOUT + MORE BASELINE COVARIATE OVERLAP
# ================================================================================

# SCENARIO 11: Severe dropout, moderate correlation, more overlap
# Key characteristics: Severe dropout but good covariate overlap and correlation
Scenario11 <- c(make_scenario(dropout = "MAR",
                              dropout_inct = c(-3.0, -2.5, 0),
                              rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68),
                              ac_meanX = c(0.5, 0.5)),
                list(truth_structure = "TVEM"))

# SCENARIO 12: Severe dropout, low correlation, more overlap
# Key characteristics: Severe dropout with good overlap but weak correlation
Scenario12 <- c(make_scenario(dropout = "MAR",
                              dropout_inct = c(-3.0, -2.5, 0),
                              rho = c(0.22, 0.18, 0.15, 0.24, 0.20, 0.25),
                              ac_meanX = c(0.5, 0.5)),
                list(truth_structure = "TVEM"))

# Rename list elements for clarity
Scenario1_TVEM <- Scenario1; names(Scenario1_TVEM) <- c("AC", "BC", "truth_structure")
Scenario2_TVEM <- Scenario2; names(Scenario2_TVEM) <- c("AC", "BC", "truth_structure")
Scenario3_TVEM <- Scenario3; names(Scenario3_TVEM) <- c("AC", "BC", "truth_structure")
Scenario4_TVEM <- Scenario4; names(Scenario4_TVEM) <- c("AC", "BC", "truth_structure")
Scenario5_TVEM <- Scenario5; names(Scenario5_TVEM) <- c("AC", "BC", "truth_structure")
Scenario6_TVEM <- Scenario6; names(Scenario6_TVEM) <- c("AC", "BC", "truth_structure")
Scenario7_TVEM <- Scenario7; names(Scenario7_TVEM) <- c("AC", "BC", "truth_structure")
Scenario8_TVEM <- Scenario8; names(Scenario8_TVEM) <- c("AC", "BC", "truth_structure")
Scenario9_TVEM <- Scenario9; names(Scenario9_TVEM) <- c("AC", "BC", "truth_structure")
Scenario10_TVEM <- Scenario10; names(Scenario10_TVEM) <- c("AC", "BC", "truth_structure")
Scenario11_TVEM <- Scenario11; names(Scenario11_TVEM) <- c("AC", "BC", "truth_structure")
Scenario12_TVEM <- Scenario12; names(Scenario12_TVEM) <- c("AC", "BC", "truth_structure")

# ================================================================================
# TIME-INVARIANT EFFECT MODIFIER (TIEM) VERSIONS
# ================================================================================
# Create TIEM versions by converting vector parameters to scalars
# Average the time-varying effects into single constant values

make_scenario_tiem <- function(base_scenario) {
        scenario <- base_scenario
        scenario$AC$beta_em <- 0.2
        scenario$AC$beta_em_trt <- 0.3
        scenario$BC$beta_em <- 0.2
        scenario$BC$beta_em_trt <- 0.3
        scenario$truth_structure <- "TIEM"
        return(scenario)
}

Scenario1_TIEM <- make_scenario_tiem(Scenario1_TVEM)
Scenario2_TIEM <- make_scenario_tiem(Scenario2_TVEM)
Scenario3_TIEM <- make_scenario_tiem(Scenario3_TVEM)
Scenario4_TIEM <- make_scenario_tiem(Scenario4_TVEM)
Scenario5_TIEM <- make_scenario_tiem(Scenario5_TVEM)
Scenario6_TIEM <- make_scenario_tiem(Scenario6_TVEM)
Scenario7_TIEM <- make_scenario_tiem(Scenario7_TVEM)
Scenario8_TIEM <- make_scenario_tiem(Scenario8_TVEM)
Scenario9_TIEM <- make_scenario_tiem(Scenario9_TVEM)
Scenario10_TIEM <- make_scenario_tiem(Scenario10_TVEM)
Scenario11_TIEM <- make_scenario_tiem(Scenario11_TVEM)
Scenario12_TIEM <- make_scenario_tiem(Scenario12_TVEM)

# Create a comprehensive master list of all scenarios
# Organized by type: 12 TVEM variants + 12 TIEM variants = 24 total scenarios
scenarios_list <- list(
        # TVEM scenarios (Time-Varying Effect Modifier)
        Scenario1_TVEM = Scenario1_TVEM,
        Scenario2_TVEM = Scenario2_TVEM,
        Scenario3_TVEM = Scenario3_TVEM,
        Scenario4_TVEM = Scenario4_TVEM,
        Scenario5_TVEM = Scenario5_TVEM,
        Scenario6_TVEM = Scenario6_TVEM,
        Scenario7_TVEM = Scenario7_TVEM,
        Scenario8_TVEM = Scenario8_TVEM,
        Scenario9_TVEM = Scenario9_TVEM,
        Scenario10_TVEM = Scenario10_TVEM,
        Scenario11_TVEM = Scenario11_TVEM,
        Scenario12_TVEM = Scenario12_TVEM,

        # TIEM scenarios (Time-Invariant Effect Modifier)
        Scenario1_TIEM = Scenario1_TIEM,
        Scenario2_TIEM = Scenario2_TIEM,
        Scenario3_TIEM = Scenario3_TIEM,
        Scenario4_TIEM = Scenario4_TIEM,
        Scenario5_TIEM = Scenario5_TIEM,
        Scenario6_TIEM = Scenario6_TIEM,
        Scenario7_TIEM = Scenario7_TIEM,
        Scenario8_TIEM = Scenario8_TIEM,
        Scenario9_TIEM = Scenario9_TIEM,
        Scenario10_TIEM = Scenario10_TIEM,
        Scenario11_TIEM = Scenario11_TIEM,
        Scenario12_TIEM = Scenario12_TIEM
)

cat("=== Scenario Definitions Loaded ===\n")
cat("Available scenarios:\n")

# Print TVEM scenarios
cat("\nTVEM Scenarios (Time-Varying Effect Modifier):\n")
tvem_scenarios <- grep("_TVEM$", names(scenarios_list), value = TRUE)
for (i in seq_along(tvem_scenarios)) {
        scenario_name <- tvem_scenarios[i]
        scenario <- scenarios_list[[scenario_name]]
        dropout_type <- scenario$AC$dropout
        ac_overlap <- scenario$AC$meanX[1]
        rho_val <- scenario$AC$rho[1]

        overlap_label <- if (ac_overlap == 0.3) "less" else "more"
        rho_label <- if (rho_val > 0.5) "high" else "low"

        cat(sprintf("  %2d. %s | Dropout: %-6s | Overlap: %-4s | Corr: %-4s\n",
                    i, scenario_name, dropout_type, overlap_label, rho_label))
}

# Print TIEM scenarios
cat("\nTIEM Scenarios (Time-Invariant Effect Modifier):\n")
tiem_scenarios <- grep("_TIEM$", names(scenarios_list), value = TRUE)
for (i in seq_along(tiem_scenarios)) {
        scenario_name <- tiem_scenarios[i]
        scenario <- scenarios_list[[scenario_name]]
        dropout_type <- scenario$AC$dropout
        ac_overlap <- scenario$AC$meanX[1]
        rho_val <- scenario$AC$rho[1]

        overlap_label <- if (ac_overlap == 0.3) "less" else "more"
        rho_label <- if (rho_val > 0.5) "high" else "low"

        cat(sprintf("  %2d. %s | Dropout: %-6s | Overlap: %-4s | Corr: %-4s\n",
                    12 + i, scenario_name, dropout_type, overlap_label, rho_label))
}

cat("\nTotal available scenarios: ", length(scenarios_list), "\n")
cat("Use scenarios_list$<scenario_name> to access individual scenarios\n\n")

# End of scenario definitions module
