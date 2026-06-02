# ================================================================================
# COVARIANCE MATRIX UTILITY FUNCTIONS
# ================================================================================
# Purpose: Functions for generating and manipulating covariance structures
#          This module handles the Unstructured (US) covariance matrix generation
# ================================================================================

# Function: gen_covmat_us
# Description: Generate an unstructured covariance matrix for residuals
#              Takes standard deviations and correlations and constructs the
#              full covariance matrix using the outer product method
#
# Parameters:
#   sd (numeric): Vector of standard deviations (length = number of time points)
#   rho (numeric): Vector of pairwise correlations in order:
#                  c(rho12, rho13, rho14, rho23, rho24, rho34)
#
# Returns:
#   Sigma (matrix): Covariance matrix (4 x 4 for typical repeated measures)
#
# Example:
#   gen_covmat_us(sd = c(1.0, 1.3, 1.7, 2.1),
#                 rho = c(0.65, 0.60, 0.58, 0.67, 0.62, 0.68))

gen_covmat_us <- function(sd, rho) {
        # Validate input dimensions
        stopifnot("sd must have length 4 (one per time point)" = length(sd) == 4)
        stopifnot("rho must have length 6 (all pairwise correlations)" = length(rho) == 6)

        # Map correlation vector to matrix positions
        # rho parameter order: rho12, rho13, rho14, rho23, rho24, rho34
        r12 <- rho[1]
        r13 <- rho[2]
        r14 <- rho[3]
        r23 <- rho[4]
        r24 <- rho[5]
        r34 <- rho[6]

        # Initialize correlation matrix with 1's on diagonal
        R <- matrix(1, length(sd), length(sd))

        # Fill in off-diagonal elements (symmetric)
        R[1, 2] <- R[2, 1] <- r12
        R[1, 3] <- R[3, 1] <- r13
        R[1, 4] <- R[4, 1] <- r14
        R[2, 3] <- R[3, 2] <- r23
        R[2, 4] <- R[4, 2] <- r24
        R[3, 4] <- R[4, 3] <- r34

        # Convert correlation matrix to covariance matrix
        # using the formula: Sigma = D * R * D, where D is diag(sd)
        # outer(sd, sd) computes sd_i * sd_j for all pairs
        Sigma <- outer(sd, sd) * R

        return(Sigma)
}

# End of covariance utility functions module
