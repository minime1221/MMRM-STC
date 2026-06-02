# ================================================================================
# SETUP AND LIBRARY LOADING
# ================================================================================
# Purpose: Load all required libraries and set global options for the simulation
# ================================================================================

# Load required packages
library(tidyverse)      # Data manipulation and visualization (ggplot2, dplyr, tidyr)
library(progress)       # Progress bar for long simulations
library(mmrm)           # Mixed model for repeated measures
# library(MASS)           # For mvrnorm (multivariate normal distribution)
library(emmeans)        # For estimated marginal means (optional, for post-hoc analysis)
library(ggplot2)        # Plotting (included in tidyverse, but explicit for clarity)
library(patchwork)      # Combining multiple plots
library(flextable)      # Publication-ready tables (for Word document export)
library(scales)         # Scaling and formatting (e.g., percent)

# Global options
options(scipen = 999)   # Disable scientific notation for cleaner output
options(width = 120)    # Set console width for better readability

# Set seed for reproducibility (can be overridden in individual simulations)
# set.seed(12345)

cat("=== Simulation Setup Complete ===\n")
cat("Libraries loaded successfully.\n")
cat("Ready to run MMRM STC/Bucher simulation study.\n")
