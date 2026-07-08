# MMRM-STC
Code for MMRM-based Simulated Treatment Comparison

This repository contains the simulation code accompanying the manuscript

> **Extending Simulated Treatment Comparison to Longitudinal Continuous Outcomes Using Mixed Models for Repeated Measures**

## Overview

This repository contains the R code used to reproduce the simulation study presented in the accompanying manuscript. The simulations evaluate an extension of Simulated Treatment Comparison (STC) for longitudinal continuous outcomes using Mixed Models for Repeated Measures (MMRM).

The simulation compares four indirect treatment comparison approaches:

- Bucher's method
- STC--ANCOVA
- STC--MMRM (TIEM)
- STC--MMRM (TVEM)

under varying levels of

- covariate overlap,
- monotone MAR dropout,
- within-subject correlation, and
- treatment--effect modifier structures (TVEM and TIEM truth).

## Requirements

The simulations were developed using **R 4.5.1**.

Required packages include

- tidyverse
- mmrm
- emmeans
- MASS
- progress

## Repository structure

```
R/
├── 01_setup.R
├── 02_covariance_utilities.R
├── 03_data_generation.R
├── 04_simulation_functions.R
├── 05_scenario_definitions.R
├── 06_main_simulation.R
└── 07_batch_simulation.R
```

### Script descriptions

| Script | Description |
|---------|-------------|
|01_setup.R|Loads required packages and global options.|
|02_covariance_utilities.R|Utility functions for constructing residual covariance matrices.|
|03_data_generation.R|Generates simulated AC and BC trial datasets under TVEM or TIEM truth, with optional monotone MAR dropout.|
|04_simulation_functions.R|Fits Bucher, STC--ANCOVA, STC--MMRM (TIEM), and STC--MMRM (TVEM) for one simulation replicate.|
|05_scenario_definitions.R|Defines the 24 simulation scenarios (12 TVEM truth and 12 TIEM truth).|
|06_main_simulation.R|Runs a single selected simulation scenario.|
|07_batch_simulation.R|Runs all simulation scenarios in batch and saves the results.|


