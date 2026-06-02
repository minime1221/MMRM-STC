# MMRM-STC
Code for MMRM-based Simulated Treatment Comparison

Description

This repository contains the simulation code used to evaluate a Mixed Model for Repeated Measures (MMRM) extension of Simulated Treatment Comparison (STC).

Requirements
R (version 4.5.1)
mmrm
emmeans
tidyverse

Repository Structure
01_setup.R: package loading and setup
02_covariance_utilities.R: covariance matrix utilities
03_data_generation.R: simulation data generation
04_simulation_functions.R: simulation functions
05_scenario_definitions.R: scenario specifications
06_main_simulation.R: single simulation execution
07_batch_simulation.R: batch simulation runs
