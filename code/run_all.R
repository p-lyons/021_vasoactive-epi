# ==============================================================================
# run_all.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Execute full site-level pipeline: setup → cohort → variables → table
# ==============================================================================

# Clear environment (optional - comment out if you want to preserve objects)
# rm(list = ls())
# gc()

# Record start time
pipeline_start = Sys.time()

message("\n", strrep("=", 70))
message("  CLIF Vasopressor Escalation Pipeline")
message("  Started: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"))
message(strrep("=", 70), "\n")

# ==============================================================================
# 00_setup.R - Load packages, set paths, connect to data
# ==============================================================================

message("\n>>> Running 00_setup.R <<<\n")
source(here::here("code/00_setup.R"))

# ==============================================================================
# 01_cohort.R - Build cohort, link hospitalizations, apply exclusions
# ==============================================================================

message("\n>>> Running 01_cohort.R <<<\n")
source(here::here("code/01_cohort.R"))

# ==============================================================================
# 02_variables.R - Define outcomes, escalation, timing groups, derived vars
# ==============================================================================

message("\n>>> Running 02_variables.R <<<\n")
source(here::here("code/02_variables.R"))

# ==============================================================================
# 03_table.R - Generate poolable Table 1 statistics
# ==============================================================================

message("\n>>> Running 03_table.R <<<\n")
source(here::here("code/03_table.R"))

# ==============================================================================
# Summary
# ==============================================================================

pipeline_end = Sys.time()
elapsed = round(difftime(pipeline_end, pipeline_start, units = "mins"), 1)

message("\n", strrep("=", 70))
message("  Pipeline Complete")
message("  Ended: ", format(pipeline_end, "%Y-%m-%d %H:%M:%S"))
message("  Elapsed: ", elapsed, " minutes")
message(strrep("=", 70))

message("\n  Outputs in upload_to_box/:")
message("    - table1_continuous_*.csv")
message("    - table1_binary_*.csv")
message("    - table1_categorical_*.csv")
message("    - table1_timing_*.csv")
message("    - table1_totals_*.csv")
message("    - flow_diagram_*.csv")
message("    - exclusion_cascade_*.csv")
message("    - qc_*.csv\n")
