# ==============================================================================
# 00_pool_load.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Coordinating Center: Load and combine site-level summary statistics
# ==============================================================================

# setup ------------------------------------------------------------------------

library(data.table)
library(tidytable)
library(collapse)
library(stringr)
library(here)

# configuration ----------------------------------------------------------------

ALLOWED_SITES = c(
  "emory",
  "hopkins",
  "ohsu",
  "rush",
  "ucmc",
  "umn",
  "upenn"
)

OUTCOME_LABELS = c(
  noesc_dead  = "No Esc + Dead/Hospice",
  esc_dead    = "Esc + Dead/Hospice",
  esc_alive   = "Esc + Alive",
  noesc_alive = "No Esc + Alive"
)

today = format(Sys.Date(), "%y%m%d")

# helper functions -------------------------------------------------------------

#' Format numbers with commas
format_n = function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' Calculate SD from pooled sum and sum of squares
calculate_sd_from_sums = function(sum_val, sumsq_val, n_val) {
  sqrt((sumsq_val - sum_val^2 / n_val) / (n_val - 1))
}

#' Calculate pooled mean from site-level sums
calculate_pooled_mean = function(sum_vec, n_vec) {
  sum(sum_vec, na.rm = TRUE) / sum(n_vec, na.rm = TRUE)
}

#' Read files matching a pattern from site folders
#' Structure: sites/{site}/*.csv or sites/{site}/upload_to_box/*.csv
read_site_files = function(main_folder, file_stem, allowed_sites = ALLOWED_SITES) {
  
  all_files = character(0)
  
  # Get site folders
  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
  site_folders = site_folders[basename(site_folders) %in% allowed_sites]
  
  for (site_folder in site_folders) {
    site_name = basename(site_folder)
    
    # Try different naming patterns and locations
    patterns = c(
      # Direct in site folder
      file.path(site_folder, paste0(file_stem, "_", site_name, ".csv")),
      file.path(site_folder, paste0(file_stem, "-", site_name, ".csv")),
      # In upload_to_box subfolder
      file.path(site_folder, "upload_to_box", paste0(file_stem, "_", site_name, ".csv")),
      file.path(site_folder, "upload_to_box", paste0(file_stem, "-", site_name, ".csv"))
    )
    
    found_file = patterns[file.exists(patterns)][1]
    
    if (!is.na(found_file)) {
      all_files = c(all_files, found_file)
    }
  }
  
  if (length(all_files) == 0) {
    warning("No files found matching pattern: ", file_stem)
    return(data.table())
  }
  
  message("  Found ", length(all_files), " files matching '", file_stem, "'")
  
  # Read and combine
  file_list = lapply(all_files, function(f) {
    dt = fread(f)
    # Ensure site column exists
    if (!"site" %in% names(dt)) {
      dt$site = str_extract(basename(f), paste(allowed_sites, collapse = "|"))
    }
    dt
  })
  
  combined = rbindlist(file_list, fill = TRUE)
  
  message("    Loaded ", format_n(nrow(combined)), " rows from ", length(file_list), " sites")
  
  return(combined)
}

# load all data ----------------------------------------------------------------

message("\n== Loading site-level data ==")

data_folder = here("sites")  # site data lives in sites/{sitename}/

## table 1 components ----------------------------------------------------------

continuous_raw  = read_site_files(data_folder, "table1_continuous")
binary_raw      = read_site_files(data_folder, "table1_binary")
categorical_raw = read_site_files(data_folder, "table1_categorical")
timing_raw      = read_site_files(data_folder, "table1_timing")
totals_raw      = read_site_files(data_folder, "table1_totals")
flow_raw        = read_site_files(data_folder, "flow_diagram")

## QC files --------------------------------------------------------------------

qc_missing_raw     = read_site_files(data_folder, "qc_missing")
qc_ranges_raw      = read_site_files(data_folder, "qc_ranges")
qc_flags_raw       = read_site_files(data_folder, "qc_flags")
qc_categories_raw  = read_site_files(data_folder, "qc_categories")
qc_diagnostics_raw = read_site_files(data_folder, "qc_diagnostics")

## Exclusion cascade (if available) --------------------------------------------

exclusion_raw = read_site_files(data_folder, "exclusion_cascade")

# pool totals ------------------------------------------------------------------

message("\n== Computing pooled Ns ==")

COHORT_N = totals_raw[, .(
  n_total    = sum(n_total,    na.rm = TRUE),
  n_patients = sum(n_patients, na.rm = TRUE)
), by = outcome_group]

setorder(COHORT_N, outcome_group)

message("  Cohort totals:")
for (i in 1:nrow(COHORT_N)) {
  message("    ", COHORT_N$outcome_group[i], ": ", 
          format_n(COHORT_N$n_total[i]), " encounters, ",
          format_n(COHORT_N$n_patients[i]), " patients")
}

SITE_N = totals_raw[, .(
  n_total = sum(n_total, na.rm = TRUE)
), by = site]

message("  Site totals: ", paste(SITE_N$site, "=", format_n(SITE_N$n_total), collapse = ", "))

# validation -------------------------------------------------------------------

message("\n== Validation ==")

# Check that all sites have all outcome groups
site_group_check = totals_raw[, .(
  n_groups = uniqueN(outcome_group)
), by = site]

if (any(site_group_check$n_groups < 3)) {
  warning("Some sites missing outcome groups:")
  print(site_group_check[n_groups < 3])
}

# Check for missing continuous variables
if (nrow(continuous_raw) > 0) {
  vars_by_site = continuous_raw[, .(vars = paste(unique(variable), collapse = ", ")), by = site]
  message("  Continuous variables by site:")
  print(vars_by_site)
}

message("\n== Data loading complete ==")
message("  Sites loaded: ", paste(unique(totals_raw$site), collapse = ", "))
message("  Total encounters: ", format_n(sum(COHORT_N$n_total)))