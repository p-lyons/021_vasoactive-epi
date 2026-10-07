# ==============================================================================
# 00_pool_load.R
# Vasopressor Escalation in Refractory Distributive Shock - CLIF Consortium
# Coordinating site: load and combine site-level summary files
#
# Site files carry a `hospital` column. Rows with hospital == "site_total" are
# whole-site summaries; all other rows are single hospitals. Every pooled
# (site-level) object below uses ONLY site_total rows, and every hospital-level
# object uses ONLY hospital rows, so nothing is counted twice.
#
# Folder layout: sites/{site}/, sites/{site}/upload_to_box/, or
#   sites/{site}/output/ (Sankey files)
# ==============================================================================

# setup ------------------------------------------------------------------------

library(data.table)
library(tidytable)
library(collapse)
library(stringr)
library(here)

# configuration ----------------------------------------------------------------

site_details  = fread(here("config", "clif_sites.csv"))
ALLOWED_SITES = tolower(site_details$site_name)

OUTCOME_LABELS = c(
  noesc_dead  = "No Esc + Dead/Hospice",
  esc_dead    = "Esc + Dead/Hospice",
  esc_alive   = "Esc + Alive",
  noesc_alive = "No Esc + Alive"
)

# rescue agents: a site with zero use across ALL its encounters is treated as
# not capturing that agent and is left out of its numerator and denominator
RESCUE_VARS = c(
  "mb_01",
  "b12_01",
  "a2_01"
)

today = format(Sys.Date(), "%y%m%d")

rm(site_details)

# helper functions -------------------------------------------------------------

#' Format numbers with commas
format_n = function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' Calculate SD from pooled sum and sum of squares
calculate_sd_from_sums = function(sum_val, sumsq_val, n_val) {
  sqrt((sumsq_val - sum_val^2 / n_val) / (n_val - 1))
}

#' Read files matching a pattern from site folders
#' Structure: sites/{site}/, sites/{site}/upload_to_box/, sites/{site}/output/
read_site_files = function(main_folder, file_stem, allowed_sites = ALLOWED_SITES) {

  all_files = character(0)

  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
  site_folders = site_folders[basename(site_folders) %in% allowed_sites]

  for (site_folder in site_folders) {
    site_name = basename(site_folder)

    patterns = c(
      file.path(site_folder,                  paste0(file_stem, "_", site_name, ".csv")),
      file.path(site_folder,                  paste0(file_stem, "-", site_name, ".csv")),
      file.path(site_folder, "upload_to_box", paste0(file_stem, "_", site_name, ".csv")),
      file.path(site_folder, "upload_to_box", paste0(file_stem, "-", site_name, ".csv")),
      file.path(site_folder, "output",        paste0(file_stem, "_", site_name, ".csv")),
      file.path(site_folder, "output",        paste0(file_stem, "-", site_name, ".csv"))
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

  file_list = lapply(all_files, function(f) {
    dt = fread(f)
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

data_folder = here("sites")

## table 1 components (hospital-stratified) ------------------------------------

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
qc_hospital_raw    = read_site_files(data_folder, "qc_hospital")

## exclusion cascade and Sankey files ------------------------------------------

exclusion_raw      = read_site_files(data_folder, "exclusion_cascade")
sankey_trans_raw   = read_site_files(data_folder, "sankey_transitions")
sankey_summary_raw = read_site_files(data_folder, "sankey_site_summary")

# hospital stratum -------------------------------------------------------------

message("\n== Checking hospital stratum ==")

hosp_tables = list(
  continuous_raw  = continuous_raw,
  binary_raw      = binary_raw,
  categorical_raw = categorical_raw,
  timing_raw      = timing_raw,
  totals_raw      = totals_raw,
  flow_raw        = flow_raw,
  qc_missing_raw  = qc_missing_raw
)

for (nm in names(hosp_tables)) {
  dt = hosp_tables[[nm]]
  if (nrow(dt) > 0 && !"hospital" %in% names(dt)) {
    stop(
      sprintf("%s has no 'hospital' column. Re-run the current 03_table.R at that site.", nm),
      call. = FALSE
    )
  }
}

rm(hosp_tables, nm, dt)

# hospital IDs can be numeric at some sites; key = site + hospital
# (one statement per table: := inside a list() would modify a copy)

continuous_raw[,  `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
binary_raw[,      `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
categorical_raw[, `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
timing_raw[,      `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
totals_raw[,      `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
flow_raw[,        `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
qc_missing_raw[,  `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]

if (nrow(qc_diagnostics_raw) > 0) {
  qc_diagnostics_raw[, hospital := as.character(hospital)]
}

if (nrow(qc_hospital_raw) > 0) {
  qc_hospital_raw[, `:=`(hospital = as.character(hospital), hosp_key = paste(site, hospital, sep = "__"))]
}

## site-level (site_total) and hospital-level splits ---------------------------

continuous_site  = continuous_raw[hospital == "site_total"]
binary_site      = binary_raw[hospital == "site_total"]
categorical_site = categorical_raw[hospital == "site_total"]
timing_site      = timing_raw[hospital == "site_total"]
totals_site      = totals_raw[hospital == "site_total"]
flow_site        = flow_raw[hospital == "site_total"]

binary_hosp = binary_raw[hospital != "site_total"]
totals_hosp = totals_raw[hospital != "site_total"]

message("  Sites: ", uniqueN(totals_site$site),
        " | Hospitals: ", uniqueN(totals_hosp$hosp_key))

# rescue agents not captured ---------------------------------------------------

message("\n== Rescue agents: sites with zero use ==")

rescue_use = binary_site[
  variable %chin% RESCUE_VARS,
  .(n_1 = sum(n_1, na.rm = TRUE)),
  by = .(site, variable)
]

NOT_CAPTURED = rescue_use[n_1 == 0, .(site, variable)]

if (nrow(NOT_CAPTURED) > 0) {
  message("  Treated as not captured (left out of numerator and denominator):")
  print(NOT_CAPTURED)
} else {
  message("  All sites recorded some use of every rescue agent")
}

binary_site = binary_site[!NOT_CAPTURED, on = .(site, variable)]
binary_hosp = binary_hosp[!NOT_CAPTURED, on = .(site, variable)]

rm(rescue_use)

# pool totals ------------------------------------------------------------------

message("\n== Computing pooled Ns ==")

COHORT_N = totals_site[, .(
  n_total    = sum(n_total,    na.rm = TRUE),
  n_patients = sum(n_patients, na.rm = TRUE)
), by = outcome_group]

setorder(COHORT_N, outcome_group)

message("  Cohort totals:")
for (i in seq_len(nrow(COHORT_N))) {
  message("    ", COHORT_N$outcome_group[i], ": ",
          format_n(COHORT_N$n_total[i]), " encounters, ",
          format_n(COHORT_N$n_patients[i]), " patients")
}

SITE_N = totals_site[, .(
  n_total = sum(n_total, na.rm = TRUE)
), by = site]

message("  Site totals: ", paste(SITE_N$site, "=", format_n(SITE_N$n_total), collapse = ", "))

# validation -------------------------------------------------------------------

message("\n== Validation ==")

## all four outcome groups at every site ---------------------------------------

site_group_check = totals_site[, .(
  n_groups = uniqueN(outcome_group)
), by = site]

if (any(site_group_check$n_groups < length(OUTCOME_LABELS))) {
  warning("Some sites missing outcome groups:")
  print(site_group_check[n_groups < length(OUTCOME_LABELS)])
}

## hospital rows sum to site_total ---------------------------------------------

hosp_sum = totals_hosp[, .(n_hosp = sum(n_total)), by = .(site, outcome_group)]
site_sum = totals_site[, .(n_site = sum(n_total)), by = .(site, outcome_group)]

hosp_check = merge(
  hosp_sum,
  site_sum,
  by  = c("site", "outcome_group"),
  all = TRUE
)

hosp_mismatch = hosp_check[is.na(n_hosp) | is.na(n_site) | n_hosp != n_site]

if (nrow(hosp_mismatch) > 0) {
  warning("Hospital rows do not sum to site_total:")
  print(hosp_mismatch)
} else {
  message("  ✅ Hospital rows sum to site_total at every site")
}

rm(hosp_sum, site_sum, hosp_check, hosp_mismatch, site_group_check)

# Sankey: pooled transitions ---------------------------------------------------

message("\n== Pooling Sankey transitions ==")

if (nrow(sankey_trans_raw) > 0) {

  SANKEY_POOLED = sankey_trans_raw[, .(
    n       = sum(n, na.rm = TRUE),
    n_sites = uniqueN(site)
  ), by = .(block_from, state_from, state_to)]

  setorder(SANKEY_POOLED, block_from, state_from, state_to)

  if (!dir.exists(here("output", "sankey"))) {
    dir.create(here("output", "sankey"), recursive = TRUE)
  }

  fwrite(SANKEY_POOLED, here("output", "sankey", paste0("sankey_transitions_pooled_", today, ".csv")))

  message("  Pooled ", format_n(nrow(SANKEY_POOLED)), " transition cells from ",
          uniqueN(sankey_trans_raw$site), " sites")

} else {
  SANKEY_POOLED = data.table()
  message("  No Sankey transition files found")
}

# Sankey: cascade must match Table 1 cascade -----------------------------------

if (nrow(sankey_summary_raw) > 0 && nrow(exclusion_raw) > 0) {

  step_cols = grep("^n_0", names(sankey_summary_raw), value = TRUE)

  sankey_cascade = melt(
    sankey_summary_raw[, c("site", step_cols), with = FALSE],
    id.vars       = "site",
    variable.name = "step",
    value.name    = "n_sankey"
  )

  sankey_cascade[, step := sub("^n_", "", as.character(step))]

  CASCADE_CHECK = merge(
    exclusion_raw[, .(site, step, n_table1 = n_remaining)],
    sankey_cascade,
    by  = c("site", "step"),
    all = TRUE
  )

  CASCADE_CHECK[, match := !is.na(n_table1) & !is.na(n_sankey) & n_table1 == n_sankey]

  if (any(!CASCADE_CHECK$match)) {
    warning("Sankey cascade differs from Table 1 cascade:")
    print(CASCADE_CHECK[match == FALSE])
  } else {
    message("  ✅ Sankey cascade matches Table 1 cascade at every site")
  }

  rm(step_cols, sankey_cascade)

} else {
  CASCADE_CHECK = data.table()
  message("  Cascade check skipped (Sankey summary or exclusion cascade missing)")
}

message("\n== Data loading complete ==")
message("  Sites loaded: ", paste(sort(unique(totals_site$site)), collapse = ", "))
