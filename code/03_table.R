# ==============================================================================
# 03_table.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Generate poolable summary statistics for Table 1
# ==============================================================================

# Requires: cohort from 02_variables.R

if (!exists("cohort") || !"outcome_group" %in% names(cohort)) {
  cohort = read_parquet(here("proj_tables", "cohort_analytic.parquet"))
}

if (!exists("site_lowercase")) {
  config         = yaml::read_yaml(here("config", "config_clif_pressors.yaml"))
  site_lowercase = config$site_lowercase
}

message(sprintf("\n== Generating Table 1 for site: %s ==", site_lowercase))
message(sprintf("  Cohort size: %d encounters", nrow(cohort)))

cohort_dt = as.data.table(cohort)

# ==============================================================================
# CONTINUOUS VARIABLES - poolable stats
# ==============================================================================

message("\n  Summarizing continuous variables...")

cont_vars = c(
  "age",
  "vw",
  "los_to_t0_d",
  "icu_los_to_t0_d",
  "ne_dose_t0",
  "vp_dose_t0",
  "max_ne_equiv_48h",
  "los_hosp_d",
  "svi_percentile",
  "adi_percentile"
)

# function: poolable stats for continuous variable
summarize_continuous = function(df, var, group_var = "outcome_group") {
  df[!is.na(get(group_var)) & get(group_var) != "other", .(
    variable  = var,
    n         = sum(!is.na(get(var))),
    n_miss    = sum(is.na(get(var))),
    sum       = sum(get(var), na.rm = TRUE),
    sumsq     = sum(get(var)^2, na.rm = TRUE),
    min       = min(get(var), na.rm = TRUE),
    max       = max(get(var), na.rm = TRUE),
    p025      = quantile(get(var), 0.025, na.rm = TRUE),
    p25       = quantile(get(var), 0.25, na.rm = TRUE),
    p50       = quantile(get(var), 0.50, na.rm = TRUE),
    p75       = quantile(get(var), 0.75, na.rm = TRUE),
    p975      = quantile(get(var), 0.975, na.rm = TRUE)
  ), by = group_var]
}

t1_continuous = lapply(cont_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_continuous(cohort_dt, v)
  } else {
    message(sprintf("    ⚠️  Variable '%s' not found", v))
    NULL
  }
}) |>
  Filter(Negate(is.null), x = _) |>
  rbindlist(use.names = TRUE, fill = TRUE)

t1_continuous$site = site_lowercase

# ==============================================================================
# CATEGORICAL VARIABLES - cell counts
# ==============================================================================

message("  Summarizing categorical variables...")

## binary 01 variables ---------------------------------------------------------

binary_vars = c(
  "female_01",
  "epi_01",
  "phenyl_01",
  "dopa_01",
  "a2_01",
  "mb_01",
  "b12_01",
  "imv_at_t0_01",
  "imv_48h_01",
  "crrt_01",
  "code_documented_01",
  "dead_01",
  "hospice_01"
)

# function: summarize binary variable (n and n with value=1)
summarize_binary = function(df, var, group_var = "outcome_group") {
  df[!is.na(get(group_var)) & get(group_var) != "other", .(
    variable = var,
    n        = sum(!is.na(get(var))),
    n_1      = sum(get(var) == 1, na.rm = TRUE)
  ), by = group_var]
}

t1_binary = lapply(binary_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_binary(cohort_dt, v)
  } else {
    message(sprintf("    ⚠️  Binary variable '%s' not found", v))
    NULL
  }
}) |>
  Filter(Negate(is.null), x = _) |>
  rbindlist(use.names = TRUE, fill = TRUE)

t1_binary$site = site_lowercase

## multi-level categorical variables -------------------------------------------

cat_vars = c(
  "race_category",
  "ethnicity_category",
  "age_cat",
  "vw_cat",
  "los_cat",
  "icu_los_cat",
  "code_status_t0"
)

# function: cell counts for categorical variable
summarize_categorical = function(df, var, group_var = "outcome_group") {
  result = df[!is.na(get(group_var)) & get(group_var) != "other", .N, by = c(group_var, var)]
  result[, variable := var]
  setnames(result, var, "category")
  result[, category := as.character(category)]
  setnames(result, "N", "n")
  result[, .(outcome_group, variable, category, n)]
}

t1_categorical = lapply(cat_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_categorical(cohort_dt, v)
  } else {
    message(sprintf("    ⚠️  Categorical variable '%s' not found", v))
    NULL
  }
}) |>
  Filter(Negate(is.null), x = _) |>
  rbindlist(use.names = TRUE, fill = TRUE)

t1_categorical$site = site_lowercase

## timing group variables (3-level: 0=none, 1=before T0, 2=T0 to endpoint) -----

timing_vars = c("imv_timing_group", "crrt_timing_group")

summarize_timing = function(df, var, group_var = "outcome_group") {
  df[!is.na(get(group_var)) & get(group_var) != "other", .(
    variable = var,
    n        = sum(!is.na(get(var))),
    n_0      = sum(get(var) == 0, na.rm = TRUE),
    n_1      = sum(get(var) == 1, na.rm = TRUE),
    n_2      = sum(get(var) == 2, na.rm = TRUE)
  ), by = group_var]
}

t1_timing = lapply(timing_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_timing(cohort_dt, v)
  } else {
    message(sprintf("    ⚠️  Timing variable '%s' not found", v))
    NULL
  }
}) |>
  Filter(Negate(is.null), x = _) |>
  rbindlist(use.names = TRUE, fill = TRUE)

t1_timing$site = site_lowercase

# ==============================================================================
# GROUP TOTALS
# ==============================================================================

message("  Computing group totals...")

message(sprintf("  outcome_group values: %s", 
                paste(unique(cohort_dt$outcome_group), collapse = ", ")))

t1_totals = cohort_dt[!is.na(outcome_group) & outcome_group != "other", .(
  n_total    = .N,
  n_patients = uniqueN(patient_id)
), by = outcome_group]

t1_totals$site = site_lowercase

# ==============================================================================
# QUALITY CONTROL
# ==============================================================================

message("\n  Quality control...")

## mask small cells (n < 5) ----------------------------------------------------

mask_small = function(dt, n_col = "n", threshold = 5) {
  dt = copy(dt)
  dt[get(n_col) > 0 & get(n_col) < threshold, (n_col) := NA_integer_]
  dt
}

n_small_binary = sum(t1_binary$n_1 > 0 & t1_binary$n_1 < 5, na.rm = TRUE)
n_small_cat    = sum(t1_categorical$n > 0 & t1_categorical$n < 5, na.rm = TRUE)

if (n_small_binary > 0 || n_small_cat > 0) {
  message(sprintf("  ⚠️  Masking %d small cells in binary, %d in categorical (n < 5)", 
                  n_small_binary, n_small_cat))
  
  t1_binary      = mask_small(t1_binary, "n_1")
  t1_categorical = mask_small(t1_categorical, "n")
}

## verify totals ---------------------------------------------------------------

if (nrow(t1_totals) == 0) {
  stop("t1_totals is empty - check outcome_group", call. = FALSE)
}

total_from_groups = sum(t1_totals$n_total)
message(sprintf("  Total across groups: %d", total_from_groups))

# ==============================================================================
# FLOW DIAGRAM
# ==============================================================================

message("  Creating flow diagram...")

flow_diagram = data.table(
  step = c(
    "Total encounters meeting T0 criteria",
    "Dead/hospice",
    "Escalated (alive)",
    "Stable (no escalation, alive)"
  ),
  n = c(
    nrow(cohort_dt[outcome_group != "other"]),
    sum(cohort_dt$outcome_group == "dead_hospice", na.rm = TRUE),
    sum(cohort_dt$outcome_group == "escalated", na.rm = TRUE),
    sum(cohort_dt$outcome_group == "stable", na.rm = TRUE)
  ),
  site = site_lowercase
)

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n== Saving outputs ==")

output_dir = here("upload_to_box")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fwrite(t1_continuous,  file.path(output_dir, sprintf("table1_continuous_%s.csv",  site_lowercase)))
fwrite(t1_binary,      file.path(output_dir, sprintf("table1_binary_%s.csv",      site_lowercase)))
fwrite(t1_categorical, file.path(output_dir, sprintf("table1_categorical_%s.csv", site_lowercase)))
fwrite(t1_timing,      file.path(output_dir, sprintf("table1_timing_%s.csv",      site_lowercase)))
fwrite(t1_totals,      file.path(output_dir, sprintf("table1_totals_%s.csv",      site_lowercase)))
fwrite(flow_diagram,   file.path(output_dir, sprintf("flow_diagram_%s.csv",       site_lowercase)))

message(sprintf("  ✅ Saved to: %s", output_dir))
message("    - table1_continuous_*.csv  (n, sum, sumsq, percentiles by group)")
message("    - table1_binary_*.csv      (n, n_1 by group)")
message("    - table1_categorical_*.csv (cell counts by group)")
message("    - table1_timing_*.csv      (n_0, n_1, n_2 by group)")
message("    - table1_totals_*.csv      (group sizes)")
message("    - flow_diagram_*.csv       (cohort flow)")

message("\n== 03_table.R complete ==")
