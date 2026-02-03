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

# ==============================================================================
# CONTINUOUS VARIABLES
# ==============================================================================

message("\n  Summarizing continuous variables...")

cont_vars = c(
  "age",
  "vw",
  "los_to_t0_d",
  "icu_los_to_t0_h",
  "ne_dose_t0",
  "vp_dose_t0",
  "max_ne_equiv_48h",
  "los_hosp_d"
)

summarize_continuous = function(df, var, group_var = "outcome_group") {
  df = as.data.table(df)
  df[!is.na(get(group_var)), .(
    variable = var,
    n        = sum(!is.na(get(var))),
    n_miss   = sum(is.na(get(var))),
    sum      = sum(get(var), na.rm = TRUE),
    sumsq    = sum(get(var)^2, na.rm = TRUE),
    min      = min(get(var), na.rm = TRUE),
    max      = max(get(var), na.rm = TRUE),
    p25      = quantile(get(var), 0.25, na.rm = TRUE),
    p50      = quantile(get(var), 0.50, na.rm = TRUE),
    p75      = quantile(get(var), 0.75, na.rm = TRUE)
  ), by = group_var]
}

t1_continuous = lapply(cont_vars, function(v) {
  if (v %in% names(cohort)) {
    summarize_continuous(cohort, v)
  } else {
    message(sprintf("    ⚠️  Variable '%s' not found", v))
    NULL
  }
}) |>
  Filter(Negate(is.null), x = _) |>
  rbindlist(use.names = TRUE, fill = TRUE)

t1_continuous$site = site_lowercase

# ==============================================================================
# CATEGORICAL VARIABLES
# ==============================================================================

message("  Summarizing categorical variables...")

cat_vars = c(
  "female_01",
  "race_category",
  "ethnicity_category",
  "age_cat",
  "vw_cat",
  "imv_at_t0_01",
  "epi_01",
  "phenyl_01",
  "dopa_01",
  "imv_48h_01",
  "dead_01",
  "hospice_01"
)

summarize_categorical = function(df, var, group_var = "outcome_group") {
  df = as.data.table(df)
  result = df[!is.na(get(group_var)), .N, by = c(group_var, var)]
  result[, variable := var]
  setnames(result, var, "category")
  result[, category := as.character(category)]
  setnames(result, "N", "n")
  result[, .(outcome_group, variable, category, n)]
}

t1_categorical = lapply(cat_vars, function(v) {
  if (v %in% names(cohort)) {
    summarize_categorical(cohort, v)
  } else {
    message(sprintf("    ⚠️  Variable '%s' not found", v))
    NULL
  }
}) |>
  Filter(Negate(is.null), x = _) |>
  rbindlist(use.names = TRUE, fill = TRUE)

t1_categorical$site = site_lowercase

# ==============================================================================
# GROUP TOTALS
# ==============================================================================

message("  Computing group totals...")

cohort_dt = as.data.table(cohort)

message(sprintf("  outcome_group values: %s", 
                paste(unique(cohort_dt$outcome_group), collapse = ", ")))

t1_totals = cohort_dt[!is.na(outcome_group), .(
  n_total    = .N,
  n_patients = uniqueN(patient_id)
), by = outcome_group]

t1_totals$site = site_lowercase

# ==============================================================================
# QUALITY CONTROL
# ==============================================================================

message("\n  Quality control...")

## mask small cells ------------------------------------------------------------

small_cat = fsubset(t1_categorical, n < 5 & n > 0)

if (nrow(small_cat) > 0) {
  message(sprintf("  ⚠️  Masking %d small cells (n < 5)", nrow(small_cat)))
  t1_categorical = ftransform(t1_categorical,
    n = fifelse(n > 0 & n < 5, NA_integer_, n)
  )
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
    "Stable/de-escalated (alive)"
  ),
  n = c(
    nrow(cohort_dt),
    sum(cohort_dt$dead_hospice_01 == 1, na.rm = TRUE),
    sum(cohort_dt$escalated_01 == 1 & cohort_dt$dead_hospice_01 == 0, na.rm = TRUE),
    sum(cohort_dt$escalated_01 == 0 & cohort_dt$dead_hospice_01 == 0, na.rm = TRUE)
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
fwrite(t1_categorical, file.path(output_dir, sprintf("table1_categorical_%s.csv", site_lowercase)))
fwrite(t1_totals,      file.path(output_dir, sprintf("table1_totals_%s.csv",      site_lowercase)))
fwrite(flow_diagram,   file.path(output_dir, sprintf("flow_diagram_%s.csv",       site_lowercase)))

message(sprintf("  ✅ Saved to: %s", output_dir))
message("    - table1_continuous_*.csv")
message("    - table1_categorical_*.csv")
message("    - table1_totals_*.csv")
message("    - flow_diagram_*.csv")

message("\n== 03_table.R complete ==")
