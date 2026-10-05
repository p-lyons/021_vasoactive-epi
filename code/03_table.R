# ==============================================================================
# 03_table.R
# Vasopressor Escalation in Refractory Distributive Shock - CLIF Consortium
# Generate poolable summary statistics for Table 1
# Stratified by hospital (column `hospital`) within site.
#   - Hospital = hospital_id_t0 from 02_variables.R (hospital at T0); switch
#     hosp_var below to "hospital_id" to use first-ADT hospital instead.
#   - Every table1/flow output also carries hospital = "site_total" rows
#     (whole-site summary). Exclude site_total rows when summing across
#     hospitals.
# Small-cell policy: counts 1-4 are masked (NA) only for sensitive demographic
#   variables (sex, race, ethnicity, language). All other counts are not PHI
#   and are left unmasked. Files are pooled at the coordinating site under the
#   consortium DUA before anything is shared.
# Dates: per-hospital study_start / study_end are reported as year-month.
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

# hospital stratum -------------------------------------------------------------

hosp_var = "hospital_id_t0"

if (!hosp_var %in% names(cohort_dt)) {
  stop(sprintf("'%s' not found in cohort. Re-run 02_variables.R.", hosp_var), call. = FALSE)
}

cohort_dt[, hospital := as.character(get(hosp_var))]
cohort_dt[is.na(hospital), hospital := "unknown"]

message(sprintf("  Hospitals: %d", uniqueN(cohort_dt$hospital)))
print(cohort_dt[, .N, by = hospital][order(hospital)])

# site-level copy for site_total rows (see header)
cohort_site_dt = copy(cohort_dt)[, hospital := "site_total"]
cohort_both_dt = rbindlist(list(cohort_dt, cohort_site_dt), use.names = TRUE)   # hospital rows + site_total rows

# sensitive variables for small-cell masking (see header) ----------------------

sensitive_binary = c(
  "female_01",
  "white_01",
  "hispanic_01",
  "english_01"
)

sensitive_categorical = c(
  "race_category",
  "ethnicity_category"
)

mask_threshold = 5L

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
  "los_from_t0_d",
  "svi_percentile",
  "adi_percentile"
)

# function: poolable stats for continuous variable
summarize_continuous = function(df, var, group_var = "outcome_group", strata = "hospital") {
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
  ), by = c(strata, group_var)]
}

t1_continuous = lapply(cont_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    rbindlist(list(
      summarize_continuous(cohort_dt, v),
      summarize_continuous(cohort_site_dt, v)
    ), use.names = TRUE)
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
  "white_01",
  "hispanic_01",
  "english_01",
  "academic_01",
  "peak_covid_01",
  "full_code_01",
  "major_procedure_01",
  "epi_01",
  "phenyl_01",
  "dopa_01",
  "a2_01",
  "mb_01",
  "b12_01",
  "imv_at_t0_01",
  "imv_48h_01",
  "crrt_01",
  "crrt_at_t0_01",
  "code_documented_01",
  "dead_01",
  "hospice_01"
)

# function: summarize binary variable (n and n with value=1)
summarize_binary = function(df, var, group_var = "outcome_group", strata = "hospital") {
  df[!is.na(get(group_var)) & get(group_var) != "other", .(
    variable = var,
    n        = sum(!is.na(get(var))),
    n_1      = sum(get(var) == 1, na.rm = TRUE)
  ), by = c(strata, group_var)]
}

t1_binary = lapply(binary_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_binary(cohort_both_dt, v)
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
summarize_categorical = function(df, var, group_var = "outcome_group", strata = "hospital") {
  result = df[!is.na(get(group_var)) & get(group_var) != "other", .N, by = c(strata, group_var, var)]
  result[, variable := var]
  setnames(result, var, "category")
  result[, category := as.character(category)]
  setnames(result, "N", "n")
  result[, .(hospital, outcome_group, variable, category, n)]
}

t1_categorical = lapply(cat_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_categorical(cohort_both_dt, v)
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

summarize_timing = function(df, var, group_var = "outcome_group", strata = "hospital") {
  df[!is.na(get(group_var)) & get(group_var) != "other", .(
    variable = var,
    n        = sum(!is.na(get(var))),
    n_0      = sum(get(var) == 0, na.rm = TRUE),
    n_1      = sum(get(var) == 1, na.rm = TRUE),
    n_2      = sum(get(var) == 2, na.rm = TRUE)
  ), by = c(strata, group_var)]
}

t1_timing = lapply(timing_vars, function(v) {
  if (v %in% names(cohort_dt)) {
    summarize_timing(cohort_both_dt, v)
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

t1_totals = cohort_both_dt[!is.na(outcome_group)
                          & outcome_group != "other", .(
  n_total    = .N,
  n_patients = uniqueN(patient_id)
), by = .(hospital, outcome_group)]

t1_totals$site = site_lowercase

# ==============================================================================
# QUALITY CONTROL
# ==============================================================================

message("\n  Quality control...")

## mask small cells for sensitive demographic variables only -------------------

n_small_binary = t1_binary[
  variable %chin% sensitive_binary & n_1 > 0 & n_1 < mask_threshold,
  .N
]

n_small_cat = t1_categorical[
  variable %chin% sensitive_categorical & n > 0 & n < mask_threshold,
  .N
]

message(sprintf("  Masking %d sensitive binary cells and %d sensitive categorical cells (n < %d)",
                n_small_binary, n_small_cat, mask_threshold))

t1_binary[
  variable %chin% sensitive_binary & n_1 > 0 & n_1 < mask_threshold,
  n_1 := NA_integer_
]

t1_categorical[
  variable %chin% sensitive_categorical & n > 0 & n < mask_threshold,
  n := NA_integer_
]

## verify totals ---------------------------------------------------------------

if (nrow(t1_totals) == 0) {
  stop("t1_totals is empty - check outcome_group", call. = FALSE)
}

total_from_groups = sum(t1_totals[hospital != "site_total"]$n_total)
message(sprintf("  Total across groups: %d", total_from_groups))

# ==============================================================================
# FLOW DIAGRAM
# ==============================================================================

message("  Creating flow diagram...")

# one block of 5 rows per hospital + site_total
flow_steps = c(
  "Total encounters meeting T0 criteria",
  "No escalation + dead/hospice",
  "Escalated + dead/hospice",
  "Escalated + alive",
  "No escalation + alive"
)

flow_diagram = cohort_both_dt[, .(
  step = flow_steps,
  n    = c(
    sum(outcome_group != "other", na.rm = TRUE),
    sum(outcome_group == "noesc_dead", na.rm = TRUE),
    sum(outcome_group == "esc_dead", na.rm = TRUE),
    sum(outcome_group == "esc_alive", na.rm = TRUE),
    sum(outcome_group == "noesc_alive", na.rm = TRUE)
  )
), by = hospital]

flow_diagram$site = site_lowercase

# ==============================================================================
# QC DIAGNOSTICS
# ==============================================================================

message("  Creating QC diagnostics...")

## missingness report ----------------------------------------------------------

all_vars = c(
  "age", "female_01", "race_category", "ethnicity_category", "vw",
  "code_status_t0", "los_to_t0_d", "icu_los_to_t0_d",
  "ne_dose_t0", "vp_dose_t0", "max_ne_equiv_48h",
  "svi_percentile", "adi_percentile", "los_hosp_d",
  "imv_dttm", "crrt_dttm", "census_block_code"
)

fn_qc_missing = function(dt, hosp_label) {
  data.table(
    hospital = hosp_label,
    variable = all_vars,
    n_total  = nrow(dt),
    n_miss   = sapply(all_vars, function(v) {
      if (v %in% names(dt)) sum(is.na(dt[[v]])) else NA_integer_
    }),
    site = site_lowercase
  )
}

qc_missing = rbindlist(c(
  lapply(sort(unique(cohort_dt$hospital)), function(h) fn_qc_missing(cohort_dt[hospital == h], h)),
  list(fn_qc_missing(cohort_dt, "site_total"))
))
qc_missing[, pct_miss := round(n_miss / n_total * 100, 1)]

## continuous variable ranges --------------------------------------------------

cont_vars_qc = c("age", "vw", "los_to_t0_d", "icu_los_to_t0_d",
                 "ne_dose_t0", "vp_dose_t0", "max_ne_equiv_48h",
                 "svi_percentile", "adi_percentile", "los_hosp_d")

qc_ranges = rbindlist(lapply(cont_vars_qc, function(v) {
  if (v %in% names(cohort_dt)) {
    data.table(
      variable = v,
      n        = sum(!is.na(cohort_dt[[v]])),
      min      = min(cohort_dt[[v]], na.rm = TRUE),
      p01      = quantile(cohort_dt[[v]], 0.01, na.rm = TRUE),
      p25      = quantile(cohort_dt[[v]], 0.25, na.rm = TRUE),
      median   = median(cohort_dt[[v]], na.rm = TRUE),
      p75      = quantile(cohort_dt[[v]], 0.75, na.rm = TRUE),
      p99      = quantile(cohort_dt[[v]], 0.99, na.rm = TRUE),
      max      = max(cohort_dt[[v]], na.rm = TRUE),
      site     = site_lowercase
    )
  }
}), fill = TRUE)

## plausibility flags ----------------------------------------------------------

qc_flags = data.table(
  check = c(
    "age_under_18",
    "age_over_110",
    "ne_dose_over_5",
    "vp_dose_over_0.1",
    "los_negative",
    "t0_before_admission",
    "t0_after_discharge",
    "endpoint_after_discharge"
  ),
  n_flagged = c(
    sum(cohort_dt$age < 18, na.rm = TRUE),
    sum(cohort_dt$age > 110, na.rm = TRUE),
    sum(cohort_dt$ne_dose_t0 > 5, na.rm = TRUE),
    sum(cohort_dt$vp_dose_t0 > 0.1, na.rm = TRUE),
    sum(cohort_dt$los_to_t0_d < 0, na.rm = TRUE),
    sum(cohort_dt$t0_dttm < cohort_dt$admission_dttm, na.rm = TRUE),
    sum(cohort_dt$t0_dttm > cohort_dt$discharge_dttm, na.rm = TRUE),
    sum(cohort_dt$endpoint_dttm > cohort_dt$discharge_dttm, na.rm = TRUE)
  ),
  site = site_lowercase
)

## categorical value inventory -------------------------------------------------

cat_vars_qc = c("race_category", "ethnicity_category", "code_status_t0",
                "discharge_category", "outcome_group")

qc_categories = rbindlist(lapply(cat_vars_qc, function(v) {
  if (v %in% names(cohort_dt)) {
    cohort_dt[, .(n = .N), by = c(v)][, .(
      variable = v,
      category = as.character(get(v)),
      n        = n,
      site     = site_lowercase
    )]
  }
}), fill = TRUE)

# same sensitive-variable masking as Table 1
qc_categories[
  variable %chin% sensitive_categorical & n > 0 & n < mask_threshold,
  n := NA_integer_
]

## site diagnostics summary ----------------------------------------------------

# Key metrics for cross-site validation; run per hospital + site_total
fn_qc_diagnostics = function(dt, hosp_label) {
  data.table(
    hospital = hosp_label,
    metric = c(
      "n_encounters",
      "n_patients",
      "study_start_ym",
      "study_end_ym",
      "pct_female",
      "median_age",
      "pct_white",
      "pct_hispanic",
      "pct_english",
      "pct_academic",
      "pct_peak_covid",
      "pct_dead_hospice",
      "pct_escalated",
      "pct_svi_linked",
      "pct_adi_linked",
      "pct_code_documented",
      "pct_major_procedure",
      "median_ne_dose_t0",
      "median_vp_dose_t0",
      "median_los_to_t0_h",
      "pct_imv_at_t0",
      "pct_crrt"
    ),
    value = c(
      nrow(dt),
      uniqueN(dt$patient_id),
      format(min(dt$t0_dttm, na.rm = TRUE), "%Y-%m"),
      format(max(dt$t0_dttm, na.rm = TRUE), "%Y-%m"),
      round(mean(dt$female_01, na.rm = TRUE) * 100, 1),
      round(median(dt$age, na.rm = TRUE), 1),
      round(mean(dt$white_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$hispanic_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$english_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$academic_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$peak_covid_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$outcome_group %in% c("noesc_dead", "esc_dead"), na.rm = TRUE) * 100, 1),
      round(mean(dt$outcome_group %in% c("esc_dead", "esc_alive"), na.rm = TRUE) * 100, 1),
      round(mean(!is.na(dt$svi_percentile)) * 100, 1),
      round(mean(!is.na(dt$adi_percentile)) * 100, 1),
      round(mean(dt$code_documented_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$major_procedure_01, na.rm = TRUE) * 100, 1),
      round(median(dt$ne_dose_t0, na.rm = TRUE), 3),
      round(median(dt$vp_dose_t0, na.rm = TRUE), 4),
      round(median(dt$los_to_t0_d, na.rm = TRUE) * 24, 1),
      round(mean(dt$imv_at_t0_01, na.rm = TRUE) * 100, 1),
      round(mean(dt$crrt_01, na.rm = TRUE) * 100, 1)
    ),
    site = site_lowercase
  )
}

qc_diagnostics = rbindlist(c(
  lapply(sort(unique(cohort_dt$hospital)), function(h) fn_qc_diagnostics(cohort_dt[hospital == h], h)),
  list(fn_qc_diagnostics(cohort_dt, "site_total"))
))

## hospital attribution --------------------------------------------------------

qc_hospital = cohort_dt[, .(
  n_encounters         = .N,
  n_first_hosp_differs = sum(as.character(hospital_id) != hospital, na.rm = TRUE),
  pct_academic         = round(mean(academic_01, na.rm = TRUE) * 100, 1)
), by = hospital][order(hospital)]
qc_hospital$site = site_lowercase

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

# QC outputs
fwrite(qc_missing,     file.path(output_dir, sprintf("qc_missing_%s.csv",     site_lowercase)))
fwrite(qc_ranges,      file.path(output_dir, sprintf("qc_ranges_%s.csv",      site_lowercase)))
fwrite(qc_flags,       file.path(output_dir, sprintf("qc_flags_%s.csv",       site_lowercase)))
fwrite(qc_categories,  file.path(output_dir, sprintf("qc_categories_%s.csv",  site_lowercase)))
fwrite(qc_diagnostics, file.path(output_dir, sprintf("qc_diagnostics_%s.csv", site_lowercase)))
fwrite(qc_hospital,    file.path(output_dir, sprintf("qc_hospital_%s.csv",    site_lowercase)))

message(sprintf("  ✅ Saved to: %s", output_dir))
message("    - all table1/flow files stratified by hospital (column: hospital)")
message("    - table1_continuous_*.csv  (n, sum, sumsq, percentiles by group)")
message("    - table1_binary_*.csv      (n, n_1 by group)")
message("    - table1_categorical_*.csv (cell counts by group)")
message("    - table1_timing_*.csv      (n_0, n_1, n_2 by group)")
message("    - table1_totals_*.csv      (group sizes)")
message("    - flow_diagram_*.csv       (cohort flow)")
message("    - qc_*.csv                 (QC diagnostics)")

message("\n== 03_table.R complete ==")
