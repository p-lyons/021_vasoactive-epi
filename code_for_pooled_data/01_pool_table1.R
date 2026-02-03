# ==============================================================================
# 01_pool_table1.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Coordinating Center: Create pooled Table 1 by outcome group
# ==============================================================================

# Requires: 00_pool_load.R

library(flextable)
library(officer)

# ensure output directory exists
if (!dir.exists(here("output", "tables"))) {
  dir.create(here("output", "tables"), recursive = TRUE)
}

# ==============================================================================
# POOL CONTINUOUS VARIABLES
# ==============================================================================

message("\n== Pooling continuous variables ==")

cont_pooled = continuous_raw[, .(
  n      = sum(n,      na.rm = TRUE),
  n_miss = sum(n_miss, na.rm = TRUE),
  sum    = sum(sum,    na.rm = TRUE),
  sumsq  = sum(sumsq,  na.rm = TRUE),
  min    = min(min,    na.rm = TRUE),
  max    = max(max,    na.rm = TRUE)
), by = .(outcome_group, variable)]

# calculate mean and SD
cont_pooled[, `:=`(
  mean = sum / n,
  sd   = calculate_sd_from_sums(sum, sumsq, n)
)]

# format as "mean (SD)"
cont_pooled[, formatted := paste0(
  round(mean, 1), " (", round(sd, 1), ")"
)]

# handle Inf/-Inf from empty groups
cont_pooled[is.infinite(min), min := NA_real_]
cont_pooled[is.infinite(max), max := NA_real_]
cont_pooled[is.nan(mean),     formatted := "—"]

message("  Pooled ", uniqueN(cont_pooled$variable), " continuous variables")

# ==============================================================================
# POOL BINARY VARIABLES
# ==============================================================================

message("  Pooling binary variables...")

binary_pooled = binary_raw[, .(
  n   = sum(n,   na.rm = TRUE),
  n_1 = sum(n_1, na.rm = TRUE)
), by = .(outcome_group, variable)]

# calculate percentage
binary_pooled[, pct := (n_1 / n) * 100]

# format as "n (%)"
binary_pooled[, formatted := paste0(format_n(n_1), " (", round(pct, 1), "%)")]

# handle NA from masking
binary_pooled[is.na(n_1), formatted := "<5"]

message("  Pooled ", uniqueN(binary_pooled$variable), " binary variables")

# ==============================================================================
# POOL CATEGORICAL VARIABLES
# ==============================================================================

message("  Pooling categorical variables...")

cat_pooled = categorical_raw[, .(
  n = sum(n, na.rm = TRUE)
), by = .(outcome_group, variable, category)]

# get group totals for percentages
cat_pooled = merge(
  cat_pooled,
  COHORT_N[, .(outcome_group, N = n_total)],
  by = "outcome_group"
)

cat_pooled[, pct := (n / N) * 100]

# format as "n (%)"
cat_pooled[, formatted := paste0(format_n(n), " (", round(pct, 1), "%)")]

# handle masked cells
cat_pooled[is.na(n), formatted := "<5"]

message("  Pooled ", uniqueN(cat_pooled$variable), " categorical variables")

# ==============================================================================
# POOL TIMING VARIABLES
# ==============================================================================

message("  Pooling timing variables...")

timing_pooled = timing_raw[, .(
  n   = sum(n,   na.rm = TRUE),
  n_0 = sum(n_0, na.rm = TRUE),
  n_1 = sum(n_1, na.rm = TRUE),
  n_2 = sum(n_2, na.rm = TRUE)
), by = .(outcome_group, variable)]

# get group totals
timing_pooled = merge(
  timing_pooled,
  COHORT_N[, .(outcome_group, N = n_total)],
  by = "outcome_group"
)

# format each level
timing_pooled[, `:=`(
  formatted_0 = paste0(format_n(n_0), " (", round(n_0/N*100, 1), "%)"),
  formatted_1 = paste0(format_n(n_1), " (", round(n_1/N*100, 1), "%)"),
  formatted_2 = paste0(format_n(n_2), " (", round(n_2/N*100, 1), "%)")
)]

message("  Pooled ", uniqueN(timing_pooled$variable), " timing variables")

# ==============================================================================
# BUILD TABLE 1
# ==============================================================================

message("\n== Building Table 1 ==")

## variable labels -------------------------------------------------------------

var_labels = data.table(
  variable = c(
    "n_total",
    "age",
    "female_01",
    "race_category",
    "ethnicity_category",
    "vw",
    "code_status_t0",
    "code_documented_01",
    "los_to_t0_d",
    "icu_los_to_t0_d",
    "ne_dose_t0",
    "vp_dose_t0",
    "max_ne_equiv_48h",
    "epi_01",
    "phenyl_01",
    "dopa_01",
    "a2_01",
    "mb_01",
    "b12_01",
    "imv_timing_group",
    "imv_at_t0_01",
    "imv_48h_01",
    "crrt_timing_group",
    "crrt_01",
    "svi_percentile",
    "adi_percentile",
    "los_hosp_d",
    "dead_01",
    "hospice_01"
  ),
  label = c(
    "N",
    "Age, years",
    "Female sex",
    "Race",
    "Ethnicity",
    "Van Walraven score",
    "Code status at T0",
    "Code status documented",
    "LOS to T0, days",
    "ICU LOS to T0, days",
    "Norepinephrine dose at T0, mcg/kg/min",
    "Vasopressin dose at T0, units/min",
    "Max NE equivalent, mcg/kg/min",
    "Epinephrine",
    "Phenylephrine",
    "Dopamine",
    "Angiotensin II",
    "Methylene blue",
    "Hydroxocobalamin (B12)",
    "IMV timing",
    "IMV at T0",
    "IMV within 48h",
    "CRRT timing",
    "CRRT",
    "SVI percentile",
    "ADI percentile",
    "Hospital LOS, days",
    "Died in hospital",
    "Discharged to hospice"
  ),
  section = c(
    "header",
    "demographics",
    "demographics",
    "demographics",
    "demographics",
    "comorbidities",
    "code_status",
    "code_status",
    "timing",
    "timing",
    "vasopressors",
    "vasopressors",
    "vasopressors",
    "escalation",
    "escalation",
    "escalation",
    "escalation",
    "escalation",
    "escalation",
    "organ_support",
    "organ_support",
    "organ_support",
    "organ_support",
    "organ_support",
    "sdoh",
    "sdoh",
    "outcomes",
    "outcomes",
    "outcomes"
  ),
  sort_order = 1:29
)

## N row -----------------------------------------------------------------------

n_row = dcast(COHORT_N, . ~ outcome_group, value.var = "n_total")
n_row[, `:=`(variable = "n_total", category = NA_character_)]
n_row[, . := NULL]

# format with commas
for (col in names(OUTCOME_LABELS)) {
  if (col %in% names(n_row)) {
    n_row[[col]] = format_n(n_row[[col]])
  }
}

## continuous variables wide ---------------------------------------------------

cont_wide = dcast(
  cont_pooled, 
  variable ~ outcome_group, 
  value.var = "formatted"
)
cont_wide[, category := NA_character_]

## binary variables wide -------------------------------------------------------

binary_wide = dcast(
  binary_pooled,
  variable ~ outcome_group,
  value.var = "formatted"
)
binary_wide[, category := NA_character_]

## categorical variables wide --------------------------------------------------

cat_wide = dcast(
  cat_pooled,
  variable + category ~ outcome_group,
  value.var = "formatted"
)

## timing variables - reshape to categorical format ----------------------------

# IMV timing
imv_timing = timing_pooled[variable == "imv_timing_group"]

imv_timing_long = rbindlist(list(
  imv_timing[, .(outcome_group, variable = "imv_timing_group", category = "None",       formatted = formatted_0)],
  imv_timing[, .(outcome_group, variable = "imv_timing_group", category = "Before T0",  formatted = formatted_1)],
  imv_timing[, .(outcome_group, variable = "imv_timing_group", category = "T0 to 48h",  formatted = formatted_2)]
))

imv_timing_wide = dcast(imv_timing_long, variable + category ~ outcome_group, value.var = "formatted")

# CRRT timing
crrt_timing = timing_pooled[variable == "crrt_timing_group"]

crrt_timing_long = rbindlist(list(
  crrt_timing[, .(outcome_group, variable = "crrt_timing_group", category = "None",       formatted = formatted_0)],
  crrt_timing[, .(outcome_group, variable = "crrt_timing_group", category = "Before T0",  formatted = formatted_1)],
  crrt_timing[, .(outcome_group, variable = "crrt_timing_group", category = "T0 to 48h",  formatted = formatted_2)]
))

crrt_timing_wide = dcast(crrt_timing_long, variable + category ~ outcome_group, value.var = "formatted")

## combine all -----------------------------------------------------------------

all_data = rbindlist(list(
  n_row,
  cont_wide,
  binary_wide,
  cat_wide,
  imv_timing_wide,
  crrt_timing_wide
), fill = TRUE)

## merge labels ----------------------------------------------------------------

all_data = merge(all_data, var_labels, by = "variable", all.x = TRUE)

## create display column -------------------------------------------------------

# Variables that are categorical (need header row + indented categories)
categorical_vars = c("race_category", "ethnicity_category", "code_status_t0", 
                     "imv_timing_group", "crrt_timing_group")

all_data[, display := fcase(
  variable == "n_total",                               "N",
  variable %in% categorical_vars & is.na(category),    label,
  variable %in% categorical_vars & !is.na(category),   paste0("    ", category),
  !is.na(label),                                       label,
  default = variable
)]

## add header rows for categorical variables -----------------------------------

header_rows = data.table(
  variable = categorical_vars,
  category = NA_character_
)
header_rows = merge(header_rows, var_labels, by = "variable")

# Check which headers are missing
existing_headers = all_data[variable %in% categorical_vars & is.na(category)]
missing_headers = setdiff(categorical_vars, existing_headers$variable)

if (length(missing_headers) > 0) {
  new_headers = header_rows[variable %in% missing_headers]
  new_headers[, display := label]
  all_data = rbindlist(list(all_data, new_headers), fill = TRUE)
}

## sort ------------------------------------------------------------------------

# sort_order already comes from the var_labels merge above
all_data[is.na(sort_order), sort_order := 999]

# For categorical vars, sort categories after header
all_data[, cat_order := fifelse(is.na(category), 0, 1)]
setorder(all_data, sort_order, cat_order, category, na.last = FALSE)

## final table -----------------------------------------------------------------

# Select and rename columns
outcome_cols = intersect(names(OUTCOME_LABELS), names(all_data))

table1 = all_data[, c("display", outcome_cols), with = FALSE]
setnames(table1, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1, "display", "Variable")

# Remove rows with all NA outcome columns (shouldn't happen but safety check)
table1 = table1[!is.na(Variable)]

message("  Table 1 created with ", nrow(table1), " rows")

# ==============================================================================
# CREATE FLEXTABLE
# ==============================================================================

message("\n== Formatting table ==")

ft1 = flextable(table1) |>
  set_header_labels(Variable = "") |>
  autofit() |>
  align(j = 2:ncol(table1), align = "center", part = "all") |>
  bold(i = 1, part = "body") |>  # N row
  hline(i = 1, border = fp_border(width = 1), part = "body") |>
  fontsize(size = 10, part = "all") |>
  padding(padding = 2, part = "all")

# Add section headers (bold the category headers)
cat_header_rows = which(table1$Variable %in% var_labels[variable %in% categorical_vars, label])
if (length(cat_header_rows) > 0) {
  ft1 = bold(ft1, i = cat_header_rows, j = 1)
}

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n== Saving outputs ==")

# Word document
save_as_docx(ft1, path = here("output", "tables", paste0("table1_characteristics_", today, ".docx")))
message("  Saved: table1_characteristics_", today, ".docx")

# CSV for further analysis
fwrite(all_data, here("output", "tables", paste0("table1_data_", today, ".csv")))
message("  Saved: table1_data_", today, ".csv")

# ==============================================================================
# POOLED FLOW DIAGRAM
# ==============================================================================

message("\n== Creating pooled flow diagram ==")

flow_pooled = flow_raw[, .(
  n = sum(n, na.rm = TRUE)
), by = step]

setorder(flow_pooled, step)

fwrite(flow_pooled, here("output", "tables", paste0("flow_diagram_pooled_", today, ".csv")))
message("  Saved: flow_diagram_pooled_", today, ".csv")

print(flow_pooled)

message("\n== Table 1 complete ==")
