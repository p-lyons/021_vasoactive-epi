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

## p-values for continuous variables (Welch's ANOVA approximation) -------------

calculate_welch_anova_p = function(means, sds, ns) {
  # Welch's ANOVA for unequal variances
  # Returns p-value for test of equal means across groups
  k = length(means)
  if (k < 2 || any(ns < 2) || any(is.na(sds)) || any(sds == 0)) return(NA_real_)
  
  weights = ns / (sds^2)
  grand_mean = sum(weights * means) / sum(weights)
  
  # F statistic numerator
  f_num = sum(weights * (means - grand_mean)^2) / (k - 1)
  
  # F statistic denominator (with Welch correction)
  lambda = sum((1 - weights/sum(weights))^2 / (ns - 1))
  f_den = 1 + (2 * (k - 2) / (k^2 - 1)) * lambda
  
  f_stat = f_num / f_den
  
  # Degrees of freedom

  df1 = k - 1
  df2 = (k^2 - 1) / (3 * lambda)
  
  p_val = pf(f_stat, df1, df2, lower.tail = FALSE)
  return(p_val)
}

cont_pvalues = cont_pooled[, .(
  p_value = calculate_welch_anova_p(mean, sd, n)
), by = variable]

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

## p-values for binary variables (chi-square test) -----------------------------

calculate_chisq_p_binary = function(n_1_vec, n_vec) {
  # Chi-square test for binary variable across groups
  n_0_vec = n_vec - n_1_vec
  if (any(is.na(n_1_vec)) || any(n_vec < 1)) return(NA_real_)
  
  obs = matrix(c(n_1_vec, n_0_vec), nrow = 2, byrow = TRUE)
  
  # Check for valid contingency table
  if (any(colSums(obs) == 0) || any(rowSums(obs) == 0)) return(NA_real_)
  
  tryCatch({
    chisq.test(obs)$p.value
  }, error = function(e) NA_real_)
}

binary_pvalues = binary_pooled[, .(
  p_value = calculate_chisq_p_binary(n_1, n)
), by = variable]

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

## p-values for categorical variables (chi-square test) ------------------------

calculate_chisq_p_cat = function(dt) {
  # Chi-square test for categorical variable across outcome groups
  # dt should have columns: outcome_group, category, n
  wide = dcast(dt, category ~ outcome_group, value.var = "n", fill = 0)
  
  # Remove category column to get matrix
  mat = as.matrix(wide[, -1, with = FALSE])
  
  # Check for valid contingency table
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  if (any(colSums(mat) == 0) || any(rowSums(mat) == 0)) return(NA_real_)
  
  tryCatch({
    chisq.test(mat)$p.value
  }, error = function(e) NA_real_)
}

cat_pvalues = cat_pooled[, .(
  p_value = calculate_chisq_p_cat(.SD)
), by = variable, .SDcols = c("outcome_group", "category", "n")]

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

# CHARACTERISTICS (top half)
char_labels = data.table(
  variable = c(
    "n_total",
    "age",
    "female_01",
    "race_category",
    "ethnicity_category",
    "vw",
    "svi_percentile",
    "adi_percentile",
    "code_status_t0",
    "los_to_t0_d",
    "icu_los_to_t0_d",
    "imv_at_t0_01",
    "crrt_at_t0_01"
  ),
  label = c(
    "N",
    "Age, years, mean (SD)",
    "Female, n (%)",
    "Race, n (%)",
    "Ethnicity, n (%)",
    "Van Walraven comorbidity score, mean (SD)",
    "Social Vulnerability Index, mean (SD)",
    "Area Deprivation Index, mean (SD)",
    "Code status at study entry, n (%)",
    "Days in hospital before study entry, mean (SD)",
    "Days in ICU before study entry, mean (SD)",
    "Invasive mechanical ventilation at study entry, n (%)",
    "Continuous renal replacement therapy at study entry, n (%)"
  ),
  sort_order = 1:13,
  section = "characteristics"
)

# OUTCOMES (bottom half)
outcome_labels = data.table(
  variable = c(
    "dead_01",
    "hospice_01",
    "los_from_t0_d",
    "max_ne_equiv_48h",
    "epi_01",
    "phenyl_01",
    "dopa_01",
    "a2_01",
    "mb_01",
    "b12_01"
  ),
  label = c(
    "Died in hospital, n (%)",
    "Discharged to hospice, n (%)",
    "Length of stay from study entry, days, mean (SD)",
    "Maximum norepinephrine equivalent dose, mcg/kg/min, mean (SD)",
    "Received epinephrine infusion, n (%)",
    "Received phenylephrine infusion, n (%)",
    "Received dopamine infusion, n (%)",
    "Received angiotensin II infusion, n (%)",
    "Received methylene blue, n (%)",
    "Received hydroxocobalamin, n (%)"
  ),
  sort_order = 101:110,
  section = "outcomes"
)

var_labels = rbindlist(list(char_labels, outcome_labels))

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

# Filter to only variables we want
cont_vars_keep = c("age", "vw", "svi_percentile", "adi_percentile", 
                   "los_to_t0_d", "icu_los_to_t0_d", "los_from_t0_d", "max_ne_equiv_48h")

cont_filtered = cont_pooled[variable %in% cont_vars_keep]

# Special formatting for vw (round to integer)
cont_filtered[variable == "vw", formatted := paste0(round(mean, 0), " (", round(sd, 0), ")")]

cont_wide = dcast(
  cont_filtered, 
  variable ~ outcome_group, 
  value.var = "formatted"
)
cont_wide[, category := NA_character_]

# Add p-values for continuous variables
cont_wide = merge(cont_wide, cont_pvalues, by = "variable", all.x = TRUE)

## binary variables wide -------------------------------------------------------

# Filter to only variables we want
binary_vars_keep = c("female_01", "imv_at_t0_01", "crrt_at_t0_01", 
                     "dead_01", "hospice_01",
                     "epi_01", "phenyl_01", "dopa_01", "a2_01", "mb_01", "b12_01")

binary_filtered = binary_pooled[variable %in% binary_vars_keep]

binary_wide = dcast(
  binary_filtered,
  variable ~ outcome_group,
  value.var = "formatted"
)
binary_wide[, category := NA_character_]

# Add p-values for binary variables
binary_wide = merge(binary_wide, binary_pvalues, by = "variable", all.x = TRUE)

## categorical variables wide --------------------------------------------------

# Only race, ethnicity, code_status_t0
cat_vars_keep = c("race_category", "ethnicity_category", "code_status_t0")

cat_filtered = cat_pooled[variable %in% cat_vars_keep]

# Ensure category is character (not factor)
cat_filtered[, category := as.character(category)]

# Clean up code status labels
cat_filtered[variable == "code_status_t0", category := fcase(
  category == "full",     "Full code",
  category == "dnr_dni",  "DNR/DNI",
  category == "partial",  "Partial limitations",
  category == "other",    "Other",
  default = category
)]

# Clean up race labels (title case)
cat_filtered[variable == "race_category", category := fcase(
  tolower(category) == "white",                                     "White",
  tolower(category) == "black or african american",                 "Black or African American",
  tolower(category) == "asian",                                     "Asian",
  tolower(category) == "american indian or alaska native",          "American Indian or Alaska Native",
  tolower(category) == "native hawaiian or other pacific islander", "Native Hawaiian or Other Pacific Islander",
  tolower(category) %in% c("other", "unknown", "other/unknown"),    "Other/Unknown",
  default = category
)]

# Clean up ethnicity labels
cat_filtered[variable == "ethnicity_category", category := fcase(
  tolower(category) == "hispanic",                        "Hispanic",
  tolower(category) %in% c("non-hispanic", "non_hispanic"), "Non-Hispanic",
  tolower(category) == "unknown",                         "Unknown",
  default = category
)]

# Re-aggregate after label cleaning (some categories may have been merged, e.g. other + unknown -> Other/Unknown)
cat_filtered = cat_filtered[, .(
  n   = sum(n, na.rm = TRUE),
  N   = first(N),
  pct = sum(n, na.rm = TRUE) / first(N) * 100
), by = .(outcome_group, variable, category)]

# Re-format after aggregation
cat_filtered[, formatted := paste0(format_n(n), " (", round(pct, 1), "%)")]
cat_filtered[is.na(n) | n == 0, formatted := "0 (0%)"]

cat_wide = dcast(
  cat_filtered,
  variable + category ~ outcome_group,
  value.var = "formatted"
)

# Ensure category stays character after dcast
cat_wide[, category := as.character(category)]

# Add p-values for categorical variables (only on header row, NA for category rows)
cat_wide = merge(cat_wide, cat_pvalues, by = "variable", all.x = TRUE)
cat_wide[!is.na(category), p_value := NA_real_]

## combine all -----------------------------------------------------------------

all_data = rbindlist(list(
  n_row,
  cont_wide,
  binary_wide,
  cat_wide
), fill = TRUE)

## merge labels ----------------------------------------------------------------

all_data = merge(all_data, var_labels, by = "variable", all.x = TRUE)

## format p-values -------------------------------------------------------------

format_pvalue = function(p) {
  ifelse(is.na(p), "", 
         ifelse(p < 0.001, "<0.001",
                ifelse(p < 0.01, sprintf("%.3f", p),
                       sprintf("%.2f", p))))
}

all_data[, p_formatted := format_pvalue(p_value)]

## create display column -------------------------------------------------------

categorical_vars = c("race_category", "ethnicity_category", "code_status_t0")

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

# Add p-values to header rows for categorical variables
header_rows = merge(header_rows, cat_pvalues, by = "variable", all.x = TRUE)
header_rows[, p_formatted := format_pvalue(p_value)]

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

## split into characteristics and outcomes -------------------------------------

outcome_cols = intersect(names(OUTCOME_LABELS), names(all_data))

# Characteristics table (sort_order < 100)
char_data = all_data[sort_order < 100]
table1_char = char_data[, c("display", outcome_cols, "p_formatted"), with = FALSE]
setnames(table1_char, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1_char, "display", "Characteristic")
setnames(table1_char, "p_formatted", "P-value")
table1_char = table1_char[!is.na(Characteristic)]

# Outcomes table (sort_order >= 100)
outcome_data = all_data[sort_order >= 100]
table1_outcomes = outcome_data[, c("display", outcome_cols, "p_formatted"), with = FALSE]
setnames(table1_outcomes, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1_outcomes, "display", "Outcome")
setnames(table1_outcomes, "p_formatted", "P-value")
table1_outcomes = table1_outcomes[!is.na(Outcome)]

message("  Characteristics table: ", nrow(table1_char), " rows")
message("  Outcomes table: ", nrow(table1_outcomes), " rows")

# Combined for full table export
table1 = all_data[, c("display", "section", outcome_cols, "p_formatted"), with = FALSE]
setnames(table1, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1, "display", "Variable")
setnames(table1, "p_formatted", "P-value")
table1 = table1[!is.na(Variable)]

message("  Combined Table 1: ", nrow(table1), " rows")

# ==============================================================================
# CREATE FLEXTABLES
# ==============================================================================

message("\n== Formatting tables ==")

## Characteristics table -------------------------------------------------------

ft_char = flextable(table1_char) |>
  set_header_labels(Characteristic = "") |>
  width(j = 1, width = 2.5) |>
  width(j = 2:4, width = 1.3) |>
  width(j = 5, width = 0.7) |>
  align(j = 2:ncol(table1_char), align = "center", part = "all") |>
  bold(i = 1, part = "body") |>  # N row
  hline(i = 1, border = fp_border(width = 1), part = "body") |>
  fontsize(size = 10, part = "all") |>
  padding(padding = 2, part = "all")

# Bold the category headers
cat_header_rows_char = which(table1_char$Characteristic %in% var_labels[variable %in% categorical_vars, label])
if (length(cat_header_rows_char) > 0) {
  ft_char = bold(ft_char, i = cat_header_rows_char, j = 1)
}

## Outcomes table --------------------------------------------------------------

ft_outcomes = flextable(table1_outcomes) |>
  set_header_labels(Outcome = "") |>
  width(j = 1, width = 2.5) |>
  width(j = 2:4, width = 1.3) |>
  width(j = 5, width = 0.7) |>
  align(j = 2:ncol(table1_outcomes), align = "center", part = "all") |>
  fontsize(size = 10, part = "all") |>
  padding(padding = 2, part = "all")

## Combined table (for single document) ----------------------------------------

# Add a section separator row
separator_row = data.table(
  Variable = "Outcomes",
  section = "separator"
)
for (col in names(OUTCOME_LABELS)) {
  nm = OUTCOME_LABELS[col]
  if (nm %in% names(table1)) {
    separator_row[[nm]] = ""
  }
}
separator_row[["P-value"]] = ""

# Find where outcomes start and insert separator
char_rows = table1[section == "characteristics"]
outcome_rows = table1[section == "outcomes"]
table1_combined = rbindlist(list(char_rows, separator_row, outcome_rows), fill = TRUE)
table1_combined[, section := NULL]

ft_combined = flextable(table1_combined) |>
  set_header_labels(Variable = "") |>
  width(j = 1, width = 2.5) |>
  width(j = 2:4, width = 1.3) |>
  width(j = 5, width = 0.7) |>
  align(j = 2:ncol(table1_combined), align = "center", part = "all") |>
  bold(i = 1, part = "body") |>  # N row
  hline(i = 1, border = fp_border(width = 1), part = "body") |>
  fontsize(size = 10, part = "all") |>
  padding(padding = 2, part = "all")

# Bold the category headers and section separator
cat_labels = var_labels[variable %in% categorical_vars, label]
bold_rows = which(table1_combined$Variable %in% c(cat_labels, "Outcomes"))
if (length(bold_rows) > 0) {
  ft_combined = bold(ft_combined, i = bold_rows, j = 1)
}

# Add line before outcomes section
outcomes_row = which(table1_combined$Variable == "Outcomes")
if (length(outcomes_row) > 0) {
  ft_combined = hline(ft_combined, i = outcomes_row - 1, border = fp_border(width = 1), part = "body")
}

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n== Saving outputs ==")

# Create landscape section properties
landscape_props = prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
)

# Characteristics table (landscape)
save_as_docx(ft_char, path = here("output", "tables", paste0("table1_characteristics_", today, ".docx")),
             pr_section = landscape_props)
message("  Saved: table1_characteristics_", today, ".docx")

# Outcomes table (landscape)
save_as_docx(ft_outcomes, path = here("output", "tables", paste0("table1_outcomes_", today, ".docx")),
             pr_section = landscape_props)
message("  Saved: table1_outcomes_", today, ".docx")

# Combined table (landscape)
save_as_docx(ft_combined, path = here("output", "tables", paste0("table1_combined_", today, ".docx")),
             pr_section = landscape_props)
message("  Saved: table1_combined_", today, ".docx")

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
