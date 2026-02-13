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

# helper functions -------------------------------------------------------------

#' Format percentage: integer if >=10%, 1 decimal if <10%
format_pct = function(x) {
  fifelse(x >= 10, sprintf("%.0f", x), sprintf("%.1f", x))
}

#' Format count with percentage
format_n_pct = function(n, pct) {
  paste0(format_n(n), " (", format_pct(pct), "%)")
}

#' Format count with denominator and percentage (for variables with missingness)
#' Shows "n/N (pct%)" format
format_n_N_pct = function(n_1, n_total, pct) {
  paste0(format_n(n_1), "/", format_n(n_total), " (", format_pct(pct), "%)")
}

#' Format p-value
format_pval = function(p) {
  fcase(
    is.na(p),    "",
    p < 0.001,   "<0.001",
    p < 0.01,    sprintf("%.3f", p),
    default      = sprintf("%.2f", p)
  )
}

#' Chi-square test from contingency table counts
#' Returns p-value for 3-group comparison
calc_chisq_pval = function(counts_by_group) {
  # counts_by_group should have columns: outcome_group, n (or n_1 for binary)
  # Returns chi-square p-value
  tryCatch({
    tab = as.matrix(dcast(counts_by_group, . ~ outcome_group, value.var = "count")[, -1, with = FALSE])
    if (any(is.na(tab)) || any(tab < 0)) return(NA_real_)
    chisq.test(tab)$p.value
  }, error = function(e) NA_real_)
}

#' ANOVA p-value from pooled summary statistics (n, sum, sumsq)
#' Uses between-group and within-group variance from summary stats
calc_anova_pval_pooled = function(dt) {
  # dt should have: outcome_group, n, sum, sumsq
  # Only compute if we have 3 groups with valid data
  dt = dt[!is.na(n) & n > 0]
  if (nrow(dt) < 2) return(NA_real_)
  
  tryCatch({
    # Calculate group means
    dt[, mean := sum / n]
    
    # Grand mean (weighted)
    N_total = sum(dt$n)
    grand_mean = sum(dt$sum) / N_total
    
    # Between-group sum of squares
    SSB = sum(dt$n * (dt$mean - grand_mean)^2)
    
    # Within-group sum of squares: sumsq - n*mean^2 for each group
    dt[, ss_within := sumsq - n * mean^2]
    SSW = sum(dt$ss_within)
    
    # Degrees of freedom
    k = nrow(dt)  # number of groups
    df_between = k - 1
    df_within = N_total - k
    
    if (df_within <= 0 || SSW <= 0) return(NA_real_)
    
    # F statistic
    MSB = SSB / df_between
    MSW = SSW / df_within
    F_stat = MSB / MSW
    
    # p-value
    pf(F_stat, df_between, df_within, lower.tail = FALSE)
  }, error = function(e) NA_real_)
}

#' ANOVA p-value from summary statistics
#' n_vec, sum_vec, sumsq_vec are vectors (one per group)
calc_anova_pval = function(n_vec, sum_vec, sumsq_vec) {
  tryCatch({
    k = length(n_vec)
    if (k < 2) return(NA_real_)
    
    # Group means
    means = sum_vec / n_vec
    
    # Overall mean
    overall_mean = sum(sum_vec) / sum(n_vec)
    
    # Between-group sum of squares
    ss_between = sum(n_vec * (means - overall_mean)^2)
    
    # Within-group sum of squares: Σ(sumsq - sum²/n)
    ss_within = sum(sumsq_vec - sum_vec^2 / n_vec)
    
    # Degrees of freedom
    df_between = k - 1
    df_within = sum(n_vec) - k
    
    if (df_within <= 0 || ss_within <= 0) return(NA_real_)
    
    # F statistic
    f_stat = (ss_between / df_between) / (ss_within / df_within)
    
    # P-value
    pf(f_stat, df_between, df_within, lower.tail = FALSE)
  }, error = function(e) NA_real_)
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

# Calculate p-values (ANOVA across 3 groups)
cont_pvals = cont_pooled[, .(
  pval = calc_anova_pval(n, sum, sumsq)
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

# calculate percentage (of non-missing)
binary_pooled[, pct := (n_1 / n) * 100]

# format as "n/N (pct%)" showing denominator
binary_pooled[, formatted := format_n_N_pct(n_1, n, pct)]

# handle NA from masking
binary_pooled[is.na(n_1), formatted := "<5"]

# Calculate p-values (chi-square for binary: success vs failure across groups)
binary_pvals = binary_pooled[, {
  # Create 2xK contingency table: (successes, failures) x groups
  n_0 = n - n_1  # failures
  counts = rbind(n_1, n_0)
  pval = tryCatch({
    if (any(is.na(counts))) NA_real_ else chisq.test(counts)$p.value
  }, error = function(e) NA_real_)
  .(pval = pval)
}, by = variable]

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

# format as "n (%)" with smart rounding
cat_pooled[, formatted := format_n_pct(n, pct)]

# handle masked cells
cat_pooled[is.na(n), formatted := "<5"]

# Calculate p-values (chi-square: categories x outcome_groups)
cat_pvals = cat_pooled[, {
  # Create contingency table: categories (rows) x outcome_groups (cols)
  wide = dcast(.SD, category ~ outcome_group, value.var = "n", fill = 0)
  mat = as.matrix(wide[, -1, with = FALSE])
  pval = tryCatch({
    if (any(is.na(mat))) NA_real_ else chisq.test(mat)$p.value
  }, error = function(e) NA_real_)
  .(pval = pval)
}, by = variable]

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
    "white_01",
    "hispanic_01",
    "english_01",
    "academic_01",
    "peak_covid_01",
    "vw",
    "svi_percentile",
    "adi_percentile",
    "full_code_01",
    "major_procedure_01",
    "los_to_t0_d",
    "icu_los_to_t0_d",
    "imv_at_t0_01",
    "crrt_at_t0_01"
  ),
  label = c(
    "N",
    "Age, years, mean (SD)",
    "Female, n/N (%)",
    "White race, n/N (%)",
    "Hispanic ethnicity, n/N (%)",
    "English language, n/N (%)",
    "Academic hospital, n/N (%)",
    "Peak COVID period, n/N (%)",
    "Van Walraven Score, mean (SD)",
    "Social Vulnerability Index percentile, mean (SD)",
    "Area Deprivation Index percentile, mean (SD)",
    "Full code at study entry, n/N (%)",
    "Major procedure within 24h of study entry, n/N (%)",
    "Days in hospital before study entry, mean (SD)",
    "Days in ICU before study entry, mean (SD)",
    "Invasive mechanical ventilation at study entry, n/N (%)",
    "Continuous renal replacement therapy at study entry, n/N (%)"
  ),
  sort_order = 1:17,
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
    "b12_01",
    "imv_48h_01"
  ),
  label = c(
    "Died in hospital, n/N (%)",
    "Discharged to hospice, n/N (%)",
    "Days in hospital after study entry, mean (SD)",
    "Maximum NE-equivalent dose, mcg/kg/min, mean (SD)",
    "Received epinephrine infusion, n/N (%)",
    "Received phenylephrine infusion, n/N (%)",
    "Received dopamine infusion, n/N (%)",
    "Received angiotensin II infusion, n/N (%)",
    "Received methylene blue, n/N (%)",
    "Received hydroxocobalamin, n/N (%)",
    "Invasive mechanical ventilation within 48h, n/N (%)"
  ),
  sort_order = 101:111,
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

# Merge p-values
cont_wide = merge(cont_wide, cont_pvals, by = "variable", all.x = TRUE)

## binary variables wide -------------------------------------------------------

# Filter to only variables we want
binary_vars_keep = c(
  # Demographics/characteristics
  "female_01", "white_01", "hispanic_01", "english_01", 
  "academic_01", "peak_covid_01", "full_code_01", "major_procedure_01",
  # Co-interventions
  "imv_at_t0_01", "imv_48h_01", "crrt_at_t0_01",
  # Outcomes
  "dead_01", "hospice_01",
  "epi_01", "phenyl_01", "dopa_01", "a2_01", "mb_01", "b12_01"
)

binary_filtered = binary_pooled[variable %in% binary_vars_keep]

binary_wide = dcast(
  binary_filtered,
  variable ~ outcome_group,
  value.var = "formatted"
)
binary_wide[, category := NA_character_]

# Merge p-values
binary_wide = merge(binary_wide, binary_pvals, by = "variable", all.x = TRUE)

## combine all -----------------------------------------------------------------

# n_row doesn't have p-value
n_row[, pval := NA_real_]

all_data = rbindlist(list(
  n_row,
  cont_wide,
  binary_wide
), fill = TRUE)

## merge labels ----------------------------------------------------------------

all_data = merge(all_data, var_labels, by = "variable", all.x = TRUE)

## create display column -------------------------------------------------------

all_data[, display := fcase(
  variable == "n_total", "N",
  !is.na(label),         label,
  default = variable
)]

## sort ------------------------------------------------------------------------

# sort_order already comes from the var_labels merge above
all_data[is.na(sort_order), sort_order := 999]
setorder(all_data, sort_order)

## split into characteristics and outcomes -------------------------------------

outcome_cols = intersect(names(OUTCOME_LABELS), names(all_data))

# Format p-values
all_data[, pval_fmt := format_pval(pval)]

# Characteristics table (sort_order < 100)
char_data = all_data[sort_order < 100]
table1_char = char_data[, c("display", outcome_cols, "pval_fmt"), with = FALSE]
setnames(table1_char, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1_char, "display", "Characteristic")
setnames(table1_char, "pval_fmt", "p-value")
table1_char = table1_char[!is.na(Characteristic)]

# Outcomes table (sort_order >= 100)
outcome_data = all_data[sort_order >= 100]
table1_outcomes = outcome_data[, c("display", outcome_cols, "pval_fmt"), with = FALSE]
setnames(table1_outcomes, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1_outcomes, "display", "Outcome")
setnames(table1_outcomes, "pval_fmt", "p-value")
table1_outcomes = table1_outcomes[!is.na(Outcome)]

message("  Characteristics table: ", nrow(table1_char), " rows")
message("  Outcomes table: ", nrow(table1_outcomes), " rows")

# Combined for full table export
table1 = all_data[, c("display", "section", outcome_cols, "pval_fmt"), with = FALSE]
setnames(table1, outcome_cols, OUTCOME_LABELS[outcome_cols])
setnames(table1, "display", "Variable")
setnames(table1, "pval_fmt", "p-value")
table1 = table1[!is.na(Variable)]

message("  Combined Table 1: ", nrow(table1), " rows")

# ==============================================================================
# CREATE FLEXTABLES
# ==============================================================================

message("\n== Formatting tables ==")

# Shorter column headers for compact display
short_headers = c(
  "No Esc + Dead/Hospice" = "No Esc\nDead/Hosp",
  "Esc + Dead/Hospice"    = "Esc\nDead/Hosp",
  "Esc + Alive"           = "Esc\nAlive",
  "No Esc + Alive"        = "No Esc\nAlive"
)

## Characteristics table -------------------------------------------------------

ft_char = flextable(table1_char) |>
  set_header_labels(Characteristic = "", `p-value` = "p") |>
  fontsize(size = 7, part = "all") |>
  fontsize(size = 7, part = "header") |>
  padding(padding.top = 1, padding.bottom = 1, padding.left = 1, padding.right = 1, part = "all") |>
  align(j = 2:ncol(table1_char), align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  bold(i = 1, part = "body") |>
  hline(i = 1, border = fp_border(width = 0.5), part = "body") |>
  width(j = 1, width = 1.9) |>
  width(j = 2:(ncol(table1_char)-1), width = 1.0) |>
  width(j = ncol(table1_char), width = 0.4) |>
  set_table_properties(layout = "fixed", width = 1)

# Apply short headers
for (old_name in names(short_headers)) {
  if (old_name %in% names(table1_char)) {
    ft_char = set_header_labels(ft_char, values = setNames(list(short_headers[old_name]), old_name))
  }
}

## Outcomes table --------------------------------------------------------------

ft_outcomes = flextable(table1_outcomes) |>
  set_header_labels(Outcome = "", `p-value` = "p") |>
  fontsize(size = 7, part = "all") |>
  fontsize(size = 7, part = "header") |>
  padding(padding.top = 1, padding.bottom = 1, padding.left = 1, padding.right = 1, part = "all") |>
  align(j = 2:ncol(table1_outcomes), align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  width(j = 1, width = 1.9) |>
  width(j = 2:(ncol(table1_outcomes)-1), width = 1.0) |>
  width(j = ncol(table1_outcomes), width = 0.4) |>
  set_table_properties(layout = "fixed", width = 1)

# Apply short headers
for (old_name in names(short_headers)) {
  if (old_name %in% names(table1_outcomes)) {
    ft_outcomes = set_header_labels(ft_outcomes, values = setNames(list(short_headers[old_name]), old_name))
  }
}

## Combined table (for single document) ----------------------------------------

# Add a section separator row
separator_row = data.table(
  Variable = "Outcomes",
  section = "separator",
  `p-value` = ""
)
for (col in names(OUTCOME_LABELS)) {
  nm = OUTCOME_LABELS[col]
  if (nm %in% names(table1)) {
    separator_row[[nm]] = ""
  }
}

# Find where outcomes start and insert separator
char_rows = table1[section == "characteristics"]
outcome_rows = table1[section == "outcomes"]
table1_combined = rbindlist(list(char_rows, separator_row, outcome_rows), fill = TRUE)
table1_combined[, section := NULL]

ft_combined = flextable(table1_combined) |>
  set_header_labels(Variable = "", `p-value` = "p") |>
  fontsize(size = 7, part = "all") |>
  fontsize(size = 7, part = "header") |>
  padding(padding.top = 1, padding.bottom = 1, padding.left = 1, padding.right = 1, part = "all") |>
  align(j = 2:ncol(table1_combined), align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  bold(i = 1, part = "body") |>
  hline(i = 1, border = fp_border(width = 0.5), part = "body") |>
  width(j = 1, width = 1.9) |>
  width(j = 2:(ncol(table1_combined)-1), width = 1.0) |>
  width(j = ncol(table1_combined), width = 0.4) |>
  set_table_properties(layout = "fixed", width = 1)

# Apply short headers
for (old_name in names(short_headers)) {
  if (old_name %in% names(table1_combined)) {
    ft_combined = set_header_labels(ft_combined, values = setNames(list(short_headers[old_name]), old_name))
  }
}

# Bold the Outcomes section separator
outcomes_row = which(table1_combined$Variable == "Outcomes")
if (length(outcomes_row) > 0) {
  ft_combined = bold(ft_combined, i = outcomes_row, j = 1)
}

# Add line before outcomes section
outcomes_row = which(table1_combined$Variable == "Outcomes")
if (length(outcomes_row) > 0) {
  ft_combined = hline(ft_combined, i = outcomes_row - 1, border = fp_border(width = 0.5), part = "body")
}

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n== Saving outputs ==")

# Page properties for landscape orientation
page_landscape = prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
)

# Page properties for portrait (narrower margins)
page_portrait = prop_section(
  page_size = page_size(orient = "portrait"),
  page_margins = page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
)

# Characteristics table (portrait)
save_as_docx(ft_char, path = here("output", "tables", paste0("table1_characteristics_", today, ".docx")),
             pr_section = page_portrait)
message("  Saved: table1_characteristics_", today, ".docx")

# Outcomes table (portrait)
save_as_docx(ft_outcomes, path = here("output", "tables", paste0("table1_outcomes_", today, ".docx")),
             pr_section = page_portrait)
message("  Saved: table1_outcomes_", today, ".docx")

# Combined table (portrait with narrow margins, or landscape if needed)
save_as_docx(ft_combined, path = here("output", "tables", paste0("table1_combined_", today, ".docx")),
             pr_section = page_portrait)
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