# ==============================================================================
# 02_pool_qc.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Coordinating Center: QC review across sites
# ==============================================================================

# Requires: 00_pool_load.R

library(flextable)
library(officer)

message("\n== QC Review ==")

# ensure output directory
if (!dir.exists(here("output", "qc"))) {
  dir.create(here("output", "qc"), recursive = TRUE)
}

# ==============================================================================
# 1. EXCLUSION CASCADE
# ==============================================================================

message("\n-- Exclusion cascade --")

if (nrow(exclusion_raw) > 0) {
  
  # Pooled cascade
  exclusion_pooled = exclusion_raw[, .(
    n_remaining = sum(n_remaining, na.rm = TRUE),
    n_excluded  = sum(n_excluded, na.rm = TRUE)
  ), by = step]
  
  exclusion_pooled[, pct_of_start := round(n_remaining / exclusion_pooled[1, n_remaining] * 100, 1)]
  
  message("  Pooled exclusion cascade:")
  print(exclusion_pooled)
  
  # By-site comparison (% remaining at each step)
  exclusion_pct = exclusion_raw[, .(
    pct_remaining = round(n_remaining / n_remaining[1] * 100, 1)
  ), by = .(site, step)]
  
  exclusion_wide = dcast(exclusion_pct, step ~ site, value.var = "pct_remaining")
  
  # Check for outlier exclusion rates
  site_cols = intersect(ALLOWED_SITES, names(exclusion_wide))
  if (length(site_cols) > 1) {
    exclusion_wide[, mean_pct := rowMeans(.SD, na.rm = TRUE), .SDcols = site_cols]
    exclusion_wide[, sd_pct   := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = site_cols]
    
    # Flag outliers (>2 SD from mean)
    for (site in site_cols) {
      outliers = exclusion_wide[abs(get(site) - mean_pct) > 2 * sd_pct & sd_pct > 5]
      if (nrow(outliers) > 0) {
        message(sprintf("  ⚠️  %s has outlier exclusion rates at: %s",
                        toupper(site), paste(outliers$step, collapse = ", ")))
      }
    }
  }
  
  fwrite(exclusion_pooled, here("output", "qc", paste0("qc_exclusion_pooled_", today, ".csv")))
  fwrite(exclusion_wide,   here("output", "qc", paste0("qc_exclusion_bysite_", today, ".csv")))
  message("  Saved exclusion cascade summaries")
  
} else {
  message("  No exclusion cascade data available")
}

# ==============================================================================
# 2. PLAUSIBILITY FLAGS
# ==============================================================================

message("\n-- Plausibility flags --")

if (nrow(qc_flags_raw) > 0) {
  
  flags_wide = dcast(qc_flags_raw, check ~ site, value.var = "n_flagged", fill = 0)
  flags_wide[, total := rowSums(.SD, na.rm = TRUE), .SDcols = ALLOWED_SITES[ALLOWED_SITES %in% names(flags_wide)]]
  
  # Flag checks with any issues
  flags_with_issues = flags_wide[total > 0]
  
  if (nrow(flags_with_issues) > 0) {
    message("  ⚠️  ISSUES FOUND:")
    print(flags_with_issues)
  } else {
    message("  ✅ No plausibility issues")
  }
  
  fwrite(flags_wide, here("output", "qc", paste0("qc_flags_summary_", today, ".csv")))
  
} else {
  message("  No QC flags data available")
}

# ==============================================================================
# 2. MISSINGNESS COMPARISON
# ==============================================================================

message("\n-- Missingness by site --")

if (nrow(qc_missing_raw) > 0) {
  
  missing_wide = dcast(qc_missing_raw, variable ~ site, value.var = "pct_miss")
  
  # Calculate mean and flag outliers (>2 SD from mean)
  site_cols = intersect(ALLOWED_SITES, names(missing_wide))
  
  if (length(site_cols) > 1) {
    missing_wide[, mean_miss := rowMeans(.SD, na.rm = TRUE), .SDcols = site_cols]
    missing_wide[, sd_miss   := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = site_cols]
    
    # Check for high missingness (>20%) or outliers
    high_miss = missing_wide[mean_miss > 20]
    
    if (nrow(high_miss) > 0) {
      message("  ⚠️  Variables with >20% mean missingness:")
      print(high_miss[, .(variable, mean_miss = round(mean_miss, 1))])
    }
    
    # Check for site outliers
    for (site in site_cols) {
      outliers = missing_wide[abs(get(site) - mean_miss) > 2 * sd_miss & sd_miss > 5]
      if (nrow(outliers) > 0) {
        message(sprintf("  ⚠️  %s has outlier missingness for: %s", 
                        toupper(site), paste(outliers$variable, collapse = ", ")))
      }
    }
  }
  
  fwrite(missing_wide, here("output", "qc", paste0("qc_missing_comparison_", today, ".csv")))
  message("  Saved missingness comparison")
  
} else {
  message("  No missingness data available")
}

# ==============================================================================
# 3. CONTINUOUS VARIABLE RANGES
# ==============================================================================

message("\n-- Continuous variable ranges --")

if (nrow(qc_ranges_raw) > 0) {
  
  # Check for impossible values
  range_checks = list(
    age              = list(min = 18, max = 120, label = "Age"),
    ne_dose_t0       = list(min = 0,  max = 10,  label = "NE dose at T0"),
    vp_dose_t0       = list(min = 0,  max = 0.2, label = "VP dose at T0"),
    max_ne_equiv_48h = list(min = 0,  max = 20,  label = "Max NE equiv"),
    los_to_t0_d      = list(min = 0,  max = 365, label = "LOS to T0"),
    svi_percentile   = list(min = 0,  max = 1,   label = "SVI"),
    adi_percentile   = list(min = 1,  max = 100, label = "ADI")
  )
  
  for (var in names(range_checks)) {
    check = range_checks[[var]]
    var_data = qc_ranges_raw[variable == var]
    
    if (nrow(var_data) > 0) {
      issues = var_data[min < check$min | max > check$max]
      if (nrow(issues) > 0) {
        message(sprintf("  ⚠️  %s out of range [%s, %s] at: %s",
                        check$label, check$min, check$max,
                        paste(issues$site, collapse = ", ")))
        message(sprintf("      Min: %s, Max: %s", 
                        paste(round(issues$min, 2), collapse = "/"),
                        paste(round(issues$max, 2), collapse = "/")))
      }
    }
  }
  
  # Create summary table with medians by site
  ranges_summary = dcast(qc_ranges_raw, variable ~ site, value.var = "median")
  fwrite(ranges_summary, here("output", "qc", paste0("qc_ranges_medians_", today, ".csv")))
  
  # Full ranges table
  fwrite(qc_ranges_raw, here("output", "qc", paste0("qc_ranges_full_", today, ".csv")))
  message("  Saved range summaries")
  
} else {
  message("  No range data available")
}

# ==============================================================================
# 4. CATEGORICAL VALUE INVENTORY
# ==============================================================================

message("\n-- Categorical value check --")

if (nrow(qc_categories_raw) > 0) {
  
  # Check for unexpected categories
  expected_categories = list(
    outcome_group = c("dead_hospice", "escalated", "stable", "other"),
    code_status_t0 = c("full", "dnr_dni", "partial", "other"),
    race_category = c("white", "black or african american", "asian", 
                      "american indian or alaska native", 
                      "native hawaiian or other pacific islander",
                      "other", "unknown", "other/unknown"),
    ethnicity_category = c("hispanic", "non-hispanic", "non_hispanic", "unknown")
  )
  
  for (var in names(expected_categories)) {
    var_data = qc_categories_raw[variable == var]
    
    if (nrow(var_data) > 0) {
      observed = unique(tolower(var_data$category))
      expected = tolower(expected_categories[[var]])
      unexpected = setdiff(observed, c(expected, NA))
      
      if (length(unexpected) > 0) {
        message(sprintf("  ⚠️  Unexpected %s values: %s", var, paste(unexpected, collapse = ", ")))
        # Show which sites have these
        for (val in unexpected) {
          sites_with = unique(var_data[tolower(category) == val, site])
          message(sprintf("      '%s' at: %s", val, paste(sites_with, collapse = ", ")))
        }
      }
    }
  }
  
  # Category distribution by site
  cat_summary = qc_categories_raw[, .(
    n_sites = uniqueN(site),
    total_n = sum(n, na.rm = TRUE)
  ), by = .(variable, category)]
  
  fwrite(cat_summary, here("output", "qc", paste0("qc_categories_summary_", today, ".csv")))
  message("  Saved category inventory")
  
} else {
  message("  No category data available")
}

# ==============================================================================
# 5. CROSS-SITE CONSISTENCY
# ==============================================================================

message("\n-- Cross-site consistency --")

## Check that totals sum correctly ---------------------------------------------

if (nrow(totals_raw) > 0 && nrow(binary_raw) > 0) {
  
  # Total N from totals file
  total_by_site = totals_raw[, .(n_from_totals = sum(n_total)), by = site]
  
  # Total N from binary file (should match)
  binary_n = binary_raw[variable == "female_01", .(n_from_binary = sum(n)), by = site]
  
  consistency = merge(total_by_site, binary_n, by = "site", all = TRUE)
  consistency[, diff := n_from_totals - n_from_binary]
  consistency[, match := fifelse(abs(diff) <= 1, "OK", "MISMATCH")]
  
  mismatches = consistency[match == "MISMATCH"]
  if (nrow(mismatches) > 0) {
    message("  ⚠️  Total N mismatch between files:")
    print(mismatches)
  } else {
    message("  ✅ Totals consistent across files")
  }
}

## Check outcome group proportions ---------------------------------------------

if (nrow(totals_raw) > 0) {
  
  outcome_props = totals_raw[, .(
    pct = n_total / sum(n_total) * 100
  ), by = .(site, outcome_group)]
  
  outcome_wide = dcast(outcome_props, site ~ outcome_group, value.var = "pct")
  
  # Flag sites with unusual distributions
  if ("dead_hospice" %in% names(outcome_wide)) {
    mean_mort = mean(outcome_wide$dead_hospice, na.rm = TRUE)
    sd_mort   = sd(outcome_wide$dead_hospice, na.rm = TRUE)
    
    outliers = outcome_wide[abs(dead_hospice - mean_mort) > 2 * sd_mort]
    if (nrow(outliers) > 0) {
      message("  ⚠️  Outlier mortality rates:")
      print(outliers[, .(site, dead_hospice = round(dead_hospice, 1))])
    }
  }
  
  fwrite(outcome_wide, here("output", "qc", paste0("qc_outcome_proportions_", today, ".csv")))
  message("  Saved outcome proportions")
}

# ==============================================================================
# 6. SITE COMPARISON TABLE
# ==============================================================================

message("\n-- Creating site comparison table --")

if (nrow(continuous_raw) > 0) {
  
  # Key variables for comparison
  key_vars = c("age", "vw", "ne_dose_t0", "los_to_t0_d")
  
  site_comparison = continuous_raw[variable %in% key_vars, .(
    mean = sum(sum, na.rm = TRUE) / sum(n, na.rm = TRUE),
    n    = sum(n, na.rm = TRUE)
  ), by = .(site, variable)]
  
  site_comparison[, formatted := paste0(round(mean, 1), " (n=", format_n(n), ")")]
  
  site_comp_wide = dcast(site_comparison, variable ~ site, value.var = "formatted")
  
  fwrite(site_comp_wide, here("output", "qc", paste0("qc_site_comparison_", today, ".csv")))
  
  # Calculate SMDs between each site and pooled
  pooled_means = continuous_raw[variable %in% key_vars, .(
    pooled_mean = sum(sum, na.rm = TRUE) / sum(n, na.rm = TRUE),
    pooled_sd   = calculate_sd_from_sums(sum(sum), sum(sumsq), sum(n))
  ), by = variable]
  
  site_means = continuous_raw[variable %in% key_vars, .(
    site_mean = sum(sum, na.rm = TRUE) / sum(n, na.rm = TRUE),
    site_sd   = calculate_sd_from_sums(sum(sum), sum(sumsq), sum(n))
  ), by = .(site, variable)]
  
  site_smds = merge(site_means, pooled_means, by = "variable")
  site_smds[, smd := abs(site_mean - pooled_mean) / sqrt((site_sd^2 + pooled_sd^2) / 2)]
  
  smd_wide = dcast(site_smds, variable ~ site, value.var = "smd")
  
  # Flag large SMDs (>0.2)
  large_smds = site_smds[smd > 0.2]
  if (nrow(large_smds) > 0) {
    message("  ⚠️  Large site-vs-pooled SMDs (>0.2):")
    print(large_smds[, .(site, variable, smd = round(smd, 3))])
  }
  
  fwrite(smd_wide, here("output", "qc", paste0("qc_site_smds_", today, ".csv")))
}

# ==============================================================================
# 7. SITE DIAGNOSTICS DASHBOARD
# ==============================================================================

message("\n-- Site diagnostics dashboard --")

if (nrow(qc_diagnostics_raw) > 0) {
  
  # Pivot to wide format for easy comparison
  diag_wide = dcast(qc_diagnostics_raw, metric ~ site, value.var = "value")
  
  # Key metrics to flag
  numeric_metrics = c("pct_dead_hospice", "pct_escalated", "pct_svi_linked", 
                      "pct_adi_linked", "pct_code_documented", "pct_female",
                      "median_age", "pct_imv_at_t0", "pct_crrt")
  
  site_cols = intersect(ALLOWED_SITES, names(diag_wide))
  
  if (length(site_cols) > 1) {
    
    # Convert to numeric for flagging
    diag_numeric = diag_wide[metric %in% numeric_metrics]
    for (col in site_cols) {
      diag_numeric[[col]] = as.numeric(diag_numeric[[col]])
    }
    
    diag_numeric[, mean_val := rowMeans(.SD, na.rm = TRUE), .SDcols = site_cols]
    diag_numeric[, sd_val   := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = site_cols]
    
    # Flag outliers
    message("  Checking for outlier sites (>2 SD from mean):")
    
    for (met in numeric_metrics) {
      row = diag_numeric[metric == met]
      if (nrow(row) > 0 && !is.na(row$sd_val) && row$sd_val > 0) {
        for (site in site_cols) {
          val = row[[site]]
          if (!is.na(val) && abs(val - row$mean_val) > 2 * row$sd_val) {
            message(sprintf("    ⚠️  %s: %s = %.1f (mean = %.1f)", 
                            toupper(site), met, val, row$mean_val))
          }
        }
      }
    }
  }
  
  # Print summary table
  message("\n  Site diagnostics summary:")
  print(diag_wide, nrow = 20)
  
  fwrite(diag_wide, here("output", "qc", paste0("qc_diagnostics_summary_", today, ".csv")))
  message("  Saved diagnostics summary")
  
} else {
  message("  No diagnostics data available")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n== QC Summary ==")
message("  Output files saved to: output/qc/")
message("    - qc_exclusion_pooled_*.csv    (pooled exclusion cascade)")
message("    - qc_exclusion_bysite_*.csv    (exclusion % by site)")
message("    - qc_flags_summary_*.csv       (plausibility flags)")
message("    - qc_missing_comparison_*.csv  (missingness by site)")
message("    - qc_ranges_medians_*.csv      (continuous var medians)")
message("    - qc_ranges_full_*.csv         (continuous var distributions)")
message("    - qc_categories_summary_*.csv  (categorical value inventory)")
message("    - qc_outcome_proportions_*.csv (outcome rates by site)")
message("    - qc_site_comparison_*.csv     (key metrics by site)")
message("    - qc_site_smds_*.csv           (site vs pooled SMDs)")
message("    - qc_diagnostics_summary_*.csv (site diagnostics dashboard)")

message("\n== QC Review complete ==")
