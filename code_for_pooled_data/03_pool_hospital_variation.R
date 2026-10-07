# ==============================================================================
# 03_pool_hospital_variation.R
# Vasopressor Escalation in Refractory Distributive Shock - CLIF Consortium
# Coordinating site: hospital-level variation in escalation and early death
#
# Unit = hospital at T0 (hosp_key = site + hospital). Uses hospital rows only.
#   - Raw ranges: hospitals with >= MIN_HOSP_N encounters.
#   - Random-intercept logistic models: ALL hospitals (shrinkage handles small
#     hospitals). Unadjusted: site files carry no patient-level covariates.
#   - Rescue agents a site does not capture are left out (see 00_pool_load.R).
# ==============================================================================

# Requires: 00_pool_load.R

library(lme4)

MIN_HOSP_N = 30L

dir_out = here("output", "hospital_variation")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

message("\n== Hospital-level variation ==")

# ==============================================================================
# STEP 1: Hospital-level counts
# ==============================================================================

## encounters, any escalation, 48h death/hospice (from totals) -----------------

hosp_totals = totals_hosp[, .(
  n_enc       = sum(n_total),
  n_escalated = sum(n_total[outcome_group %chin% c("esc_dead", "esc_alive")]),
  n_dead_48h  = sum(n_total[outcome_group %chin% c("esc_dead", "noesc_dead")])
), by = .(hosp_key, site, hospital)]

hosp_totals = hosp_totals[hospital != "unknown"]

## hospital type ---------------------------------------------------------------

hosp_type = binary_hosp[
  variable == "academic_01",
  .(pct_academic = round(sum(n_1, na.rm = TRUE) / sum(n) * 100, 1)),
  by = hosp_key
]

hosp_totals = merge(
  hosp_totals,
  hosp_type,
  by    = "hosp_key",
  all.x = TRUE
)

## long format: one row per hospital x outcome ---------------------------------

agent_vars = c(
  "epi_01",
  "phenyl_01",
  "dopa_01",
  "a2_01",
  "mb_01",
  "b12_01"
)

hosp_agents = binary_hosp[
  variable %chin% agent_vars,
  .(n = sum(n), n_1 = sum(n_1)),
  by = .(hosp_key, outcome = variable)
]

hosp_any = hosp_totals[, .(hosp_key, outcome = "any_escalation",  n = n_enc, n_1 = n_escalated)]
hosp_dead = hosp_totals[, .(hosp_key, outcome = "dead_hospice_48h", n = n_enc, n_1 = n_dead_48h)]

hosp_long = rbindlist(
  list(
    hosp_any,
    hosp_agents,
    hosp_dead
  ),
  use.names = TRUE
)

hosp_long = merge(
  hosp_long,
  hosp_totals[, .(hosp_key, site, hospital, n_enc, pct_academic)],
  by = "hosp_key"
)

hosp_long[, pct_raw := n_1 / n * 100]
hosp_long[, eligible := n_enc >= MIN_HOSP_N]

message(sprintf("  Hospitals: %d total, %d with >= %d encounters",
                uniqueN(hosp_long$hosp_key),
                uniqueN(hosp_long[eligible == TRUE]$hosp_key),
                MIN_HOSP_N))

rm(hosp_type, hosp_agents, hosp_any, hosp_dead)

# ==============================================================================
# STEP 2: Raw ranges across eligible hospitals
# ==============================================================================

outcome_order = c(
  "any_escalation",
  agent_vars,
  "dead_hospice_48h"
)

range_summary = hosp_long[
  eligible == TRUE,
  .(
    n_hospitals = .N,
    pooled_pct  = sum(n_1) / sum(n) * 100,
    min_pct     = min(pct_raw),
    p25_pct     = quantile(pct_raw, 0.25),
    median_pct  = median(pct_raw),
    p75_pct     = quantile(pct_raw, 0.75),
    max_pct     = max(pct_raw)
  ),
  by = outcome
]

range_summary[, outcome := factor(outcome, levels = outcome_order)]
setorder(range_summary, outcome)

message("\n  Raw hospital ranges (hospitals with >= ", MIN_HOSP_N, " encounters):")
print(range_summary, digits = 3)

# ==============================================================================
# STEP 3: Random-intercept models (all hospitals)
# ==============================================================================
# logit(p_hospital) = b0 + u_hospital, u ~ N(0, tau^2)
# ICC (latent scale) = tau^2 / (tau^2 + pi^2 / 3)
# MOR = exp(sqrt(2 * tau^2) * qnorm(0.75))

model_rows  = list()
shrunk_rows = list()

for (oc in outcome_order) {

  dt = hosp_long[outcome == oc & n > 0]

  if (nrow(dt) < 2 || sum(dt$n_1) == 0) {
    message(sprintf("  Skipping %s: fewer than 2 hospitals or no events", oc))
    next
  }

  fit = tryCatch(
    glmer(
      cbind(n_1, n - n_1) ~ 1 + (1 | hosp_key),
      data    = dt,
      family  = binomial,
      control = glmerControl(optimizer = "bobyqa")
    ),
    error = function(e) {
      message(sprintf("  Model failed for %s: %s", oc, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) next

  tau2 = as.numeric(VarCorr(fit)$hosp_key)
  b0   = fixef(fit)[["(Intercept)"]]

  model_rows[[oc]] = data.table(
    outcome           = oc,
    n_hospitals       = nrow(dt),
    n_encounters      = sum(dt$n),
    n_events          = sum(dt$n_1),
    pct_mean_hospital = plogis(b0) * 100,
    tau2              = tau2,
    icc               = tau2 / (tau2 + pi^2 / 3),
    mor               = exp(sqrt(2 * tau2) * qnorm(0.75)),
    singular          = isSingular(fit)
  )

  re = ranef(fit)$hosp_key

  shrunk_rows[[oc]] = data.table(
    outcome    = oc,
    hosp_key   = rownames(re),
    pct_shrunk = plogis(b0 + re[["(Intercept)"]]) * 100
  )
}

model_summary = rbindlist(model_rows)
hosp_shrunk   = rbindlist(shrunk_rows)

rm(model_rows, shrunk_rows, oc, dt, fit, tau2, b0, re)

message("\n  Random-intercept models (all hospitals):")
print(model_summary, digits = 3)

# ==============================================================================
# STEP 4: Hospital-level table and save
# ==============================================================================

hosp_table = merge(
  hosp_long,
  hosp_shrunk,
  by    = c("outcome", "hosp_key"),
  all.x = TRUE
)

hosp_table[, outcome := factor(outcome, levels = outcome_order)]
setorder(hosp_table, outcome, -n_enc)

fwrite(hosp_table,    file.path(dir_out, paste0("hospital_level_",   today, ".csv")))
fwrite(range_summary, file.path(dir_out, paste0("hospital_ranges_",  today, ".csv")))
fwrite(model_summary, file.path(dir_out, paste0("hospital_models_",  today, ".csv")))

message("\n  Saved to: ", dir_out)
message("    - hospital_level_*.csv   (raw and shrunken % per hospital and outcome)")
message("    - hospital_ranges_*.csv  (raw ranges, hospitals >= ", MIN_HOSP_N, " encounters)")
message("    - hospital_models_*.csv  (ICC and MOR, all hospitals)")
message("\n  NOTE: hospital IDs are site-internal; anonymize before sharing outside the team.")

message("\n== Hospital variation complete ==")
