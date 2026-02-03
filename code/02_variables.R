# ==============================================================================
# 02_variables.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Define: 48h window, escalation, outcome groups
# ==============================================================================

# Requires: cohort, vasoactive_doses from 01_cohort.R

if (!exists("cohort")) {
  cohort            = read_parquet(here("proj_tables", "cohort.parquet"))
  hid_jid_crosswalk = read_parquet(here("proj_tables", "hid_jid_crosswalk.parquet"))
}

if (nrow(cohort) == 0) {
  stop("Cohort is empty. Check 01_cohort.R output.", call. = FALSE)
}

message(sprintf("\n== Preparing outcome variables for %d encounters ==", nrow(cohort)))

cohort_jids = funique(cohort$joined_hosp_id)
cohort_hids = funique(hid_jid_crosswalk$hospitalization_id)

# ==============================================================================
# STEP 1: Define 48h outcome window
# ==============================================================================

message("  Defining 48h outcome window...")

endpoint_hours = 48L

cohort = ftransform(cohort, 
  endpoint_dttm = t0_dttm + lubridate::dhours(endpoint_hours)
)

date_frame = fselect(cohort, joined_hosp_id, t0_dttm, endpoint_dttm, admission_dttm, discharge_dttm)

# ==============================================================================
# STEP 2: Identify escalation and extract doses
# ==============================================================================

message("  Loading vasoactive doses...")

va = read_parquet(here("proj_tables", "vasoactive_doses.parquet"))
setDT(va)

## escalation = epi, phenyl, or dopa after T0 within 48h -----------------------

message("  Identifying escalation events...")

escalation = 
  join(va, date_frame, how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm > t0_dttm & admin_dttm <= endpoint_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    epi_01    = as.integer(any(!is.na(epi_dose) & epi_dose > 0)),
    phenyl_01 = as.integer(any(!is.na(phenyl_dose) & phenyl_dose > 0)),
    dopa_01   = as.integer(any(!is.na(dopa_dose) & dopa_dose > 0))
  )

escalation = ftransform(escalation,
  escalated_01 = as.integer(epi_01 == 1 | phenyl_01 == 1 | dopa_01 == 1)
)

cohort = join(cohort, escalation, how = "left", multiple = FALSE)

# fill NA with 0
for (v in c("epi_01", "phenyl_01", "dopa_01", "escalated_01")) {
  if (!v %in% names(cohort)) {
    cohort[[v]] = rep(0L, nrow(cohort))
  } else {
    cohort[[v]] = fifelse(is.na(cohort[[v]]), 0L, as.integer(cohort[[v]]))
  }
}

message(sprintf("    Escalated: %d (%.1f%%)", 
                sum(cohort$escalated_01), 100 * mean(cohort$escalated_01)))

## doses at T0 -----------------------------------------------------------------

message("  Extracting doses at T0...")

t0_doses = 
  join(va, date_frame, how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm == t0_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    ne_dose_t0 = fmax(ne_dose, na.rm = TRUE),
    vp_dose_t0 = fmax(vp_dose, na.rm = TRUE)
  )

t0_doses[is.infinite(ne_dose_t0), ne_dose_t0 := NA_real_]
t0_doses[is.infinite(vp_dose_t0), vp_dose_t0 := NA_real_]

cohort = join(cohort, t0_doses, how = "left", multiple = FALSE)

## max NE equivalent within 48h ------------------------------------------------

message("  Calculating max NE equivalent...")

# NE equivalent: NE 1:1, Epi 1:1, VP 2.5 per 0.01 units/min, Phenyl 0.1:1, Dopa 0.01:1

dose_in_window = 
  join(va, date_frame, how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm >= t0_dttm & admin_dttm <= endpoint_dttm)

dose_in_window = ftransform(dose_in_window,
  ne_equiv = 
    fifelse(is.na(ne_dose),     0, ne_dose) +
    fifelse(is.na(epi_dose),    0, epi_dose) +
    fifelse(is.na(vp_dose),     0, 2.5 * vp_dose) +
    fifelse(is.na(phenyl_dose), 0, 0.1 * phenyl_dose) +
    fifelse(is.na(dopa_dose),   0, 0.01 * dopa_dose)
)

max_ne_equiv = 
  fgroup_by(dose_in_window, joined_hosp_id) |>
  fsummarize(max_ne_equiv_48h = fmax(ne_equiv, na.rm = TRUE))

max_ne_equiv[is.infinite(max_ne_equiv_48h), max_ne_equiv_48h := NA_real_]

cohort = join(cohort, max_ne_equiv, how = "left", multiple = FALSE)

rm(va, escalation, t0_doses, dose_in_window, max_ne_equiv, date_frame)
gc()

# ==============================================================================
# STEP 3: ICU and IMV timing relative to T0
# ==============================================================================

message("  Adding ICU timing...")

icu_times = 
  dplyr::filter(data_list$adt, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(tolower(location_category) == "icu") |>
  dplyr::select(hospitalization_id, in_dttm) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  fselect(joined_hosp_id, in_dttm) |>
  distinct()

icu_times = 
  join(icu_times, fselect(cohort, joined_hosp_id, admission_dttm, discharge_dttm), 
       how = "inner", multiple = TRUE) |>
  fsubset(in_dttm >= admission_dttm & in_dttm <= discharge_dttm) |>
  roworder(in_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(first_icu_dttm = ffirst(in_dttm))

cohort = join(cohort, icu_times, how = "left", multiple = FALSE)

cohort = ftransform(cohort,
  icu_los_to_t0_h = as.numeric(difftime(t0_dttm, first_icu_dttm, units = "hours"))
)

rm(icu_times)

## IMV flags -------------------------------------------------------------------

message("  Adding IMV flags...")

# IMV at T0 (within 6h before T0)
cohort = ftransform(cohort,
  imv_at_t0_01 = fifelse(
    !is.na(imv_dttm) & imv_dttm <= t0_dttm & imv_dttm >= (t0_dttm - lubridate::dhours(6)),
    1L, 0L
  )
)

# IMV within 48h of T0
cohort = ftransform(cohort,
  imv_48h_01 = fifelse(
    !is.na(imv_dttm) & imv_dttm <= endpoint_dttm,
    1L, 0L
  )
)

# ==============================================================================
# STEP 4: Derived variables
# ==============================================================================

message("  Creating derived variables...")

## LOS to T0 -------------------------------------------------------------------

cohort = ftransform(cohort,
  los_to_t0_d = as.numeric(difftime(t0_dttm, admission_dttm, units = "hours")) / 24
)

## age categories --------------------------------------------------------------

cohort = ftransform(cohort,
  age_cat = cut(age, 
                breaks = c(18, 40, 50, 60, 70, 80, Inf), 
                labels = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+"), 
                right = FALSE)
)

## Elixhauser categories -------------------------------------------------------

cohort = ftransform(cohort,
  vw_cat = cut(vw, 
               breaks = c(-Inf, 1, 5, 10, 15, Inf), 
               labels = c("<=0", "1-4", "5-9", "10-14", ">=15"), 
               right = FALSE)
)

# ==============================================================================
# STEP 5: Outcome groups
# ==============================================================================

message("  Finalizing outcome groups...")

## death/hospice flag ----------------------------------------------------------

cohort = ftransform(cohort,
  dead_hospice_01 = fifelse(dead_01 == 1 | hospice_01 == 1, 1L, 0L)
)

cohort$dead_hospice_01 = fifelse(is.na(cohort$dead_hospice_01), 0L, cohort$dead_hospice_01)
cohort$escalated_01    = fifelse(is.na(cohort$escalated_01), 0L, cohort$escalated_01)

## three-level outcome ---------------------------------------------------------

# Group 1: Died/hospice
# Group 2: Escalated (alive)
# Group 3: Stable/de-escalated (alive, no escalation)

cohort = ftransform(cohort,
  outcome_group = case_when(
    dead_hospice_01 == 1                     ~ "dead_hospice",
    escalated_01 == 1 & dead_hospice_01 == 0 ~ "escalated",
    escalated_01 == 0 & dead_hospice_01 == 0 ~ "stable_deescalated",
    TRUE                                     ~ "other"
  )
)

cohort = ftransform(cohort,
  outcome_group = factor(outcome_group, 
                         levels = c("dead_hospice", "escalated", "stable_deescalated", "other"))
)

message("\n  Outcome distribution:")
print(table(cohort$outcome_group, useNA = "ifany"))

# ==============================================================================
# STEP 6: Save analytic cohort
# ==============================================================================

message("\n== Saving analytic cohort ==")

write_parquet(cohort, here("proj_tables", "cohort_analytic.parquet"))

message(sprintf("  ✅ Saved cohort_analytic.parquet: %d encounters", nrow(cohort)))
message("\n== 02_variables.R complete ==")

# go to 03
