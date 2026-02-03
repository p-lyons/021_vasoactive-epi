# ==============================================================================
# 02_variables.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Define: 48h window, escalation, timing groups, SVI/ADI, outcome groups
# ==============================================================================

# Requires: cohort, vasoactive_doses, data_list from 01_cohort.R

if (!exists("cohort")) {
  cohort            = read_parquet(here("proj_tables", "cohort.parquet"))
  hid_jid_crosswalk = read_parquet(here("proj_tables", "hid_jid_crosswalk.parquet"))
}

if (!exists("hid_jid_crosswalk")) {
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
# STEP 2: Identify escalation from vasoactive_doses (continuous drugs)
# ==============================================================================

message("  Loading vasoactive doses...")

va = read_parquet(here("proj_tables", "vasoactive_doses.parquet"))
setDT(va)

## escalation = epi, phenyl, dopa, or a2 after T0 within 48h -------------------

message("  Identifying escalation events (continuous drugs)...")

escalation_cont = 
  
  join(va, date_frame, how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm > t0_dttm & admin_dttm <= endpoint_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    epi_01    = as.integer(any(!is.na(epi_dose) & epi_dose > 0)),
    phenyl_01 = as.integer(any(!is.na(phenyl_dose) & phenyl_dose > 0)),
    dopa_01   = as.integer(any(!is.na(dopa_dose) & dopa_dose > 0)),
    a2_01     = as.integer(any(!is.na(a2_dose) & a2_dose > 0))
  )

cohort = join(cohort, escalation_cont, how = "left", multiple = FALSE)

# fill NA with 0
for (v in c("epi_01", "phenyl_01", "dopa_01", "a2_01")) {
  if (!v %in% names(cohort)) {
    cohort[[v]] = rep(0L, nrow(cohort))
  } else {
    cohort[[v]] = fifelse(is.na(cohort[[v]]), 0L, as.integer(cohort[[v]]))
  }
}

rm(escalation_cont)

# ==============================================================================
# STEP 3: Methylene blue and B12 from medication_admin_intermittent
# ==============================================================================

message("  Extracting methylene blue and B12 from intermittent meds...")

## check if medication_admin_intermittent is available -------------------------

if ("medication_admin_intermittent" %in% names(data_list)) {
  
  ## methylene blue -------------------------------------------------------------
  
  mb_raw = tryCatch({
    dplyr::filter(data_list$medication_admin_intermittent, 
                  tolower(med_category) == "methylene_blue") |>
      dplyr::filter(hospitalization_id %in% cohort_hids) |>
      dplyr::select(hospitalization_id, admin_dttm, med_dose, med_dose_unit) |>
      dplyr::collect() |>
      distinct()
  }, error = function(e) {
    message("    ⚠️  Could not extract methylene blue: ", e$message)
    data.table()
  })
  
  if (nrow(mb_raw) > 0) {
    mb = 
      join(mb_raw, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
      join(date_frame, how = "inner", multiple = TRUE) |>
      fsubset(admin_dttm > t0_dttm & admin_dttm <= endpoint_dttm) |>
      fsubset(!is.na(med_dose) & med_dose > 0) |>
      fgroup_by(joined_hosp_id) |>
      fsummarize(mb_01 = 1L)
    
    cohort = join(cohort, mb, how = "left", multiple = FALSE)
    rm(mb)
    message(sprintf("    MB administrations found: %d encounters", sum(!is.na(cohort$mb_01) & cohort$mb_01 == 1)))
  }
  
  rm(mb_raw)
  
  ## hydroxocobalamin (B12) -----------------------------------------------------
  
  b12_raw = tryCatch({
    dplyr::filter(data_list$medication_admin_intermittent, 
                  tolower(med_category) == "hydroxocobalamin") |>
      dplyr::filter(hospitalization_id %in% cohort_hids) |>
      dplyr::select(hospitalization_id, admin_dttm, med_dose, med_dose_unit) |>
      dplyr::collect() |>
      distinct()
  }, error = function(e) {
    message("    ⚠️  Could not extract hydroxocobalamin: ", e$message)
    data.table()
  })
  
  if (nrow(b12_raw) > 0) {
    b12 = 
      join(b12_raw, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
      join(date_frame, how = "inner", multiple = TRUE) |>
      fsubset(admin_dttm > t0_dttm & admin_dttm <= endpoint_dttm) |>
      fsubset(!is.na(med_dose) & med_dose > 0) |>
      fgroup_by(joined_hosp_id) |>
      fsummarize(b12_01 = 1L)
    
    cohort = join(cohort, b12, how = "left", multiple = FALSE)
    rm(b12)
    message(sprintf("    B12 administrations found: %d encounters", sum(!is.na(cohort$b12_01) & cohort$b12_01 == 1)))
  }
  
  rm(b12_raw)
  
} else {
  message("    ⚠️  medication_admin_intermittent not in data_list - skipping MB/B12")
}

# fill NA with 0 for MB and B12
for (v in c("mb_01", "b12_01")) {
  if (!v %in% names(cohort)) {
    cohort[[v]] = rep(0L, nrow(cohort))
  } else {
    cohort[[v]] = fifelse(is.na(cohort[[v]]), 0L, as.integer(cohort[[v]]))
  }
}

## escalated flag = any of the escalation drugs --------------------------------

cohort = ftransform(cohort,
                    escalated_01 = as.integer(
                      epi_01 == 1 | phenyl_01 == 1 | dopa_01 == 1 | a2_01 == 1 | mb_01 == 1 | b12_01 == 1
                    )
)

message(sprintf("    Escalated (any drug): %d (%.1f%%)", 
                sum(cohort$escalated_01), 100 * mean(cohort$escalated_01)))

# ==============================================================================
# STEP 4: Doses at T0 and max NE equivalent
# ==============================================================================

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

# NE equivalent: NE 1:1, Epi 1:1, VP 2.5 per 0.01 units/min, Phenyl 0.1:1, A2 0.01:1

dose_in_window = 
  join(va, date_frame, how = "inner", multiple = TRUE) |>
  fsubset(admin_dttm >= t0_dttm & admin_dttm <= endpoint_dttm)

dose_in_window = ftransform(dose_in_window,
                            ne_equiv = 
                              fifelse(is.na(ne_dose),     0, ne_dose) +
                              fifelse(is.na(epi_dose),    0, epi_dose) +
                              fifelse(is.na(vp_dose),     0, 2.5 * vp_dose) +
                              fifelse(is.na(phenyl_dose), 0, 0.1 * phenyl_dose) +
                              fifelse(is.na(a2_dose),     0, 0.01 * a2_dose)
)

max_ne_equiv = 
  fgroup_by(dose_in_window, joined_hosp_id) |>
  fsummarize(max_ne_equiv_48h = fmax(ne_equiv, na.rm = TRUE))

max_ne_equiv[is.infinite(max_ne_equiv_48h), max_ne_equiv_48h := NA_real_]

cohort = join(cohort, max_ne_equiv, how = "left", multiple = FALSE)

rm(va, t0_doses, dose_in_window, max_ne_equiv, date_frame)
gc()

# ==============================================================================
# STEP 5: ICU timing and timing groups
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

## IMV timing group ------------------------------------------------------------

message("  Adding IMV timing group...")

# 0 = no IMV
# 1 = IMV before T0
# 2 = IMV T0 to endpoint

cohort = ftransform(cohort,
                    imv_timing_group = case_when(
                      is.na(imv_dttm)                                    ~ 0L,
                      imv_dttm < t0_dttm                                 ~ 1L,
                      imv_dttm >= t0_dttm & imv_dttm <= endpoint_dttm   ~ 2L,
                      TRUE                                              ~ 0L
                    )
)

# also create simple binary flags
cohort = ftransform(cohort,
                    imv_at_t0_01 = fifelse(imv_timing_group == 1L, 1L, 0L),
                    imv_48h_01   = fifelse(imv_timing_group %in% c(1L, 2L), 1L, 0L)
)

## CRRT timing group -----------------------------------------------------------

message("  Adding CRRT timing group...")

# 0 = no CRRT
# 1 = CRRT before T0
# 2 = CRRT T0 to endpoint

cohort = ftransform(cohort,
                    crrt_timing_group = case_when(
                      is.na(crrt_dttm)                                    ~ 0L,
                      crrt_dttm < t0_dttm                                 ~ 1L,
                      crrt_dttm >= t0_dttm & crrt_dttm <= endpoint_dttm   ~ 2L,
                      TRUE                                                ~ 0L
                    )
)

cohort = ftransform(cohort,
                    crrt_01       = fifelse(crrt_timing_group %in% c(1L, 2L), 1L, 0L),
                    crrt_at_t0_01 = fifelse(crrt_timing_group == 1L, 1L, 0L)
)

# ==============================================================================
# STEP 6: Code status at T0
# ==============================================================================

message("  Adding code status...")

cohort_pats = funique(cohort$patient_id)

code_status_raw = 
  dplyr::filter(data_list$code_status, patient_id %in% cohort_pats) |>
  dplyr::select(patient_id, start_dttm, code_status_category) |>
  dplyr::collect() |>
  distinct()

# standardize code status categories
code_status_raw = ftransform(code_status_raw,
                             code_status_category = tolower(code_status_category) |>
                               str_replace_all("/", "_") |>
                               str_replace_all(" ", "_") |>
                               str_replace_all("dnar", "dnr")
)

# collapse to: full, dnr_dni, partial
# - full = full code
# - dnr_dni = DNR + DNI combined
# - partial = DNR only, DNI only, or other limitations
code_status_raw = ftransform(code_status_raw,
                             code_status_category = case_when(
                               code_status_category %in% c("full", "full_code")
                               ~ "full",
                               code_status_category %in% c("dnr_dni", "dnr_and_dni", "dni_dnr")
                               ~ "dnr_dni",
                               code_status_category %in% c("dnr", "dni", "partial", "special", "special_partial", "limited")
                               ~ "partial",
                               TRUE ~ "other"
                             )
)

# join to cohort to get admission and T0 times
code_status = 
  join(code_status_raw, 
       fselect(cohort, patient_id, joined_hosp_id, admission_dttm, t0_dttm), 
       how = "inner", multiple = TRUE) |>
  fsubset(start_dttm >= admission_dttm & start_dttm <= t0_dttm) |>
  ftransform(code_status_category = tolower(code_status_category))

## Step 1: Initial code status = LAST code status in hours 0-12 ----------------
## (avoids placeholder orders that quickly get changed)

initial_window = 
  ftransform(code_status,
             hours_from_admit = as.numeric(difftime(start_dttm, admission_dttm, units = "hours"))
  ) |>
  fsubset(hours_from_admit >= 0 & hours_from_admit <= 12) |>
  roworder(start_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(initial_code_status = flast(code_status_category))

## Step 2: Most recent code status at T0 (LOCF) --------------------------------

code_at_t0 = 
  code_status |>
  fsubset(start_dttm <= t0_dttm) |>
  roworder(start_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(code_status_t0 = flast(code_status_category))

## Step 3: Combine -------------------------------------------------------------
## Use code_status_t0 if available, else initial_code_status, else presume_full

cohort = join(cohort, initial_window, how = "left", multiple = FALSE)
cohort = join(cohort, code_at_t0, how = "left", multiple = FALSE)

# flag: was code status documented (vs presumed full due to missing)?
# true if ANY code status was documented by T0
cohort = ftransform(cohort,
                    code_documented_01 = fifelse(!is.na(code_status_t0) | !is.na(initial_code_status), 1L, 0L)
)

# if no code status in first 12h, presume full
cohort = ftransform(cohort,
                    initial_code_status = fifelse(is.na(initial_code_status), "full", initial_code_status)
)

# code status at T0: use most recent if available, otherwise fall back to initial
cohort = ftransform(cohort,
                    code_status_t0 = case_when(
                      !is.na(code_status_t0)      ~ code_status_t0,
                      !is.na(initial_code_status) ~ initial_code_status,
                      TRUE                        ~ "full"
                    )
)

message(sprintf("    Code status distribution at T0:"))
print(table(cohort$code_status_t0, useNA = "ifany"))

rm(code_status_raw, code_status, initial_window, code_at_t0)
gc()

# ==============================================================================
# STEP 7: SVI and ADI (from config files)
# ==============================================================================

message("  Adding SVI/ADI...")

## SVI (Social Vulnerability Index) --------------------------------------------

svi_file = here("config", "svi_2020.parquet")

if (file.exists(svi_file)) {
  svi_list = read_parquet(svi_file)
  
  # extract census tract from census_block_code (first 11 digits)
  cohort = ftransform(cohort,
                      census_tract = substr(census_block_code, 1, 11)
  )
  
  cohort = join(cohort, svi_list, on = "census_tract", how = "left", multiple = FALSE)
  
  message(sprintf("    SVI linked: %d of %d encounters (%.1f%%)",
                  sum(!is.na(cohort$svi_percentile)),
                  nrow(cohort),
                  100 * mean(!is.na(cohort$svi_percentile))))
} else {
  message("    ⚠️  SVI file not found: config/svi_2020.parquet")
  cohort$svi_percentile = NA_real_
}

## ADI (Area Deprivation Index) ------------------------------------------------

adi_file = here("config", "adi_2020.parquet")

if (file.exists(adi_file)) {
  adi_list = read_parquet(adi_file)
  
  # extract census block group from census_block_code (first 12 digits)
  # only if not already present
  if (!"census_block_group_code" %in% names(cohort) || all(is.na(cohort$census_block_group_code))) {
    cohort = ftransform(cohort,
                        census_block_group_code = substr(census_block_code, 1, 12)
    )
  }
  
  cohort = join(cohort, adi_list, on = "census_block_group_code", how = "left", multiple = FALSE)
  
  message(sprintf("    ADI linked: %d of %d encounters (%.1f%%)",
                  sum(!is.na(cohort$adi_percentile)),
                  nrow(cohort),
                  100 * mean(!is.na(cohort$adi_percentile))))
} else {
  message("    ⚠️  ADI file not found: config/adi_2020.parquet")
  cohort$adi_percentile = NA_real_
}

# ==============================================================================
# STEP 8: Derived variables
# ==============================================================================

message("  Creating derived variables...")

## LOS to T0 -------------------------------------------------------------------

cohort = ftransform(cohort,
                    los_to_t0_d = as.numeric(difftime(t0_dttm, admission_dttm, units = "hours")) / 24
)

## LOS from T0 to discharge ----------------------------------------------------

cohort = ftransform(cohort,
                    los_from_t0_d = as.numeric(difftime(discharge_dttm, t0_dttm, units = "hours")) / 24
)

## ICU LOS to T0 (in days for table) -------------------------------------------

cohort = ftransform(cohort,
                    icu_los_to_t0_d = icu_los_to_t0_h / 24
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

## LOS categories (days) -------------------------------------------------------

cohort = ftransform(cohort,
                    los_cat = cut(los_to_t0_d, 
                                  breaks = c(0, 2, 4, 7, 14, Inf), 
                                  labels = c("0-47h", "48h-96h", "96h-1wk", "1wk-2wk", "2wk+"), 
                                  right = FALSE)
)

## ICU LOS categories (days) ---------------------------------------------------

cohort = ftransform(cohort,
                    icu_los_cat = cut(icu_los_to_t0_d, 
                                      breaks = c(-Inf, 0, 2, 4, 7, 14, Inf), 
                                      labels = c("pre-ICU", "0-47h", "48h-96h", "96h-1wk", "1wk-2wk", "2wk+"), 
                                      right = FALSE)
)

# ==============================================================================
# STEP 9: Outcome groups (3 groups)
# ==============================================================================

message("  Finalizing outcome groups...")

## death/hospice flag ----------------------------------------------------------

cohort = ftransform(cohort,
                    dead_hospice_01 = fifelse(dead_01 == 1 | hospice_01 == 1, 1L, 0L)
)

cohort$dead_hospice_01 = fifelse(is.na(cohort$dead_hospice_01), 0L, cohort$dead_hospice_01)
cohort$escalated_01    = fifelse(is.na(cohort$escalated_01), 0L, cohort$escalated_01)

## three-level outcome ---------------------------------------------------------

# Group 1: Died/hospice (regardless of escalation)
# Group 2: Escalated + alive
# Group 3: No escalation + alive

cohort = ftransform(cohort,
                    outcome_group = case_when(
                      dead_hospice_01 == 1                     ~ "dead_hospice",
                      escalated_01 == 1 & dead_hospice_01 == 0 ~ "escalated",
                      escalated_01 == 0 & dead_hospice_01 == 0 ~ "stable",
                      TRUE                                     ~ "other"
                    )
)

cohort = ftransform(cohort,
                    outcome_group = factor(outcome_group, 
                                           levels = c("dead_hospice", "escalated", "stable", "other"))
)

message("\n  Outcome distribution:")
print(table(cohort$outcome_group, useNA = "ifany"))

# ==============================================================================
# STEP 10: Save analytic cohort
# ==============================================================================

message("\n== Saving analytic cohort ==")

write_parquet(cohort, here("proj_tables", "cohort_analytic.parquet"))

message(sprintf("  ✅ Saved cohort_analytic.parquet: %d encounters", nrow(cohort)))
message("\n== 02_variables.R complete ==")

# go to 03
