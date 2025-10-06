

# Requires data_list to be loaded/validated from 00_*

# resources for RAM heavy wrangling --------------------------------------------

## cores & RAM (reuse from 00 if available) ------------------------------------

os_type   = if (exists("os_type"))   os_type   else Sys.info()[["sysname"]]
all_cores = if (exists("all_cores")) all_cores else {
  x = parallel::detectCores(logical = TRUE); if (is.na(x)) 1L else as.integer(x)
}

if (!exists("avail_ram_gb") || !is.finite(avail_ram_gb)) {
  get_ram_gb = function() {
    tryCatch({
      if (Sys.info()[["sysname"]] == "Darwin") {
        bytes = suppressWarnings(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)))
        if (length(bytes) > 0 && !is.na(bytes)) bytes / 1024^3 else NA_real_
      } else if (file.exists("/proc/meminfo")) {
        kb = suppressWarnings(as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)))
        if (length(kb) > 0 && !is.na(kb)) kb / 1024^2 else NA_real_
      } else if (requireNamespace("ps", quietly = TRUE)) {
        ps::ps_system_memory()[["available"]] / 1024^3
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  avail_ram_gb = get_ram_gb()
}

## choose threads (2 GB per thread; leave 1 core free) -------------------------

reserve_cores  = 1L
gb_per_thread  = 2.0
max_by_cores   = max(1L, all_cores - reserve_cores)
max_by_memory  = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads      = as.integer(max(1L, min(max_by_cores, max_by_memory)))
n_math_threads = as.integer(max(1L, min(n_threads, 8L)))

## apply thread settings -------------------------------------------------------

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

message(
  sprintf("01 resources | OS=%s | Cores=%d | Threads=%d | MathThreads=%d | Avail RAM≈%s GB (rule: 2 GB/core)",
          os_type, all_cores, n_threads, n_math_threads,
          ifelse(is.finite(avail_ram_gb), round(avail_ram_gb, 1), "NA"))
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# cohort identification --------------------------------------------------------

## start by linking contiguous hospitalizations --------------------------------

### encounters with age >= 18 and dates within study window --------------------

hosp_blocks = 
  dplyr::filter(data_list$hospitalization, age_at_admission >= 18) |>
  dplyr::filter(admission_dttm >= start_date & admission_dttm <= end_date) |> 
  dplyr::filter(admission_dttm < discharge_dttm & !is.na(discharge_dttm)) |>
  dplyr::select(patient_id, hospitalization_id, admission_dttm, discharge_dttm) |>
  dplyr::collect() |>
  roworder(patient_id, admission_dttm)

### use data.table to find joined hospitalizations with <= 6h gaps -------------

link_hours = 6L
linked     = as.data.table(hosp_blocks)
setorder(linked, patient_id, admission_dttm)

#### calculate gaps between encounters
linked[, next_admit := shift(admission_dttm, type = "lead"), by = patient_id]
linked[, next_gap   := as.numeric(difftime(next_admit, discharge_dttm, units = "hours"))]
linked[, prev_dc    := shift(discharge_dttm, type = "lag"), by = patient_id]  
linked[, prev_gap   := as.numeric(difftime(admission_dttm, prev_dc, units = "hours"))]

#### mark encounters that should be linked
linked[, link_flag := (next_gap < link_hours | prev_gap < link_hours)]
linked[is.na(link_flag), link_flag := FALSE]

#### create unique joined hospitalization ID, new group whenever gap > link_hours from  previous discharge
linked[, new_group := is.na(prev_gap) | prev_gap >= link_hours]
linked[, joined_hosp_id := .GRP, by = .(patient_id, cumsum(new_group))]

#### create hid_jid_crosswalk --------------------------------------------------
hid_jid_crosswalk = select(linked, ends_with("id")) |> as_tidytable()

## hospital ward admissions ----------------------------------------------------

### stay requires ed or icu ----------------------------------------------------

inpatient_hids = 
  dplyr::filter(data_list$adt, tolower(location_category) %in% c("ed", "icu")) |>
  dplyr::select(hospitalization_id) |>
  dplyr::collect() |>
  funique() |>
  tibble::deframe()

inpatient_jids = 
  fsubset(hid_jid_crosswalk, hospitalization_id %in% inpatient_hids) |>
  select(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

### don't want to include psych ------------------------------------------------

drop_ob = 
  dplyr::filter(data_list$adt, tolower(location_category) %in% c("psych", "rehab")) |>
  dplyr::select(hospitalization_id) |>
  dplyr::collect() |>
  funique() |>
  tibble::deframe()

drop_ob_jids = 
  fsubset(hid_jid_crosswalk, hospitalization_id %in% drop_ob) |>
  select(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

linked = fsubset(linked,  joined_hosp_id %in% inpatient_jids) 
linked = fsubset(linked, !joined_hosp_id %in% drop_ob_jids)
linked = select(linked, ends_with("id"), ends_with("dttm"))

### an admission requires at least 1 full set of vital signs -------------------

has_vital_signs = 
  dplyr::filter(data_list$vitals, hospitalization_id %in% inpatient_hids) |>
  dplyr::filter(vital_category %in% req_vitals) |>
  dplyr::select(hospitalization_id, vital_category, recorded_dttm) |>
  dplyr::collect() |>
  funique() 

has_vital_signs = 
  join(has_vital_signs, linked, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= admission_dttm & recorded_dttm <= recorded_dttm) |>
  fselect(joined_hosp_id, vital_category) |>
  funique() |>
  summarize(n = n(), .by = joined_hosp_id) |>
  fsubset(n == length(req_vitals)) |>
  pull(joined_hosp_id) 

linked            = fsubset(linked, joined_hosp_id %in% has_vital_signs) 
hid_jid_crosswalk = select(linked, ends_with("id"))
cohort_hids       = funique(hid_jid_crosswalk$hospitalization_id)
cohort_pats       = funique(hid_jid_crosswalk$patient_id)

### clean up helpers -----------------------------------------------------------

rm(inpatient_hids, inpatient_jids, drop_ob, drop_ob_jids)
rm(has_vital_signs, hosp_blocks, link_hours)
gc()

## assemble cohort data frame --------------------------------------------------

#### pull additional data for cohort filtering and final variables
cohort_data = 
  dplyr::filter(data_list$hospitalization, hospitalization_id %in% cohort_hids) |>
  dplyr::select(ends_with("id"), age_at_admission, discharge_category) |> 
  dplyr::collect()

hid_dups_source =
  fcount(cohort_data, hospitalization_id) |>
  fsubset(N > 1) 

if (nrow(hid_dups_source) > 0) {
  stop(
    sprintf("Source has duplicate hospitalization_id: %s",
            paste(head(hid_dups_source$hospitalization_id, 50), collapse = ", ")),
    call. = FALSE
  )
}

#### create final cohort - 1 row per joined_hosp_id
cohort = 
  join(linked, cohort_data, how = "left", multiple = T) |>
  roworder(admission_dttm) |>
  fgroup_by(patient_id, joined_hosp_id) |>
  fsummarize(
    age                = ffirst(age_at_admission),
    admission_dttm     = ffirst(admission_dttm),
    discharge_dttm     = flast(discharge_dttm),
    discharge_category = flast(discharge_category)
  )

#### clean up temporary variables
rm(linked, cohort_data); gc()

## quality control -------------------------------------------------------------

### check for duplicates -------------------------------------------------------

dupes = cohort |> janitor::get_dupes(patient_id, admission_dttm)

if (nrow(dupes) > 0) {
  dup_ids = funique(dupes$joined_hosp_id)
  stop(
    sprintf("Found %d duplicate joined_hosp_id(s): %s", length(dup_ids), paste(dup_ids, collapse = ", ")),
    call. = FALSE
  )
}

message("✅ No duplicate joined_hosp_id found.")

### YODO (you only die once) ---------------------------------------------------

#### death duplicates ----------------------------------------------------------

dup_deaths = 
  fsubset(cohort, discharge_category == "Expired") |>
  roworder(admission_dttm, discharge_dttm) |>
  fgroup_by(patient_id) |>
  fmutate(one = 1L) |>
  fmutate(n_deaths = fsum(one), counter = fcumsum(one)) |>
  fungroup() |>
  fsubset(n_deaths > 1) 

if (nrow(dup_deaths) > 0) {
  
  encdrop = 
    fsubset(dup_deaths, counter == 1) |>
    select(joined_hosp_id) |>
    tibble::deframe()
  
  cohort  = fsubset(cohort, !joined_hosp_id %in% encdrop)
  n_dupes = sum(dup_deaths$counter > 1)
  n_pats  = length(unique(dup_deaths$patient_id))
  
  message(
    sprintf("Removed %d duplicate deaths for %d patients: %s",
            n_dupes, n_pats, paste(unique(dup_deaths$patient_id), collapse = ", "))
  )
}

#### readmissions following death ----------------------------------------------

death_times = 
  fsubset(cohort, tolower(discharge_category) == "expired") |>
  roworder(discharge_dttm) |>
  fgroup_by(patient_id) |>
  fsummarise(death_instant = ffirst(discharge_dttm))

post_death_admissions = 
  join(cohort, death_times, how = "inner", multiple = T) |>
  fsubset(admission_dttm >= death_instant) 

if (nrow(post_death_admissions) > 0) {
  cohort = fsubset(cohort, !joined_hosp_id %in% post_death_admissions$joined_hosp_id)
  
  message(
    sprintf("Removed %d post-death encs for %d patients (? organ donors): %s",
            nrow(post_death_admissions),
            length(unique(post_death_admissions$patient_id)),
            paste(unique(post_death_admissions$patient_id), collapse = ", "))
  )
}

message("✅ Cleaned duplicate deaths and post-death encounters.")

rm(dupes, dup_deaths, death_times, post_death_admissions, start_date, end_date)
rm(encdrop, file_type, n_dupes, n_pats); gc()

cohort_pats = funique(cohort$patient_id)
cohort_jids = funique(cohort$joined_hosp_id)
cohort_hids = funique(hid_jid_crosswalk$hospitalization_id)
date_frame  = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))

## vasopressors ----------------------------------------------------------------

## vasoactives -----------------------------------------------------------------

vlist = c(
  "norepinephrine", 
  "vasopressin"
)

va = 
  dplyr::filter(data_list$medication_admin_continuous, med_category %in% vlist) |>
  dplyr::filter(hospitalization_id %in% cohort_hids) |>
  dplyr::select(
    hospitalization_id, 
    admin_dttm, 
    med_category, 
    strength, 
    med_dose,
    med_dose_unit,
    mar_action_category
  ) |>
  dplyr::collect() |>
  funique()

va = 
  join(va, hid_jid_crosswalk, how = "inner", multiple = T) |>
  join(date_frame,            how = "inner", multiple = T) |>
  fsubset(admin_dttm >= admission_dttm) |>
  fsubset(admin_dttm <= discharge_dttm) |>
  fsubset() |> # make sure med is GOING
  ftransform(admin_dttm = lubridate::floor_date(admin_dttm, unit = "hours")) |>
  ftransform() |> # update dose units if need be
  fgroup_by(joined_hosp_id, admin_dttm, med_category) |>
  fsummarize(med_dose = fmax(med_dose)) |>
  pivot_wider(names_from = med_category, values_from = med_dose)

va = fsubset(va, !is.na(vasopressin))
va = fsubset(va, norepinephrine >= 0.2)

cohort_jids       = funique(va$joined_hosp_id)
cohort            = fsubset(cohort, joined_hosp_id %in% cohort_jids)
hid_jid_crosswalk = fsubset(hid_jid_crosswalk, joined_hosp_id %in% cohort_jids)
cohort_hids       = funique(hid_jid_crosswalk$hospitalization_id)
cohort_pats       = funique(cohort$patient_id)
date_frame        = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))

## time starts at the moment vasopressor requirements are met ------------------

time_0 = 
  roworder(va, admin_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(time_0 = ffirst(admin_dttm))

cohort = join(cohort, time_0, how = "left", multiple = F)

# prepare additional cohort details --------------------------------------------

## patient demographics --------------------------------------------------------

### pull demographics from arrow table -----------------------------------------

cohort_demographics = 
  dplyr::filter(data_list$patient, patient_id %in% cohort_pats) |>
  dplyr::select(patient_id, death_dttm, ends_with("category")) |>
  dplyr::collect() |>
  funique()

pt_dups = 
  fcount(cohort_demographics, patient_id, name = "n") |>
  fsubset(n > 1) 

if (nrow(pt_dups) > 0) {
  stop(
    sprintf("Duplicate patient_id(s): %s", paste(pt_dups$patient_id, collapse = ", ")), call. = F
  )
}

### add demographics to cohort df ----------------------------------------------

cohort = 
  join(cohort, cohort_demographics, how = "left", multiple = F) |>
  fmutate(age        = if_else(age > 90, 90.9, age)) |>
  fmutate(female_01  = if_else(tolower(sex_category) == "female", 1L, 0L)) |>
  fmutate(dead_01    = if_else(tolower(discharge_category) == "expired", 1L, 0L)) |>
  fmutate(hospice_01 = if_else(tolower(discharge_category) == "hospice", 1L, 0L)) |>
  fmutate(los_hosp_d = as.numeric(difftime(discharge_dttm, admission_dttm), "hours")/24) |>
  mutate(across(
    .cols = where(is.character) & !any_of("patient_id"),
    .fns  = ~tolower(.x)
  )) |>
  fselect(-sex_category, -discharge_category)

rm(pt_dups, cohort_demographics); gc()

## admission code status -------------------------------------------------------

codes = 
  data_list[["code_status"]] |>
  dplyr::select(patient_id, start_dttm, code_status_category) |>
  dplyr::collect() |>
  funique()

codes = 
  join(cohort, codes, how = "left", multiple = T) |>
  fsubset(start_dttm >= admission_dttm - lubridate::ddays(1)) |>
  fsubset(start_dttm <= discharge_dttm) |>
  roworder(start_dttm) 

events = select(codes, patient_id, event_dttm = start_dttm, event = code_status_category)

codes = 
  fgroup_by(codes, joined_hosp_id) |>
  fsummarize(initial_code_status = ffirst(code_status_category)) |>
  fmutate(initial_code_status = if_else(tolower(initial_code_status) == "dnr", "Special/Partial", initial_code_status))

cohort = join(cohort, codes, how = "left", multiple = F)

rm(codes); gc()

# outcomes ---------------------------------------------------------------------

## death -----------------------------------------------------------------------

death = 
  fsubset(cohort, dead_01 == 1) |>
  select(joined_hosp_id, event_dttm = discharge_dttm) |>
  fmutate(event = "death") 

## hospice ---------------------------------------------------------------------

hospice = 
  fsubset(cohort, hospice_01 == 1) |>
  select(joined_hosp_id, event_dttm = discharge_dttm) |>
  fmutate(event = "hospice") 

## vasopressors ----------------------------------------------------------------

va_list = c(
  "norepinephrine",
  "vasopressin",
  "phenylephrine",
  "epinephrine",
  "dopamine",
  "angiotensin_2",
  "dobutamine",
  "milrinone"
)

meds = 
  dplyr::filter(data_list$medication_admin_continuous, tolower(med_category) %in% va_list) |>
  dplyr::select(hospitalization_id, admin_dttm, ends_with("category")) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(joined_hosp_id, admin_dttm, med_category, med_dose, mar_action_category) |>
  funique()

meds = 
  join(meds, date_frame, how = "inner", multiple = T) |>
  fsubset(admin_dttm >= admission_dttm) |>
  fsubset(admin_dttm <= discharge_dttm) |>
  fmutate(
    med_category = tolower(med_category),
    event = case_when(
      med_category == "dobutamine" ~ "inotrope",
      med_category == "milrinone"  ~ "inotrope",
      TRUE                         ~ "vasopressor"
    )
  )|>
  select(joined_hosp_id, event_dttm = admin_dttm, event, med_category) |> # needs adjusting to include dosage
  funique()

rm(va_list); gc()

### respiratory support --------------------------------------------------------

resp_list = c(
  "imv",
  "nippv",
  "high flow nc"
)

resp = 
  dplyr::filter(data_list$respiratory_support, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(tolower(device_category) %in% resp_list) |>
  dplyr::select(hospitalization_id, recorded_dttm, device_category) |>
  dplyr::collect() |>
  join(hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(joined_hosp_id, recorded_dttm, device_category) |>
  funique()

resp = 
  join(resp, date_frame, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= admission_dttm) |>
  fsubset(recorded_dttm <= discharge_dttm) |>
  fmutate(event = tolower(device_category)) |>
  select(joined_hosp_id, event_dttm = recorded_dttm, event) |>
  funique()

present_events  = sort(funique(tolower(resp$event)))
missing_events  = setdiff(resp_list, present_events)

if (length(missing_events) > 0) {
  stop(
    sprintf("Resp need at least one %s, but %s not found.",
            if (length(missing_events) == 1) "event" else "events",
            paste(missing_events, collapse = ", ")
    ),
    call. = FALSE
  )
}

rm(resp_list, present_events, missing_events)
gc()

### combine and save -----------------------------------------------------------

## outcomes data frame ---------------------------------------------------------

df_outcomes = rowbind(events, death, hospice, meds, resp) 

fwrite(df_outcomes, here("proj_tables", "outcome_times.csv"))

rm(df_outcomes, death, hospice, events); gc()

rowbind(meds, resp, fill = T) |> 
  write_parquet(here("proj_tables", "careprocess.parquet"))

## cohort (1 row per encounter) ------------------------------------------------

cohort = 
  funique(cohort) |>
  fmutate(wicu_01 = if_else(joined_hosp_id %in% ward_icu_tx, 1L, 0L, 0L)) |>
  fmutate(icu_01  = if_else(joined_hosp_id %in% icu_encs,    1L, 0L, 0L)) |>
  fmutate(imv_01  = if_else(joined_hosp_id %in% imv_encs,    1L, 0L, 0L)) |>
  fmutate(va_01   = if_else(joined_hosp_id %in% va_encs,     1L, 0L, 0L)) |>
  select(
    patient_id, 
    joined_hosp_id, 
    admission_dttm, 
    discharge_dttm,
    first_ward_dttm,
    age, 
    race_category, 
    ethnicity_category, 
    ends_with("01"),
    initial_code_status,
    los_hosp_d
  ) |>
  mutate(across(
    .cols = ends_with("category"),
    .fns  = ~if_else(is.na(.x), "unknown", tolower(.x))
  )) |> 
  mutate(across(
    .cols = ends_with("01"),
    .fns  = ~if_else(is.na(.x), 0L, .x)
  ))

## add elixhauser --------------------------------------------------------------

### vector common codes to avoid for RAM's sake (to repl w/ keep list) ---------

unused_vect = c(
  "Z79.899",
  "E78.5",
  "Z87.891",
  "K21.9",
  "Z20.822",
  "F17.210",
  "Z79.01",
  "N17.9"
)

### only relevant codes from relevant encounters -------------------------------

elix = 
  dplyr::filter(data_list$hospital_diagnosis, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(!toupper(diagnosis_code) %in% unused_vect) |>
  dplyr::select(hospitalization_id, poa_present, diagnosis_code) |>
  dplyr::collect() |>
  funique()

### assign elixhauser diagnosis dummies ----------------------------------------

elix = 
  comorbidity::comorbidity(
    elix, 
    id      = "hospitalization_id", 
    code    = "diagnosis_code", 
    map     = "elixhauser_icd10_quan", 
    assign0 = T
  )

### link to joined-hosp-id -----------------------------------------------------

elix = 
  join(elix, hid_jid_crosswalk, how = "left", multiple = T) |>
  fselect(-patient_id, -hospitalization_id) |>
  fgroup_by(joined_hosp_id) |>
  fmax()

### calculate vw scores --------------------------------------------------------

vw        = comorbidity::score(elix, weights = "vw", assign0 = T)
elix      = cbind(elix, vw = vw) |> fselect(joined_hosp_id, vw)
cohort    = join(cohort, elix, how = "left", multiple = F)
cohort$vw = if_else(is.na(cohort$vw), 0L, cohort$vw)

rm(elix, vw); gc()

## add hospital_id -------------------------------------------------------------

hospital = 
  dplyr::select(data_list$adt, hospitalization_id, in_dttm, hospital_id) |>
  dplyr::filter(hospitalization_id %in% cohort_hids) |>
  dplyr::collect()

hospital = 
  join(hospital, hid_jid_crosswalk, how = "left", multiple = T) |>
  roworder(in_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(hospital_id = ffirst(hospital_id))

cohort = join(cohort, hospital, how = "left", multiple = F)

rm(hospital); gc()

## sanity check before saving --------------------------------------------------

props = 
  tidytable(
    icu     = fmean(cohort$icu_01,     na.rm=TRUE),
    dead    = fmean(cohort$dead_01,    na.rm=TRUE),
    hospice = fmean(cohort$hospice_01, na.rm=TRUE),
    imv     = fmean(cohort$imv_01,     na.rm=TRUE),
    va      = fmean(cohort$va_01,      na.rm=TRUE)
  )

if (
  props$icu     > 0.50 | props$icu == 0     |
  props$dead    > 0.20 | props$dead == 0    |
  props$hospice > 0.10 | props$hospice == 0 |
  props$imv     > 0.40 | props$imv == 0     |
  props$va      > 0.40 | props$va == 0) {
  stop("Sanity check failed: Outcome distribution out of expected range.")
}

write_parquet(cohort,            here("proj_tables", "cohort.parquet"))
write_parquet(hid_jid_crosswalk, here("proj_tables", "hid_jid_crosswalk.parquet"))

# go to 02
