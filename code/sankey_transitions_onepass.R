# ==============================================================================
# sankey_transitions_onepass.R
# Vasopressor Escalation in Refractory Distributive Shock - CLIF Consortium
# SITE-SIDE EXPORT: CLIF 2.1 tables -> block-to-block state transition counts
#
# WHAT THIS SCRIPT DOES (one pass, no other pipeline script needed):
#   1. Builds the primary cohort exactly as 01_cohort.R STEPS 1-6 do
#      (linkage, ED/ICU, psych/rehab, disposition, YODO, uninterpretable MAR,
#      T0, third-line-before-T0 exclusion). Exclusion step names match the
#      01_cohort.R cascade.
#   2. Assigns the Sec 6.3.1 state at each 6h block boundary (boundary
#      snapshot + within-block look-back). State assignment exists only in
#      this script.
#   3. Counts transitions (block_from, state_from, state_to) across encounters.
#   4. Counts complete paths: state at 0, 12, 24, 36 h + outcome at 48 h.
#      Each 12-h state = the more severe of the two 6-h block states in the
#      preceding 12 h (D if dead/hospice by that hour). Patients discharged
#      alive stay W. 48-h outcome: dead_48h, alive_48h_died_later,
#      alive_48h_survived (hospice discharge counts as death throughout).
#
# OUTPUTS (under upload_to_box/, with the Table 1 files) -- aggregate counts only:
#   sankey_transitions_<site>.csv   block_from | state_from | state_to | n | site
#   sankey_paths_<site>.csv         state_0h | state_12h | state_24h | state_36h | outcome_48h | n | site
#   sankey_site_summary_<site>.csv  cohort size, exclusions, settings
#   No patient-level rows, identifiers, or timestamps leave this script.
#
# SHARING: counts are UNSUPPRESSED. Transition counts contain no demographic
#   variables, and all files are pooled at the coordinating site under the
#   consortium DUA before anything is shared.
#
# DELIBERATE DEVIATIONS FROM PIPELINE CONVENTIONS:
#   - Standalone: bootstraps its own packages and reads config itself rather
#     than chaining off 01's session, so a site can run it with nothing else.
#   - Cohort logic is a PORT of 01_cohort.R, not a source() of it. If the
#     cohort definition changes in 01, mirror the change here.
#   - Drops everything that does not affect cohort membership or state:
#     comorbidity, demographics, SVI/ADI, IMV/CRRT, NEE, inotropes.
#   - Encounters with ANY unresolved boundary (assign_state() returns NA) are
#     dropped from the counts and tallied in the site summary. With the
#     revised S2a rule, every active-agent combination resolves, so this count
#     should be 0.
#
# CONFIG: config/config_clif_pressors.yaml (same file as 01_cohort.R). Uses
#   clif_data_location, file_type, site_lowercase, time_zone.
# ==============================================================================

# packages ---------------------------------------------------------------------

pkgs = c("data.table", "arrow", "dplyr", "yaml", "here")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try(install.packages(p, dependencies = TRUE), silent = TRUE)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

rm(pkgs, p)

# config -----------------------------------------------------------------------

config_path = file.path(here::here(), "config", "config_clif_pressors.yaml")

if (!file.exists(config_path)) {
  stop(sprintf("Config file not found: %s", config_path), call. = FALSE)
}

config = yaml::read_yaml(config_path)

fn_cfg = function(k) {
  v = config[[k]]
  if (is.null(v) || identical(v, "")) {
    stop(sprintf("config is missing '%s'", k), call. = FALSE)
  }
  v
}

tables_location  = normalizePath(path.expand(fn_cfg("clif_data_location")), mustWork = FALSE)
file_type        = tolower(fn_cfg("file_type"))
site_lowercase   = tolower(fn_cfg("site_lowercase"))
site_tz          = fn_cfg("time_zone")

if (!site_tz %in% OlsonNames()) {
  stop(
    sprintf("config 'time_zone' is not a valid IANA time zone: '%s'", site_tz),
    call. = FALSE
  )
}

# same folder as the Table 1 files, so sites send one folder
dir_out = here::here("upload_to_box")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

# open CLIF tables (root + 2 levels, case-insensitive) -------------------------

fn_open_clif = function(nm) {
  ext  = if (file_type == "csv") "csv" else "parquet"
  lvl0 = list.files(tables_location, full.names = TRUE)
  lvl1 = unlist(lapply(lvl0[dir.exists(lvl0)], list.files, full.names = TRUE))
  ents = c(lvl0, lvl1)
  stem = sub(paste0("\\.", ext, "$"), "", basename(ents), ignore.case = TRUE)
  hit  = ents[grepl(sprintf("^(clif_)?%s$", nm), stem, ignore.case = TRUE)]
  if (length(hit) == 0) {
    stop(sprintf("CLIF table '%s' not found under %s", nm, tables_location), call. = FALSE)
  }
  hit = hit[order(lengths(strsplit(hit, .Platform$file.sep)))][1]
  if (file_type == "csv") arrow::read_csv_arrow(hit) else arrow::open_dataset(hit)
}

hospitalization             = fn_open_clif("hospitalization")
adt                         = fn_open_clif("adt")
vitals                      = fn_open_clif("vitals")
medication_admin_continuous = fn_open_clif("medication_admin_continuous")

# ------------------------------------------------------------------------------
# Constants (must match 01_cohort.R)
# ------------------------------------------------------------------------------

start_date           = as.POSIXct("2016-01-01 00:00:00", tz = site_tz)
end_date             = as.POSIXct("2024-12-31 23:59:59", tz = site_tz)
cf_window_h          = 4      # carry-forward, Sec 5.4.2
ne_threshold         = 0.2    # NE mcg/kg/min for T0, Sec 5.4
link_hours           = 6L     # contiguous-encounter gap, Sec 5.3
third_line_grace_min = 0      # 0 = primary; 60 = sensitivity arm

# implausible-dose limits after unit conversion (records above are deleted
# before carry-forward, so the previous valid dose carries forward)
dose_limits = c(
  norepinephrine = 5,      # mcg/kg/min
  epinephrine    = 5,      # mcg/kg/min
  vasopressin    = 0.5,    # units/min
  phenylephrine  = 10,     # mcg/kg/min
  dopamine       = 50,     # mcg/kg/min
  angiotensin    = 200     # ng/kg/min
)
block_hours          = 6L     # Sec 6.2, fixed by protocol
window_hours         = 48L    # tunable (48 or 72)
path_hours           = c(12L, 24L, 36L)   # path states after T0 (0 h added below)
path_outcome_hours   = 48L                # path outcome time point

scoped_pressors = c(
  "norepinephrine",
  "vasopressin",
  "epinephrine",
  "phenylephrine",
  "dopamine",
  "angiotensin"
)

ang2_aliases = c(
  "angiotensin_ii",
  "angiotensin ii",
  "angiotensin_2",
  "angiotensin2",
  "ang2",
  "angii"
)

t0_block_agents   = c("epinephrine", "phenylephrine", "dopamine", "angiotensin")
adrenergic_agents = c("norepinephrine", "epinephrine", "phenylephrine", "dopamine")
adr_escalation    = c("epinephrine", "phenylephrine", "dopamine")
state_levels      = c("W", "S0", "S1", "S2a", "S2b", "S3", "D")
n_blocks          = as.integer(window_hours / block_hours)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

cascade = data.table(step = character(), n_remaining = integer())

fn_log = function(step, n) {
  cascade <<- rbind(cascade, data.table(step = step, n_remaining = as.integer(n)))
  message(sprintf("  [%s] %d", step, n))
}

# 4h carry-forward for ONE agent at arbitrary evaluation points.
# pts: joined_hosp_id, eval_dttm (unique). Returns pts + <agent>_dose/_active,
# keyed by (joined_hosp_id, eval_dttm) so callers MERGE rather than cbind.
fn_active_at = function(rec, pts, agent) {
  dcol = paste0(agent, "_dose")
  acol = paste0(agent, "_active")
  out  = copy(pts)
  
  if (nrow(rec) == 0L) {
    return(out[, (dcol) := NA_real_][, (acol) := FALSE][])
  }
  
  r = copy(rec)[, rec_dttm := admin_dttm]
  j = r[
    out,
    on   = .(joined_hosp_id, admin_dttm = eval_dttm),
    roll = TRUE
  ]
  setnames(j, "admin_dttm", "eval_dttm")
  
  j[, gap_h := as.numeric(difftime(eval_dttm, rec_dttm, units = "hours"))]
  j[, (acol) := !is.na(dose) & !is.na(gap_h) & gap_h <= cf_window_h &
                dose > 0 & (is.na(is_stop) | !is_stop)]
  j[, (dcol) := fifelse(get(acol), dose, NA_real_)]
  
  j[, c("joined_hosp_id", "eval_dttm", dcol, acol), with = FALSE]
}

fn_snapshot = function(mac, pts) {
  out = copy(pts)
  
  for (a in scoped_pressors) {
    rec = mac[agent == a, .(joined_hosp_id, admin_dttm, dose = med_dose, is_stop)]
    out = merge(
      out,
      fn_active_at(rec, unique(pts[, .(joined_hosp_id, eval_dttm)]), a),
      by    = c("joined_hosp_id", "eval_dttm"),
      all.x = TRUE
    )
  }
  
  for (a in scoped_pressors) {
    ac = paste0(a, "_active")
    out[is.na(get(ac)), (ac) := FALSE]
  }
  
  out
}

# Sec 6.3.1 hierarchy. S2a (protocol v8): two or more adrenergic agents active
# at once, checked after S3 and S2b. With this rule every active-agent
# combination maps to a state.
assign_state = function(dt) {
  act = function(a) dt[[paste0(a, "_active")]]
  
  n_total     = Reduce(`+`, lapply(scoped_pressors,   function(a) as.integer(act(a))))
  n_adr       = Reduce(`+`, lapply(adrenergic_agents, function(a) as.integer(act(a))))
  has_adr_esc = Reduce(`|`, lapply(adr_escalation,    function(a) act(a)))
  ang2        = act("angiotensin")
  vp          = act("vasopressin")
  is_dead     = !is.na(dt$d_instant) & dt$d_instant <= dt$eval_dttm
  
  fcase(
    is_dead,                                 "D",
    n_total == 0L,                           "W",
    (n_total >= 4L) | (has_adr_esc & ang2),  "S3",
    ang2 & (n_total >= 2L),                  "S2b",
    n_adr >= 2L,                             "S2a",
    vp & (n_adr == 1L) & (n_total == 2L),    "S0",
    n_total == 1L,                           "S1",
    default = NA_character_
  )
}

# ==============================================================================
# STEP 1: Adult inpatients in study period; link contiguous stays (<=6h gap)
# ==============================================================================

message("\n== Cohort ==")

hosp_raw = 
  dplyr::filter(hospitalization, age_at_admission >= 18) |>
  dplyr::filter(admission_dttm >= start_date & admission_dttm <= end_date) |>
  dplyr::filter(admission_dttm < discharge_dttm & !is.na(discharge_dttm)) |>
  dplyr::select(
    patient_id,
    hospitalization_id,
    admission_dttm,
    discharge_dttm,
    age_at_admission,
    discharge_category
  ) |>
  dplyr::collect()

hosp = as.data.table(hosp_raw)
rm(hosp_raw)

if (anyDuplicated(hosp$hospitalization_id)) {
  stop("Source has duplicate hospitalization_id.", call. = FALSE)
}

setorder(hosp, patient_id, admission_dttm)

hosp[, prev_gap := as.numeric(difftime(admission_dttm, shift(discharge_dttm), units = "hours")), by = patient_id]
hosp[, joined_hosp_id := .GRP, by = .(patient_id, cumsum(is.na(prev_gap) | prev_gap >= link_hours))]

fn_log("01_adult_inpatients_study_period", uniqueN(hosp$joined_hosp_id))

# ==============================================================================
# STEP 2: Require ED or ICU; exclude psych/rehab
# ==============================================================================

adt_loc = 
  dplyr::select(adt, hospitalization_id, location_category) |>
  dplyr::collect() |>
  as.data.table()

adt_loc[, loc := tolower(location_category)]

keep_jid = hosp[
  hospitalization_id %in% adt_loc[loc %chin% c("ed", "icu"), hospitalization_id],
  unique(joined_hosp_id)
]

hosp = hosp[joined_hosp_id %in% keep_jid]
fn_log("02_required_ed_or_icu", uniqueN(hosp$joined_hosp_id))

drop_jid = hosp[
  hospitalization_id %in% adt_loc[loc %chin% c("psych", "rehab"), hospitalization_id],
  unique(joined_hosp_id)
]

hosp = hosp[!joined_hosp_id %in% drop_jid]
fn_log("03_excluded_psych_rehab", uniqueN(hosp$joined_hosp_id))

rm(adt_loc, keep_jid, drop_jid)

# ==============================================================================
# STEP 3: Encounter frame; disposition; duplicate deaths; post-death encounters
# ==============================================================================
# first/last use the first/last NON-missing value, matching collapse::ffirst/
# flast defaults used in 01_cohort.R.

fn_first = function(x) {
  x = x[!is.na(x)]
  if (length(x)) x[1] else x[NA_integer_]
}

fn_last = function(x) {
  x = x[!is.na(x)]
  if (length(x)) x[length(x)] else x[NA_integer_]
}

setorder(hosp, admission_dttm)

cohort = hosp[
  ,
  .(
    admission_dttm     = fn_first(admission_dttm),
    discharge_dttm     = fn_last(discharge_dttm),
    discharge_category = fn_last(discharge_category)
  ),
  by = .(patient_id, joined_hosp_id)
]

cohort[, dispo := tolower(trimws(as.character(discharge_category)))]
cohort = cohort[!is.na(dispo) & dispo != ""]
fn_log("04_missing_discharge_category", nrow(cohort))

setorder(cohort, admission_dttm, discharge_dttm)

dd = cohort[dispo == "expired"][, `:=`(n_deaths = .N, counter = seq_len(.N)), by = patient_id]
cohort = cohort[!joined_hosp_id %in% dd[n_deaths > 1 & counter == 1L, joined_hosp_id]]

death_times = cohort[dispo == "expired", .(death_instant = min(discharge_dttm)), by = patient_id]
post_death  = merge(cohort, death_times, by = "patient_id")[admission_dttm >= death_instant, joined_hosp_id]
cohort      = cohort[!joined_hosp_id %in% post_death]
fn_log("04b_yodo_cleanup", nrow(cohort))

hid_jid = hosp[joined_hosp_id %in% cohort$joined_hosp_id, .(hospitalization_id, joined_hosp_id)]

rm(dd, death_times, post_death)

# ==============================================================================
# STEP 4: First in-encounter weight (dose-unit correction)
# ==============================================================================

w_raw = 
  dplyr::filter(vitals, vital_category == "weight_kg") |>
  dplyr::filter(hospitalization_id %in% !!hid_jid$hospitalization_id) |>
  dplyr::select(hospitalization_id, recorded_dttm, weight_kg = vital_value) |>
  dplyr::collect()

w = unique(as.data.table(w_raw))
rm(w_raw)

w = merge(w, hid_jid, by = "hospitalization_id")
w = merge(w, cohort[, .(joined_hosp_id, admission_dttm, discharge_dttm)], by = "joined_hosp_id")
w = w[recorded_dttm >= admission_dttm & recorded_dttm <= discharge_dttm]

setorder(w, recorded_dttm)
w = w[, .(weight_kg = weight_kg[1L]), by = joined_hosp_id]

# ==============================================================================
# STEP 5: Continuous vasopressor MAR: extract, normalize, clean
# ==============================================================================

raw_match = c(scoped_pressors, ang2_aliases)

mac_cols = c(
  "hospitalization_id",
  "admin_dttm",
  "med_category",
  "med_dose",
  "med_dose_unit",
  "mar_action_category"
)

mac = tryCatch(
  {
    dplyr::filter(medication_admin_continuous, tolower(med_category) %in% raw_match) |>
      dplyr::select(dplyr::all_of(mac_cols)) |>
      dplyr::collect()
  },
  error = function(e) {
    # list<string> med_category at some sites has no Arrow string kernel
    message("  NOTE: lazy med_category filter failed; collecting then flattening.")
    raw = dplyr::select(medication_admin_continuous, dplyr::all_of(mac_cols)) |>
      dplyr::collect()
    if (is.list(raw$med_category)) {
      raw$med_category = vapply(
        raw$med_category,
        function(z) if (length(z) == 0 || is.null(z)) NA_character_ else as.character(z[[1]]),
        character(1)
      )
    }
    raw[tolower(raw$med_category) %in% raw_match, ]
  }
)

mac = unique(as.data.table(mac))
mac[, agent := tolower(trimws(as.character(med_category)))]
mac[agent %chin% ang2_aliases, agent := "angiotensin"]
mac = mac[agent %chin% scoped_pressors]

mac = merge(mac, hid_jid, by = "hospitalization_id")
mac = merge(mac, cohort[, .(joined_hosp_id, admission_dttm, discharge_dttm)], by = "joined_hosp_id")
mac = mac[admin_dttm >= admission_dttm & admin_dttm <= discharge_dttm]

## weight and unit correction (mirrors 01_cohort.R) ----------------------------
## output units: catecholamines mcg/kg/min, vasopressin units/min,
## angiotensin ng/kg/min

mac = merge(mac, w, by = "joined_hosp_id", all.x = TRUE)
mac[is.na(weight_kg), weight_kg := 70]
mac[, med_dose_unit := tolower(trimws(med_dose_unit))]

mac[, dose_mult := fcase(
  med_dose_unit %chin% c("mcg/min", "ng/min"), 1 / weight_kg,
  med_dose_unit == "units/kg/min",             weight_kg,
  med_dose_unit == "units/hr",                 1 / 60,
  default = 1
)]

mac[, med_dose := med_dose * dose_mult]

## stop actions and one record per agent-timestamp -----------------------------

mac[, is_stop := !is.na(mar_action_category) & tolower(mar_action_category) == "stopped"]
mac[is_stop == TRUE, med_dose := 0]

## implausible doses: delete record (never a stop record; stops are 0) ---------

mac[, dose_limit := dose_limits[agent]]

implausible = mac[
  !is.na(med_dose) & med_dose > dose_limit,
  .(n_records = .N, n_encounters = uniqueN(joined_hosp_id)),
  by = agent
]

n_implausible = sum(implausible$n_records)

message(sprintf("  Implausible dose records deleted: %d", n_implausible))
if (n_implausible > 0) print(implausible)

mac = mac[is.na(med_dose) | med_dose <= dose_limit]
mac[, dose_limit := NULL]

mac = mac[
  ,
  .(
    med_dose = suppressWarnings(max(med_dose, na.rm = TRUE)),
    is_stop  = all(is_stop)
  ),
  by = .(joined_hosp_id, agent, admin_dttm)
]

mac[!is.finite(med_dose), med_dose := NA_real_]

## exclude encounters with no usable NE or VP dose -----------------------------

bad_mar = mac[
  agent %chin% c("norepinephrine", "vasopressin"),
  .(usable = any(is.finite(med_dose) & med_dose > 0)),
  by = joined_hosp_id
][usable == FALSE, joined_hosp_id]

cohort = cohort[!joined_hosp_id %in% bad_mar]
mac    = mac[joined_hosp_id %in% cohort$joined_hosp_id]
fn_log("05_uninterpretable_mar", nrow(cohort))

rm(bad_mar)

# ==============================================================================
# STEP 6: T0 (NE >= 0.2 + VP, no other pressor) + third-line-before-T0
# ==============================================================================

grid = unique(mac[, .(joined_hosp_id, eval_dttm = admin_dttm)])
wide = fn_snapshot(mac, grid)

wide[, any_block := Reduce(`|`, lapply(t0_block_agents, function(a) get(paste0(a, "_active"))))]
wide[, t0_ok := norepinephrine_active & norepinephrine_dose >= ne_threshold &
                vasopressin_active & !any_block]

setorder(wide, joined_hosp_id, eval_dttm)

t0_tab = wide[t0_ok == TRUE, .(t0_dttm = eval_dttm[1L]), by = joined_hosp_id]
cohort = merge(cohort, t0_tab, by = "joined_hosp_id")
fn_log("06_met_t0_criteria", nrow(cohort))

wide = merge(wide, t0_tab, by = "joined_hosp_id")

pre_t0_block = wide[
  eval_dttm < t0_dttm - third_line_grace_min * 60 & any_block,
  unique(joined_hosp_id)
]

cohort = cohort[!joined_hosp_id %in% pre_t0_block]
mac    = mac[joined_hosp_id %in% cohort$joined_hosp_id]
fn_log("07_third_line_prior_to_t0", nrow(cohort))

rm(grid, wide, t0_tab, pre_t0_block, w, hosp)

# ==============================================================================
# STEP 7: Boundary snapshot states
# ==============================================================================

message("\n== Block states ==")

cohort[, d_instant := fifelse(
  dispo %chin% c("expired", "hospice"),
  discharge_dttm,
  as.POSIXct(NA, tz = site_tz)
)]

bnd = cohort[, .(block_idx = 0:n_blocks), by = .(joined_hosp_id, t0_dttm)]
bnd[, eval_dttm := t0_dttm + block_idx * block_hours * 3600]

states = fn_snapshot(mac, bnd)
states = merge(states, cohort[, .(joined_hosp_id, d_instant)], by = "joined_hosp_id")
states[, state := assign_state(states)]

# ==============================================================================
# STEP 8: Within-block look-back
# ==============================================================================
# Block k covers (boundary_{k-1}, boundary_k]. Evaluate the snapshot at every
# in-window MAR timestamp, take the highest-severity non-D state in the block,
# and upgrade the boundary state if it is higher (non-D, k > 0 only).

state_rank = c(W = 0L, S1 = 1L, S0 = 2L, S2a = 3L, S2b = 4L, S3 = 5L)

evt = merge(
  unique(mac[, .(joined_hosp_id, eval_dttm = admin_dttm)]),
  cohort[, .(joined_hosp_id, t0_dttm, d_instant)],
  by = "joined_hosp_id"
)

evt[, hrs := as.numeric(difftime(eval_dttm, t0_dttm, units = "hours"))]
evt = evt[hrs > 0 & hrs <= window_hours]
evt[, block_idx := as.integer(ceiling(hrs / block_hours))]

n_upgrade = 0L

if (nrow(evt) > 0L) {
  
  ev = fn_snapshot(mac, evt)   # merged by (encounter, timestamp)
  ev[, rank_evt := state_rank[assign_state(ev)]]
  
  blk_max = ev[
    !is.na(rank_evt),
    .(rank_within = max(rank_evt)),
    by = .(joined_hosp_id, block_idx)
  ]
  
  states = merge(
    states,
    blk_max,
    by    = c("joined_hosp_id", "block_idx"),
    all.x = TRUE
  )
  
  up = states$block_idx > 0L &
    !is.na(states$state) &
    states$state != "D" &
    !is.na(states$rank_within) &
    states$rank_within > state_rank[states$state]
  
  up[is.na(up)] = FALSE
  n_upgrade     = sum(up)
  
  states[up, state := names(state_rank)[match(rank_within, state_rank)]]
  states[, rank_within := NULL]
}

message(sprintf("  %d boundary state(s) upgraded by within-block look-back", n_upgrade))

n_not_s0 = states[block_idx == 0L & (is.na(state) | state != "S0"), .N]

if (n_not_s0 > 0) {
  warning(sprintf("%d encounter(s) not S0 at T0; inspect before sharing.", n_not_s0), call. = FALSE)
}

# ==============================================================================
# STEP 9: Transition counts
# ==============================================================================

message("\n== Transitions ==")

unres = states[is.na(state), unique(joined_hosp_id)]
st    = states[!joined_hosp_id %in% unres]
setorder(st, joined_hosp_id, block_idx)

trans_enc = st[
  ,
  .(
    block_from = block_idx[-.N],
    state_from = state[-.N],
    state_to   = state[-1L]
  ),
  by = joined_hosp_id
]

trans = trans_enc[, .(n = .N), by = .(block_from, state_from, state_to)]
setorder(trans, block_from, state_from, state_to)
trans[, site := site_lowercase]

rm(trans_enc)

# conservation: inflow to (k, s) must equal outflow from (k, s) for 0 < k < n_blocks
inflow  = trans[, .(n_in  = sum(n)), by = .(k = block_from + 1L, s = state_to)]
outflow = trans[, .(n_out = sum(n)), by = .(k = block_from,      s = state_from)]

cons = merge(inflow, outflow, by = c("k", "s"), all = TRUE)[k > 0L & k < n_blocks]
cons[is.na(n_in),  n_in  := 0L]
cons[is.na(n_out), n_out := 0L]

if (cons[n_in != n_out, .N] > 0) {
  stop("Flow conservation failed; inspect `cons`.", call. = FALSE)
}

# ==============================================================================
# STEP 10: Complete paths (0, 12, 24, 36 h states + 48 h outcome)
# ==============================================================================

message("\n== Paths ==")

st[, rank := state_rank[state]]   # D -> NA

## state at T0 -----------------------------------------------------------------

path_parts = list(
  st[block_idx == 0L, .(joined_hosp_id, col = "state_0h", state)]
)

## 12-h states: more severe of the two 6-h block states; D if dead by then ------

for (h in path_hours) {
  
  k = as.integer(h / block_hours)
  
  win = st[block_idx %in% c(k - 1L, k)]
  
  part = win[, .(
    is_dead  = any(block_idx == k & state == "D"),
    max_rank = if (all(is.na(rank))) NA_integer_ else max(rank, na.rm = TRUE)
  ), by = joined_hosp_id]
  
  part[, state := fifelse(
    is_dead,
    "D",
    names(state_rank)[match(max_rank, state_rank)]
  )]
  
  path_parts[[length(path_parts) + 1L]] = part[, .(joined_hosp_id, col = paste0("state_", h, "h"), state)]
}

rm(h, k, win, part)

## 48-h outcome ----------------------------------------------------------------

outcome_part = cohort[
  joined_hosp_id %in% st$joined_hosp_id,
  .(
    joined_hosp_id,
    col   = "outcome_48h",
    state = fcase(
      !is.na(d_instant) & d_instant <= t0_dttm + path_outcome_hours * 3600, "dead_48h",
      !is.na(d_instant),                                                    "alive_48h_died_later",
      default = "alive_48h_survived"
    )
  )
]

path_parts[[length(path_parts) + 1L]] = outcome_part

path_long = rbindlist(path_parts, use.names = TRUE)

path_wide = dcast(
  path_long,
  joined_hosp_id ~ col,
  value.var = "state"
)

path_cols = c(
  "state_0h",
  paste0("state_", path_hours, "h"),
  "outcome_48h"
)

paths = path_wide[, .(n = .N), by = path_cols]
setorderv(paths, "n", order = -1L)
paths[, site := site_lowercase]

if (sum(paths$n) != uniqueN(st$joined_hosp_id) || anyNA(paths[, ..path_cols])) {
  stop("Path counts do not cover every encounter exactly once; inspect `path_wide`.", call. = FALSE)
}

message(sprintf("  %d encounters on %d unique paths", sum(paths$n), nrow(paths)))

st[, rank := NULL]
rm(path_parts, outcome_part, path_long, path_wide)

# ==============================================================================
# STEP 11: Site summary and save
# ==============================================================================

summary_tab = data.table(
  site                  = site_lowercase,
  n_cohort              = nrow(cohort),
  n_unresolved_excluded = length(unres),
  n_in_transitions      = uniqueN(st$joined_hosp_id),
  n_implausible_doses   = n_implausible,
  n_lookback_upgrades   = n_upgrade,
  n_unique_paths        = nrow(paths),
  window_hours          = window_hours,
  block_hours           = block_hours,
  third_line_grace_min  = third_line_grace_min,
  run_date              = as.character(Sys.Date())
)

cascade_wide = as.data.table(
  as.list(setNames(cascade$n_remaining, paste0("n_", cascade$step)))
)

summary_tab = cbind(summary_tab, cascade_wide)

fwrite(trans,       file.path(dir_out, sprintf("sankey_transitions_%s.csv",  site_lowercase)))
fwrite(summary_tab, file.path(dir_out, sprintf("sankey_site_summary_%s.csv", site_lowercase)))
fwrite(paths,       file.path(dir_out, sprintf("sankey_paths_%s.csv",        site_lowercase)))

message(sprintf("  %d encounters -> %d transition cells (%d unresolved excluded)",
                summary_tab$n_in_transitions, nrow(trans), length(unres)))
message("\n== sankey_transitions_onepass.R complete ==")

# end sankey_transitions_onepass.R ---------------------------------------------
