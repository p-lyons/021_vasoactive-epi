

# t02. characteristics/outcomes by cancer status (0 = none, 1 = cancer) --------

## prepare table component = continuous variables ------------------------------

t2_cont = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fgroup_by(ca_01) |>
  fsummarize(
    n         = fnobs(joined_hosp_id),
    age_sum   = fsum(age),
    age_sumsq = fsum(age^2),
    age_p025  = fquantile(age, 0.025),
    age_p975  = fquantile(age, 0.975),
    vw_sum    = fsum(vw),
    vw_sumsq  = fsum(vw^2),
    los_sum   = fsum(los_hosp_d),
    los_sumsq = fsum(los_hosp_d^2),
    los_p025  = fquantile(los_hosp_d, 0.025),
    los_p975  = fquantile(los_hosp_d, 0.975)
  ) |>
  ftransform(var_type = "continuous") |>
  ftransform(site     = paste0(site_lowercase))

## prepare table component = categorized continuous variables ------------------

### age ------------------------------------------------------------------------

age_breaks = c( 18,      40,      50,      60,      70,      80,  Inf)
age_labs   = c("18_39", "40_49", "50_59", "60_69", "70_79", "80_plus")

ages_cat = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(ca_01, age) |>
  fmutate(a = cut(age, breaks = age_breaks, labels = age_labs, right = F)) |>
  fgroup_by(ca_01, a) |>
  fnobs() |>
  select(ca_01, age_cat = a, n = age) |>
  ftransform(age_cat = paste0("age_", age_cat)) |>
  ftransform(var = "age", category = str_remove(age_cat, "age_")) |>
  fselect(ca_01, var, category, n) 

### elixhauser -----------------------------------------------------------------

elix_breaks = c(-Inf,  0, 4,  9, 14,  Inf)
elix_labs   = c("<= 0", "1-4", "5-9", "10-14", ">= 15")

elix_cat = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(ca_01, vw) |>
  fmutate(a = cut(vw, breaks = elix_breaks, labels = elix_labs, right = F)) |>
  fgroup_by(ca_01, a) |>
  fnobs() |>
  select(ca_01, elix_cat = a, n = vw) |>
  ftransform(elix_cat = paste0("vw_", elix_cat)) |>
  ftransform(var = "vw", category = str_remove(elix_cat, "vw_")) |>
  fselect(ca_01, var, category, n) 

### los (days) -----------------------------------------------------------------

l_breaks = c( 0,        2,         4,         7,         14, Inf)
l_labs   = c("0-47h", "48h_96h", "96h_1wk", "1wk_2wk", "2wk_plus")

los_cat = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(ca_01, los_hosp_d) |>
  fmutate(los_cat = cut(los_hosp_d, breaks = l_breaks, labels = l_labs, right = F)) |>
  fgroup_by(ca_01, los_cat) |>
  fnobs() |>
  select(ca_01, los_cat, n = los_hosp_d) |>
  ftransform(var = "los", category = los_cat) |>
  fselect(ca_01, var, category, n) 

## prepare table component = categorical variables -----------------------------

t2_cat = 
  fsubset(cohort, ed_admit_01 == 1) |>
  select(-ends_with("id"), -ends_with("dttm"),  -age, -vw, -los_hosp_d) |>
  pivot_longer(-ca_01, names_to = "var", values_to = "val") |>
  fsubset(!is.na(val)) |>
  fmutate(n = val) |>
  fgroup_by(ca_01, var, val) |>
  fnobs() |>
  ftransform(var      = str_remove(var, "_01")) |>
  ftransform(category = tolower(str_replace_all(as.character(val), "-", "_"))) |>
  fselect(ca_01, var, category, n) 

## quality control -------------------------------------------------------------

### check for small n in cells -------------------------------------------------

ages_cat$site = site_lowercase
elix_cat$site = site_lowercase
los_cat$site  = site_lowercase
t2_cat$site   = site_lowercase

all_cat = rowbind(ages_cat, elix_cat, los_cat, t2_cat)
small_c = fsubset(all_cat, n < 5 & n > 0)

if (nrow(small_c) > 0) {
  msg = paste0(
    "ERROR: ", nrow(small_c), " cells have n < 5:\n",
    capture.output(print(small_c |> fselect(var, category, ca_01, n))) |> 
      paste(collapse = "\n")
  )
  stop(msg)
}

### check for sample size mismatch in continuous table -------------------------

if (sum(t2_cont$n) != nrow(cohort[ed_admit_01 == 1])) {
  stop("ERROR: Sample size mismatch! Sum of t2_cont$n != nrow(df)")
}

## export table 2 --------------------------------------------------------------

fwrite(all_cat, here("proj_output", paste0("table_02_cat_",  site_lowercase, ".csv")))
fwrite(t2_cont, here("proj_output", paste0("table_02_cont_", site_lowercase, ".csv")))

# final cleanup ----------------------------------------------------------------

keep = c(
  "data_list",
  "site_lowercase",
  "cohort",
  "cohort_hids",
  "cohort_jids",
  "cohort_pats",
  "hid_jid_crosswalk",
  "ward_times",
  "req_vitals",
  "req_labs"
)

rm(list = setdiff(ls(), keep)); gc()
