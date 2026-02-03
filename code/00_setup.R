# ==============================================================================
# 00_setup.R
# Vasopressor Escalation in Septic Shock - CLIF Consortium
# Setup: packages, environment, data loading, validation
# ==============================================================================

# packages ---------------------------------------------------------------------

packages_to_install = c(
  "comorbidity",
  "data.table",
  "tidyverse",
  "tidytable",
  "collapse",
  "janitor",
  "arrow",
  "yaml",
  "here",
  "tableone",
  "smd"
)

packages_to_load = c(
  "data.table",
  "tidytable",
  "collapse",
  "stringr",
  "arrow",
  "here"
)

fn_install_if_missing = function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try(install.packages(p, dependencies = TRUE), silent = TRUE)
  }
}

fn_load_quiet = function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

invisible(lapply(packages_to_install, fn_install_if_missing))
invisible(lapply(packages_to_load, fn_load_quiet))
options(dplyr.summarise.inform = FALSE)

rm(packages_to_install, packages_to_load, fn_install_if_missing, fn_load_quiet); gc()

# environment ------------------------------------------------------------------

## threads and RAM -------------------------------------------------------------

os_type   = Sys.info()[["sysname"]]
all_cores = parallel::detectCores(logical = TRUE)
all_cores = if (is.na(all_cores)) 1L else as.integer(all_cores)

get_ram_gb = function() {
  tryCatch({
    if (os_type == "Darwin") {
      bytes = suppressWarnings(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)))
      if (length(bytes) > 0 && !is.na(bytes)) bytes / 1024^3 else NA_real_
    } else {
      if (file.exists("/proc/meminfo")) {
        kb = suppressWarnings(
          as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
        )
        if (length(kb) > 0 && !is.na(kb)) kb / 1024^2 else NA_real_
      } else if (requireNamespace("ps", quietly = TRUE)) {
        ps::ps_system_memory()[["available"]] / 1024^3
      } else {
        NA_real_
      }
    }
  }, error = function(e) NA_real_)
}

avail_ram_gb = get_ram_gb()

## thread settings -------------------------------------------------------------

reserve_cores  = 1L
gb_per_thread  = 0.50
max_by_cores   = max(1L, all_cores - reserve_cores)
max_by_memory  = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads      = as.integer(max(1L, min(max_by_cores, max_by_memory, 8L)))

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

message(
  sprintf(
    "Environment OK | OS=%s | Cores=%d | Threads=%d | RAM≈%s GB",
    os_type, all_cores, n_threads,
    ifelse(is.finite(avail_ram_gb), round(avail_ram_gb, 1), "NA")
  )
)

# site configuration -----------------------------------------------------------

config        = yaml::read_yaml(here("config", "config_clif_pressors.yaml"))
site_details  = fread(here("config", "clif_sites.csv"))
allowed_sites = tolower(site_details$site_name)
allowed_files = c("parquet", "csv", "fst")

site_lowercase   = tolower(config$site_lowercase)
file_type        = tolower(config$file_type)
tables_location  = config$clif_data_location
project_location = config$project_location

## validate site name ----------------------------------------------------------

if (!(site_lowercase %in% allowed_sites)) {
  stop(
    paste0("Invalid site '", site_lowercase, "'. Expected one of: ", 
           paste(allowed_sites, collapse = ", ")),
    call. = FALSE
  )
}

## validate file type ----------------------------------------------------------

if (!(file_type %in% allowed_files)) {
  stop(
    paste0("Invalid file type '", file_type, "'. Expected one of: ", 
           paste(allowed_files, collapse = ", ")),
    call. = FALSE
  )
}

## create output directories ---------------------------------------------------

if (!dir.exists(paste0(project_location, "/proj_tables"))) {
  dir.create(paste0(project_location, "/proj_tables"), recursive = TRUE)
}

if (!dir.exists(paste0(project_location, "/upload_to_box"))) {
  dir.create(paste0(project_location, "/upload_to_box"), recursive = TRUE)
}

## study dates -----------------------------------------------------------------

start_date = as.POSIXct("2016-01-01", tz = "UTC")
end_date   = as.POSIXct("2024-12-31", tz = "UTC")
today      = format(Sys.Date(), "%y%m%d")

## required vital signs for cohort inclusion -----------------------------------

req_vitals = c(
  "heart_rate",
  "sbp",
  "dbp",
  "resp_rate",
  "spo2"
)

# data loading -----------------------------------------------------------------

## required tables for this project --------------------------------------------

required_tables = c(
  "patient",
  "hospital_diagnosis",
  "hospitalization",
  "adt",
  "vitals",
  "labs",
  "medication_admin_continuous",
  "medication_admin_intermittent",
  "respiratory_support",
  "crrt_therapy",
  "code_status"
)

## check available tables ------------------------------------------------------

clif_table_filenames = list.files(
  path       = tables_location,
  pattern    = paste0("^clif_.*\\.", file_type, "$"),
  full.names = TRUE
)

clif_table_basenames = basename(clif_table_filenames) |>
  str_remove(paste0("\\.", file_type, "$")) |>
  str_remove("^clif_") |>
  str_remove("(_\\d{4}(_\\d{4})*)$")

table_file_map     = setNames(clif_table_filenames, clif_table_basenames)
missing_tables     = setdiff(required_tables, clif_table_basenames)
required_filenames = table_file_map[required_tables]

if (length(missing_tables) > 0) {
  stop(paste("Missing required tables:", paste(missing_tables, collapse = ", ")), call. = FALSE)
} else {
  message("All required tables are present.")
}

## load tables as arrow datasets -----------------------------------------------

if (file_type == "parquet") {
  data_list = lapply(required_filenames, open_dataset)
} else if (file_type == "csv") {
  data_list = lapply(required_filenames, \(f) read_csv_arrow(f))
} else if (file_type == "fst") {
  data_list = lapply(required_filenames, \(f) {
    tmp = fst::read.fst(f, as.data.table = TRUE)
    arrow_table(tmp)
  })
} else {
  stop("Unsupported file format", call. = FALSE)
}

names(data_list) = names(required_filenames)

rm(site_details, allowed_files, allowed_sites, missing_tables)
rm(clif_table_basenames, clif_table_filenames, required_filenames, table_file_map)
gc()

# validation -------------------------------------------------------------------

## case-insensitive column resolver --------------------------------------------

.resolve = function(tbl, v) {
  nm = names(tbl)
  nm[match(tolower(v), tolower(nm))]
}

## validate single table -------------------------------------------------------

validate_table = function(tbl, table_name, req_vars = NULL, req_values = list()) {
  
  if (is.null(req_vars)) req_vars = character(0)
  
  problems     = character()
  actual_names = tolower(names(tbl))
  req_lower    = tolower(req_vars)
  missing_vars = req_vars[!req_lower %in% actual_names]
  
  if (length(missing_vars)) {
    problems = c(problems, sprintf("Missing required vars: %s", paste(missing_vars, collapse = ", ")))
  }
  
  for (var in names(req_values)) {
    resolved_var = .resolve(tbl, var)
    if (is.na(resolved_var)) {
      problems = c(problems, sprintf("Missing '%s' needed for value checks.", var))
      next
    }
    
    if (inherits(tbl, "FileSystemDataset")) {
      present_vals = character(0)
      tryCatch({
        value_counts = dplyr::select(tbl, !!rlang::sym(resolved_var)) |>
          dplyr::group_by(!!rlang::sym(resolved_var)) |>
          dplyr::summarize(n = dplyr::n()) |>
          dplyr::arrange(dplyr::desc(n)) |>
          dplyr::collect()
        if (nrow(value_counts) > 0) {
          present_vals = na.omit(as.character(value_counts[[resolved_var]]))
        }
      }, error = function(e) {})
      
      if (length(present_vals) == 0) {
        tryCatch({
          sample_data = dplyr::select(tbl, !!rlang::sym(resolved_var)) |>
            utils::head(1000) |>
            dplyr::collect()
          if (nrow(sample_data) > 0) {
            present_vals = na.omit(unique(as.character(sample_data[[resolved_var]])))
          }
        }, error = function(e) {})
      }
    } else if (!is.data.frame(tbl) && inherits(tbl, "Table")) {
      present_vals = character(0)
      tryCatch({
        distinct_vals = tbl |>
          dplyr::select(!!rlang::sym(resolved_var)) |>
          dplyr::distinct() |>
          dplyr::collect()
        if (nrow(distinct_vals) > 0) {
          present_vals = na.omit(as.character(distinct_vals[[resolved_var]]))
        }
      }, error = function(e) {})
    } else {
      present_vals = na.omit(unique(as.character(tbl[[resolved_var]])))
    }
    
    present_vals  = tolower(trimws(as.character(present_vals)))
    expected_vals = unique(tolower(trimws(req_values[[var]])))
    missing_vals  = setdiff(expected_vals, present_vals)
    
    if (length(missing_vals)) {
      problems = c(
        problems,
        sprintf("Variable '%s' is missing expected values: %s", var, paste(missing_vals, collapse = ", "))
      )
    }
  }
  
  if (length(problems)) {
    return(sprintf("Table '%s':\n- %s", table_name, paste(problems, collapse = "\n- ")))
  }
  
  invisible(NULL)
}

## validate all tables ---------------------------------------------------------

validate_all_tables = function(data_list, validation_specs) {
  all_problems = character()
  
  for (spec in validation_specs) {
    tbl_name = spec$table_name
    if (!tbl_name %in% names(data_list)) {
      all_problems = c(all_problems, sprintf("Table '%s' is missing entirely.", tbl_name))
      next
    }
    
    tbl = data_list[[tbl_name]]
    
    tbl_problems = validate_table(
      tbl        = tbl,
      table_name = tbl_name,
      req_vars   = spec$req_vars,
      req_values = spec$req_values
    )
    
    if (!is.null(tbl_problems)) all_problems = c(all_problems, tbl_problems)
  }
  
  if (length(all_problems)) {
    stop("Validation errors found:\n", paste(all_problems, collapse = "\n\n"), call. = FALSE)
  }
  
  message("✅ Validation passed: all required tables present with needed values.")
}

## validation specifications ---------------------------------------------------

patient_spec = list(
  table_name = "patient",
  req_vars   = c("patient_id", "race_category", "ethnicity_category", "sex_category"),
  req_values = list(
    sex_category       = c("Female", "Male"),
    race_category      = c("White", "Black or African American", "Asian"),
    ethnicity_category = c("Hispanic", "Non-Hispanic")
  )
)

hosp_spec = list(
  table_name = "hospitalization",
  req_vars   = c(
    "patient_id",
    "hospitalization_id",
    "age_at_admission",
    "admission_dttm",
    "discharge_dttm",
    "discharge_category"
  ),
  req_values = list(discharge_category = c("Hospice", "Expired"))
)

adt_spec = list(
  table_name = "adt",
  req_vars   = c("hospitalization_id", "hospital_id", "location_category", "in_dttm", "out_dttm"),
  req_values = list(location_category = c("icu"))
)

dx_spec = list(
  table_name = "hospital_diagnosis",
  req_vars   = c("hospitalization_id", "diagnosis_code", "diagnosis_code_format"),
  req_values = list(diagnosis_code_format = c("ICD10CM"))
)

med_spec_c = list(
  table_name = "medication_admin_continuous",
  req_vars   = c("hospitalization_id", "admin_dttm", "med_category", "med_dose", "med_dose_unit"),
  req_values = list(
    med_category = c("norepinephrine", "vasopressin", "epinephrine", "phenylephrine")
  )
)

med_spec_i = list(
  table_name = "medication_admin_intermittent",
  req_vars   = c("hospitalization_id", "admin_dttm", "med_category", "med_dose", "med_dose_unit")
)

resp_spec = list(
  table_name = "respiratory_support",
  req_vars   = c("hospitalization_id", "device_category", "recorded_dttm"),
  req_values = list(device_category = c("IMV"))
)

code_spec = list(
  table_name = "code_status",
  req_vars   = c("patient_id", "start_dttm", "code_status_category")
)

crrt_spec = list(
  table_name = "crrt_therapy",
  req_vars   = c("hospitalization_id", "recorded_dttm")
)
 
validation_specs = list(
  patient_spec,
  hosp_spec,
  adt_spec,
  dx_spec,
  med_spec_c,
  med_spec_i,
  resp_spec,
  code_spec,
  crrt_spec
)

## run validation --------------------------------------------------------------

validate_all_tables(data_list, validation_specs)

rm(patient_spec, hosp_spec, adt_spec, dx_spec, resp_spec, validation_specs)
rm(crrt_spec, code_spec, med_spec_c, med_spec_i); gc()

message("\n✅ 00_setup.R complete. Proceed to 01_cohort.R")

# end 00_setup.R ---------------------------------------------------------------
