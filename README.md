# CLIF Vasopressor Escalation Study

## Site-Level Analysis Pipeline

This pipeline identifies patients with refractory septic shock receiving norepinephrine ≥0.2 mcg/kg/min plus vasopressin, and characterizes outcomes at 48 hours based on vasopressor escalation and mortality.

---

## Table of Contents

1. [Overview](#overview)
2. [Cohort Definition](#cohort-definition)
3. [Outcome Groups](#outcome-groups)
4. [Prerequisites](#prerequisites)
5. [Setup Instructions](#setup-instructions)
6. [Running the Pipeline](#running-the-pipeline)
7. [Output Files](#output-files)
8. [Script Descriptions](#script-descriptions)
9. [Variables](#variables)
10. [Troubleshooting](#troubleshooting)

---

## Overview

This study examines vasopressor escalation patterns in patients with septic shock who are already receiving high-dose norepinephrine plus vasopressin. The pipeline:

- Links contiguous hospitalizations (≤6 hour gaps)
- Identifies T0: first time with NE ≥0.2 mcg/kg/min AND vasopressin concurrently (no other vasopressors running)
- Tracks escalation to additional agents within 48 hours
- Generates poolable summary statistics for multi-site analysis

---

## Cohort Definition

**Inclusion criteria:**
- Age ≥18 years at admission
- Admission between 2016-01-01 and 2024-12-31
- Had ED or ICU stay during hospitalization
- Received norepinephrine ≥0.2 mcg/kg/min AND vasopressin concurrently
- At T0, no other vasopressors were running (epinephrine, phenylephrine, dopamine, angiotensin II)

**Exclusion criteria:**
- Psychiatric or rehabilitation unit stay
- Missing discharge disposition

**T0 definition:**
First timepoint where:
- Norepinephrine dose ≥0.2 mcg/kg/min (carry-forward)
- Vasopressin is running (any dose, carry-forward)
- No epinephrine, phenylephrine, dopamine, or angiotensin II running

---

## Outcome Groups

Patients are classified into 4 mutually exclusive groups based on status at 48 hours post-T0:

| Group | Label | Definition |
|-------|-------|------------|
| 0 | No Esc + Dead/Hospice | No escalation, died or discharged to hospice within 48h |
| 1 | Esc + Dead/Hospice | Escalated, died or discharged to hospice within 48h |
| 2 | Esc + Alive | Escalated, alive at 48h |
| 3 | No Esc + Alive | No escalation, alive at 48h |

**Escalation** is defined as receipt of any of the following after T0 and within 48 hours:
- Epinephrine (continuous infusion)
- Phenylephrine (continuous infusion)
- Dopamine (continuous infusion)
- Angiotensin II (continuous infusion)
- Methylene blue (intermittent)
- Hydroxocobalamin (intermittent)

---

## Prerequisites

### Required CLIF Tables

Your site must have the following CLIF tables in parquet, CSV, or FST format:

| Table | Required Columns |
|-------|------------------|
| `clif_patient` | patient_id, race_category, ethnicity_category, sex_category |
| `clif_hospitalization` | patient_id, hospitalization_id, age_at_admission, admission_dttm, discharge_dttm, discharge_category |
| `clif_adt` | hospitalization_id, hospital_id, location_category, in_dttm, out_dttm |
| `clif_hospital_diagnosis` | hospitalization_id, diagnosis_code, diagnosis_code_format |
| `clif_medication_admin_continuous` | hospitalization_id, admin_dttm, med_category, med_dose, med_dose_unit |
| `clif_medication_admin_intermittent` | hospitalization_id, admin_dttm, med_category, med_dose, med_dose_unit |
| `clif_respiratory_support` | hospitalization_id, device_category, recorded_dttm |
| `clif_crrt_therapy` | hospitalization_id, recorded_dttm |
| `clif_code_status` | patient_id, start_dttm, code_status_category |
| `clif_vitals` | hospitalization_id, recorded_dttm, vital_category, vital_value |
| `clif_labs` | hospitalization_id, lab_result_dttm, lab_category, lab_value |

### Required med_category Values

In `clif_medication_admin_continuous`:
- `norepinephrine`
- `vasopressin`
- `epinephrine`
- `phenylephrine`
- `dopamine`
- `angiotensin` (angiotensin II)

In `clif_medication_admin_intermittent`:
- `methylene_blue`
- `hydroxocobalamin`

### Required Config Files

Place these in `config/`:

1. **config_clif_pressors.yaml** - Site-specific configuration
2. **clif_sites.csv** - List of valid site names
3. **svi_2020.parquet** - Social Vulnerability Index data (optional)
4. **adi_2020.parquet** - Area Deprivation Index data (optional)

---

## Setup Instructions

### 1. Clone or Download the Repository

```bash
git clone <repository_url>
cd clif_pressors
```

### 2. Create Configuration File

Create `config/config_clif_pressors.yaml`:

```yaml
# Site identifier (lowercase, must match clif_sites.csv)
site_lowercase: "ohsu"

# File format of your CLIF tables: "parquet", "csv", or "fst"
file_type: "parquet"

# Path to folder containing CLIF tables
clif_data_location: "/path/to/clif_tables"

# Path to this project folder
project_location: "/path/to/clif_pressors"
```

### 3. Verify clif_sites.csv

Ensure `config/clif_sites.csv` contains your site name:

```csv
site_name
emory
jhu
northwestern
ohsu
...
```

### 4. Install R Packages

The pipeline will auto-install missing packages, but you can pre-install:

```r
install.packages(c(
  "data.table", "tidyverse", "tidytable", "collapse",
  "arrow", "here", "comorbidity", "yaml", "janitor",
  "tableone", "smd"
))
```

---

## Running the Pipeline

### Option 1: Run All Scripts

```r
source("run_all.R")
```

This executes all scripts in sequence and reports elapsed time.

### Option 2: Run Scripts Individually

```r
source("00_setup.R")    # Load packages, validate data
source("01_cohort.R")   # Build cohort, apply exclusions
source("02_variables.R") # Define outcomes and variables
source("03_table.R")    # Generate summary statistics
```

---

## Output Files

All files for upload are saved to `upload_to_box/`:

### Table 1 Components (for pooling)

| File | Description |
|------|-------------|
| `table1_continuous_{site}.csv` | Continuous variable summary stats (n, sum, sumsq, percentiles) |
| `table1_binary_{site}.csv` | Binary variable counts (n, n_1) |
| `table1_categorical_{site}.csv` | Categorical variable counts by category |
| `table1_timing_{site}.csv` | Timing group distributions (before T0, after T0, none) |
| `table1_totals_{site}.csv` | Total N per outcome group |

### Flow Diagram

| File | Description |
|------|-------------|
| `flow_diagram_{site}.csv` | Counts by outcome group |
| `exclusion_cascade_{site}.csv` | Step-by-step exclusion counts |

### QC Files

| File | Description |
|------|-------------|
| `qc_missingness_{site}.csv` | Missing data rates by variable |
| `qc_ranges_{site}.csv` | Min/max/median for continuous variables |
| `qc_categories_{site}.csv` | Category frequencies |
| `qc_diagnostics_{site}.csv` | Key site metrics for validation |

### Intermediate Files (proj_tables/)

These are saved locally for debugging but do **not** need to be uploaded:

- `cohort.parquet` - Final cohort with all variables
- `hid_jid_crosswalk.parquet` - Hospitalization ID linkage
- `vasoactive_doses.parquet` - Time-series of vasopressor doses

---

## Script Descriptions

### 00_setup.R
- Loads and installs required packages
- Configures parallel processing (threads based on available RAM)
- Reads site configuration from YAML
- Loads CLIF tables as Arrow datasets
- Validates required tables and columns exist

### 01_cohort.R
- Links contiguous hospitalizations (≤6 hour gaps)
- Applies inclusion/exclusion criteria
- Identifies T0 for each encounter
- Builds vasoactive_doses table with carry-forward logic
- Adds demographics, Elixhauser comorbidities, IMV/CRRT times
- Saves exclusion cascade

### 02_variables.R
- Defines 48-hour outcome window
- Identifies escalation events (continuous and intermittent drugs)
- Calculates maximum NE-equivalent dose (with dose caps)
- Determines IMV and CRRT timing relative to T0
- Adds code status, SVI, ADI
- Assigns 4-way outcome groups

### 03_table.R
- Generates poolable summary statistics by outcome group
- Continuous variables: n, sum, sum of squares, percentiles
- Binary variables: n, n with value=1
- Categorical variables: counts per category
- Creates QC diagnostics and flow diagram

---

## Variables

### Norepinephrine Equivalent Calculation

```
NE-equiv = NE + Epi + (2.5 × VP) + (0.1 × Phenyl) + (0.01 × Dopa) + (0.01 × A2)
```

**Dose caps applied before calculation:**

| Drug | Cap | Units |
|------|-----|-------|
| Norepinephrine | 5 | mcg/kg/min |
| Epinephrine | 5 | mcg/kg/min |
| Vasopressin | 0.1 | units/min |
| Phenylephrine | 2 | mcg/kg/min |
| Dopamine | 20 | mcg/kg/min |
| Angiotensin II | 80 | ng/kg/min |

### Table 1 Variables

**Characteristics:**
- Age (years)
- Female sex
- White race (vs non-white; unknown coded as NA)
- Hispanic ethnicity (vs non-Hispanic; unknown coded as NA)
- Van Walraven comorbidity score
- Social Vulnerability Index percentile
- Area Deprivation Index percentile
- Full code status at T0
- Days in hospital before T0
- Days in ICU before T0
- Invasive mechanical ventilation at T0
- CRRT at T0

**Outcomes:**
- In-hospital mortality
- Discharge to hospice
- Length of stay after T0
- Maximum NE-equivalent dose in 48h
- Receipt of epinephrine, phenylephrine, dopamine, angiotensin II
- Receipt of methylene blue, hydroxocobalamin

---

## Troubleshooting

### "Missing required tables"
Ensure all CLIF tables are in the folder specified by `clif_data_location` and follow the naming convention `clif_{table_name}.parquet`.

### "Invalid site"
Your `site_lowercase` in the config file must match a row in `config/clif_sites.csv`.

### Validation errors for missing values
Some categories (e.g., race, med_category) must contain expected values. Check that your CLIF tables use standard CLIF category names.

### Empty cohort
- Verify your site has patients receiving norepinephrine AND vasopressin concurrently
- Check that med_dose_unit is standardized (mcg/kg/min for pressors, units/min for vasopressin)
- Review vasoactive_doses.parquet for dose distributions

### Memory errors
The pipeline auto-configures threads based on available RAM. If you still encounter memory issues, you can manually reduce threads in 00_setup.R:
```r
n_threads = 2L  # Override automatic detection
```

---

## Contact

For questions about the pipeline, contact the coordinating center.

For site-specific data issues, contact your local CLIF data team.
