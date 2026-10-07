# CLIF Vasopressor Escalation Study

Vasopressor escalation in refractory distributive shock: a federated, multicenter CLIF cohort.

This repository holds the site-level pipeline that each CLIF site runs locally, and the coordinating-site scripts that pool the site outputs. No patient-level data leave a site.

---

## Contents

1. [Repository layout](#repository-layout)
2. [Cohort definition](#cohort-definition)
3. [Outcome groups](#outcome-groups)
4. [Prerequisites](#prerequisites)
5. [Setup](#setup)
6. [Running the pipeline](#running-the-pipeline)
7. [Output files](#output-files)
8. [Small-cell policy](#small-cell-policy)
9. [Variables](#variables)
10. [Troubleshooting](#troubleshooting)

---

## Repository layout

```
code/                          site-level pipeline (run at each site)
  00_setup.R                   packages, config, CLIF table loading, validation
  01_cohort.R                  cohort, T0, exclusion cascade, vasoactive doses
  02_variables.R               escalation, outcomes, hospital at T0, covariates
  03_table.R                   poolable Table 1 and QC summaries, by hospital
  run_all.R                    runs 00 -> 03 in order
  sankey_transitions_onepass.R standalone: state-transition counts for the Sankey figure
code_for_pooled_data/          coordinating site only: pools site outputs (run in order)
  00_pool_load.R               loads sites/{site}/ files; site_total vs hospital rows; Sankey pooling
  01_pool_table1.R             pooled Table 1 and flow diagram
  02_pool_qc.R                 cross-site QC, hospital attribution, Sankey cascade check
  03_pool_hospital_variation.R hospital ranges (>= 30 encounters) and random-intercept ICC/MOR
config/
  config_clif_pressors_EXAMPLE.yaml   template for your site config
  clif_sites.csv               valid site names and time zones
  svi_2020.parquet             Social Vulnerability Index (optional)
  adi_2020.parquet             Area Deprivation Index (optional)
```

---

## Cohort definition

**Inclusion criteria**

- Age 18 years or older at admission
- Admission between 2016-01-01 and 2024-12-31
- ED or ICU stay during the hospitalization
- Norepinephrine >= 0.2 mcg/kg/min and vasopressin active at the same time, with no other continuous vasopressor active (T0)

**Exclusion criteria** (in cascade order, as written to `exclusion_cascade_{site}.csv`)

| Step | Excludes |
|------|----------|
| `02_required_ed_or_icu` | No ED or ICU stay |
| `03_excluded_psych_rehab` | Psychiatric or rehabilitation unit stay |
| `04_missing_discharge_category` | Missing discharge disposition |
| `04b_yodo_cleanup` | Duplicate death records and encounters after death |
| `05_uninterpretable_mar` | No usable (finite, > 0) norepinephrine or vasopressin dose |
| `06_met_t0_criteria` | Never met T0 criteria |
| `07_third_line_prior_to_t0` | Any other continuous vasopressor active before T0 |

Hospitalizations with gaps of less than 6 hours are linked into one encounter.

**Dose activity rule.** A continuous dose is active only if it is > 0, it is not a "stopped" action, and it was charted within the previous 4 hours (carry-forward window). Doses are converted to mcg/kg/min (catecholamines), units/min (vasopressin), and ng/kg/min (angiotensin II). The first in-encounter weight is used; 70 kg if none is charted.

**T0** is the first MAR timestamp at which norepinephrine >= 0.2 mcg/kg/min and vasopressin are active, and epinephrine, phenylephrine, dopamine, and angiotensin II are not.

The cohort logic in `01_cohort.R` and `sankey_transitions_onepass.R` must stay identical. Change both together.

---

## Outcome groups

Each encounter is classified at 48 hours after T0:

| Group | Label | Definition |
|-------|-------|------------|
| `noesc_dead` | No escalation + dead/hospice | No escalation; died or discharged to hospice within 48 h |
| `esc_dead` | Escalation + dead/hospice | Escalated; died or discharged to hospice within 48 h |
| `esc_alive` | Escalation + alive | Escalated; alive and not in hospice at 48 h |
| `noesc_alive` | No escalation + alive | No escalation; alive and not in hospice at 48 h |

**Escalation** is any of the following after T0 and within 48 hours:

- Continuous epinephrine, phenylephrine, dopamine, or angiotensin II (active dose)
- Intermittent methylene blue or hydroxocobalamin

---

## Prerequisites

### Required CLIF tables

CLIF 2.1 tables in parquet, CSV, or FST format, named `clif_{table}.{ext}`:

| Table | Columns used |
|-------|--------------|
| `patient` | patient_id, sex_category, race_category, ethnicity_category, language_category |
| `hospitalization` | patient_id, hospitalization_id, age_at_admission, admission_dttm, discharge_dttm, discharge_category, census_block_code, census_block_group_code |
| `adt` | hospitalization_id, hospital_id, hospital_type, location_category, in_dttm, out_dttm |
| `vitals` | hospitalization_id, recorded_dttm, vital_category, vital_value |
| `labs` | loaded and validated only |
| `hospital_diagnosis` | hospitalization_id, diagnosis_code, diagnosis_code_format, poa_present |
| `medication_admin_continuous` | hospitalization_id, admin_dttm, med_category, med_dose, med_dose_unit, mar_action_category |
| `medication_admin_intermittent` | hospitalization_id, admin_dttm, med_category, med_dose, med_dose_unit |
| `respiratory_support` | hospitalization_id, device_category, recorded_dttm |
| `crrt_therapy` | hospitalization_id, recorded_dttm |
| `code_status` | patient_id, start_dttm, code_status_category |
| `patient_procedures` | hospitalization_id, procedure_code, procedure_code_format, procedure_billed_dttm |

### Required `med_category` values

- `medication_admin_continuous`: `norepinephrine`, `vasopressin`, `epinephrine`, `phenylephrine`, `dopamine`, `angiotensin`
- `medication_admin_intermittent`: `methylene_blue`, `hydroxocobalamin`

### Optional files in `config/`

- `svi_2020.parquet`, `adi_2020.parquet`: neighborhood indices. SVI/ADI are set to NA if absent or if `census_block_code` is unusable.
- `PClassR_v2026-1.csv`: AHRQ procedure classes for the major-procedure flag. If absent, `major_procedure_01` is 0 for all encounters.

---

## Setup

1. Clone the repository.

   ```bash
   git clone https://github.com/p-lyons/021_vasoactive-epi.git
   ```

2. Copy the config template and edit it for your site. `config/config_clif_pressors.yaml` is git-ignored.

   ```bash
   cp config/config_clif_pressors_EXAMPLE.yaml config/config_clif_pressors.yaml
   ```

3. Confirm that your `site_lowercase` value is in `config/clif_sites.csv`.

4. Open `021_vasoactive-epi.Rproj` in RStudio. Run `renv::restore()` to install the package versions in `renv.lock`, or let `00_setup.R` install missing packages.

---

## Running the pipeline

Run the full site pipeline from the project root:

```r
source(here::here("code", "run_all.R"))
```

Then run the Sankey export. It is standalone (it rebuilds the cohort itself), so it can run in a fresh session:

```r
source(here::here("code", "sankey_transitions_onepass.R"))
```

**Check after both runs:** the exclusion cascade counts in `upload_to_box/exclusion_cascade_{site}.csv` must match the `n_*` columns in `output/sankey_site_summary_{site}.csv`.

---

## Output files

### `upload_to_box/` (send to the coordinating site)

All `table1_*` and `flow_diagram` files are stratified by hospital (`hospital` = hospital at T0) and also carry `hospital = "site_total"` rows. Exclude `site_total` rows when summing across hospitals.

| File | Contents |
|------|----------|
| `table1_continuous_{site}.csv` | n, missing, sum, sum of squares, min, max, percentiles by hospital and outcome group |
| `table1_binary_{site}.csv` | n and n_1 by hospital and outcome group |
| `table1_categorical_{site}.csv` | Cell counts by hospital, outcome group, and category |
| `table1_timing_{site}.csv` | IMV and CRRT timing groups (none / before T0 / T0 to 48 h) |
| `table1_totals_{site}.csv` | Encounters and patients by hospital and outcome group |
| `flow_diagram_{site}.csv` | Outcome-group counts by hospital |
| `exclusion_cascade_{site}.csv` | Exclusion counts by step |
| `qc_missing_{site}.csv` | Missingness by variable and hospital |
| `qc_ranges_{site}.csv` | Distribution of continuous variables (site level) |
| `qc_flags_{site}.csv` | Plausibility flag counts (site level) |
| `qc_categories_{site}.csv` | Category frequencies (site level) |
| `qc_diagnostics_{site}.csv` | Key metrics by hospital; study period as year-month |
| `qc_hospital_{site}.csv` | Hospital attribution checks |

### `output/` (send to the coordinating site)

| File | Contents |
|------|----------|
| `sankey_transitions_{site}.csv` | Counts of 6-hour block-to-block state transitions over 48 h |
| `sankey_site_summary_{site}.csv` | Cohort size, exclusion cascade, settings |

### `proj_tables/` (local only; never upload or commit)

Patient-level intermediates: `cohort.parquet`, `cohort_analytic.parquet`, `hid_jid_crosswalk.parquet`, `vasoactive_doses.parquet`, `exclusion_cascade.csv`.

---

## Small-cell policy

- Counts of 1–4 are masked (NA) at the site only for sensitive demographic variables: sex, race, ethnicity, and language.
- All other counts are not PHI and are not masked.
- Site files are pooled at the coordinating site under the consortium data use agreement before anything is shared.
- Hospital-level variation ranges are reported only for hospitals with at least 30 T0-eligible encounters.

---

## Variables

### Norepinephrine-equivalent dose

```
NEE = NE + Epi + (2.5 x VP) + (0.1 x Phenyl) + (0.01 x Dopa) + (0.01 x A2)
```

Doses are capped before the calculation:

| Drug | Cap | Units |
|------|-----|-------|
| Norepinephrine | 5 | mcg/kg/min |
| Epinephrine | 5 | mcg/kg/min |
| Vasopressin | 0.1 | units/min |
| Phenylephrine | 2 | mcg/kg/min |
| Dopamine | 20 | mcg/kg/min |
| Angiotensin II | 80 | ng/kg/min |

### Hospital at T0

`hospital_id_t0` is the `hospital_id` of the last ADT row at or before T0 (else the first ADT row of the encounter). `academic_01` comes from the `hospital_type` of the same ADT row.

### Sankey states (Protocol Sec 6.3.1)

States are assigned at each 6-hour block boundary from T0 to 48 h, with a within-block look-back that keeps the most severe state reached in the block.

| State | Definition (active agents) |
|-------|----------------------------|
| D | Died or discharged to hospice |
| W | No vasopressor |
| S1 | One agent |
| S0 | Vasopressin + one adrenergic agent |
| S2a | Two or more adrenergic agents |
| S2b | Angiotensin II + norepinephrine and/or vasopressin (three or fewer agents) |
| S3 | Four or more agents, or angiotensin II + epinephrine, phenylephrine, or dopamine |

---

## Troubleshooting

**"Missing required tables"**: confirm that all tables are in `clif_data_location` and are named `clif_{table}.{ext}`.

**"Invalid site"**: `site_lowercase` must match a row in `config/clif_sites.csv`.

**Empty cohort**: check `med_dose_unit` values and the dose distributions that `01_cohort.R` prints after unit correction.

**Memory errors**: reduce threads manually in `00_setup.R` (for example, `n_threads = 2L`).

---

## Contact

Pipeline questions: the coordinating site (OHSU). Site data questions: your local CLIF data team.
