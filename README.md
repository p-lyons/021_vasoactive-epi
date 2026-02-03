# Coordinating Center - Vasopressor Escalation Study

## Folder Structure

```
coordinating_center/
├── 00_pool_load.R       # Load site data, compute pooled Ns
├── 01_pool_table1.R     # Create pooled Table 1
├── 02_pool_qc.R         # QC review across sites
├── sites/               # Site-level outputs (gitignored)
│   ├── emory/
│   ├── hopkins/
│   ├── ohsu/
│   ├── rush/
│   ├── ucmc/
│   ├── umn/
│   └── upenn/
└── output/              # Generated outputs (gitignored)
    ├── tables/
    └── qc/
```

## Setup

1. Create `sites/` folder
2. For each site, create subfolder (e.g., `sites/ohsu/`)
3. Copy site's `upload_to_box/*.csv` files into their folder

## Expected Files Per Site

**Table 1 components:**
- `table1_continuous_{site}.csv`
- `table1_binary_{site}.csv`
- `table1_categorical_{site}.csv`
- `table1_timing_{site}.csv`
- `table1_totals_{site}.csv`
- `flow_diagram_{site}.csv`

**QC files:**
- `qc_missing_{site}.csv`
- `qc_ranges_{site}.csv`
- `qc_flags_{site}.csv`
- `qc_categories_{site}.csv`
- `qc_diagnostics_{site}.csv`
- `exclusion_cascade_{site}.csv`

## Running

```r
source("00_pool_load.R")
source("01_pool_table1.R")  
source("02_pool_qc.R")
```

## Outputs

**Tables (`output/tables/`):**
- `table1_characteristics_YYMMDD.docx` - Final Table 1
- `table1_data_YYMMDD.csv` - Underlying data
- `flow_diagram_pooled_YYMMDD.csv` - Pooled flow

**QC (`output/qc/`):**
- `qc_exclusion_pooled_*.csv` - Pooled exclusion cascade
- `qc_exclusion_bysite_*.csv` - Exclusion % by site
- `qc_flags_summary_*.csv` - Plausibility flags
- `qc_missing_comparison_*.csv` - Missingness by site
- `qc_ranges_*.csv` - Continuous variable distributions
- `qc_categories_*.csv` - Categorical value inventory
- `qc_outcome_proportions_*.csv` - Outcome rates by site
- `qc_site_comparison_*.csv` - Key metrics comparison
- `qc_site_smds_*.csv` - Site vs pooled SMDs
- `qc_diagnostics_summary_*.csv` - Site diagnostics dashboard
