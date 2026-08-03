# PLT Trajectory Analysis

Reproducible R analysis for platelet trajectories in acute pancreatitis.

## Input data

The script loads four files (edit paths at the top of `PLT_analysis.R` if needed):

- Local, TRACE and MIMIC-IV cohort files (CSV, GBK encoding)
- `training_data.xlsx` - random forest training data

Cohort files are located automatically by pattern under `data_dir`.

## Run

```r
Rscript PLT_analysis.R
```

Missing R packages are installed automatically.

## Analyses

- Cohort description
- Cumulative platelet exposure (CumPLT)
- GBTM trajectory modeling
- Cox models for mortality (in-hospital / 28-day / 90-day)
- Restricted cubic splines
- Blood culture associations
- Mediation analysis
- Double machine learning (exploratory)
- Random forest classifier
- Sensitivity analyses (K-means, multiple imputation, E-values)

Trajectory class labels are taken from the cohort files (frozen GBTM output).
Some exploratory estimates may differ slightly from the manuscript because the
archived files are processed versions of the source data.
