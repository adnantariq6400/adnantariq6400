## Project classification
- **Domain**: Epidemiology / clinical biomarkers (not bioinformatics; not metagenomics).
- **Main purpose**: Detect, for each biomarker, (1) existence of a risk-minimizing optimal subrange within the population "normal" range and (2) setpoint-like intra-individual stability; evaluate whether personalized reference intervals improve risk prediction.
  - Clarification for interviews: While bioinformatics can include biomarker analytics, this project focuses on population-based clinical lab measures and outcomes modeling (epidemiology/biostatistics, precision medicine), not metagenomic sequencing or omics pipelines.

## Data and cohorts
- BigQuery-accessible datasets: UK Biobank, ARIC, MESA, CARDIA, NHANES (cross-sectional).
- Repeated biomarker measures (except NHANES), longitudinal outcomes, demographics and covariates.

## Harmonization strategy
- Create cohort-specific extraction views to a harmonized schema:
  - Longitudinal table: `participant_id`, `visit`, `time` (years since baseline), `age`, `sex`, `biomarker` (one analyte per run).
  - Survival table: `participant_id`, `time` (follow-up), `event` (0/1), baseline covariates (`age0`, `sex`), and optional adjustments.
- Use unit normalization (e.g., mmol/L to mg/dL if needed) with explicit mapping per cohort.
- Enforce schemas with pandera validations.

## Modeling plan

### 1) Optimal ranges
- Fit Cox proportional hazards model with cubic regression spline for the biomarker (restricted to standard normal interval or full range with later restriction).
- Generate relative hazard curve across biomarker values, centered at median hazard.
- Define "practical optimal band" as the maximal contiguous interval with relative hazard ≤ (1+δ)×min hazard (δ ≈ 2–5%).
- Validate across cohorts via meta-analysis of curve minima and band overlap; subgroup analyses by age/sex and baseline disease.

### 2) Setpoint detection
- Estimate individual setpoints using random-intercept mixed models on repeated measures. Report variance components and ICC to quantify intra- vs inter-individual variability.
- Compare prediction using (a) population deviation (value − population mean) vs (b) personalized deviation (value − individual setpoint) in Cox models; compare concordance index and calibration.
- Optionally define personalized reference intervals: setpoint ± z×within-person SD; evaluate association of out-of-interval status with outcomes.

### 3) Minimal number of tests
- For k=1..K, estimate setpoint by mean of first k visits and compare to "gold" (all visits) using MAE/RMSE. Optionally, model error vs k via asymptotic curve to recommend minimum k achieving target MAE.
- Sensitivity analyses: increased measurement error, drift, missingness; simulate to assess robustness of k.

## Confounding, measurement error, generalizability
- Adjust for demographics, baseline comorbidities, seasonality (where relevant), and lab batch if available.
- Address regression dilution bias via mixed models and simulation-extrapolation (SIMEX) if strong error suspected.
- Use cross-cohort validation and transportability checks; quantify heterogeneity.

## Workflow (high level)
1. Extract & harmonize from BigQuery to CSV/parquet per biomarker.
2. Validate schemas and units.
3. Optimal range modeling → risk curve & band detection.
4. Setpoint estimation → variance components and personalized intervals; risk model comparison.
5. Minimal repeats evaluation and sensitivity analyses.
6. Reporting: tables, plots, and summary JSON.

### Workflow diagram (ASCII)
```
BigQuery cohort tables
      │
      ├─ Cohort-specific SQL views (per biomarker)
      │
      ├─ Harmonize columns & units → validated longitudinal/survival tables
      │
      ├─ Analysis 1: Optimal range (Cox + splines)
      │        └─ Risk curve + optimal band
      │
      ├─ Analysis 2: Setpoints (mixed models)
      │        └─ Variance components + personalized deviation models
      │
      ├─ Analysis 3: Minimal repeats (k=1..K)
      │        └─ MAE/RMSE vs k + simulations
      │
      └─ Outputs: CSV tables, PNG plots, summary
```

## Outputs
- Risk curve with candidate optimal band overlay.
- Trajectories with setpoint overlay; distribution of setpoints.
- Error vs number of repeats curve.

## BigQuery extraction & harmonization (example)
Example longitudinal extraction (pseudo-SQL):
```sql
SELECT
  person_id        AS participant_id,
  exam_no          AS visit,
  TIMESTAMP_DIFF(visit_date, baseline_date, DAY)/365.25 AS time,
  age_at_visit     AS age,
  IF(sex='Male',1,0) AS sex,
  sodium_mmol_L    AS biomarker
FROM `project.dataset.lab_table`
WHERE NOT SAFE_CONVERT_FLOAT64(sodium_mmol_L) IS NULL;
```
Then apply a dataset-specific mapping to harmonized columns and units.

## How to run (demo)
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m biomarker_pipeline.cli demo --output-dir outputs/demo
```

To run on harmonized CSVs:
```bash
python -m biomarker_pipeline.cli run --config configs/demo.yaml --output-dir outputs/run1
```

## Tools and libraries
- Python: lifelines (Cox), statsmodels (mixed models), patsy (splines), pandas, numpy, seaborn/matplotlib; BigQuery via `google-cloud-bigquery`.
 - CLI: Typer. Validation: pandera. Metrics: lifelines C-index, custom MAE/RMSE. Optional: scikit-learn for additional metrics and CV.

## Assumptions and limitations
- Cox PH holds; competing risks and time-varying exposures are not fully modeled in the demo.
- Last-observation value used for outcome modeling in demo; real analyses should consider time-updated exposures with time-varying Cox models.
- Measurement error modeled simply; lab drift/batch should be handled when available.
 - NHANES is cross-sectional; use logistic or Poisson models for associations (no survival), primarily for setpoint estimation and range description, not risk over time.
