# Biomarker Optimal Ranges & Setpoints Pipeline

A reproducible pipeline to analyze biomarker optimal ranges and individual setpoints across large cohorts. Includes a synthetic demo run with end-to-end outputs.

## Quick start

```bash
# 1) Create and activate a virtual environment (recommended)
python -m venv .venv && source .venv/bin/activate

# 2) Install dependencies
pip install -r requirements.txt

# 3) Run the synthetic demo (no external data required)
python -m biomarker_pipeline.cli demo --output-dir outputs/demo

# Outputs: CSVs and PNGs under outputs/demo
```

## CLI usage

```bash
python -m biomarker_pipeline.cli --help
python -m biomarker_pipeline.cli demo --output-dir outputs/demo
python -m biomarker_pipeline.cli run --config configs/demo.yaml --output-dir outputs/run1
```

## Contents
- `biomarker_pipeline/`: Python package with CLI, pipeline, analysis, and viz
- `configs/`: YAML configs for demo and templates for real data
- `docs/design.md`: Design document with methods, workflow and example visuals
- `outputs/`: Default output directory (created on first run)

## BigQuery (optional)
Set the following environment variables if connecting to BigQuery:

```
GOOGLE_CLOUD_PROJECT=your-project-id
GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json
```

Then implement dataset-specific queries in `biomarker_pipeline/data_access/bigquery_client.py` and mappings in `biomarker_pipeline/harmonization/`.
