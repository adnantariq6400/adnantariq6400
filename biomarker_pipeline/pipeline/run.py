from __future__ import annotations

import pathlib
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import yaml

from biomarker_pipeline.analysis.optimal_range import fit_optimal_range_models, detect_optimal_band
from biomarker_pipeline.analysis.setpoint import estimate_setpoints, evaluate_personalized_intervals
from biomarker_pipeline.analysis.repeats import evaluate_min_repeats
from biomarker_pipeline.viz.plots import plot_risk_curve, plot_setpoint_variability, plot_repeats_accuracy


@dataclass
class Config:
    # Input CSVs (harmonized) or BigQuery settings can be plugged in here later
    longitudinal_csv: Optional[str]
    survival_csv: Optional[str]
    biomarker_col: str
    time_col: str
    event_col: str
    id_col: str
    adjust_cols: List[str]
    normal_range: Tuple[float, float]


def load_config(path: str) -> Config:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return Config(
        longitudinal_csv=cfg.get("longitudinal_csv"),
        survival_csv=cfg.get("survival_csv"),
        biomarker_col=cfg["biomarker_col"],
        time_col=cfg["time_col"],
        event_col=cfg["event_col"],
        id_col=cfg["id_col"],
        adjust_cols=cfg.get("adjust_cols", []),
        normal_range=tuple(cfg.get("normal_range", [None, None])),
    )


def run_pipeline(config_path: str, output_dir: str) -> None:
    out = pathlib.Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    cfg = load_config(config_path)

    assert cfg.longitudinal_csv and cfg.survival_csv, "Provide harmonized input CSVs in config."

    longitudinal = pd.read_csv(cfg.longitudinal_csv)
    survival = pd.read_csv(cfg.survival_csv)

    or_models, risk_df = fit_optimal_range_models(
        longitudinal=longitudinal,
        survival=survival,
        biomarker_col=cfg.biomarker_col,
        time_col=cfg.time_col,
        event_col=cfg.event_col,
        id_col=cfg.id_col,
        adjust_cols=cfg.adjust_cols,
        normal_range=cfg.normal_range if all(cfg.normal_range) else None,
    )
    band = detect_optimal_band(
        risk_df,
        normal_range=cfg.normal_range if all(cfg.normal_range) else None,
        relative_margin=0.02,
    )
    risk_df.to_csv(out / "risk_curve.csv", index=False)
    if band:
        with open(out / "optimal_band.txt", "w") as f:
            f.write(f"Optimal band: {band[0]:.2f} to {band[1]:.2f}\n")

    fig1 = plot_risk_curve(
        risk_df,
        normal_range=cfg.normal_range if all(cfg.normal_range) else None,
        optimal_band=band,
        title="Risk curve",
    )
    fig1.savefig(out / "risk_curve.png", dpi=150, bbox_inches="tight")

    setpoint_df, variance_components = estimate_setpoints(
        longitudinal,
        id_col=cfg.id_col,
        value_col=cfg.biomarker_col,
    )
    perf_df = evaluate_personalized_intervals(
        longitudinal=longitudinal,
        survival=survival,
        setpoint_df=setpoint_df,
        id_col=cfg.id_col,
        value_col=cfg.biomarker_col,
        time_col=cfg.time_col,
        event_col=cfg.event_col,
        adjust_cols=cfg.adjust_cols,
    )

    setpoint_df.to_csv(out / "setpoints.csv", index=False)
    pd.DataFrame([variance_components]).to_csv(out / "variance_components.csv", index=False)
    perf_df.to_csv(out / "personalized_vs_population.csv", index=False)

    fig2 = plot_setpoint_variability(
        longitudinal, setpoint_df, id_col=cfg.id_col, value_col=cfg.biomarker_col
    )
    fig2.savefig(out / "setpoint_variability.png", dpi=150, bbox_inches="tight")

    repeats_df = evaluate_min_repeats(
        longitudinal=longitudinal,
        id_col=cfg.id_col,
        value_col=cfg.biomarker_col,
        max_repeats=6,
    )
    repeats_df.to_csv(out / "min_repeats_accuracy.csv", index=False)

    fig3 = plot_repeats_accuracy(repeats_df)
    fig3.savefig(out / "min_repeats_accuracy.png", dpi=150, bbox_inches="tight")
