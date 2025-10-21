from __future__ import annotations

import pathlib
import typing as t

import numpy as np
import pandas as pd

from biomarker_pipeline.analysis.optimal_range import fit_optimal_range_models, detect_optimal_band
from biomarker_pipeline.analysis.setpoint import estimate_setpoints, evaluate_personalized_intervals
from biomarker_pipeline.analysis.repeats import evaluate_min_repeats
from biomarker_pipeline.viz.plots import (
    plot_risk_curve,
    plot_setpoint_variability,
    plot_repeats_accuracy,
)


def simulate_synthetic_cohort(
    n_participants: int = 5000,
    n_visits: int = 4,
    seed: int | None = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)

    # Baseline covariates
    participant_id = np.arange(n_participants)
    age0 = rng.normal(55, 8, size=n_participants)
    sex = rng.integers(0, 2, size=n_participants)  # 0 F, 1 M

    # Subject-specific setpoints for biomarker X
    # Population mean 100, SD 15; individual setpoint SD 8 around population mean
    setpoint = rng.normal(100, 8, size=n_participants)

    # Longitudinal visits with measurement error and small drift
    rows = []
    for v in range(n_visits):
        time_years = v * 2.0
        drift = rng.normal(0, 0.5, size=n_participants) * (v)
        measurement_error = rng.normal(0, 5, size=n_participants)
        value = setpoint + drift + measurement_error
        age = age0 + time_years
        rows.append(
            pd.DataFrame(
                {
                    "participant_id": participant_id,
                    "visit": v,
                    "time": time_years,
                    "age": age,
                    "sex": sex,
                    "biomarker": value,
                }
            )
        )
    longitudinal = pd.concat(rows, ignore_index=True)

    # Survival outcome simulated from hazard depending on biomarker (U-shaped risk)
    # True optimal around 100 with increased risk at high/low
    baseline_hazard = 0.02
    beta_quadratic = 0.002  # U-shape around 100

    # Use last observed biomarker as covariate for event generation
    last_by_id = longitudinal.sort_values(["participant_id", "time"]).groupby("participant_id").tail(1)
    biom = last_by_id["biomarker"].to_numpy()

    # Individual hazard multiplier
    hazard = baseline_hazard * np.exp(beta_quadratic * (biom - 100.0) ** 2 + 0.02 * (sex) + 0.01 * (age0 - 55))
    follow_up = 10.0
    # Exponential time-to-event with censoring at follow_up
    u = rng.uniform(size=n_participants)
    tte = -np.log(1 - u) / hazard
    event = (tte <= follow_up).astype(int)
    time = np.minimum(tte, follow_up)

    survival = pd.DataFrame(
        {
            "participant_id": participant_id,
            "time": time,
            "event": event,
            "age0": age0,
            "sex": sex,
        }
    )

    return longitudinal, survival


def run_synthetic_demo(output_dir: str, seed: int = 42) -> None:
    out = pathlib.Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    longitudinal, survival = simulate_synthetic_cohort(seed=seed)

    longitudinal.to_csv(out / "longitudinal.csv", index=False)
    survival.to_csv(out / "survival.csv", index=False)

    # Optimal range modeling (restricted to a plausible normal range, e.g., 80-120)
    or_models, risk_df = fit_optimal_range_models(
        longitudinal=longitudinal,
        survival=survival,
        biomarker_col="biomarker",
        time_col="time",
        event_col="event",
        id_col="participant_id",
        adjust_cols=["sex", "age0"],
        normal_range=(80, 120),
    )
    band = detect_optimal_band(risk_df, normal_range=(80, 120), relative_margin=0.02)

    risk_df.to_csv(out / "risk_curve.csv", index=False)
    if band:
        with open(out / "optimal_band.txt", "w") as f:
            f.write(f"Optimal band: {band[0]:.2f} to {band[1]:.2f}\n")

    fig1 = plot_risk_curve(risk_df, normal_range=(80, 120), optimal_band=band, title="Synthetic biomarker risk curve")
    fig1.savefig(out / "risk_curve.png", dpi=150, bbox_inches="tight")

    # Setpoint estimation and personalized intervals
    setpoint_df, variance_components = estimate_setpoints(longitudinal, id_col="participant_id", value_col="biomarker")
    perf_df = evaluate_personalized_intervals(
        longitudinal=longitudinal,
        survival=survival,
        setpoint_df=setpoint_df,
        id_col="participant_id",
        value_col="biomarker",
        time_col="time",
        event_col="event",
        adjust_cols=["sex", "age0"],
    )

    setpoint_df.to_csv(out / "setpoints.csv", index=False)
    pd.DataFrame([variance_components]).to_csv(out / "variance_components.csv", index=False)
    perf_df.to_csv(out / "personalized_vs_population.csv", index=False)

    fig2 = plot_setpoint_variability(longitudinal, setpoint_df, id_col="participant_id", value_col="biomarker")
    fig2.savefig(out / "setpoint_variability.png", dpi=150, bbox_inches="tight")

    # Minimal repeats evaluation
    repeats_df = evaluate_min_repeats(
        longitudinal=longitudinal,
        id_col="participant_id",
        value_col="biomarker",
        max_repeats=6,
    )
    repeats_df.to_csv(out / "min_repeats_accuracy.csv", index=False)

    fig3 = plot_repeats_accuracy(repeats_df)
    fig3.savefig(out / "min_repeats_accuracy.png", dpi=150, bbox_inches="tight")
