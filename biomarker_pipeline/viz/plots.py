from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


sns.set_context("talk")


def plot_risk_curve(
    risk_df: pd.DataFrame,
    normal_range: Optional[Tuple[float, float]] = None,
    optimal_band: Optional[Tuple[float, float]] = None,
    title: str = "Risk curve",
):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(risk_df["biomarker"], risk_df["relative_risk"], color="#1f77b4", lw=2)
    ax.set_xlabel("Biomarker")
    ax.set_ylabel("Relative hazard")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    if normal_range is not None:
        ax.axvspan(normal_range[0], normal_range[1], color="#2ca02c", alpha=0.08, label="Normal range")
    if optimal_band is not None:
        ax.axvspan(optimal_band[0], optimal_band[1], color="#ff7f0e", alpha=0.15, label="Optimal band")

    ax.legend(frameon=False)
    return fig


def plot_setpoint_variability(
    longitudinal: pd.DataFrame,
    setpoint_df: pd.DataFrame,
    id_col: str,
    value_col: str,
):
    # Sample 50 participants for visualization
    sample_ids = (
        longitudinal[id_col]
        .drop_duplicates()
        .sample(n=min(50, longitudinal[id_col].nunique()), random_state=1)
        .tolist()
    )
    sub = longitudinal[longitudinal[id_col].isin(sample_ids)].copy()
    sub = sub.merge(setpoint_df[[id_col, "setpoint_estimate"]], on=id_col, how="left")

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.lineplot(data=sub, x="time", y=value_col, hue=id_col, alpha=0.4, legend=False, ax=ax)
    sns.lineplot(data=sub, x="time", y="setpoint_estimate", color="black", lw=2, ax=ax)
    ax.set_title("Within-person trajectories and setpoint estimate")
    ax.set_xlabel("Time")
    ax.set_ylabel("Biomarker")
    ax.grid(True, alpha=0.3)
    return fig


def plot_repeats_accuracy(repeats_df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.lineplot(data=repeats_df, x="k", y="mae", marker="o", label="MAE", ax=ax)
    sns.lineplot(data=repeats_df, x="k", y="rmse", marker="o", label="RMSE", ax=ax)
    ax.set_xlabel("Number of repeated measurements (k)")
    ax.set_ylabel("Error")
    ax.set_title("Accuracy vs number of repeats")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    return fig
