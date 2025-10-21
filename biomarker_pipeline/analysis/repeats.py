from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd


def _estimate_setpoint_from_k(values: np.ndarray) -> float:
    # Simple estimator: average of first k values
    return float(values.mean())


def evaluate_min_repeats(
    longitudinal: pd.DataFrame,
    id_col: str,
    value_col: str,
    max_repeats: int = 6,
) -> pd.DataFrame:
    """
    For each k=1..max_repeats, estimate setpoint as mean of first k visits, then
    evaluate accuracy vs "gold" estimate (mean of all visits).

    Returns wide and long forms with MAE and RMSE.
    """
    df = longitudinal[[id_col, value_col, "visit"]].dropna().copy()
    df = df.sort_values([id_col, "visit"])  # rely on numeric visit order

    gold = df.groupby(id_col)[value_col].mean().rename("setpoint_gold")
    results = []

    for k in range(1, max_repeats + 1):
        est_k = (
            df.groupby(id_col)
            .head(k)
            .groupby(id_col)[value_col]
            .mean()
            .rename("setpoint_k")
        )
        merged = pd.concat([gold, est_k], axis=1).dropna()
        err = merged["setpoint_k"] - merged["setpoint_gold"]
        mae = float(err.abs().mean())
        rmse = float(np.sqrt((err ** 2).mean()))
        results.append({"k": k, "mae": mae, "rmse": rmse})

    return pd.DataFrame(results)
