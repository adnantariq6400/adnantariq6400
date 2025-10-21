from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from typing import List, Optional, Tuple, Dict, Any


def _prepare_survival_df(
    longitudinal: pd.DataFrame,
    survival: pd.DataFrame,
    biomarker_col: str,
    id_col: str,
    time_col: str,
    adjust_cols: List[str],
) -> pd.DataFrame:
    # Use last available biomarker per participant for baseline risk model.
    last_values = (
        longitudinal.sort_values([id_col, time_col])
        .groupby(id_col)
        .tail(1)[[id_col, biomarker_col]]
        .rename(columns={biomarker_col: f"{biomarker_col}_last"})
    )
    df = survival.merge(last_values, on=id_col, how="left")

    # Add adjustment covariates from survival (assumed merged/available)
    covs = []
    for c in adjust_cols:
        if c in df.columns:
            covs.append(c)
    return df


def fit_optimal_range_models(
    longitudinal: pd.DataFrame,
    survival: pd.DataFrame,
    biomarker_col: str,
    time_col: str,
    event_col: str,
    id_col: str,
    adjust_cols: List[str],
    normal_range: Optional[Tuple[float, float]] = None,
    spline_df: int = 4,
) -> tuple[Dict[str, Any], pd.DataFrame]:
    """
    Fit CoxPH with spline basis for biomarker to estimate non-linear risk.

    Returns:
      models: dict with fitted Cox model and grid
      risk_df: DataFrame with relative risk across biomarker grid
    """
    df = _prepare_survival_df(
        longitudinal, survival, biomarker_col, id_col, time_col, adjust_cols
    )

    biom = f"{biomarker_col}_last"
    use_df = df[[biom, event_col, "time"] + [c for c in adjust_cols if c in df.columns]].dropna()

    # Optionally restrict to normal range for modeling or for risk curve
    if normal_range is not None:
        lo, hi = normal_range
        use_df = use_df[(use_df[biom] >= lo) & (use_df[biom] <= hi)]

    # Build spline basis via patsy
    spline_terms = dmatrix(f"bs({biom}, df={spline_df}, include_intercept=False)", data=use_df, return_type="dataframe")
    design = pd.concat([use_df[[event_col, "time"] + [c for c in adjust_cols if c in use_df.columns]], spline_terms], axis=1)

    cph = CoxPHFitter()
    cph.fit(design, duration_col="time", event_col=event_col)

    # Risk prediction across grid
    grid = np.linspace(use_df[biom].quantile(0.01), use_df[biom].quantile(0.99), 200)
    spline_grid = dmatrix(f"bs(x, df={spline_df}, include_intercept=False)", data={"x": grid}, return_type="dataframe")

    # Center risk at median to get relative hazard
    # Create a dataframe of zeros for covariates other than spline to isolate biomarker effect
    zero_cov = {c: 0.0 for c in design.columns if c not in [event_col, "time"] and not c.startswith("bs(")}
    pred_df = pd.concat([spline_grid, pd.DataFrame([zero_cov] * len(spline_grid))], axis=1)
    log_hazard = cph.predict_log_partial_hazard(pred_df).to_numpy().ravel()
    rel_risk = np.exp(log_hazard - np.median(log_hazard))

    risk_df = pd.DataFrame({"biomarker": grid, "relative_risk": rel_risk})
    if normal_range is not None:
        risk_df = risk_df[(risk_df["biomarker"] >= lo) & (risk_df["biomarker"] <= hi)]

    models = {"cox": cph, "grid": grid}
    return models, risk_df


def detect_optimal_band(
    risk_df: pd.DataFrame,
    normal_range: Optional[Tuple[float, float]] = None,
    relative_margin: float = 0.02,
) -> Optional[Tuple[float, float]]:
    """
    Detect a practical optimal band as the maximal contiguous interval where risk <= (1 + margin) * min risk.
    """
    if risk_df.empty:
        return None
    rr = risk_df["relative_risk"].to_numpy()
    x = risk_df["biomarker"].to_numpy()

    min_rr = rr.min()
    threshold = (1.0 + relative_margin) * min_rr
    mask = rr <= threshold

    if not mask.any():
        return None

    # Find longest contiguous segment
    best_len = 0
    best_interval = None
    start = None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        elif (not m) and start is not None:
            if i - start > best_len:
                best_len = i - start
                best_interval = (x[start], x[i - 1])
            start = None
    if start is not None:
        if len(mask) - start > best_len:
            best_interval = (x[start], x[len(mask) - 1])

    # Intersect with normal range if provided
    if best_interval and normal_range is not None:
        lo, hi = normal_range
        band_lo = max(best_interval[0], lo)
        band_hi = min(best_interval[1], hi)
        if band_lo < band_hi:
            return (float(band_lo), float(band_hi))
        return None

    return best_interval if best_interval else None
