from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from lifelines import CoxPHFitter


def estimate_setpoints(
    longitudinal: pd.DataFrame,
    id_col: str,
    value_col: str,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Estimate individual setpoints via random-intercept mixed model and extract variance components.

    Returns:
      setpoint_df: participant_id, setpoint_estimate, n_obs
      variance_components: dict with between_subject_var and within_subject_var
    """
    # Mixed model with random intercept per individual (no fixed covariates for simplicity)
    # Note: statsmodels MixedLM requires an explicit intercept and groups.
    df = longitudinal[[id_col, value_col]].dropna().copy()
    df["Intercept"] = 1.0

    # MixedLM: value ~ 1 + (1 | id)
    model = sm.MixedLM(endog=df[value_col], exog=df[["Intercept"]], groups=df[id_col])
    result = model.fit(method="lbfgs", reml=True, disp=False)

    # Empirical Bayes posterior means for random intercepts (BLUPs)
    re = result.random_effects  # dict id -> array([b_i])
    overall_mean = float(result.fe_params["Intercept"])  # population mean

    setpoints = []
    for pid, b in re.items():
        setpoints.append((pid, overall_mean + float(b["Intercept"])) )

    setpoint_df = (
        pd.DataFrame(setpoints, columns=[id_col, "setpoint_estimate"])
        .merge(df.groupby(id_col).size().rename("n_obs"), on=id_col, how="left")
        .reset_index(drop=True)
    )

    var_between = float(result.cov_re.iloc[0, 0])
    var_within = float(result.scale)
    variance_components = {
        "between_subject_var": var_between,
        "within_subject_var": var_within,
        "icc": var_between / (var_between + var_within) if (var_between + var_within) > 0 else np.nan,
        "overall_mean": overall_mean,
    }

    return setpoint_df, variance_components


def evaluate_personalized_intervals(
    longitudinal: pd.DataFrame,
    survival: pd.DataFrame,
    setpoint_df: pd.DataFrame,
    id_col: str,
    value_col: str,
    time_col: str,
    event_col: str,
    adjust_cols: List[str],
) -> pd.DataFrame:
    """
    Compare risk prediction using population-normalized deviation vs personalized (setpoint) deviation.

    Builds two Cox models with cubic spline of deviation and compares concordance index.
    """
    # Merge last observed value and setpoint
    last_val = (
        longitudinal.sort_values([id_col, time_col]).groupby(id_col).tail(1)[[id_col, value_col]].rename(columns={value_col: "value_last"})
    )
    df = survival.merge(last_val, on=id_col, how="left").merge(setpoint_df, on=id_col, how="left")

    # Population mean deviation
    pop_mean = df["value_last"].mean()
    df["dev_population"] = df["value_last"] - pop_mean

    # Personalized deviation from setpoint
    df["dev_personal"] = df["value_last"] - df["setpoint_estimate"]

    # Build and fit Cox models (simple linear+quadratic terms to capture U-shape)
    covs = [c for c in adjust_cols if c in df.columns]

    def fit_and_score(feature: str) -> float:
        use = df[["time", event_col, feature] + covs].dropna()
        # Quadratic expansion
        use[f"{feature}_2"] = use[feature] ** 2
        cph = CoxPHFitter()
        cph.fit(use, duration_col="time", event_col=event_col)
        return float(cph.concordance_index_)

    c_idx_pop = fit_and_score("dev_population")
    c_idx_per = fit_and_score("dev_personal")

    return pd.DataFrame(
        {
            "model": ["population_deviation", "personal_deviation"],
            "c_index": [c_idx_pop, c_idx_per],
        }
    )
