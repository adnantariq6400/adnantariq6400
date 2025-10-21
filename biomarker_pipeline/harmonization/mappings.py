from __future__ import annotations

from typing import Dict, Optional, Tuple

import pandas as pd


def apply_column_mapping(
    df: pd.DataFrame,
    mapping: Dict[str, str],
    units_scale: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """
    Rename dataset-specific columns to harmonized names and apply unit scaling.

    Parameters
    - mapping: dict like {"eid": "participant_id", "exam": "visit", "value": "biomarker"}
    - units_scale: dict like {"biomarker": 0.001} to convert to target units
    """
    out = df.rename(columns=mapping).copy()
    if units_scale:
        for col, factor in units_scale.items():
            if col in out.columns:
                out[col] = out[col] * factor
    return out
