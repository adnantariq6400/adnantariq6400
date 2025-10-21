from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import pandera as pa
from pandera import Column, DataFrameSchema


# Harmonized longitudinal schema
longitudinal_schema = DataFrameSchema(
    {
        "participant_id": Column(int),
        "visit": Column(int),
        "time": Column(float),
        "age": Column(float),
        "sex": Column(int),
        "biomarker": Column(float),
    },
    coerce=True,
)

# Harmonized survival schema
survival_schema = DataFrameSchema(
    {
        "participant_id": Column(int),
        "time": Column(float),
        "event": Column(int),
        "age0": Column(float),
        "sex": Column(int),
    },
    coerce=True,
)
