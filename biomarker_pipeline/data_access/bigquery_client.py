from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

try:
    from google.cloud import bigquery  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    bigquery = None


@dataclass
class BigQueryConfig:
    project: str
    location: Optional[str] = None


class BigQueryClient:
    def __init__(self, config: BigQueryConfig):
        if bigquery is None:
            raise RuntimeError("google-cloud-bigquery not installed. Install optional deps.")
        self.client = bigquery.Client(project=config.project, location=config.location)

    def read_sql(self, sql: str) -> pd.DataFrame:
        job = self.client.query(sql)
        return job.result().to_dataframe(create_bqstorage_client=True)
