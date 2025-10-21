import os
import pathlib
import typing as t

import typer

from biomarker_pipeline.pipeline.run import run_pipeline
from biomarker_pipeline.pipeline.synthetic import run_synthetic_demo

app = typer.Typer(add_completion=False, help="Biomarker Optimal Range & Setpoint pipeline")


@app.command()
def demo(
    output_dir: str = typer.Option("outputs/demo", "--output-dir", help="Output directory"),
    seed: int = typer.Option(42, help="Random seed for reproducibility"),
):
    """Run synthetic end-to-end demo: simulate, analyze, and plot."""
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    run_synthetic_demo(output_dir=output_dir, seed=seed)


@app.command()
def run(
    config: str = typer.Option(..., "--config", help="YAML config path"),
    output_dir: str = typer.Option("outputs/run", "--output-dir", help="Output directory"),
):
    """Run pipeline using a YAML config (real data or demo)."""
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    run_pipeline(config_path=config, output_dir=output_dir)


if __name__ == "__main__":
    app() 
