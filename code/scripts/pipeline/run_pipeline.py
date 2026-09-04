from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from common import FINAL_GPKG, REPO_ROOT, get_logger

log = get_logger("pipeline")

HERE = Path(__file__).parent
STEPS = [
    ("fetch_boroughs.py",       "boundaries"),
    ("fetch_crimes.py",         "crime records"),
    ("aggregate_temperature.py","temperature"),
    ("build_gpkg.py",           "assemble gpkg"),
]

R_SCRIPT = REPO_ROOT / "code" / "model" / "GTWR_crime.r"


def _run_py(script: Path) -> None:
    log.info("→ %s", script.name)
    p = subprocess.run([sys.executable, str(script)], cwd=HERE)
    if p.returncode:
        raise SystemExit(f"{script.name} failed with exit {p.returncode}")


def _run_r() -> None:
    log.info("→ Rscript %s", R_SCRIPT)
    env = {
        "CRIME_DATA_GPKG": str(FINAL_GPKG),
        "GTWR_MINIMAL": "1",
        "GTWR_SAMPLE_N": os.environ.get("GTWR_SAMPLE_N", "800"),
    }
    full_env = {**os.environ, **env}
    if "GWMODEL_REPO" not in full_env:
        candidate = REPO_ROOT.parent / "GWmodel"
        if (candidate / "R" / "gtwr.R").exists():
            full_env["GWMODEL_REPO"] = str(candidate)
    p = subprocess.run(["Rscript", str(R_SCRIPT)], cwd=REPO_ROOT, env=full_env)
    if p.returncode:
        raise SystemExit(f"R script failed with exit {p.returncode}")


def main() -> None:
    for name, label in STEPS:
        log.info("=== %s (%s) ===", name, label)
        _run_py(HERE / name)
    log.info("=== GTWR fit ===")
    _run_r()
    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()
