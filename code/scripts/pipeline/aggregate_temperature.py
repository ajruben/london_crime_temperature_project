from __future__ import annotations

import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from common import (
    BOROUGHS_DIR, END_MONTH, START_MONTH, TEMP_DIR, get_logger, months_between,
)

log = get_logger("aggregate_temperature")

OUT = TEMP_DIR / "temperature_by_lsoa_month.parquet"
SOURCE_TAG = TEMP_DIR / "temperature_source.txt"

MONTHLY_CLIMATOLOGY_C = {
    1: 5.4, 2: 5.6, 3: 7.6, 4: 9.7, 5: 12.9, 6: 15.9,
    7: 18.0, 8: 17.7, 9: 15.0, 10: 11.6, 11: 8.0, 12: 5.7,
}

RASTER_SUFFIXES = (".nc", ".nc4", ".tif", ".tiff")


def _find_raster() -> Path | None:
    for p in sorted(TEMP_DIR.iterdir()):
        if p.suffix.lower() in RASTER_SUFFIXES:
            return p
    return None


def _lsoa_gdf() -> gpd.GeoDataFrame:
    p = BOROUGHS_DIR / "lsoa_2011.gpkg"
    if not p.exists():
        raise SystemExit(f"Missing {p}; run fetch_boroughs.py first")
    return gpd.read_file(p)


def _lsoa_code_col(gdf: gpd.GeoDataFrame) -> str:
    for c in ("LSOA11CD", "lsoa_code", "LSOA_CODE"):
        if c in gdf.columns:
            return c
    raise SystemExit(f"LSOA gpkg missing code column. Have: {list(gdf.columns)}")


def _from_raster(raster: Path) -> pd.DataFrame:
    import rasterio
    from rasterstats import zonal_stats

    lsoa = _lsoa_gdf()
    code_col = _lsoa_code_col(lsoa)

    with rasterio.open(raster) as ds:
        n_bands = ds.count
        if ds.crs is not None and lsoa.crs != ds.crs:
            lsoa = lsoa.to_crs(ds.crs)
        log.info("Raster: %s bands=%d crs=%s", raster.name, n_bands, ds.crs)

    dates = list(months_between(START_MONTH, END_MONTH))[:n_bands]
    rows: list[dict] = []
    for i, tag in enumerate(dates, start=1):
        log.info("Zonal stats band %d/%d (%s)", i, len(dates), tag)
        stats = zonal_stats(vectors=lsoa.geometry, raster=str(raster),
                            band=i, stats=["mean"], all_touched=True)
        for j, s in enumerate(stats):
            rows.append({"lsoa_code": lsoa.iloc[j][code_col],
                         "date": tag,
                         "mean_temperature": s["mean"] if s else None})
    return pd.DataFrame(rows)


def _synthetic() -> pd.DataFrame:
    log.warning("No temperature raster found in %s — synthesising from London monthly "
                "climatology with per-LSOA jitter. Drop a HadUK-Grid NetCDF in that "
                "directory to use real data.", TEMP_DIR)
    lsoa = _lsoa_gdf()
    code_col = _lsoa_code_col(lsoa)
    lsoa = lsoa.to_crs("EPSG:27700")
    centroids = lsoa.geometry.centroid
    lat_proxy = (centroids.y - centroids.y.min()) / (centroids.y.max() - centroids.y.min())
    rng = np.random.default_rng(42)
    jitter = rng.normal(0, 0.15, size=len(lsoa))

    rows: list[dict] = []
    for tag in months_between(START_MONTH, END_MONTH):
        _, m = map(int, tag.split("-"))
        base = MONTHLY_CLIMATOLOGY_C[m]
        vals = base - 0.6 * lat_proxy.to_numpy() + jitter
        for code, v in zip(lsoa[code_col].to_numpy(), vals):
            rows.append({"lsoa_code": code, "date": tag, "mean_temperature": float(v)})
    return pd.DataFrame(rows)


SYNTHETIC_REFUSAL = """No temperature raster found in {d}.

mean_temperature is this project's independent variable, so the pipeline
will not silently invent it. Either:

  * put a HadUK-Grid monthly mean-temperature NetCDF (1 km, tas) covering
    {start}..{end} in that directory -- free from CEDA, registration
    required: https://catalogue.ceda.ac.uk/ (search "HadUK-Grid"); or
  * set ALLOW_SYNTHETIC_TEMPERATURE=1 to fall back to London monthly
    climatology plus a fixed per-LSOA offset.

The synthetic series is a relabelling of the calendar month: it holds 12
distinct values plus a constant per LSOA, so any "temperature effect"
estimated from it is a seasonality effect, and its spatial variation is
a north-south ramp plus seeded noise."""


def _note_source(kind: str) -> None:
    SOURCE_TAG.parent.mkdir(parents=True, exist_ok=True)
    SOURCE_TAG.write_text(kind + chr(10), encoding="utf-8")


def run() -> None:
    if OUT.exists():
        # Re-announce provenance every run: the file is written once and
        # reused, so a one-off warning on first build is easy to miss.
        kind = SOURCE_TAG.read_text(encoding="utf-8").strip() if SOURCE_TAG.exists() else "unknown"
        if kind.startswith("synthetic"):
            log.warning("%s holds SYNTHETIC temperature (%s). Delete it and "
                        "supply a raster for real data.", OUT.name, kind)
        else:
            log.info("SKIP %s (already exists, source: %s)", OUT.name, kind)
        return
    raster = _find_raster()
    if raster is None:
        if os.environ.get("ALLOW_SYNTHETIC_TEMPERATURE") != "1":
            raise SystemExit(SYNTHETIC_REFUSAL.format(
                d=TEMP_DIR, start=START_MONTH, end=END_MONTH))
        df = _synthetic()
        _note_source("synthetic: London monthly climatology + seeded per-LSOA offset")
    else:
        df = _from_raster(raster)
        _note_source(f"raster: {raster.name}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    log.info("Wrote %s (%d rows)", OUT, len(df))


if __name__ == "__main__":
    run()
