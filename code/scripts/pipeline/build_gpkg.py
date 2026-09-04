from __future__ import annotations

import re

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from common import (
    API_CRS, BOROUGHS_DIR, CRIME_DIR, FINAL_GPKG, TARGET_CRS, TEMP_DIR, get_logger,
)

log = get_logger("build_gpkg")


def _load_crimes() -> pd.DataFrame:
    files = sorted(CRIME_DIR.glob("crimes_*.parquet"))
    if not files:
        raise SystemExit(f"No crime parquets in {CRIME_DIR}; run fetch_crimes.py first")
    log.info("Loading %d crime files", len(files))
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def _load_lsoa() -> gpd.GeoDataFrame:
    p = BOROUGHS_DIR / "lsoa_2011.gpkg"
    if not p.exists():
        raise SystemExit(f"Missing {p}; run fetch_boroughs.py first")
    return gpd.read_file(p).to_crs(TARGET_CRS)


def _load_temperature() -> pd.DataFrame:
    p = TEMP_DIR / "temperature_by_lsoa_month.parquet"
    if not p.exists():
        raise SystemExit(f"Missing {p}; run aggregate_temperature.py first")
    return pd.read_parquet(p)


def _lsoa_code_col(gdf) -> str:
    for c in ("LSOA11CD", "lsoa_code", "LSOA_CODE"):
        if c in gdf.columns:
            return c
    raise SystemExit(f"LSOA missing code column. Have: {list(gdf.columns)}")


def run() -> None:
    crimes = _load_crimes()
    lsoa = _load_lsoa()
    code_col = _lsoa_code_col(lsoa)

    lon = crimes["location.longitude"].astype(float)
    lat = crimes["location.latitude"].astype(float)
    pts = gpd.GeoDataFrame(
        crimes.assign(geometry=[Point(x, y) for x, y in zip(lon, lat)]),
        geometry="geometry", crs=API_CRS,
    ).to_crs(TARGET_CRS)

    log.info("Spatial join %d crimes → %d LSOAs", len(pts), len(lsoa))
    joined = gpd.sjoin(pts, lsoa[[code_col, "geometry"]],
                       how="inner", predicate="within")
    joined = joined.rename(columns={code_col: "lsoa_code"})
    joined["year"] = joined["month"].str.slice(0, 4).astype(int)
    joined["month_num"] = joined["month"].str.slice(5, 7).astype(int)

    keys = ["lsoa_code", "year", "month_num"]
    counts = (joined.groupby(keys, as_index=False)
                    .size().rename(columns={"size": "crime_count"}))
    counts["log_crime_count"] = np.log1p(counts["crime_count"])

    # Per-category counts. Categories respond to temperature with opposite
    # signs -- bicycle theft and anti-social behaviour rise with it while
    # burglary and theft-from-the-person fall -- so a single total cancels
    # real signal. Keep the total for continuity and add one column per
    # category; an absent category in a present LSOA-month is a true zero.
    cats = joined.assign(category=joined["category"].fillna("unknown"))
    wide = (cats.groupby(keys + ["category"], as_index=False).size()
                .pivot_table(index=keys, columns="category", values="size",
                             fill_value=0)
                .reset_index())
    wide.columns.name = None
    cat_src = [c for c in wide.columns if c not in keys]
    ren = {c: "n_" + re.sub(r"[^0-9a-zA-Z]+", "_", c).strip("_").lower()
           for c in cat_src}
    wide = wide.rename(columns=ren)
    cat_cols = sorted(ren.values())
    counts = counts.merge(wide, on=keys, how="left")
    counts[cat_cols] = counts[cat_cols].fillna(0).astype(int)
    log.info("Aggregated: %d LSOA×month rows, %d crime categories",
             len(counts), len(cat_cols))

    temp = _load_temperature()
    temp[["year", "month_num"]] = temp["date"].str.split("-", expand=True).astype(int)
    merged = counts.merge(
        temp[["lsoa_code", "year", "month_num", "mean_temperature"]],
        on=["lsoa_code", "year", "month_num"], how="left",
    )

    # Residential population ships with the LSOA boundary file and was being
    # discarded. POPDEN is renamed to the name the model already expects.
    pop_map = {"USUALRES": "population", "POPDEN": "population_density",
               "HHOLDS": "households"}
    have_pop = {k: v for k, v in pop_map.items() if k in lsoa.columns}
    lsoa_out = (lsoa.rename(columns={code_col: "lsoa_code", **have_pop})
                    [["lsoa_code", "geometry"] + list(have_pop.values())])
    merged = merged.merge(lsoa_out, on="lsoa_code", how="left")
    merged = merged.rename(columns={"month_num": "month"})

    keep = (["lsoa_code", "year", "month", "crime_count", "log_crime_count",
             "mean_temperature"] + list(have_pop.values()) + cat_cols
            + ["geometry"])
    gdf = gpd.GeoDataFrame(merged[keep], geometry="geometry", crs=TARGET_CRS)

    if FINAL_GPKG.exists():
        FINAL_GPKG.unlink()
    gdf.to_file(FINAL_GPKG, driver="GPKG")
    log.info("Wrote %s (%d rows, %d cols)", FINAL_GPKG, len(gdf), len(gdf.columns))

    n_nan_temp = int(gdf["mean_temperature"].isna().sum())
    if n_nan_temp:
        log.warning("%d rows have NA mean_temperature", n_nan_temp)


if __name__ == "__main__":
    run()
