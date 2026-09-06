from __future__ import annotations

import hashlib
import json
from pathlib import Path

import geopandas as gpd
import pandas as pd
import requests
from shapely.geometry import MultiPolygon, Polygon

from common import (
    API_CRS, BOROUGHS_DIR, CACHE_DIR, CRIME_DIR, END_MONTH, POLICE_RATE,
    RAW_DIR, START_MONTH, get_logger, months_between,
)

BULK_DIR = RAW_DIR / "crime_bulk"

log = get_logger("fetch_crimes")

POLY_ENDPOINT = "https://data.police.uk/api/crimes-street/all-crime"
MAX_VERTICES = 100


def _cache_path(key: str) -> Path:
    return CACHE_DIR / f"{hashlib.sha1(key.encode()).hexdigest()}.json"


def _post(params: dict) -> list:
    key = f"{POLY_ENDPOINT}?{sorted(params.items())}"
    cp = _cache_path(key)
    if cp.exists():
        return json.loads(cp.read_text())
    POLICE_RATE.acquire()
    r = requests.post(POLY_ENDPOINT, data=params, timeout=60)
    if r.status_code == 503:
        log.warning("503 for date=%s (poly too large, > 10k crimes)", params.get("date"))
        return []
    r.raise_for_status()
    data = r.json()
    cp.write_text(json.dumps(data))
    return data


def _simplify(geom, max_vertices: int) -> Polygon:
    if isinstance(geom, MultiPolygon):
        geom = max(geom.geoms, key=lambda g: g.area)
    tol = 10.0
    while len(geom.exterior.coords) > max_vertices and tol < 10_000:
        geom = geom.simplify(tol, preserve_topology=True)
        tol *= 2
    return geom


def _poly_str(geom_wgs84: Polygon) -> str:
    return ":".join(f"{y},{x}" for x, y in geom_wgs84.exterior.coords)


def _fetch_month(month: str, boroughs_wgs84: gpd.GeoDataFrame) -> list[dict]:
    crimes: list[dict] = []
    for _, row in boroughs_wgs84.iterrows():
        name = row.get("NAME") or row.get("name") or "borough"
        poly = _poly_str(_simplify(row.geometry, MAX_VERTICES))
        data = _post({"poly": poly, "date": month})
        for c in data:
            c["_borough"] = name
        crimes.extend(data)
        log.info("  %s %s: %d", month, name, len(data))
    return crimes


def _ingest_bulk(month: str) -> int:
    subdir = BULK_DIR / month
    if not subdir.is_dir():
        return 0
    csvs = sorted(subdir.glob("*-street.csv")) or sorted(subdir.glob("*.csv"))
    if not csvs:
        return 0
    frames = []
    for c in csvs:
        df = pd.read_csv(c, low_memory=False)
        rename = {
            "Longitude": "location.longitude",
            "Latitude":  "location.latitude",
            "Month":     "month",
            "Crime type":"category",
            "Location":  "location.street.name",
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
        keep = [c for c in ("location.longitude", "location.latitude", "month",
                            "category", "location.street.name") if c in df.columns]
        frames.append(df[keep])
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.dropna(subset=["location.longitude", "location.latitude"])
    out = CRIME_DIR / f"crimes_{month}.parquet"
    combined.to_parquet(out, index=False)
    log.info("BULK %s: wrote %s (%d rows from %d files)",
             month, out.name, len(combined), len(csvs))
    return len(combined)


def run() -> None:
    boroughs_gpkg = BOROUGHS_DIR / "boroughs.gpkg"
    if not boroughs_gpkg.exists():
        raise SystemExit(f"Missing {boroughs_gpkg}; run fetch_boroughs.py first")
    boroughs = gpd.read_file(boroughs_gpkg).to_crs(API_CRS)
    log.info("Boroughs: %d rows", len(boroughs))
    if BULK_DIR.is_dir():
        log.info("Bulk-archive dir present at %s — will use per-month subdirs where available", BULK_DIR)

    for month in months_between(START_MONTH, END_MONTH):
        out = CRIME_DIR / f"crimes_{month}.parquet"
        if out.exists():
            log.info("SKIP %s (already ingested)", out.name)
            continue
        if _ingest_bulk(month):
            continue
        log.info("Month %s (live API)", month)
        crimes = _fetch_month(month, boroughs)
        if not crimes:
            log.warning("No crimes for %s — skip", month)
            continue
        pd.json_normalize(crimes).to_parquet(out, index=False)
        log.info("Wrote %s (%d rows)", out, len(crimes))


if __name__ == "__main__":
    run()
