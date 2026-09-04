from __future__ import annotations

import io
import zipfile
from pathlib import Path

import geopandas as gpd
import requests

from common import BOROUGHS_DIR, TARGET_CRS, get_logger

log = get_logger("fetch_boroughs")

ZIP_URL = (
    "https://data.london.gov.uk/download/statistical-gis-boundary-files-london/"
    "9ba8c833-6370-4b11-abdc-314aa020d5e0/statistical-gis-boundaries-london.zip"
)


def _download_zip() -> bytes:
    zip_cache = BOROUGHS_DIR / "statistical-gis-boundaries-london.zip"
    if zip_cache.exists():
        log.info("Zip already cached at %s", zip_cache)
        return zip_cache.read_bytes()
    log.info("GET %s", ZIP_URL)
    r = requests.get(ZIP_URL, timeout=180)
    r.raise_for_status()
    zip_cache.write_bytes(r.content)
    return r.content


def _extract_layer(zip_bytes: bytes, name_hint: str, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        members = [n for n in zf.namelist() if name_hint in n.lower() and not n.endswith("/")]
        if not members:
            raise RuntimeError(f"No members matching '{name_hint}' in zip")
        for m in members:
            (dest / Path(m).name).write_bytes(zf.read(m))
    return next(dest.glob("*.shp"))


def _to_gpkg(shp: Path, out: Path) -> None:
    gdf = gpd.read_file(shp).to_crs(TARGET_CRS)
    log.info("%s: %d rows, columns=%s", out.name, len(gdf), list(gdf.columns))
    if out.exists():
        out.unlink()
    gdf.to_file(out, driver="GPKG")
    log.info("Wrote %s", out)


def run() -> None:
    boroughs_gpkg = BOROUGHS_DIR / "boroughs.gpkg"
    lsoa_gpkg = BOROUGHS_DIR / "lsoa_2011.gpkg"
    if boroughs_gpkg.exists() and lsoa_gpkg.exists():
        log.info("Both boundary gpkgs already exist — skip")
        return
    zbytes = _download_zip()
    if not boroughs_gpkg.exists():
        borough_shp = _extract_layer(zbytes, "london_borough", BOROUGHS_DIR / "borough_raw")
        _to_gpkg(borough_shp, boroughs_gpkg)
    if not lsoa_gpkg.exists():
        lsoa_shp = _extract_layer(zbytes, "lsoa_2011", BOROUGHS_DIR / "lsoa_raw")
        _to_gpkg(lsoa_shp, lsoa_gpkg)


if __name__ == "__main__":
    run()
