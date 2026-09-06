from __future__ import annotations

from pathlib import Path
import logging
import time

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = REPO_ROOT / "cache"
BOROUGHS_DIR = RAW_DIR / "boroughs"
CRIME_DIR = RAW_DIR / "crime"
TEMP_DIR = RAW_DIR / "temperature"
FINAL_GPKG = DATA_DIR / "finished_data.gpkg"

TARGET_CRS = "EPSG:27700"
API_CRS = "EPSG:4326"

START_MONTH = "2024-01"
END_MONTH = "2024-12"

for d in (DATA_DIR, RAW_DIR, CACHE_DIR, BOROUGHS_DIR, CRIME_DIR, TEMP_DIR):
    d.mkdir(parents=True, exist_ok=True)


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(name)


class RateLimiter:
    def __init__(self, max_calls: int, per_seconds: float):
        self.min_interval = per_seconds / max_calls
        self._last = 0.0

    def acquire(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()


POLICE_RATE = RateLimiter(max_calls=14, per_seconds=1.0)


def months_between(start: str, end: str):
    y, m = map(int, start.split("-"))
    ey, em = map(int, end.split("-"))
    while (y, m) <= (ey, em):
        yield f"{y:04d}-{m:02d}"
        m += 1
        if m > 12:
            m = 1
            y += 1
