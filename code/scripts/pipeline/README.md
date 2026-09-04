# Data pipeline

End-to-end ETL that produces `data/finished_data.gpkg` and runs the R GTWR fit.
No CLI flags, no threading, no fluff — one command replicates the run.

## Layout

```
code/scripts/pipeline/
├── common.py                # paths, CRS, rate limiter, month iterator, date window
├── fetch_boroughs.py        # London borough + LSOA boundaries (data.london.gov.uk)
├── fetch_crimes.py          # UK police crimes-street/all-crime, cached under cache/
├── aggregate_temperature.py # HadUK-Grid raster → LSOA×month (falls back to synthetic)
├── build_gpkg.py            # join everything → data/finished_data.gpkg
└── run_pipeline.py          # runs all four + Rscript code/model/GTWR_crime.r
```

## Install

Python:
```bash
python -m pip install requests geopandas pandas pyarrow rasterio rasterstats shapely
```

R (from the sibling patched GWmodel checkout):
```bash
R -e 'install.packages(c("sf","sp","spdep","spatialreg","tidyverse","tmap","mapview","car","cowplot","leafsync","leaflet.extras2"))'
```

## Run

```bash
cd code/scripts/pipeline
python run_pipeline.py
```

That executes:
1. `fetch_boroughs.py` — downloads the London statistical GIS boundary zip (cached), writes `data/raw/boroughs/{boroughs,lsoa_2011}.gpkg`.
2. `fetch_crimes.py` — iterates months from `START_MONTH` to `END_MONTH` (in `common.py`, defaults `2011-01`..`2013-12`), one POST per (borough, month) to `data.police.uk/api/crimes-street/all-crime`. Rate-limited at 14 req/s. Cached under `cache/*.json` — re-runs are free.
3. `aggregate_temperature.py` — if a NetCDF/GeoTIFF is in `data/raw/temperature/`, does zonal stats per LSOA per band. Otherwise falls back to London monthly climatology + tiny per-LSOA jitter (logged warning).
4. `build_gpkg.py` — spatial-joins crime points to LSOAs, aggregates counts per LSOA × month, merges temperature, writes `data/finished_data.gpkg`.
5. `Rscript code/model/GTWR_crime.r` — with `GTWR_MINIMAL=1` (fits `crime_count ~ mean_temperature`) and `GWMODEL_REPO` auto-pointed at the sibling `GWmodel` checkout. Writes fitted model + diagnostics under `results/`.

## Historic crime data (pre-2023)

`data.police.uk`'s live API only serves the last ~3 years. For older data, drop the extracted monthly folders from the [Police bulk archive](https://data.police.uk/data/archive/) into `data/raw/crime_bulk/`. Layout:

```
data/raw/crime_bulk/
├── 2011-01/
│   ├── 2011-01-metropolitan-street.csv
│   └── 2011-01-city-of-london-street.csv
├── 2011-02/
│   └── …
└── …
```

`fetch_crimes.py` uses those in preference to the live API when a month's subdirectory exists. Any bulk CSV with columns `Longitude`, `Latitude`, `Month`, `Crime type` will do — force names don't matter.

## Real temperature data

Drop a HadUK-Grid monthly mean temperature raster (NetCDF) into `data/raw/temperature/`. Get it from [CEDA — HadUK-Grid 1km gridded temperature](https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb) (free CEDA registration needed). Any file with suffix `.nc`, `.nc4`, `.tif`, or `.tiff` in that directory will be auto-detected and used instead of the synthetic fallback.

## Sources

- Boundaries: [data.london.gov.uk statistical GIS boundaries](https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london)
- Crimes: [data.police.uk](https://data.police.uk/) — [API docs](https://data.police.uk/docs/), [rate limits](https://data.police.uk/docs/api-call-limits/), [crimes-street/all-crime](https://data.police.uk/docs/method/crimes-street/)
- Temperature: [Met Office HadUK-Grid via CEDA](https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb) (or [Met Office DataHub](https://www.metoffice.gov.uk/services/data/met-office-weather-datahub) if you have an API key)

## Notes

- Changing the date window: edit `START_MONTH` / `END_MONTH` in `common.py`. Cached responses from previous windows stay valid — only new months hit the API.
- The full R model uses 16 covariates. Only `mean_temperature` and `crime_count` come out of this pipeline; the others land as NA. `GTWR_MINIMAL=1` (set automatically by `run_pipeline.py`) makes the R script fit `crime_count ~ mean_temperature` only. Drop `GTWR_MINIMAL` if you have augmented the GPKG with the other covariates.
