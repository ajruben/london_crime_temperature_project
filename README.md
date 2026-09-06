rubenswarts.nl
https://www.rubenswarts.nl/projects/2

# London crime & temperature

Monthly recorded crime for London LSOAs joined to gridded temperature, fitted
with a geographically & temporally weighted regression (GTWR).

## Setup from a clean machine

Both repos must sit **side by side** — the pipeline finds the patched GWmodel
fork as a sibling directory:

```
parent/
  GWmodel/                        # branch: vectorize-gtwr
  london_crime_temperature_project/
```

```bash
git clone -b vectorize-gtwr https://github.com/ajruben/GWmodel.git
git clone -b apply-vectorized-gtwr https://github.com/ajruben/london_crime_temperature_project.git
cd london_crime_temperature_project

python -m pip install -r requirements.txt          # Python 3.11+
Rscript code/scripts/setup/install_r_deps.R        # CRAN deps + GWmodel
python code/scripts/pipeline/run_pipeline.py
```

**Requirements:** Python 3.11+, R >= 4.5, and a matching Rtools on Windows
(Rtools45 covers R 4.5 and 4.6) — GWmodel contains C++ and will not install
without a toolchain. `winget install RProject.Rtools` is enough.

First run fetches boundaries, crime and temperature; responses are cached under
`cache/` keyed by request, so re-runs resume rather than refetch. The crime API
(`data.police.uk`) returns intermittent 502s — just re-run, the cache means it
picks up where it stopped.

**The crime API serves a rolling 36-month window** (as of 2026-09: 2023-08 to
2026-07). `START_MONTH`/`END_MONTH` in `common.py` currently ask for 2024, so
that becomes unfetchable from around 2027-01. The per-request cache under
`cache/` is what makes an old window reproducible — but only the boundary
downloads are committed, not the crime responses. Archive
`data/raw/crime/*.parquet` if you need this exact panel later.

## Temperature data: read this first

`aggregate_temperature.py` uses a real raster if one is present in
`data/raw/temperature/` and otherwise **synthesises** the series. The
synthetic fallback is London monthly climatology minus a north-south ramp
plus a seeded per-LSOA offset -- i.e. 12 distinct values and a constant per
LSOA, correlating 0.999 with the hardcoded climatology.

A "temperature effect" estimated from the synthetic series is a **seasonality**
effect wearing degC units, and its spatial variation is a ramp plus
`default_rng(42)` noise, so a geographically weighted temperature coefficient
fits that noise.

The pipeline now refuses to synthesise unless you pass
`ALLOW_SYNTHETIC_TEMPERATURE=1`, records provenance in
`data/raw/temperature/temperature_source.txt`, and warns on every run that
reuses a synthetic file.

For real data use HadUK-Grid monthly mean temperature (1 km, `tas`) from CEDA
(free, registration required: <https://catalogue.ceda.ac.uk/>, search
"HadUK-Grid"); drop the NetCDF in `data/raw/temperature/` and delete
`temperature_by_lsoa_month.parquet` to force a rebuild.

## Pipeline

`run_pipeline.py` runs four Python steps then the R fit:

| Step | Output |
|---|---|
| `fetch_boroughs.py` | borough + LSOA boundaries |
| `fetch_crimes.py` | `data/raw/crime/crimes_<month>.parquet` |
| `aggregate_temperature.py` | temperature per LSOA-month |
| `build_gpkg.py` | `data/finished_data.gpkg` |
| `code/model/GTWR_crime.r` | `results/gtwr_fit_<stamp>.rds` + diagnostics |

## The modelling data

`finished_data.gpkg` is one row per LSOA-month with:

- `crime_count`, `log_crime_count` — total recorded crime
- `n_<type>` — a count per crime type (12 of them: `n_burglary`,
  `n_bicycle_theft`, `n_violent_crime`, …). Categories non-zero in under
  15% of LSOA-months are dropped at build time — below that there is no
  within-LSOA variation left to fit.
- `mean_temperature`
- `population`, `population_density`, `households` (2011 census, from the LSOA
  boundary file)

**Model crime types separately.** They respond to temperature with opposite
signs — bicycle theft rises about +3.4%/degC and anti-social behaviour +2.5%,
while burglary falls -1.4% and theft-from-the-person -2.4%. Summing them into
one total cancels most of the signal.

## Running the model

Environment variables:

| Variable | Default | Meaning |
|---|---|---|
| `GWMODEL_REPO` | `../GWmodel` | patched GWmodel checkout |
| `CRIME_DATA_GPKG` | `<repo>/data/finished_data.gpkg` | input data |
| `GTWR_RESPONSE` | `crime_count` | which count to model, e.g. `n_burglary` |
| `GTWR_SAMPLE_N` | unset = all rows; `run_pipeline.py` passes `800` | subsample size |
| `GTWR_MINIMAL` | `0` | `1` = temperature only; `0` adds any of the covariate list present in the data (currently `population_density`) |
| `GTWR_COVARIATES` | *(unset)* | explicit comma-separated covariates, overriding `GTWR_MINIMAL` |

```bash
export GWMODEL_REPO=/path/to/GWmodel
GTWR_RESPONSE=n_burglary Rscript code/model/GTWR_crime.r
```

Results are written to `results/` timestamped per run. The fit uses all but one
physical core.

## Scale

The full panel is 34,013 LSOA-months. A full-data fit takes ~13 min and peaks
around 52 GB of RAM on 15 cores; bandwidth selection is ~85% of that time (it
is still a serial loop in `bw.gtwr`). `GTWR_SAMPLE_N` subsamples for quicker
iteration — n=14,000 runs in ~5 min under 21 GB.

## Interpreting the fit

Two things to keep in mind before reading much into a high GTWR R-squared:

- 85% of log-crime variance is *which LSOA*, and the local intercept absorbs
  it. Plain LSOA means score R-squared 0.85 with no temperature term at all,
  which is higher than the full GTWR fit. Benchmark against that, not against
  the global regression.
- `mean_temperature` is 99.8% *month* — within any month London spans only
  about 1.3 degC. So it is close to a seasonal indicator, and the temperature
  effect cannot be separated from daylight, school holidays or tourism with a
  single year of monthly data.

---

## Authorship

This project is a mix of my own work and AI-assisted work. The split:

**Mine, without AI**

- Identifying and gathering the data sources.
- Extraction, cleaning and reshaping.
- Exploratory analysis.
- Identifying the GTWR performance bottleneck (2025) and the original refactor
  from nested loops to vectorised arithmetic — see the `GWmodel` fork.

**AI-assisted (Claude Code), under my direction and review**

- Making the pipeline reproducible: pinned dependencies, the R setup script,
  the documented environment, and this README.
- Re-implementing the vectorised GTWR in R, plus a number of additional fixes
  found along the way — most significantly **blocking** the distance matrix
  construction and the related memory work, which is what lets the full
  34,013-row panel be fitted at all. Details in the `GWmodel` fork.
- Splitting recorded crime into per-type counts, carrying population through
  from the boundary file, and making the temperature step refuse to substitute
  a synthetic series without being told to.

**Why the vectorisation was redone with AI rather than restored**

1. I no longer had the code for the original GTWR fix.
2. AI has become good enough that, with my oversight and my familiarity with
   this project, I was comfortable having it redo the work.

The data work and the bottleneck analysis are mine and predate the AI work; the
blocking that made the full panel fit is not.
