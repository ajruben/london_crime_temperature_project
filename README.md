# London crime & temperature

[Project page](https://www.rubenswarts.nl/projects/2.html) ·
[rubenswarts.nl](https://rubenswarts.nl)

Monthly recorded crime for London LSOAs, joined to gridded temperature and
modelled with geographically and temporally weighted regression (GTWR).

This project is set up based on the joint work of my and my peers during a school project. we did the research and data work without AI. 
In 2025, I found a performance bottleneck in `GWmodel`'s GTWR implementation and replaced its nested loops with vectorised arithmetic. 
I later lost that code while moving between machines. 
In this repo, I used AI to help rebuild the fix in R, make the work reproducible and benchmark it.
Running the model on the data we had was not possible before the improvements.

## Authorship and use of AI

### My work, without AI

- Finding and assessing the data sources.
- Extracting the source data.
- Cleaning, reshaping and joining the data.
- Exploring the data and developing the modelling approach.
- Finding the GTWR performance bottleneck in 2025.
- Designing and implementing the original vectorisation of the nested GTWR
  loops.

Some of this work is in `code/data_munching_and_exploration/`.

### Work I used AI for

- Turning my data work into a reproducible pipeline.
- Rewriting my lost vectorisation fix in R, based on my description of how it
  worked.
- Finding more performance improvements. The main one was processing the
  distance matrix in **blocks**, along with related memory fixes. This made it
  possible to fit all 34,013 rows.

I directed and checked the AI-assisted work. AI did not choose the data, do the
original extraction or exploration, find the original GTWR bottleneck, or
design the original vectorisation. Blocking and the related extra fixes came
from the AI-assisted rewrite. See the
[`GWmodel` fork](https://github.com/ajruben/GWmodel/tree/vectorize-gtwr) for the
implementation.

## Why AI rewrote the vectorisation

I lost the code for my original GTWR fix while moving between machines. When I
returned to the project, my R was rusty, so I used AI to rewrite the fix from
my description of it. During the rewrite, it also found the blocking and memory
improvements.

## Scale

The full panel is 34,013 LSOA-months. A full-data fit takes ~13 min and peaks
around 52 GB of RAM on 15 cores; bandwidth selection is ~85% of that time (it
is still a serial loop in `bw.gtwr`). `GTWR_SAMPLE_N` subsamples for quicker
iteration. A sample of 14,000 rows runs in ~5 min using under 21 GB.

These timings are for the rebuilt version of my original vectorisation,
including the extra AI-assisted fixes.

## Reproducing the analysis

Both repos must sit **side by side**. The pipeline finds the patched GWmodel
fork in the sibling directory:

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
(Rtools45 covers R 4.5 and 4.6). GWmodel contains C++ and will not install
without a toolchain. `winget install RProject.Rtools` is enough.

First run fetches boundaries, crime and temperature; responses are cached under
`cache/` keyed by request, so re-runs resume rather than refetch. The crime API
(`data.police.uk`) returns intermittent 502s. Re-run it and the cache lets it
pick up where it stopped.

**The crime API serves a rolling 36-month window** (as of 2026-09: 2023-08 to
2026-07). `START_MONTH`/`END_MONTH` in `common.py` currently ask for 2024, so
that becomes unavailable from around 2027-01. After that, archived data must be
downloaded separately.

## Reproducible pipeline

`run_pipeline.py` runs four Python steps then the R fit:

| Step | Output |
|---|---|
| `fetch_boroughs.py` | borough + LSOA boundaries |
| `fetch_crimes.py` | `data/raw/crime/crimes_<month>.parquet` |
| `aggregate_temperature.py` | temperature per LSOA-month |
| `build_gpkg.py` | `data/finished_data.gpkg` |
| `code/model/GTWR_crime.r` | `results/gtwr_fit_<stamp>.rds` + diagnostics |

I used AI to set up this pipeline. The data choices, transformations and
modelling decisions came from my earlier work.

## Modelling data

`finished_data.gpkg` is one row per LSOA-month with:

- `crime_count`, `log_crime_count`: total recorded crime
- `n_<type>`: a count per crime type (12 of them: `n_burglary`,
  `n_bicycle_theft`, `n_violent_crime`, …). Categories non-zero in under
  15% of LSOA-months are dropped at build time. Below that there is no
  within-LSOA variation left to fit.
- `mean_temperature`
- `population`, `population_density`, `households` (2011 census, from the LSOA
  boundary file)

**Model crime types separately.** They respond to temperature with opposite
signs. Bicycle theft rises about +3.4%/degC and anti-social behaviour +2.5%,
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

