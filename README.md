rubenswarts.nl
https://www.rubenswarts.nl/projects/2

## Reproducing the GTWR fit

The GTWR model script (`code/model/GTWR_crime.r`) is wired to a patched local copy of the `GWmodel` R package (vectorized `st.dist`, multi-core calibration, guarded AICc). Two environment variables control where it looks:

- `GWMODEL_REPO` — path to the patched GWmodel checkout (defaults to a sibling directory `../GWmodel` relative to this repo).
- `CRIME_DATA_GPKG` — absolute path to `finished_data.gpkg` (defaults to `<repo>/data/finished_data.gpkg`).

```bash
export GWMODEL_REPO=/path/to/GWmodel
export CRIME_DATA_GPKG=/path/to/finished_data.gpkg
Rscript code/model/GTWR_crime.r
```

The script writes the fitted model and diagnostics under `<repo>/results/` timestamped per run. Use all-but-one physical cores by default; override with `parallel::detectCores` semantics inside the script if needed.
