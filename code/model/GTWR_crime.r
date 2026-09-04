.required_pkgs <- c("GWmodel", "sf", "sp", "spdep", "spatialreg",
                    "tmap", "mapview", "car", "RColorBrewer",
                    "cowplot", "leafsync", "leaflet.extras2", "tidyverse")
.missing <- setdiff(.required_pkgs, rownames(installed.packages()))
if (length(.missing)) install.packages(.missing, repos = "https://cloud.r-project.org")
invisible(lapply(.required_pkgs, function(p) suppressPackageStartupMessages(
  library(p, character.only = TRUE))))

resolve_gwmodel_repo <- function() {
  env <- Sys.getenv("GWMODEL_REPO", unset = "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env))
  here <- normalizePath(getwd())
  candidates <- c(
    file.path(dirname(dirname(dirname(here))), "GWmodel"),
    file.path(dirname(here), "GWmodel"),
    file.path(here, "..", "GWmodel"),
    file.path(here, "..", "..", "..", "GWmodel")
  )
  for (p in candidates) if (dir.exists(p) && file.exists(file.path(p, "R", "gtwr.R"))) return(normalizePath(p))
  stop("Could not locate the patched GWmodel repo. Set GWMODEL_REPO=/path/to/GWmodel or run from a sibling of the GWmodel checkout.")
}
gwmodel_repo <- resolve_gwmodel_repo()
message("Sourcing patched GWmodel R files from: ", gwmodel_repo)
sys.source(file.path(gwmodel_repo, "R", "gtwr_parallel.R"), envir = globalenv())
sys.source(file.path(gwmodel_repo, "R", "gtwr.R"),          envir = globalenv())

resolve_data_path <- function() {
  env <- Sys.getenv("CRIME_DATA_GPKG", unset = "")
  if (nzchar(env) && file.exists(env)) return(env)
  here <- normalizePath(getwd())
  candidates <- c(
    file.path(here, "data", "finished_data.gpkg"),   # run from the repo root
    file.path(dirname(dirname(here)), "data",  "finished_data.gpkg"),
    file.path(dirname(dirname(here)), "cache", "finished_data.gpkg"),
    file.path(here, "finished_data.gpkg"),
    file.path(here, "..", "..", "data", "finished_data.gpkg")
  )
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop("Could not locate finished_data.gpkg. Set CRIME_DATA_GPKG=/absolute/path/to/finished_data.gpkg or place it under <repo>/data/.")
}
data_path <- resolve_data_path()
message("Reading crime data from: ", data_path)

crimedata <- st_read(data_path)

# GTWR_RESPONSE picks which count to model: the total ("crime_count", the
# default) or a single crime type such as "n_burglary". Types matter here --
# they respond to temperature with opposite signs, so the total cancels
# signal. Everything is modelled on the log1p scale, matching the
# precomputed log_crime_count.
.response <- Sys.getenv("GTWR_RESPONSE", unset = "crime_count")
if (identical(.response, "crime_count")) {
  if ("log_crime_count" %in% names(crimedata)) {
    if ("crime_count" %in% names(crimedata))
      crimedata <- crimedata %>% select(-crime_count)
    crimedata <- crimedata %>% rename(crime_count = log_crime_count)
  }
} else {
  if (!(.response %in% names(crimedata)))
    stop(sprintf("GTWR_RESPONSE='%s' is not a column. Available: %s",
                 .response,
                 paste(grep("^n_", names(crimedata), value = TRUE), collapse = ", ")))
  message("Response: log1p(", .response, ")")
  crimedata$crime_count <- log1p(as.numeric(crimedata[[.response]]))
}
.sample_n <- suppressWarnings(as.integer(Sys.getenv("GTWR_SAMPLE_N", unset = "")))
if (!is.na(.sample_n) && .sample_n > 0 && .sample_n < nrow(crimedata)) {
  message(sprintf("GTWR_SAMPLE_N=%d — subsampling from %d rows", .sample_n, nrow(crimedata)))
  set.seed(1L)
  crimedata <- crimedata[sample(nrow(crimedata), .sample_n), ]
}

crimedata$year_month <- factor(
  paste(crimedata$year, crimedata$month, sep = "-"),
  levels = unique(paste(crimedata$year, crimedata$month, sep = "-"))
)
crimedata_spt <- as_Spatial(crimedata)

crimedata_spt$time_var    <- as.numeric(crimedata_spt$year_month)
crimedata_spt$crime_count <- as.numeric(crimedata_spt$crime_count)
crimedata_spt$lsoa_code   <- as.factor(crimedata_spt$lsoa_code)
.numeric_cols <- c("mean_temperature", "percent_youth", "population_density",
                   "ethnic_diversity_index_11", "inf_change_1month", "inf_base2015",
                   "Avg_Daylight_Hours", "Vacation_Days", "unemployment_unadj",
                   "median_age", "per_unemployment_lsoa_2011", "per_social_grade_DE",
                   "parks_count_within", "bars_pubs_count_within",
                   "sports_count_within", "pct_park_area_within")
for (.col in intersect(.numeric_cols, names(crimedata_spt@data))) {
  crimedata_spt@data[[.col]] <- as.numeric(crimedata_spt@data[[.col]])
}

time_effect       <- factor(crimedata_spt$time_var)
individual_effect <- factor(crimedata_spt$lsoa_code)

vars_full <- c("mean_temperature",
               "percent_youth",
               "population_density",
               "ethnic_diversity_index_11",
               "inf_change_1month",
               "inf_base2015",
               "Avg_Daylight_Hours",
               "Vacation_Days",
               "unemployment_unadj",
               "median_age",
               "per_unemployment_lsoa_2011",
               "per_social_grade_DE",
               "parks_count_within",
               "bars_pubs_count_within",
               "sports_count_within",
               "pct_park_area_within")
minimal <- isTRUE(Sys.getenv("GTWR_MINIMAL", unset = "0") == "1")
vars <- if (minimal) c("mean_temperature") else vars_full
available <- intersect(vars, colnames(crimedata_spt@data))
missing_vars <- setdiff(vars, available)
if (length(missing_vars)) {
  message("Dropping covariates not present in data: ",
          paste(missing_vars, collapse = ", "))
}
if (length(available) == 0L) stop("No usable covariates in the input data")
eq_crime <- reformulate(available, response = "crime_count")
print(eq_crime)

use_adaptive_bw <- TRUE
kernel_function <- "bisquare"
lambda_val      <- 0.5
ksi_val         <- 0
p_val           <- 2
theta_val       <- 0
is_longlat      <- FALSE

n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
message("Using cores = ", n_cores)

subset_size          <- 25000
total_rows           <- nrow(crimedata_spt)
random_indices       <- sample(seq_len(total_rows), size = min(subset_size, total_rows), replace = FALSE)
crimedata_spt_subset <- crimedata_spt[random_indices, ]
time_tags_obs_sub    <- crimedata_spt_subset$time_var
point_coordinates    <- coordinates(crimedata_spt_subset)

spatial_distance_matrix <- GWmodel::gw.dist(
  dp.locat = point_coordinates,
  rp.locat = point_coordinates,
  p        = p_val,
  theta    = theta_val,
  longlat  = is_longlat
)

message("Starting bandwidth selection at: ", Sys.time())
optimal_bw_info <- tryCatch({
  GWmodel::bw.gtwr(
    formula   = eq_crime,
    data      = crimedata_spt_subset,
    obs.tv    = time_tags_obs_sub,
    kernel    = kernel_function,
    adaptive  = use_adaptive_bw,
    p         = p_val,
    theta     = theta_val,
    longlat   = is_longlat,
    lamda     = lambda_val,
    ksi       = ksi_val,
    approach  = "AIC",
    st.dMat   = spatial_distance_matrix,
    verbose   = TRUE
  )
}, error = function(e) {
  message("Error during bandwidth selection: ", conditionMessage(e))
  NULL
})
message("Bandwidth selection done at: ", Sys.time())
message("Selected bandwidth: ", optimal_bw_info)

# The bandwidth-search distance matrix is dead once the bandwidth is chosen
# and is 8*min(subset_size, n)^2 bytes -- 5 GB at 25000 points. Holding it
# through the full fit is pure headroom lost.
rm(spatial_distance_matrix, crimedata_spt_subset)
invisible(gc())

dp.locat          <- coordinates(crimedata_spt)
regression.points <- crimedata_spt
obs.tv            <- crimedata_spt$time_var
reg.tv            <- obs.tv

st.dMat <- st.dist(dp.locat, dp.locat, obs.tv, reg.tv, focus = 0,
                   p = p_val, theta = theta_val, longlat = is_longlat,
                   lamda = lambda_val, t.units = "auto", ksi = ksi_val)

st.bw <- if (!is.null(optimal_bw_info)) optimal_bw_info else {
  message("Falling back to a heuristic bandwidth: 5% of n.")
  max(20L, as.integer(nrow(crimedata_spt) * 0.05))
}
message("Fitting GTWR at bandwidth ", st.bw, " on n = ", nrow(crimedata_spt),
        " with ", n_cores, " cores")

modelgtwr <- gtwr(
  eq_crime,
  data     = crimedata_spt,
  obs.tv   = obs.tv,
  st.bw    = st.bw,
  kernel   = kernel_function,
  adaptive = use_adaptive_bw,
  p        = p_val,
  theta    = theta_val,
  longlat  = is_longlat,
  lamda    = lambda_val,
  t.units  = "auto",
  ksi      = ksi_val,
  st.dMat  = st.dMat,
  cores    = n_cores,
  verbose  = TRUE
)

message("GTWR fit done at: ", Sys.time())
print(modelgtwr)

results_dir <- file.path(getwd(), "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

model_path <- file.path(results_dir, sprintf("gtwr_fit_%s.rds", stamp))
saveRDS(modelgtwr, model_path)
message("Saved model to: ", model_path)

diag_path <- file.path(results_dir, sprintf("gtwr_diagnostics_%s.txt", stamp))
sink(diag_path)
cat("GTWR diagnostics — ", stamp, "\n", sep = "")
cat("bandwidth = ",       st.bw,   "\n")
cat("adaptive  = ",       use_adaptive_bw, "\n")
cat("kernel    = ",       kernel_function, "\n")
cat("lamda     = ",       lambda_val, "\n")
cat("cores     = ",       n_cores, "\n")
cat("n         = ",       nrow(crimedata_spt), "\n\n")
cat("GTW.diagnostic:\n")
print(modelgtwr$GTW.diagnostic)
cat("\nGTW.arguments:\n")
print(modelgtwr$GTW.arguments)
sink()
message("Saved diagnostics to: ", diag_path)

X <- model.matrix(eq_crime, data = crimedata_spt@data)

colnames(modelgtwr$SDF@data)[colnames(modelgtwr$SDF@data) == "Intercept"] <- "(Intercept)"
common_columns <- intersect(colnames(X), colnames(modelgtwr$SDF@data))
coefs <- as.matrix(modelgtwr$SDF@data[, common_columns])
predicted_values <- rowSums(X[, common_columns, drop = FALSE] * coefs)
residuals <- crimedata_spt$crime_count - predicted_values
RSS <- sum(residuals^2)
message("RSS = ", RSS, "   log(RSS) = ", log(RSS))
