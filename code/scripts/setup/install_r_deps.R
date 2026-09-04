# Installs the R side of the pipeline, including the patched GWmodel fork.
#
#   Rscript code/scripts/setup/install_r_deps.R
#
# Requires R >= 4.5 and a matching Rtools (Rtools45 covers R 4.5 and 4.6);
# GWmodel contains C++ and will not install without a toolchain.
#
# Note there is deliberately no maptools here. Upstream GWmodel 2.2-0 depended
# on it, but it was archived from CRAN in 2023 and no longer compiles on
# R >= 4.5 (it still uses the removed Calloc/Free macros). This fork drops the
# dependency and writes shapefiles with sf instead.

pkgs <- c(
  "Rcpp", "RcppArmadillo", "sp", "sf", "spdep", "spatialreg", "robustbase",
  "spacetime", "FNN", "tidyverse", "tmap", "mapview", "car", "RColorBrewer",
  "cowplot", "leafsync", "leaflet.extras2"
)

# Versions this pipeline was last verified against, for reference when a
# future CRAN release changes behaviour.
verified <- c(
  sf = "1.1.2", sp = "2.2.3", spdep = "1.4.2", spatialreg = "1.4.3",
  tidyverse = "2.0.0", tmap = "4.4.1", mapview = "2.11.4", car = "3.1.5",
  RColorBrewer = "1.1.3", cowplot = "1.2.0", leafsync = "0.1.0",
  leaflet.extras2 = "1.3.2", Rcpp = "1.1.2", RcppArmadillo = "15.4.2.1",
  robustbase = "0.99.7", spacetime = "1.3.4", FNN = "1.1.4.1"
)

lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) lib <- .libPaths()[1]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN dependencies already present.")
}

# GWmodel is not installed from CRAN: the pipeline needs the patched fork,
# which must sit next to this repo (or be pointed at by GWMODEL_REPO).
repo <- Sys.getenv("GWMODEL_REPO", unset = "")
if (!nzchar(repo)) repo <- file.path(dirname(getwd()), "GWmodel")
if (!file.exists(file.path(repo, "DESCRIPTION")))
  stop("Could not find the GWmodel checkout at '", repo, "'. Clone it beside ",
       "this repo, or set GWMODEL_REPO to its path.")
message("Installing GWmodel from: ", repo)
install.packages(repo, lib = lib, repos = NULL, type = "source")

still <- setdiff(c(pkgs, "GWmodel"), rownames(installed.packages()))
if (length(still)) stop("Failed to install: ", paste(still, collapse = ", "))

now <- vapply(names(verified), function(p) as.character(packageVersion(p)), "")
drift <- names(verified)[now != verified]
if (length(drift))
  message("Note - versions differ from the verified set: ",
          paste(sprintf("%s %s (verified %s)", drift, now[drift], verified[drift]),
                collapse = "; "))
message("R environment ready.")
