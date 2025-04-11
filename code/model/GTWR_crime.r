#install and load libraries
if (!require("easypackages")) install.packages("easypackages")
easypackages::packages ("GWmodel","sf", "sp", "spdep", "spatialreg", "GWmodel", "tmap", "mapview", "car", "RColorBrewer", 
                        "cowplot", "leafsync", "leaflet.extras2", "mapview", "tidyverse")

# read and filter data
crimedata <- st_read("C:\\Users\\Ruben\\Documents\\Project_crime\\code\\model\\finished_data.gpkg")
crimedata <- crimedata %>% filter(year < 2014, year>2010)
crimedata <- crimedata %>% rename(crime_count = log_crime_count)

# create year_month
crimedata$year_month <- factor(
  paste(crimedata$year, crimedata$month, sep = "-"),
  levels = unique(paste(crimedata$year, crimedata$month, sep = "-"))
)
# load spatial dataframe
crimedata_spt <- as_Spatial(crimedata)

#create time_var (index/coordianate for each time) and reformat vars
crimedata_spt$time_var <- as.numeric(crimedata_spt$year_month)
crimedata_spt$crime_count <- as.numeric(crimedata_spt$crime_count)
crimedata_spt$mean_temperature <- as.numeric(crimedata_spt$mean_temperature)
crimedata_spt$lsoa_code <- as.factor(crimedata_spt$lsoa_code)
crimedata_spt$time_var <- as.numeric(crimedata_spt$time_var)
crimedata_spt$percent_youth <- as.numeric(crimedata_spt$percent_youth)
crimedata_spt$population_density <- as.numeric(crimedata_spt$population_density)

#effects
time_effect <- factor(crimedata_spt$time_var)
individual_effect <-factor(crimedata_spt$lsoa_code)

# equation for GTWR
vars <- c("mean_temperature",
          #"factor(lsoa_code)",
          #"factor(time_var)",
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
eq_crime <- reformulate(vars, response = "crime_count")
print(eq_crime)

#bandwidth selection           #--Default is:
use_adaptive_bw <- TRUE        # False
kernel_function <- "bisquare"  # bisquare
lambda_val <- 0.5              # 0.05, chose it to balance time and space
ksi_val <- 0                   # 0
p_val <- 2                     # euclidean distance will suffice for London
theta_val <- 0                 # 0
is_longlat <- FALSE            # epsg = 27700 (GB), is projected

#subset for bandwidht selection
subset_size <- 25000
total_rows <- nrow(crimedata_spt)
random_indices <- sample(1:total_rows, size = subset_size, replace = FALSE)
crimedata_spt_subset <- crimedata_spt[random_indices, ]
time_tags_obs_sub <- crimedata_spt_subset$time_var                          #obs.tv: a vector of time tags for each observation
total_rows <- nrow(crimedata_spt)
point_coordinates <- coordinates(crimedata_spt_subset)
spatial_distance_matrix <- GWmodel::gw.dist(
  dp.locat = point_coordinates,
  rp.locat = point_coordinates, # Or omit rp.locat
  p = p_val,                       # Euclidean distance
  theta = theta_val,                   # No rotation
  longlat = is_longlat             
)


print(paste("Starting bandwidth selection at:", Sys.time()))
optimal_bw_info <- tryCatch({
                      GWmodel::bw.gtwr(
                        formula = eq_crime,
                        data = crimedata_spt_subset,     #regression.points will be SpatialPointsDataFrame(coords=rp.locat, data=data) or regression.points <-SpatialPolygonsDataFrame(Sr=polygons, data=data,match.ID=F)
                        obs.tv = time_tags_obs_sub, #reg will be overwritten with obs.tv if no regression poitns given
                        kernel = kernel_function,
                        adaptive = use_adaptive_bw,
                        p = p_val,
                        theta = theta_val,
                        longlat = is_longlat,
                        lamda = lambda_val, 
                        ksi = ksi_val,
                        approach = "AIC", 
                        st.dMat = spatial_distance_matrix
                        # Other parameters like search range (bw.sel.range) might be useful
                      )
}, error = function(e) {
  print(paste("Error during bandwidth selection:", e$message))
  return(NULL) # Return NULL or some indicator of failure
})

# retrieve coords and time
dp.locat <- coordinates(crimedata_spt)  #get centroid coords for each polygon
regression.points <- crimedata_spt      # Regressiepunten = dataset zelf
obs.tv <- crimedata_spt$time_var        # Tijdstempels van de observaties
reg.tv <- obs.tv                        # Tijdstempels voor regressiepunten (zelfde als obs.tv)

s.dMat <- gw.dist(dp.locat = dp.locat, rp.locat = rp.locat, longlat = FALSE)
t.dMat <- as.matrix(dist(obs.tv))

st.dMat <- st.dist(dp.locat, rp.locat, obs.tv, reg.tv,focus=0, p=2, 
        theta=0, longlat=F,lamda=0.05,t.units = "auto",
        ksi=0, s.dMat,t.dMat)

# bandwidth selection 
bw <- bw.gtwr(eq_crime, crimedata_spt, crimedata_spt@data$time_var, approach = 'aic', kernel = 'gaussian', adaptive = TRUE, st.dMat = st.dMat)

bw
# 37 of 568 observations for small dataset
st.bw <- bw  # De eerder berekende bandbreedte

modelgtwr <- gtwr(eq_crime, 
                  data = crimedata_spt, 
                  regression.points = regression.points, 
                  obs.tv = obs.tv, 
                  reg.tv = reg.tv, 
                  st.bw = st.bw, 
                  kernel = "bisquare",
                  adaptive = TRUE, 
                  p = 2, 
                  theta = 0, 
                  longlat = TRUE, 
                  lamda = 0.05, 
                  t.units = "auto", 
                  ksi = 0,
                  st.dMat = st.dMat)

modelgtwr
names(modelgtwr$SDF@data)

## Model does not output AIC, RMSE nor significance of the coefficients.

X <- model.matrix(~ mean_temperature + factor(lsoa_code) + factor(time_var) + percent_youth + population_density, data = crimedata_spt@data)
dim(X)
head(X)

# Explicitly match columns names
colnames(modelgtwr$SDF@data)[colnames(modelgtwr$SDF@data) == "Intercept"] <- "(Intercept)"
common_columns <- intersect(colnames(X), colnames(modelgtwr$SDF@data))
print(common_columns)
coefs <- as.matrix(modelgtwr$SDF@data[, common_columns])
dim(t(coefs))
head(coefs)

#No predict function in GTWR, so we need to calculate the predicted values manually
# Assuming you have the model matrix X and coefficients (coefs)
# coefs should be a matrix with the same number of rows as the number of observations

predicted_values <- rowSums(X * coefs)  #  # Multiply X with the coefficients to get predicted values

head(predicted_values)
length(predicted_values)

head(crimedata_spt$crime_count)
head(predicted_values)

residuals <- crimedata_spt$crime_count - predicted_values  # Replace dependent_variable with your actual dependent variable
RSS <- sum(residuals^2)
print(log(RSS))


import geopandas as gpd
import matplotlib.pyplot as plt




