if (!require("easypackages")) install.packages("easypackages")

#install and load libraries
easypackages::packages ("GWmodel","sf", "sp", "spdep", "spatialreg", "GWmodel", "tmap", "mapview", "car", "RColorBrewer", 
                        "cowplot", "leafsync", "leaflet.extras2", "mapview", "tidyverse", "dplyr")


#read data
crimedata <- st_read("finished_data.gpkg")
crimedata <- crimedata %>% filter(year < 2020)
crimedata_lsoa <- st_read("crime_data_refined.gpkg", layer = "LSAO") # read the shapefile with the LSOA polygons


dim(crimedata)
names(crimedata)

#Columns that we will use for GWR
# "lsoa_code"
# "log_crime_count"
#"mean_temperature" 
#"max_mean_temperature"
# "percent_youth" 
#"population_density"
# "ethnic_diversity_index_11"
# "unemployment_unadj"
# "median_age"
# "per_unemployment_lsoa_2011"
# "per_social_grade_DE"
# "parks_count_within"
# "bars_pubs_count_within" 
# "sports_count_within"
# "pct_park_area_within"
# "geom"


geom_lsoa <- crimedata %>%
  group_by(lsoa_code) %>%
  slice(1) %>%
  select(lsoa_code, geom)

## Remove geometry and summarise the data
crimedata_summary <- crimedata %>%
  st_drop_geometry() %>%  # Drop geometry to avoid conflict
  group_by(lsoa_code) %>%
  summarise(
    log_crime_count = median(log_crime_count, na.rm = TRUE),
    mean_temperature = median(mean_temperature, na.rm = TRUE),
    max_mean_temperature = median(max_mean_temperature, na.rm = TRUE),
    percent_youth = median(percent_youth, na.rm = TRUE),
    population_density = median(population_density, na.rm = TRUE),
    ethnic_diversity_index_11 = median(ethnic_diversity_index_11, na.rm = TRUE),
    unemployment_unadj = median(unemployment_unadj, na.rm = TRUE),
    median_age = median(median_age, na.rm = TRUE),
    per_unemployment_lsoa_2011 = first(per_unemployment_lsoa_2011),
    per_social_grade_DE = first(per_social_grade_DE),
    parks_count_within = first(parks_count_within),
    bars_pubs_count_within = first(bars_pubs_count_within),
    sports_count_within = first(sports_count_within),
    pct_park_area_within = first(pct_park_area_within)
  )

# Join back the geometry
crimedata_summary_sf <- left_join(crimedata_summary, geom_lsoa, by = "lsoa_code") %>%
  st_as_sf()

#convert sf object into a sp spatial object. As the GWmodel package works on sp objects.
crimedata_sp <- as_Spatial(crimedata_summary_sf)


# Define model formula
eq_crime <- log_crime_count ~ max_mean_temperature + percent_youth + population_density + ethnic_diversity_index_11 + unemployment_unadj + median_age + per_unemployment_lsoa_2011 + per_social_grade_DE + parks_count_within + bars_pubs_count_within + sports_count_within + pct_park_area_within


#fit the adaptive kernel 
abw <- bw.gwr(eq_crime,
             approach = "AIC", #specified by CV for cross-validation approach or by AIC corrected (AICc), we used AIC 
             adaptive = TRUE,
             kernel="gaussian", #this can be different function e.g., bisquare,exponential, depend on the prior understanding or choice of the modeler
             data=crimedata_sp) #give the sp data created earlier

#abw <-3641 #Gaussian 3641

#fitting the model with gwr.basic function
a.gwr <- gwr.basic(eq_crime, #the equation
             adaptive = TRUE,
             kernel="gaussian", #indicate the Kernel again
             bw = abw, #give the optimal bandwidth we found in the last stage
             data=crimedata_sp) 

#print the model result
a.gwr

##CHECK MORAN I FOR RESIDUALS
residual <- a.gwr$SDF$residual #this will give the names of the variables in the model

#Spatial weight matrix
# Create a spatial weights matrix using Queen's contiguity
crimedata_nbq <- poly2nb(crimedata_lsoa, queen=TRUE) #Queen’s Contiguity neighborhood
summary(crimedata_nbq)
crimedata_nbq_w <- nb2listw(crimedata_nbq, style="W") #Queen’s neighborhood wights
summary(crimedata_nbq_w)

# GLobal Moran's I
#Check spatial autocorrelation of the dependent variable
mc_global <- moran.mc(residual, crimedata_nbq_w, 2999, alternative="greater") #here 2999 is the simulation number, if taking too long you can also use 999
#plot the  Moran’s I
plot(mc_global)
mc_global

#extract statistics and p-value
moran_I <- mc_global$statistic
p_value <- mc_global$p.value
cat("Moran's I:", moran_I, "\n")
cat("p-value:", p_value, "\n")

######### PLOTS ##########

library(tmap)
tmap_mode("view")

#get the SDF out of modle object as an sf object using st_as_sf function
agwr_sf = st_as_sf(a.gwr$SDF)


# determine which are significant for the greenspace variable
Crimevalaw = agwr_sf %>% dplyr::select(all_of("max_mean_temperature_TV")) %>% st_drop_geometry()
agwrCrimesig = Crimevalaw < -1.96 | Crimevalaw > 1.96 #this will create a new layer polygons with these value ranges
#plot the map with the significant values
#Go for tmap viewer mode for dynamically explore the map
tmap_mode("plot")
tm_shape(agwr_sf) +
  tm_shape(agwr_sf[agwrCrimesig,]) +
  tm_fill("max_mean_temperature", palette = "-RdYlBu", auto.palette.mapping = TRUE) + 
  tm_borders(col = "red", lwd = 0.01) +
  tm_layout(
    title = "Significant Areas: Max Mean Temperature Coefficient",
    legend.position = c("left", "bottom")
  ) +
  tm_legend(title = "Temperature Range")

