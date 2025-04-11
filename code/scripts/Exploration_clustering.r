if (!require("easypackages")) install.packages("easypackages")

#install and load libraries
easypackages::packages ("pathwork", "GWmodel","sf", "sp", "spdep", "spatialreg", "GWmodel", "tmap", "mapview", "car", "RColorBrewer", 
                        "cowplot", "leafsync", "leaflet.extras2", "mapview", "tidyverse")

#read data
crimedata <- st_read("finished_data.gpkg")
crimedata_lsoa <- st_read("crime_data_refined.gpkg", layer = "LSAO") # read the shapefile with the LSOA polygons
crimedata <- crimedata %>% filter(year < 2020)


## PLOTS
## Plot log crime_count, max_mean_temperature for a specific month
selected_year <- 2013
selected_month <- 7
# Filter the dataset
crime_data_filtered <- crimedata %>%
  filter(year == selected_year, month == selected_month)

## Crime count map
crime_plot <- ggplot(crime_data_filtered) +
  geom_sf(aes(fill = log_crime_count), color = "black", size = 0.05) +
  scale_fill_viridis_c(name = "Log Crime Count") +
  labs(title = paste("Log Crime Count per LSOA - Month", selected_month),
       subtitle = "Log-transformed ASB count") +
  theme_minimal()

# Mean temperature map
temp_plot <- ggplot(crime_data_filtered) +
  geom_sf(aes(fill = max_mean_temperature), color = "black", size = 0.05) +
  scale_fill_viridis_c(name = "Mean Temperature (°C)") +
  labs(title = paste("Mean Temperature per LSOA - Month", selected_month),
       subtitle = "Max monthly average temperature") +
  theme_minimal()

# Combine the plots side by side
crime_plot + temp_plot



###### GLOBAL MORAN I #########

# Create a spatial weights matrix using Queen's contiguity
crimedata_nbq <- poly2nb(crimedata_lsoa, queen=TRUE) #Queen’s Contiguity neighborhood
summary(crimedata_nbq)
crimedata_nbq_w <- nb2listw(crimedata_nbq, style="W") #Queen’s neighborhood wights
summary(crimedata_nbq_w)


# GLobal Moran's I
#Check spatial autocorrelation of the dependent variable
mc_global <- moran.mc(crime_data_filtered$log_crime_count, crimedata_nbq_w, 2999, alternative="greater") #here 2999 is the simulation number, if taking too long you can also use 999
mc_global
#plot the  Moran’s I
plot(mc_global)


#extract statistics and p-value
moran_I <- mc_global$statistic
p_value <- mc_global$p.value
cat("Moran's I:", moran_I, "\n")
cat("p-value:", p_value, "\n")


# Compute Moran I for each month
#Create unique time_var
crimedata$year_month <- as.factor(paste(crimedata$year, crimedata$month, sep = "-")) # create a year_month variable
unique_time_values <- unique(crimedata$year_month)
crimedata$time_var <- match(crimedata$year_month, unique_time_values)

length(crimedata$time_var)

# Create dataframe for moran_results
moran_results <- data.frame(time_var = integer(), moran_I = numeric(), p_value = numeric())

# Get Global Moran I for each time_var
for (t in unique(crimedata$time_var)) {
  crimedata_filtered <- crimedata %>% filter(time_var == t)

  # Calculate Moran’s I
  mc_global <- moran.mc(crimedata_filtered$max_mean_temperature, 
                          crimedata_nbq_w, 
                          999, alternative = "greater")

  # Store in dataframe
  moran_results <- rbind(moran_results, 
                           data.frame(time_var = t, 
                                      moran_I = mc_global$statistic, 
                                      p_value = mc_global$p.value))
} 


############# LOCAL MORAN I #############

# Local Moran's I for specific month and year
gm_mortality_LISA <- localmoran(crime_data_filtered$log_crime_count, crimedata_nbq_w)  #using Queen's contiguity 
summary(gm_mortality_LISA)

# to visualize this statistic the relevant information needs to be extracted
# extract local Moran's I values and attache them to our sf object 
crime_data_filtered$gm_mortality_LISA <- gm_mortality_LISA[,1] 
# extract p-values
crime_data_filtered$gm_mortality_LISA_p <- gm_mortality_LISA[,5] 

#Here we can map the local Moran's I with t-map, and show which areas have significant clusters
map_LISA <- tm_shape(crime_data_filtered) + 
  tm_polygons(col= "gm_mortality_LISA", title= "Local Moran’s I", midpoint=0,
              palette = "RdYlBu", breaks= c(-10, -5, 0, 5, 10, 20)) 
map_LISA_p <- tm_shape(crime_data_filtered) + 
  tm_polygons(col= "gm_mortality_LISA_p", title= "p-values",
              breaks= c(0, 0.01, 0.05, 1), palette = "Reds") 

tmap_arrange(map_LISA, map_LISA_p)


#Plot LISA
library(tmap)

# Set tmap to plotting mode
tmap_mode("plot")

# Get the bounding box manually and check the coordinates
bbox_coords <- st_bbox(crime_data_filtered)

# Now generate the map with the specified projection
map_LISA <- tm_shape(crime_data_filtered, bbox = st_bbox(crime_data_filtered)) + 
  tm_polygons(col = "gm_mortality_LISA", 
              title = "Local Moran’s I", 
              palette = "RdYlBu", 
              style = "fixed", 
              breaks = c(-10, -5, 0, 5, 10, 20), 
              midpoint = 0,
              border.col = "gray70", 
              border.alpha = 0.3) +
  tm_layout(title = "Local Moran’s I for Crime Counts", 
            legend.outside = TRUE)

# Improved significance map
map_LISA_p <- tm_shape(crime_data_filtered, bbox = st_bbox(crime_data_filtered)) + 
  tm_polygons(col = "gm_mortality_LISA_p", 
              title = "p-values (Significance)", 
              palette = "Reds", 
              style = "fixed", 
              breaks = c(0, 0.01, 0.05, 1), 
              border.col = "gray70", 
              border.alpha = 0.3) +
  tm_layout(title = "Significance of Local Moran’s I", 
            legend.outside = TRUE)


