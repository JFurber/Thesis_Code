library(sp)
library(rgdal)
library(sf)

#LATLONG TO UTM
## shortcuts
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
utm <- "+init=epsg:32630"

## Create coordinates variable
# Combine the columns longitude and latitude into 1 column
coords <- cbind(Easting = as.numeric(as.character(Ireland$Long)),
                Northing = as.numeric(as.character(Ireland$Lat)))

## Create the SpatialPointsDataFrame - this identifies what type of coordinates they are
# CRS(ukgrid), i.e. Cornwall Data
# CRS(latlong), for coordinates in lat/long coordinates
xy_SP <- SpatialPointsDataFrame(coords, #column made above
                                 data = Ireland , #choose data file
                                 proj4string = CRS(latlong))

## Convert the column to UTM coordinates
xy_SP <- spTransform(xy_SP, CRS(utm))

# create new column for the Easting/Northing Coordinates
xy_SP@data$Easting <- coordinates(xy_SP)[, 1]
xy_SP@data$Northing <- coordinates(xy_SP)[, 2]

#Find minimum of columns (Not neccessary)
easting_min <- min(xy_SP$Easting)
northing_min <- min(xy_SP$Northing)
#Create a new column with translated coordinates to 0,0
xy_SP@data$ZeroEasting <- (xy_SP$Easting - easting_min)
xy_SP@data$ZeroNorthing <- xy_SP$Northing - northing_min

#Then save the new data set xy_SP as a CSV filem with the now UTM coordinates at the end of the file
write.csv(xy_SP,"Ireland_UTM.csv", row.names = TRUE)


#-------------------------
#UTM TO LATLONG
## shortcuts
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
utm <- "+init=epsg:32630"

## Create coordinates variable
# Combine the columns longitude and latitude into 1 column
coords <- cbind(Easting = as.numeric(as.character(Traj$V1)),
                Northing = as.numeric(as.character(Traj$V2)))

## Create the SpatialPointsDataFrame - this identifies what type of coordinates they are
# CRS(ukgrid), i.e. Cornwall Data
# CRS(latlong), for coordinates in lat/long coordinates
xy_SP <- SpatialPointsDataFrame(coords, #column made above
                                data = Traj , #choose data file
                                proj4string = CRS(utm))

## Convert the column to UTM coordinates
xy_SP <- spTransform(xy_SP, CRS(latlong))

# create new column for the Easting/Northing Coordinates
xy_SP@data$Longitude <- coordinates(xy_SP)[, 1]
xy_SP@data$Latitude <- coordinates(xy_SP)[, 2]

#Find minimum of columns (Not necessary)
easting_min <- min(xy_SP$Easting)
northing_min <- min(xy_SP$Northing)
#Create a new column with translated coordinates to 0,0
xy_SP@data$ZeroEasting <- (xy_SP$Easting - easting_min)
xy_SP@data$ZeroNorthing <- xy_SP$Northing - northing_min

#Then save the new data set xy_SP as a CSV filem with the now UTM coordinates at the end of the file
write.csv(xy_SP,"Traj_LatLong.csv", row.names = TRUE)

