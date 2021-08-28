#2. Reprojecting BA shapefile
#This file reprojects the DHS shapefile to the LCC coordinate system to match the NARR data. It also saves some meta data from the shapefile (BA_SF.tsv)


rm(list = ls())
#A: Adjust BA SF to match projection of NARR Data

library(rgdal)
library(raster)
library(ncdf4)  #package update

#read in BA shapefile
sf <- readOGR("BA_SF/Control_Areas.shp")
#remove HI and AK
sf <- sf[!sf$STATE %in% c("HI", "AK"), ]
plot(sf) #original projection

#read in temp data
temp.brick <- brick("NARR_data/air.2m.1995.nc")

plot(temp.brick$X1995.05.05.21.31.26)

#convert shapefile projection to match temp data
sf <- spTransform(sf, crs(temp.brick)) #match CRS to temp data
sf@bbox <- matrix(c(4e6, 1.5e6, 9e6, 5.5e6), ncol = 2) #extend extent to make sure we get all border cells
writeOGR(sf, "BA_SF", "Control_Areas_LCC", driver = "ESRI Shapefile") #this gives warnings on the data but not used so I think it's ok
sf <- readOGR("BA_SF/Control_Areas_LCC.shp")
fwrite(sf@data, "BA_SF/sf_meta.tsv", sep = "\t")
plot(sf)
plot(temp.brick$X1995.05.05.21.31.26, add = T)
plot(sf, add = T)
