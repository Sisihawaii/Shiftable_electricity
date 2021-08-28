# visualization for overall US
# 8/18/2020
# variable in overlapped area as mean of all overlapped area
# ggplot with sf for map; viridis for color scale

rm(list = ls())

library(tidyverse)
library(rgdal)
library(data.table)
library(ggplot2)
library(rgeos)
library(maptools)
library(sp)
library(viridis)
library(sf)
library(raster)

# load data
peakstats_all <- readRDS("peakstatss.RDS")
peakstats_e <- peakstats_all[[1]]
peakstats_w <- peakstats_all[[2]]
peakstats_ercot <- peakstats_all[[3]]
summary(peakstats_e)

# avg peak reduction per BA
avgsf_all <- list()
for (i in 1:3){
  avgsf_all[[i]] <- peakstats_all[[i]] %>% 
    group_by(sf_id) %>% 
    summarise(peakred_avg = mean(peakred), baseinc_avg = mean(baseinc),
              SDred_avg = mean(sdred, na.rm = TRUE), CVred_avg = mean(cvred, na.rm = TRUE),
              flatshare_p = round(sum(flat)/length(flat) * 100))
  avgsf_all[[i]] <- as.data.table(avgsf_all[[i]])
}
avgsf_all[[2]]$sf_id <- 15 + avgsf_all[[2]]$sf_id     # stack id numbers 1-31
avgsf_all[[3]]$sf_id <- 31
data_us <- rbindlist(avgsf_all)

###############################################################
############### plot the mean values for overlaps #############
###############################################################

# load shapefile
east15_sf <- readOGR("BA_SF/aggregated/east/east_15.shp")
west15_sf <- readOGR("BA_SF/aggregated/west/west_15.shp")
ercot_sf <- readOGR("BA_SF/aggregated/ercot/ercot.shp")
ercot_sf <- ercot_sf[,1]
colnames(ercot_sf@data) <- "ID"

# change all ID from level to numeric and stack id number 1-31
east15_sf@data$ID <- as.numeric(levels(east15_sf@data$ID))[east15_sf@data$ID]
west15_sf@data$ID <- as.numeric(levels(west15_sf@data$ID))[west15_sf@data$ID] + 15
ercot_sf@data$ID <- as.numeric(levels(ercot_sf@data$ID))[ercot_sf@data$ID] + 16

# combine all the shapefile
us_sf <- rbind(east15_sf, west15_sf, ercot_sf)
# plot(us_sf)
# assign data to shapefile - all area
us_sf@data[2:6] <- data.frame(data_us[, 2:6])

##=========================================================##
# set up function to find overlap for two polygons
overlap2 <- function(sf, data, n){
  # sf: shapefiles
  # data: variables to fill in overlapped area
  # n: number of polygons in the shapefile
  spdf <- NULL
  for (i in 1:n){
    for (j in 1:n){
      if (j > i){
        itc <- gIntersection(sf[i,], sf[j,], drop_lower_td = TRUE)
        if (is_empty(itc)) {
          next
        } else {
          itc_rows <- rbind(data[i,], data[j,])
          df <- data.frame(t(colMeans(itc_rows[,2:6])))
          spdf0 <- SpatialPolygonsDataFrame(itc, df)
          if (is_empty(spdf)){
            spdf <- spdf0
          } else {
            spdf <- rbind(spdf, spdf0)
          }
        }
      }
    }
  }  
  return(spdf)
}
##==========================================================##

# use the function calculate overlapped area variable mean
overlap_all <- overlap2(us_sf, data_us, 31)
overlap_all@data$ID <- 31 + seq.int(nrow(overlap_all@data))     # add row number as ID, starting from 32
overlap_all <- cbind(overlap_all[,6], overlap_all[,1:5])          # rearranging 
saveRDS(overlap_all, "overlap2_new2.RDS")

######################################################

overlap_all <- readRDS("overlap2_new2.RDS")

# combine overlapped and original data together
comb_all <- rbind(us_sf, overlap_all)
comb_all_sf <- st_as_sf(comb_all)
saveRDS(comb_all, "comb_new2.RDS")

### plot all the area with ggplot
plot_all <- ggplot(data = comb_all_sf) +
  geom_sf(aes(fill = SDred_avg), colour = NA) +
  theme_void() +
  scale_fill_viridis_c(trans = "sqrt", limits = c(0.65, 1), 
                       breaks = c(0.75, 0.85, 0.95)) +
  labs(title = "SD Reduction for Continental US \
  (alpha = 1)", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.text  = element_text(size = 7)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 6))


