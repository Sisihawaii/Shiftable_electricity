#2: Merge BAs
#This file consolidates the BAs in the DHS shapefile into aggregated BAs. 
#Aggregation process:
#goal - consolidate the many BAs in the DHS shapefile into large, contiguous areas that presumably have substantial transmission connectivity and can plausibly be modeled as a single unit. Run this process for East and West, leave ERCOT as one.
#A: choose the number of aggregated BAs. This is arbitrary at the moment and could in theory be improved with a more rigourous ML process but currently the aggregation and downstream data work is slow.
#B: For each BA i calculate the area that is overlapping with each BA j
#C: Select the i-j pair that has the largest area overlapping. Break ties by choosing the smallest i. 
#D: combine BA i and BA j into a new BA. Note in the data the IDs are set to i.j (eg merge 2 and 3 the new ID is 2.3)
#E: if the resulting number of BAs is higher than your threshold return to step B.


rm(list = ls())
gc()

library(rgdal)
library(raster)
library(data.table)
library(ggplot2)
library(rgeos)
library(maptools)
library(parallel)
library(RColorBrewer)
library(gridExtra)

sf <- readOGR("BA_SF/Control_Areas_LCC.shp")
meta <- fread("eia_sf_match.tsv")

#create SF for eastern IC
sf.e <- sf[sf$OBJECTID %in% meta[IC == "East", id], ]
sf.e <- gSimplify(sf.e, tol = 0.00001)
sf.e <- gBuffer(sf.e, byid=TRUE, width=0)
sum(gIsValid(sf.e, byid=TRUE)==FALSE)
sf.e <- SpatialPolygonsDataFrame(sf.e, sf[sf$OBJECTID %in% meta[IC == "East", id], ]@data)
#create SF for western IC
sf.w <- sf[sf$OBJECTID %in% meta[IC == "West", id], ]
sf.w <- gSimplify(sf.w, tol = 0.00001)
sf.w <- gBuffer(sf.w, byid=TRUE, width=0)
sum(gIsValid(sf.w, byid=TRUE)==FALSE)
sf.w <- SpatialPolygonsDataFrame(sf.w, sf[sf$OBJECTID %in% meta[IC == "West", id], ]@data)
#create SF for ERCOT
sf.t <- sf[sf$OBJECTID %in% meta[IC == "ERCOT", id], ]
sf.t <- gSimplify(sf.t, tol = 0.00001)
sf.t <- gBuffer(sf.t, byid=TRUE, width=0)
sum(gIsValid(sf.t, byid=TRUE)==FALSE)
sf.t <- SpatialPolygonsDataFrame(sf.t, sf[sf$OBJECTID %in% meta[IC == "ERCOT", id], ]@data)
calc.intarea.ij <- function(sf.i, sf.j){ 
  #calculate area intersected by shape i and j
  #return NA if no intersection
  #return 0 if only borders touch
  #return area if overlap
  if (gIntersects(sf.i, sf.j)){
    sf.ij <- gIntersection(sf.i, sf.j, drop_lower_td = T)
    if (is.null(sf.ij)){
      return(data.table(rn = rownames(sf.j@data), area.ij = 0))  
    } else {
      return(data.table(rn = rownames(sf.j@data), area.ij = area(sf.ij)))  
    }
  } else{
    return(data.table(rn = rownames(sf.j@data), area.ij = NA))
  }
}

merge.smallestba.v2 <- function(sf.base, max.cl = 32){ #adjusted to try to deal with overlap
  area.dt <- as.data.table(area(sf.base))
  area.dt[, rn := rownames(sf.base@data)]
  setnames(area.dt, c("area", "rn"))
  smallest.ba <- area.dt[area == min(area)][1, rn]
  sf.i <- sf.base[rownames(sf.base@data) == smallest.ba, ]
  sf.js <- sf.base[rownames(sf.base@data) != smallest.ba, ]
  sf.js <- split(sf.js, rownames(sf.js@data))
  ncl <- min(max.cl, length(sf.js))
  cl <- makeCluster(ncl)
  clusterEvalQ(cl, library(rgeos))
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(raster))
  overlap.dt <- parLapply(cl, sf.js, calc.intarea.ij, sf.i = sf.i)
  stopCluster(cl)
  overlap.dt <- rbindlist(overlap.dt)
  overlap.dt <- overlap.dt[area.dt, on = "rn"]
  join.ba <- overlap.dt[area.ij == max(area.ij, na.rm = T)][area == min(area, na.rm = T)][1, rn] #select largest overlapping area with smallest ba
  merge.id <- rownames(sf.base@data)
  merge.id[merge.id %in% c(smallest.ba, join.ba)] <- paste(smallest.ba, join.ba, sep = ".")
  sf.base.m <- unionSpatialPolygons(sf.base, merge.id)
  pid <- sapply(slot(sf.base.m, "polygons"), function(x) slot(x, "ID"))
  p.df <- data.frame( ID=1:length(sf.base.m), row.names = pid)
  return(SpatialPolygonsDataFrame(sf.base.m, p.df))
}

# plot(sf.w, col = rainbow(length(sf.w), alpha = 0.5))
# plot(sf.e, col = rainbow(length(sf.e), alpha = 0.5))
# plot(sf.t, col = rainbow(length(sf.t), alpha = 0.5))

#collapse western IC
t0 <- proc.time()
sf.w.15 <- copy(sf.w)
while(length(sf.w.15) > 15){
  print(length(sf.w.15))
  sf.w.15 <- merge.smallestba.v2(sf.w.15)
}
t1 <- proc.time() #~4.5 mins @ 32
sf.w.10 <- copy(sf.w.15)
while(length(sf.w.10) > 10){
  print(length(sf.w.10))
  sf.w.10 <- merge.smallestba.v2(sf.w.10)
}
t2 <- proc.time() #~0.8 mins @ 32
sf.w.5 <- copy(sf.w.10)
while(length(sf.w.5) > 5){
  print(length(sf.w.5))
  sf.w.5 <- merge.smallestba.v2(sf.w.5)
}
t3 <- proc.time() #~0.7 mins @ 32

#collapse Eastern IC
t0 <- proc.time()
sf.e.15 <- copy(sf.e)
while(length(sf.e.15) > 15){
  print(length(sf.e.15))
  sf.e.15 <- merge.smallestba.v2(sf.e.15)
}
t1 <- proc.time() #~13 mins @ 32
sf.e.10 <- copy(sf.e.15)
while(length(sf.e.10) > 10){
  print(length(sf.e.10))
  sf.e.10 <- merge.smallestba.v2(sf.e.10)
}
t2 <- proc.time() #~ 5mins @ 32
sf.e.5 <- copy(sf.e.10)
while(length(sf.e.5) > 5){
  print(length(sf.e.5))
  sf.e.5 <- merge.smallestba.v2(sf.e.5)
}
t3 <- proc.time() #~5 mins @32

#plot SFs
plot.bamap <- function(sf){
  sf.df <- fortify(sf)
  ggplot(sf.df) + 
    aes(long,lat,group=group,fill=id) + 
    geom_polygon(alpha = 0.5) +
    geom_path(color="black") +
    coord_equal() +
    scale_fill_discrete(rainbow(length(sf)), guide = F) +
    theme_void()
}

p.w.all <- plot.bamap(sf.w)
p.w.15 <- plot.bamap(sf.w.15)
p.w.10 <- plot.bamap(sf.w.10)
p.w.5 <- plot.bamap(sf.w.5)
pdf("figs/data/ba_west.pdf", width = 11, height = 8.5)
grid.arrange(p.w.all, p.w.15, p.w.10, p.w.5)
dev.off()

p.e.all <- plot.bamap(sf.e)
p.e.15 <- plot.bamap(sf.e.15)
p.e.10 <- plot.bamap(sf.e.10)
p.e.5 <- plot.bamap(sf.e.5)
pdf("figs/data/ba_east.pdf", width = 11, height = 8.5)
grid.arrange(p.e.all, p.e.15, p.e.10, p.e.5)
dev.off()

  
#save SFs
writeOGR(sf.w.15, "BA_SF/aggregated/west", "west_15", driver = "ESRI Shapefile")
writeOGR(sf.w.10, "BA_SF/aggregated/west", "west_10", driver = "ESRI Shapefile")
writeOGR(sf.w.5, "BA_SF/aggregated/west", "west_5", driver = "ESRI Shapefile")
writeOGR(sf.e.15, "BA_SF/aggregated/east", "east_15", driver = "ESRI Shapefile")
writeOGR(sf.e.10, "BA_SF/aggregated/east", "east_10", driver = "ESRI Shapefile")
writeOGR(sf.e.5, "BA_SF/aggregated/east", "east_5", driver = "ESRI Shapefile")
writeOGR(sf.t, "BA_SF/aggregated/ercot", "ercot", driver = "ESRI Shapefile")

#save objectIDS
east.ids <- data.table(rowid = as.integer(rownames(sf.e@data)), sf_id = sf.e@data$OBJECTID)
east15.ids <- data.table(rowid = rownames(sf.e.15@data), sf_id = sf.e.15@data$ID)
east10.ids <- data.table(rowid = rownames(sf.e.10@data), sf_id = sf.e.10@data$ID)
east5.ids <- data.table(rowid = rownames(sf.e.5@data), sf_id = sf.e.5@data$ID)
saveRDS(east.ids, "BA_SF/ids/east_all.RDS")
saveRDS(east15.ids, "BA_SF/ids/east_15.RDS")
saveRDS(east10.ids, "BA_SF/ids/east_10.RDS")
saveRDS(east5.ids, "BA_SF/ids/east_5.RDS")
west.ids <- data.table(rowid = as.integer(rownames(sf.w@data)), sf_id = sf.w@data$OBJECTID)
west15.ids <- data.table(rowid = rownames(sf.w.15@data), sf_id = sf.w.15@data$ID)
west10.ids <- data.table(rowid = rownames(sf.w.10@data), sf_id = sf.w.10@data$ID)
west5.ids <- data.table(rowid = rownames(sf.w.5@data), sf_id = sf.w.5@data$ID)
saveRDS(west.ids, "BA_SF/ids/west_all.RDS")
saveRDS(west15.ids, "BA_SF/ids/west_15.RDS")
saveRDS(west10.ids, "BA_SF/ids/west_10.RDS")
saveRDS(west5.ids, "BA_SF/ids/west_5.RDS")
