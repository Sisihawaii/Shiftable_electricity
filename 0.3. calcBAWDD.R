#3. Calculate population weighted DDs for each consolidated BA
#This files calculates weighted degree hour measures for each BA. Here BA refers to the consolidated BAs calculated in mergeBAs.R
#process
#A: calculate population weight for each cell in the BA. Note that these weights change for each set of shapefiles (15 vs 10 vs 5). If Aij is the share of the area of cell i that is located in BA j (this will be 1 for cells entirely in the BA, 0 for cells outside the BA, and 0-1 where the BA borders intersect part of the cell) and Pi is the population living in cell i. Assume that population is uniformly distributed in the cell (not a realistic assumption but simplifies things greatly) then the population of the cell that lives in BA is Pij = Aij x Pi. The weight given to cell i for calculating weather in BA j is Pij / sum_i(Pij). This process is somewhat slow (on the order of several minutes). Note that these weights will change for each SF set because the BA boundries change.
#B: Interpolate 8x daily NARR data to hourly using spline. Calculate CDH and HDH in each cell/hour
#C: Merge the two datasets and take weighted average of DH in each BA-hour


rm(list = ls())
gc()

library(rgdal)
library(raster)
library(ncdf4)
library(data.table)
library(fasttime)
library(ggplot2)
library(zoo)
library(parallel)

do.fn <- T

if (do.fn){
  
  extract.narr.at.cells <- function(narr.brick, cells){
    narr.values <- setDT(extract(narr.brick, cells, df = T))[, cell := as.integer(cells)][, ID := NULL] #note that ID here is jsut rowname from extract
    setkey(narr.values, cell)
    narr.values <- melt(narr.values, id.vars = c("cell"), variable.name = "date_time")
    narr.values[, date_time_f := as.POSIXct(round(fastPOSIXct(substring(date_time, 2), tz = "UTC"), "hours")), by = date_time]
    narr.values[, date_time := NULL]
    setnames(narr.values, "date_time_f", "date_time")
    return(narr.values)
  }
  
  interpolate.narr <- function(dt){
    #note that this modifies in place. Do not call dt <- interpolate.narr(dt)
    dt[, ID.cell := paste(ID, cell, sep = "."), by = list(ID, cell)]
    setkey(dt, ID.cell, date_time)
    
    #fill in NA values
    dt <- dt[CJ(unique(ID.cell), seq(min(date_time), max(date_time), by = "hour"))]
    dt[, `:=`(ID = NULL, cell = NULL)]
    dt[, c("ID", "cell") := tstrsplit(ID.cell, ".", fixed = T, type.convert = T), by = ID.cell]
    dt[, ID.cell := NULL]
    setkey(dt, ID, cell, date_time)
    
    #interpolate NAs with spline
    dt[, int.value := as.numeric(na.spline(ts(value))), by = list(ID, cell)]
    return(dt)
  }
  
  interpolate.narr.by.IDcell <- function(dt){
    #note this is designed to work only on dt that contains a single ID.cell
    setkey(dt, cell, date_time)
    
    #fill in NA values
    dt <- dt[CJ(unique(cell), seq(min(date_time), max(date_time), by = "hour"))]
    id <- dt[, unique(ID)]
    dt[, ID := id[!is.na(id)]]
    setkey(dt, ID, cell, date_time)
    
    #interpolate NAs with spline
    dt[, int.value := as.numeric(na.spline(ts(value))), by = list(ID, cell)]
    return(dt)
  }
  
  
  interpolate.narr.by.ID <- function(dt){
    #note this is designed to work only on dt that contains a single ID
    setkey(dt, cell, date_time)
    
    #fill in NA values
    dt <- dt[CJ(unique(cell), seq(min(date_time), max(date_time), by = "hour"))]
    id <- dt[, unique(ID)]
    dt[, ID := id[!is.na(id)]]
    setkey(dt, ID, cell, date_time)
    
    #interpolate NAs with spline
    dt[, value := as.numeric(na.spline(ts(value))), by = cell]
    
    #calc cdh and hdh
    dt[, value := value - 273.15]
    dt[, cdh_18 := pmax(value - 18, 0)]
    dt[, hdh_18 := pmax(18 - value, 0)]
    return(dt)
  }
  
  interpolate.narr.by.cell <- function(cell.dt){
    #fill in NA values
    cell.dt <- cell.dt[CJ(unique(cell), seq(min(date_time), max(date_time), by = "hour"))]
    setkey(cell.dt, cell, date_time)
    #interpolate missing with spline
    cell.dt[, value := as.numeric(na.spline(ts(value)))]
    #calc cdh and hdh
    cell.dt[, value := value - 273.15] #this is to convert from K (narr default) to C
    cell.dt[, cdh_18 := pmax(value - 18, 0)]
    cell.dt[, hdh_18 := pmax(18 - value, 0)]
    return(cell.dt)
  }
  
  calc.pop.weight <- function(sf, pop){
    #sf shoudl be a single shape file for one BA
    #pop should be the full population grid
    pop.dt <- setDT(extract(pop, sf, df = T, cellnumbers = T, weights = T, normalizeWeights = F))[, ID := sf$ID]
    pop.dt <- pop.dt[uspop10 > 0] #remove cells with zero population and zero weight
    pop.dt[, `:=`(ID = as.integer(ID, cell = as.integer(cell)))]
    pop.dt[, w.pop := uspop10 * weight] #this is the weighted population: number of people in cell * share of cell that is in sf. Assume that population is uniform over the cell
    pop.dt[, w := w.pop / sum(w.pop), by = ID]
    pop.dt[, `:=`(uspop10 = NULL, weight = NULL, w.pop = NULL)] #remove unused columns to save memory
    setkey(pop.dt, ID, cell)
    return(pop.dt)
  }
  
  calc.pop.dt <- function(sf.list, pop, ncl){
    #calc pop weight
    cl <- makeCluster(ncl)
    clusterEvalQ(cl, {
      library(data.table)
      library(raster)
    })
    pop.dt <- parLapplyLB(cl, sf.list, calc.pop.weight, pop = pop)
    stopCluster(cl)
    return(rbindlist(pop.dt))
  }
  
  calc.wdh.by.ID <- function(temp.dt, full.pop.dt){
    setkey(temp.dt, ID, cell, date_time)
    temp.dt[full.pop.dt, w := i.w]
    return(temp.dt[, list(cdh_18 = sum(cdh_18 * w), hdh_18 = sum(hdh_18 * w)), by = list(ID, date_time)])
  }
  
  calc.wdh.by.ID.v2 <- function(ba.pop.dt, full.temp.dt){
    #pop.dt should be the population weights for all of the cells in BA i
    #full.temp.dt should be the full temperature history for all cells in US
    ba.temp <- full.temp.dt[ba.pop.dt, nomatch = 0] ## mistake here! temp.dt doesn't exist!!!
    setkey(ba.temp, date_time)
    
    wdh <- ba.temp[, list(temp = sum(value * w), hdh_18 = sum(hdh_18 * w), cdh_18 = sum(cdh_18 * w)), by = list(sf, ID, date_time)]
    setkey(wdh, date_time)
    return(wdh)
  }
  
}

#read in BA shapefile
sf <- readOGR("BA_SF/Control_Areas_LCC.shp")
sf@bbox <- matrix(c(4e6, 1.5e6, 9e6, 5.5e6), ncol = 2) #extend extent to make sure we get all border cells
sf.w.15 <- readOGR("BA_SF/aggregated/west/west_15.shp")
sf.w.10 <- readOGR("BA_SF/aggregated/west/west_10.shp")
sf.w.5 <- readOGR("BA_SF/aggregated/west/west_5.shp")
sf.e.15 <- readOGR("BA_SF/aggregated/east/east_15.shp")
sf.e.10 <- readOGR("BA_SF/aggregated/east/east_10.shp")
sf.e.5 <- readOGR("BA_SF/aggregated/east/east_5.shp")
sf.t <- readOGR("BA_SF/aggregated/ercot/ercot.shp")
#split shapefiles
sf.w.15.l <- split(sf.w.15, sf.w.15$ID)
sf.w.10.l <- split(sf.w.10, sf.w.10$ID)
sf.w.5.l <- split(sf.w.5, sf.w.5$ID)
sf.e.15.l <- split(sf.e.15, sf.e.15$ID)
sf.e.10.l <- split(sf.e.10, sf.e.10$ID)
sf.e.5.l <- split(sf.e.5, sf.e.5$ID)
sf.list <-list(
  "w.15" = sf.w.15.l,
  "w.10" = sf.w.10.l,
  "w.5" = sf.w.5.l,
  "e.15" = sf.e.15.l,
  "e.10" = sf.e.10.l,
  "e.5" = sf.e.5.l
)

#read in temp brick
all.yrs <- 2015:2018   # select the range of years
all.temp <- lapply(all.yrs, function(yr) brick(paste0("NARR_data/air.2m.", yr, ".nc"))) #this gives a warning that I don't understand
names(all.temp) <- all.yrs

#crop temp bricks
t0 <- proc.time()
cl <- makeCluster(4)
  clusterEvalQ(cl, {
    library(raster)
  })
  all.temp<- parLapply(cl, all.temp, crop, y = sf)
stopCluster(cl)
t1 <- proc.time() #~15 sec

#read in pop data for 2010
pop <- raster("usgrid_data_2010/geotiff/uspop10.tif")
pop <- aggregate(pop, fact=32, fun = sum)
pop <- projectRaster(pop, to = all.temp[[1]]) #match CRS to temp data and match resolution

#calc pop weight
t0 <- proc.time()
big.cl <- makeCluster(6)
clusterEvalQ(big.cl, {
  library(data.table)
  library(raster)
  library(parallel)
})
clusterExport(big.cl, "calc.pop.weight")
pop.dt.list <- parLapply(big.cl, sf.list, calc.pop.dt, pop = pop, ncl = 5)
stopCluster(big.cl)
t1 <- proc.time() #~8 mins @ 6x5, probably faster at 6x6
lapply(pop.dt.list, function(dt) setkey(dt, ID, cell))
lapply(1:length(pop.dt.list), function(k) pop.dt.list[[k]][, sf := names(pop.dt.list)[[k]]])
pop.dt <- rbindlist(pop.dt.list)
pop.dt[, sf := as.factor(sf)]
pop.dt[, cell := as.integer(cell)]
setkey(pop.dt, sf, ID, cell)
saveRDS(pop.dt, "popweight.RDS")
pop.dt <- readRDS("popweight.RDS")
pop.cells <- pop.dt[, unique(cell)]

#extract temp values in all pop cells
t0 <- proc.time()
cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(raster)
  library(data.table)
  library(fasttime)
})
temp.dt <- parLapply(cl, all.temp, extract.narr.at.cells, cells = pop.cells) 
stopCluster(cl)
temp.dt <- rbindlist(temp.dt)
t1 <- proc.time() #~15 sec @ 4 cores
setkey(temp.dt, cell, date_time)
temp.dt <- split(temp.dt, by = "cell")

#interpolate missing temp in each cell and calc DH
t0 <- proc.time()
cl <- makeCluster(32)
clusterEvalQ(cl, {
  library(data.table)
  library(zoo)
})
temp.dt <- parLapply(cl, temp.dt, interpolate.narr.by.cell) # no error - depends on pc
stopCluster(cl)
temp.dt <- rbindlist(temp.dt) # need to run this one again
setkey(temp.dt, cell, date_time)
t1 <- proc.time() #~1.3 mins @ 32 cores
#saveRDS(temp.dt, "temp_by_cell.RDS") #this takes very long to save
temp.dt1 <- readRDS("temp_by_cell.RDS")

gc()

#merge in temp values to pop
pop.dt <- split(pop.dt, by = c("sf", "ID"), drop = T)
lapply(pop.dt, function(dt) setkey(dt, cell))
setkey(temp.dt, cell)

gc()
t0 <- proc.time()
cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(data.table)
})
wdh.dt <- parLapplyLB(cl, pop.dt, calc.wdh.by.ID.v2, full.temp.dt = temp.dt)
stopCluster(cl)
wdh.dt <- rbindlist(wdh.dt)
setkey(wdh.dt, sf, ID, date_time)
t1 <- proc.time()

gc()
t0 <- proc.time()
wdh.dt <- lapply(pop.dt, calc.wdh.by.ID.v2, full.temp.dt = temp.dt)
wdh.dt <- rbindlist(wdh.dt)
setkey(wdh.dt, sf, ID, date_time)
t1 <- proc.time() #5.3 mins @ lapply

#save
saveRDS(wdh.dt, "wdh_by_ba.RDS")
