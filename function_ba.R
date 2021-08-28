# Function ba
# all functions for demand and weather data cleaning/aggregation

# read eba.txt
eba.to.dt <- function(eba_line){
  eba_json <- fromJSON(eba_line)
  if(is.null(eba_json$series_id)){
    return(NULL)
  } else {
    meta <- data.table(
      series_id = eba_json$series_id,
      name = eba_json$name,
      units = eba_json$units,
      f = eba_json$f,
      description = eba_json$description,
      start = eba_json$start,
      end = eba_json$end,
      last_updated = eba_json$last_updated,
      geoset_id = eba_json$geoset_id
    )
    data <- data.table(eba_json$data)
    setnames(data, c("date_time", "value"))
    data[, date_time_f := as.POSIXct(date_time, format = "%Y%m%dT%HZ", tz = "UTC")]
    data[, date_time := NULL]
    setnames(data, "date_time_f", "date_time")
    data[, series_id := meta[, series_id]]
    return(list(
      "meta" = meta,
      "data" = data
    ))
  }
}


###########################################################
####### demand cleaning: missing value and outliers #######
###########################################################

# function to impute missing or outlying values
#    -Requires data for *one* balancing authority/year and vector of data to impute
#    -Ouliers:  missing (zeros) and unusually low values indicated by miss.obs
impute.ba = function(ba.dat, miss.obs){
  ba.dat$time = 1:(nrow(ba.dat))
  ba.dat$weekdays = weekdays(ba.dat$date_time)
  fit.dat = ba.dat[-miss.obs,]
  fit =  lm( value ~ ns(time, df=5)*ns(hour, df=7)*as.factor(weekdays), data=fit.dat)
  ba.dat$r = rep(NA, nrow(ba.dat))
  ba.dat$r[-miss.obs] = resid(fit)
  fc.dat = ba.dat[miss.obs,]
  pred = predict(fit, newdata=fc.dat)
  r.pred = approx( ba.dat$time[-miss.obs], ba.dat$r[-miss.obs], xout= ba.dat$time[miss.obs])
  ba.dat$value[miss.obs] = pred + r.pred$y
  ba.dat$Flag.Interpolated <- 0
  ba.dat$Flag.Interpolated[miss.obs] = rep(1, length(miss.obs))
  ba.dat
}

# function to find outliers and call interpolation function
checkOutage = function(ba1, data){
  dat = data[ba == ba1, ]
  dat[is.na(value), value := 0]
  min.v = min(dat$value)
  q.2   = quantile(dat$value, 0.02)
  if ( min.v < q.2/2 ) {
    missing = which(dat$value < q.2/1.2)
    dat = impute.ba(dat, missing)
  }
  dat
}


################################
###### weather functions #######
################################

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
