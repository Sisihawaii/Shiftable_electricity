#1: Download
#this file downloads the data needed and includes links to data that needs to be downoaded manually
#notes: NARR temp data is ~6.5 GB

#A: Download Air Temperature
dir.create("NARR_data/") #this will throw warning if it already exists, continue through
noaa_url <- "ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/"
for (yr in 1995:2019){
  print(yr)
  f_url <- paste0(noaa_url, "air.2m.", yr, ".nc")
  f_dest <- paste0("NARR_data/air.2m.", yr, ".nc")
  download.file(f_url, f_dest, mode = "wb")
}

#B: Download EIA 930
dir.create("EIA930/")
eia_url <- "http://api.eia.gov/bulk/EBA.zip"
f_dest <- paste0("EIA930/EBA.zip")
download.file(eia_url, f_dest)
unzip(f_dest, exdir = "EIA930/")

#C: Download shapefile
dir.create("BA_SF/")
#download via web @ https://hifld-geoplatform.opendata.arcgis.com/datasets/control-areas and unzip

#D: Download Pop Grid
dir.create("usgrid_data_2010/")
#download via web @ http://sedac.ciesin.columbia.edu/data/set/usgrid-summary-file1-2010 and unzip

#E: Download EIA Meta
#download via web @ https://www.eia.gov/realtime_grid/#acronyms and adjust to fit into CSV
#save as EIA930/EIA_meta.tsv
