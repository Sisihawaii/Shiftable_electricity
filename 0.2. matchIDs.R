#2: Matching IDs
#This file creates a crosswalk between the IDs in the DHS shapefile and the balancing authorities in the EIA930 data.
#The process first matches by name which is successful for ~50 BAs when the DHS shapefile and EIA930 have the same name. Then goes through and manually adjusts the DHS names to match EIA930 names for remaining BAs. 
#two unmatched shapes: New Brunswick (Canada) and gridforce south (Unknown?)
#the crosswalk is saved as eia_sf_match.tsv

rm(list = ls())
gc()

library(data.table)

eia_meta <- fread("EIA930/EIA_meta.csv")
setnames(eia_meta, c("ba", "name", "timezone", "region"))
sf_meta <- fread("BA_SF/sf_meta.tsv")[, list(OBJECTID, NAME)]
setnames(sf_meta, c("id", "name"))

eia_meta[, name := gsub(" ", "", tolower(name))]
eia_meta <- eia_meta[timezone != ""] #remove mexico and canada
sf_meta[, name := gsub(" ", "", tolower(name))]

setkey(eia_meta, name)
setkey(sf_meta, name)

eia_meta[, .N]
sf_meta[, .N]

matched <- eia_meta[sf_meta, nomatch = 0] #match 50 bas on first run

eia_meta <- eia_meta[!(ba %in% matched[, ba])]
sf_meta <- sf_meta[!(id %in% matched[, id])]

sf_meta[name == "cityofhomestead", name := "homestead,cityof"]
sf_meta[name == "cityoftallahassee", name := "tallahassee,cityof"]
sf_meta[name == "floridapower&lightcompany", name := "floridapower&lightco."]
#sf_meta[name == "gridforcesouth", name := ""] not matched
sf_meta[name == "isonewenglandinc.", name := "newenglandiso"]
sf_meta[name == "louisvillegasandelectriccompanyandkentuckyutilities", name := "lg&eandkuservicescompanyasagentforlouisvillegasandelectriccompanyandkentuckyutilitiescompany"]
sf_meta[name == "midcontinentindependenttransmissionsystemoperator,inc..", name := "midcontinentindependentsystemoperator,inc."]
sf_meta[name == "northwesternenergy(nwmt)", name := "northwesterncorporation"]
sf_meta[name == "pacificorp-east", name := "pacificorpeast"]
sf_meta[name == "pacificorp-west", name := "pacificorpwest"]
#sf_meta[name == "progressenergycarolinas-eastandwest", name := "pacificorpwest"] not matched
sf_meta[name == "progressenergyflorida", name := "dukeenergyflorida,inc."]
sf_meta[name == "pugetsoundenergy", name := "pugetsoundenergy,inc."]
sf_meta[name == "saltriverproject", name := "saltriverprojectagriculturalimprovementandpowerdistrict"]
sf_meta[name == "tucsonelectricpowercompany", name := "tucsonelectricpower"]
sf_meta[name == "westernareapoweradministrationugpwest", name := "westernareapoweradministration-uppergreatplainswest"]
setkey(sf_meta, name)

#!
#NOTE
#!
#important kludge: combine duke progress east and west into a single BA in EIA 930 data
eia_meta[name == "dukeenergyprogresseast" | name == "dukeenergyprogresswest", name := "progressenergycarolinas-eastandwest"]
setkey(eia_meta, name)

matched2 <- eia_meta[sf_meta, nomatch = 0]


eia_meta <- eia_meta[!(ba %in% matched2[, ba])]
sf_meta <- sf_meta[!(id %in% matched2[, id])]
#two unmatched shapes: New Brunswick (Canada) and gridforce south (Unknown?)

matched <- rbind(matched, matched2)
eia_meta <- fread("EIA930/EIA_meta.csv")
setnames(eia_meta, c("ba", "eia_name", "timezone", "region"))
eia_meta <- eia_meta[, list(ba, eia_name)]
setkey(eia_meta, ba)
setkey(matched, ba)
matched <- matched[eia_meta, nomatch = 0]
matched[, name := NULL]
matched[, IC := "East"]
matched[region %in% c("California", "Northwest", "Southwest"), IC := "West"]
matched[region == "Electric Reliability Council of Texas, Inc. ", IC := "ERCOT"]
fwrite(matched, file = "eia_sf_match.tsv", sep = "\t")
