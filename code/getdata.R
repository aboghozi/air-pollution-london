##*********************************
##
## Get and Investigate KCL Data
##
##*********************************

## clear workspace
rm(list=ls())

## setwd
setwd("~/Documents/MIT/Spring2019/6.862/project/code/air-pollution-monitoring/")

#require(devtools)
#install_github('davidcarslaw/openair')

library(openair)
library(maptools)
library(rgeos)
library(tidyverse)
library(rgdal)
library(raster)
library(ggplot2)

##===============================================================================
## read in data
## metadata about monitoring programs
#aurn_meta <- importMeta(source = "aurn", all = TRUE)
#saqn_meta <- importMeta(source = "saqn", all = TRUE)
#kcl_meta <- importMeta(source = "kcl", all=TRUE)
kcl_meta <- readRDS("data/kcl_meta_data.rds")

## read in london boundaries
ldnb <- readOGR(dsn = "data/statistical-gis-boundaries-london/ESRI/",
                layer="London_Borough_Excluding_MHW") %>% spTransform(CRS("+proj=longlat +datum=WGS84"))

##===============================================================================
## subset to london

## take out NAs in lat/long for sites
kcl_meta <- kcl_meta[which(!is.na(kcl_meta$latitude)),]
## transform to spatial data object
coords <- cbind(longitude=kcl_meta$longitude, latitude=kcl_meta$latitude)
kcl_meta_SP <- SpatialPointsDataFrame(coords, data=kcl_meta, proj4string = CRS("+proj=longlat +datum=WGS84"))
## subset to sensors within the London boundaries
kcl_meta_SP <- kcl_meta_SP[ldnb,]
## convert back to data frame
kcl_meta2 <- as.data.frame(kcl_meta_SP)
## save sites
sites <- kcl_meta2$code
## access data
#kcl_data <- importKCL(site=sites, year=1990:2019, meta=TRUE, met=TRUE)
kcl_data <- readRDS("data/kcl_london_data.rds")

## create year and month-year columns
#kcl_data$year <- format(as.Date(kcl_data$date, format="%Y/%m/%d"),"%Y")
#kcl_data$month_year <- format(as.Date(kcl_data$date, format="%Y/%m/%d"),"%m-%Y")

##===============================================================================
## plot sites on top of boroughs
pdf("figures/kcl_sites_london.pdf")
plot(ldnb)
points(kcl_meta_SP, pch=16, col=rgb(0,0,0,0.8))
title(main = list("KCL Sites in Greater London Area", cex=0.8))
title(sub = list(paste("Number of Sites: ", length(unique(kcl_meta2$code)), sep=""), cex=0.8))
dev.off()

## subsets based on pollutant
pollutants <- c('nox', 'no2', 'so2', 'co', 'pm10', 'pm25', 'o3')
for (pollutant in pollutants){
  sub_pollutant <- kcl_data[!is.na(kcl_data[names(kcl_data) == pollutant]),]
  sub_sites <- kcl_meta2[which(kcl_meta2$code %in% unique(sub_pollutant$code)),]
  ## transform to spatial data object
  coords <- cbind(longitude=sub_sites$longitude, latitude=sub_sites$latitude)
  sub_sites_SP <- SpatialPointsDataFrame(coords, data=sub_sites, proj4string = CRS("+proj=longlat +datum=WGS84"))
  pdf(paste("figures/site_map_", pollutant, ".pdf", sep=""))
  plot(ldnb)
  points(sub_sites_SP, pch=19, col=rgb(0,0,0,0.8))
  title(main = list(paste("Map of Greater London with ", pollutant, " Sites", sep=""), cex=0.8))
  title(sub = list(paste("Number of Sites: ", length(unique(sub_sites$code)), sep=""), cex=0.8))
  dev.off()
  
  dta_table <- as.data.frame(table(sub_pollutant$code, sub_pollutant$year))
  names(dta_table) <- c('code', 'year', 'freq')
  write.csv(dta_table, paste("data/site_year_", pollutant, ".csv", sep=""), row.names=F)
  
  saveRDS(sub_pollutant, paste("data/kcl_data_", pollutant, ".rds", sep=""))
}

dta <- data.frame(year=character(), num_sites=character(), pollutant=character(), stringsAsFactors=FALSE)
for (pollutant in pollutants){
  dta.sub <- read.csv(paste("data/site_year_", pollutant, ".csv", sep=""), as.is=T, header=T, row.names=NULL)
  print(paste(pollutant, sum(dta.sub$freq), sep=" "))
  dta.sub <- as.data.frame(table(dta.sub$year[dta.sub$freq > 0]))
  names(dta.sub) <- c('year', 'num_sites')
  dta.sub$pollutant <- pollutant
  dta <- rbind(dta, dta.sub)
}
dta$year <- as.numeric(as.character(dta$year))

## plot all lines
pdf("figures/number_of_sites_per_pollutant_per_year.pdf", height=4, width=6)
ggplot(data=dta, aes(x=year, y=num_sites)) + geom_line(aes(colour=pollutant))+
  labs(x="year",y="number of sites")
dev.off()

calendarPlot(sub_pollutant, pollutant="o3")

##===============================================================================
## write out data
#saveRDS(kcl_data, "data/kcl_london_data.rds")
#saveRDS(kcl_meta, "data/kcl_meta_data.rds")
