##*********************************
##
## Prep Data for Modeling
##
##*********************************

## clear workspace
rm(list=ls())

## setwd
setwd("~/Documents/MIT/Spring2019/6.862/project/air-pollution-london/")

library(openair)
library(TSA)
library(forecast)
library(ggplot2)
library(maptools)
library(rgeos)
library(tidyverse)
library(rgdal)
library(raster)

## turn on and off plotting
plotting <- 1

##===============================================================================
## read in munged data
#dta_pm10_daily <- readRDS("data/kcl_pm10_daily_agg.rds")
dta_nox_daily <- readRDS("data/kcl_nox_daily_agg.rds")

#dta_pm10_monthly <- readRDS("data/kcl_pm10_monthly_agg.rds")
dta_nox_monthly <- readRDS("data/kcl_nox_monthly_agg.rds")

if (plotting == 1){
  ## read in london boundaries
  ldnb <- readOGR(dsn = "data/statistical-gis-boundaries-london/ESRI/",
                  layer="London_Borough_Excluding_MHW") %>% spTransform(CRS("+proj=longlat +datum=WGS84"))
}

##===============================================================================

## subset data to more relevant years
sub_month <- dta_nox_monthly[which(dta_nox_monthly$year >= 2000),]
sub <- dta_nox_daily[which(dta_nox_daily$year >= 2000),]
## isolate just month
sub$month <- format(as.Date(sub$day_month_year, format="%d-%m-%Y"),"%m")
sub_month$month <- vapply(strsplit(sub_month$month_year, "-"), '[', 1, FUN.VALUE=character(1))
## create indicator for winter (Sept, Oct, Nov, Dec, Jan, Feb, March)
sub$winter <- ifelse(sub$month %in% c('09', '10', '11', '12', '01', '02', '03'), TRUE, FALSE)

## t-test of significance
t.test(sub$nox~sub$winter)

## subset to winter and not winter
sub_winter <- sub[which(sub$winter == TRUE),]
sub_nowinter <- sub[which(sub$winter == FALSE),]

## collapse each point over time (winter/not winter)
sub_win_collapse <- aggregate(list(nox=sub_winter$nox),
                         by=list(site=sub_winter$site, code=sub_winter$code,
                                 latitude=sub_winter$latitude, longitude=sub_winter$longitude,
                                 site_type=sub_winter$site_type),
                         FUN=mean, na.rm=T)

sub_nowin_collapse <- aggregate(list(nox=sub_nowinter$nox),
                         by=list(site=sub_nowinter$site, code=sub_nowinter$code,
                                 latitude=sub_nowinter$latitude, longitude=sub_nowinter$longitude,
                                 site_type=sub_nowinter$site_type),
                         FUN=mean, na.rm=T)

sub_combined <- rbind(sub_win_collapse, sub_nowin_collapse)

## aggregate each point over time (winter/not winter)
sub_win_agg <- aggregate(list(nox=sub_winter$nox),
                         by=list(site=sub_winter$site, code=sub_winter$code,
                                 latitude=sub_winter$latitude, longitude=sub_winter$longitude,
                                 site_type=sub_winter$site_type, year=sub_winter$year),
                         FUN=mean, na.rm=T)

sub_nowin_agg <- aggregate(list(nox=sub_nowinter$nox),
                                by=list(site=sub_nowinter$site, code=sub_nowinter$code,
                                        latitude=sub_nowinter$latitude, longitude=sub_nowinter$longitude,
                                        site_type=sub_nowinter$site_type, year=sub_nowinter$year),
                                FUN=mean, na.rm=T)


## aggregate each point monthly (2000-2019)
sub_monthly <- sub_month[,c('nox', 'month', 'year', 'site', 'code',
                              'latitude', 'longitude', 'site_type')]

##===============================================================================
if (plotting == 1) {
  ## plot both periods
  ## winter
  ## transform to spatial data object
  coords <- cbind(longitude=sub_win_collapse$longitude, latitude=sub_win_collapse$latitude)
  sub_win_sp <- SpatialPointsDataFrame(coords, data=sub_win_collapse, proj4string = CRS("+proj=longlat +datum=WGS84"))
  ## assign color based on nox value
  val <- sub_win_collapse$nox
  valcol <- (val + abs(min(sub_combined$nox)))/max(sub_combined$nox + abs(min(sub_combined$nox)))
  pdf("figures/avg_kcl_sites_london_winter.pdf")
  plot(ldnb)
  points(sub_win_sp, pch=as.integer(sub_win_sp$site_type), col=rgb(1,0,0,valcol))
  title(main = list("Averaged KCL Sites in Greater London Area (Winter)", cex=0.8))
  title(sub = list(paste("Number of Sites: ", length(unique(sub_win_collapse$code)), sep=""), cex=0.8))
  legend("topleft", inset=0.53, legend=unique(sub_win_collapse$site_type),
         pch=as.integer(unique(sub_win_collapse$site_type)), col=rgb(1,0,0,0.8),
         cex=0.45, y.intersp=0.10, x.intersp=0.15, bty='n')
  dev.off()
  
  ## not winter
  ## transform to spatial data object
  coords <- cbind(longitude=sub_nowin_collapse$longitude, latitude=sub_nowin_collapse$latitude)
  sub_nowin_sp <- SpatialPointsDataFrame(coords, data=sub_nowin_collapse, proj4string = CRS("+proj=longlat +datum=WGS84"))
  ## assign color based on nox value
  val <- sub_nowin_collapse$nox
  valcol <- (val + abs(min(sub_combined$nox)))/max(sub_combined$nox + abs(min(sub_combined$nox)))
  pdf("figures/avg_kcl_sites_london_no_winter.pdf")
  plot(ldnb)
  points(sub_nowin_sp, pch=as.integer(sub_nowin_collapse$site_type), col=rgb(1,0,0,valcol))
  title(main = list("Averaged KCL Sites in Greater London Area (Not Winter)", cex=0.8))
  title(sub = list(paste("Number of Sites: ", length(unique(sub_nowin_collapse$code)), sep=""), cex=0.8))
  legend("topleft", inset=0.53, legend=unique(sub_nowin_collapse$site_type),
         pch=as.integer(unique(sub_nowin_collapse$site_type)), col=rgb(1,0,0,0.8),
         cex=0.45, y.intersp=0.10, x.intersp=0.15, bty='n')
  dev.off()
}


##===============================================================================
## write out data
write.csv(sub_win_collapse, "data/kcl_london_model_data_winter_collapsed.csv", row.names=F)
write.csv(sub_nowin_collapse, "data/kcl_london_model_data_nowinter_collapsed.csv", row.names=F)
write.csv(sub_win_agg, "data/kcl_london_model_data_winter_agg_time.csv", row.names=F)
write.csv(sub_nowin_agg, "data/kcl_london_model_data_nowinter_agg_time.csv", row.names=F)
write.csv(sub_monthly, "data/kcl_london_model_data_monthly.csv", row.names=F)
