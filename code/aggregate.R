##*********************************
##
## Condense Time-Series KCL Data
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
## read in raw data (per pollutant and subset to London)
kcl_meta <- readRDS("data/kcl_meta_data.rds")
kcl_nox <- readRDS("data/kcl_data_nox.rds")
#kcl_pm10 <- readRDS("data/kcl_data_pm10.rds")

## add day field (taking out hourly stamp)
#kcl_pm10$day_month_year <- format(as.Date(kcl_pm10$date, format="%Y/%m/%d"),"%d-%m-%Y")
kcl_nox$day_month_year <- format(as.Date(kcl_nox$date, format="%Y/%m/%d"),"%d-%m-%Y")
#kcl_pm10$id <- 1
kcl_nox$id <- 1

##===============================================================================
## Aggregate over year-month

#dta <- kcl_pm10
dta <- kcl_nox

dta_agg <- aggregate(list(nox=dta$nox, no2=dta$no2, so2=dta$so2,
                          co=dta$co, o3=dta$o3, pm2.5=dta$pm25,
                          pm10=dta$pm10, pm10_raw=dta$pm10_raw,
                          v10=dta$v10, v2.5=dta$v2.5,
                          nv10=dta$nv10, nv2.5=dta$nv2.5,
                          ws=dta$ws, wd=dta$wd, solar=dta$solar,
                          rain=dta$rain, temp=dta$temp,
                          bp=dta$bp, rhum=dta$rhum),
                     by=list(site=dta$site, code=dta$code,
                             month_year=dta$month_year, year=dta$year,
                             latitude=dta$latitude, longitude=dta$longitude,
                             site_type=dta$site.type),
                     FUN=mean, na.rm=T)
## add counts of observations per day and per month
dta_agg_num_obs <- aggregate(list(num_observations=dta$id),
                             by=list(code=dta$code, month_year=dta$month_year),
                             FUN=sum, na.rm=T)
## add counts of days
dta_agg_num_days <- as.data.frame(table(dta$code, dta$month_year, dta$day_month_year))
names(dta_agg_num_days) <- c('code', 'month_year', 'day_month_year', 'Freq')
dta_agg_num_days <- dta_agg_num_days[dta_agg_num_days$Freq > 0,]
dta_agg_num_days$id <- 1
dta_agg_num_days <- aggregate(list(num_days=dta_agg_num_days$id),
                             by=list(code=dta_agg_num_days$code, month_year=dta_agg_num_days$month_year),
                             FUN=sum, na.rm=T)
## combine
dta_agg <- merge(dta_agg, dta_agg_num_obs, by=c("code", "month-year"), all.x=TRUE)
dta_agg <- merge(dta_agg, dta_agg_num_days, by=c("code", "month_year"), all.x=TRUE)

## write out data
#saveRDS(dta_agg, "data/kcl_pm10_monthly_agg.rds")
#saveRDS(dta_agg, "data/kcl_nox_monthly_agg.rds")

##===============================================================================
## Aggregate over day-year-month

dta <- kcl_pm10
#dta <- kcl_nox

dta_agg <- aggregate(list(pm10=dta$pm10, pm10_raw=dta$pm10_raw,
                          nox=dta$nox, no2=dta$no2, so2=dta$so2,
                          co=dta$co, o3=dta$o3, pm2.5=dta$pm25,
                          v10=dta$v10, v2.5=dta$v2.5,
                          nv10=dta$nv10, nv2.5=dta$nv2.5,
                          ws=dta$ws, wd=dta$wd, solar=dta$solar,
                          rain=dta$rain, temp=dta$temp,
                          bp=dta$bp, rhum=dta$rhum),
                     by=list(site=dta$site, code=dta$code,
                             day_month_year=dta$day_month_year,
                             month_year=dta$month_year, year=dta$year,
                             latitude=dta$latitude, longitude=dta$longitude,
                             site_type=dta$site.type),
                     FUN=mean, na.rm=T)

## add counts of observations per day and per month
dta_agg_num_obs <- aggregate(list(num_observations=dta$id),
                             by=list(code=dta$code, day_month_year=dta$day_month_year),
                             FUN=sum, na.rm=T)
## combine
dta_agg <- merge(dta_agg, dta_agg_num_obs, by=c("code", "day_month_year"), all.x=TRUE)

## reorder
dta_agg <- dta_agg[order(dta_agg$site, dta_agg$year, dta_agg$month_year, dta_agg$day_month_year),]

## write out data
#saveRDS(dta_agg, "data/kcl_pm10_daily_agg.rds")
#saveRDS(dta_agg, "data/kcl_nox_daily_agg.rds")
