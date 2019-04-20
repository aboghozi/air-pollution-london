##*********************************
##
## Visualize Time-Series KCL Data
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

##===============================================================================
## read in munged data
dta_pm10_daily <- readRDS("data/kcl_pm10_daily_agg.rds")
dta_nox_daily <- readRDS("data/kcl_nox_daily_agg.rds")

dta_pm10_monthly <- readRDS("data/kcl_pm10_monthly_agg.rds")
dta_nox_monthly <- readRDS("data/kcl_nox_monthly_agg.rds")

##===============================================================================

## subset data to more relevant years
sub <- dta_nox_monthly[which(dta_nox_monthly$year >= 2000),]
sub$full_date <- paste("01-", sub$month_year, sep="")
sub$full_date <- format(as.Date(sub$full_date, format="%d-%m-%Y"),"%Y-%m-%d")
sub$full_date <- as.Date(sub$full_date, format="%Y-%m-%d")
sub$month <- format(as.Date(sub$full_date, format="%d-%m-%Y"),"%m")


## FULL DATASET
## collapse data over months for each site
sub_agg <- aggregate(sub$nox, by=list(sub$month, sub$code), FUN=mean, na.rm=T)
names(sub_agg) <- c('month', 'site', 'nox')
## number of years
sub$counter <- 1
sub_agg_count <- aggregate(sub$counter, by=list(sub$month, sub$code), FUN=sum, na.rm=T)
names(sub_agg_count) <- c('month', 'site', 'num_years')
## merge
sub_agg <- merge(sub_agg, sub_agg_count, by=c('month', 'site'), all.x=T)

sub_agg$month <- as.numeric(as.character(sub_agg$month))
sub_agg$nox <- as.numeric(as.character(sub_agg$nox))

## take out observations with low number of years
sub_agg <- sub_agg[which(sub_agg$num_years > 4),]

length(unique(sub_agg$site))

## plot all lines
pdf("figures/nox_sites_aggregated.pdf", height=5, width=8)
ggplot(data=sub_agg, aes(x=month, y=nox, colour=factor(site))) +
  geom_line() + labs(x="month",y="nox") + theme(legend.position="none")
dev.off()


## find the site codes with the most data
dta_site_freq <- as.data.frame(table(sub$code))
names(dta_site_freq) <- c('site', 'Freq')
dta_site_freq <- dta_site_freq[order(dta_site_freq$Freq, decreasing=TRUE),]
## subset to the top 75% sites to sample from
dta_site_freq_sub <- dta_site_freq[dta_site_freq$Freq >= quantile(dta_site_freq$Freq, 0.75),]
## randomly sample sites
set.seed(9011)
index <- sample(nrow(dta_site_freq_sub), 5, replace=FALSE)
sub_sites <- as.character(dta_site_freq_sub$site[index])



## ALL SITES -- Time Series

## take out LB4
#sub <- sub[which(sub$code != 'LB4'),]

## collapse data over sites
sub_agg_ts <- aggregate(sub$nox, by=list(sub$full_date), FUN=mean, na.rm=T)
names(sub_agg_ts) <- c('full_date', 'nox')
## number of sites
length(unique(sub_agg_count$site))
sub_agg_count_s <- aggregate(sub$counter, by=list(sub$full_date), FUN=sum, na.rm=T)
names(sub_agg_count_s) <- c('full_date','num_sites')
## merge
sub_agg_ts <- merge(sub_agg_ts, sub_agg_count_s, by=c('full_date'), all.x=T)

dates <- seq(from=min(sub$full_date), to=max(sub$full_date), by='month')
dta <- data.frame(full_date=dates)
dta <- merge(dta, sub_agg_ts[,c('full_date', 'nox')], by='full_date', all.x=T)

dta_ts <- ts(dta$nox,
             start=c(as.numeric(format(min(dta_site$full_date), "%Y")),
                     as.numeric(format(min(dta_site$full_date), "%m"))),
             end=c(as.numeric(format(max(dta_site$full_date), "%Y")),
                   as.numeric(format(max(dta_site$full_date), "%m"))), frequency=12)

## plot time series
pdf(paste("figures/time_series_all.pdf", sep=""), height=4, width=7)
plot(dta_ts, ylab="nox", xlab="year-month")
dev.off()

decomposedRes <- decompose(dta_ts, type="additive") 

pdf(paste("figures/decomposed_time_series_all.pdf", sep=""), height=6, width=9)
plot(decomposedRes)
dev.off()

dta$month <- format(as.Date(dta$full_date, format="%d-%m-%Y"),"%m")
lm.1 = lm(nox ~ factor(month), data = dta)
summary(lm.1)

## take out LB4



## PER SITE
site <- "KC1"
#for (site in sub_sites){
  #num_obvs <- dta_site_freq_sub$Freq[index]
  dta_site <- sub[which(sub$code == site),]
  dates <- seq(from=min(dta_site$full_date), to=max(dta_site$full_date), by='month')
  dta <- data.frame(full_date=dates)
  dta <- merge(dta, dta_site[,c('full_date', 'nox')], by='full_date', all.x=T)
  ## convert to time-series object
  dta_ts <- ts(dta$nox,
               start=c(as.numeric(format(min(dta_site$full_date), "%Y")),
                       as.numeric(format(min(dta_site$full_date), "%m"))),
               end=c(as.numeric(format(max(dta_site$full_date), "%Y")),
                       as.numeric(format(max(dta_site$full_date), "%m"))), frequency=12)
  ## plot time series
  pdf(paste("figures/time_series_", site, ".pdf", sep=""), height=4, width=7)
  plot(dta_ts, ylab="nox", xlab="year-month")
  dev.off()
  
  #dta_years <- ts(dta$nox, start=c(2010,1), end=c(as.numeric(format(max(dta_site$full_date), "%Y")),
  #                                         as.numeric(format(max(dta_site$full_date), "%m"))), frequency=12)
  #plot(dta_years)
  
  decomposedRes <- decompose(dta_ts, type="additive") 
  
  pdf(paste("figures/decomposed_time_series_", site, ".pdf", sep=""), height=4, width=7)
  plot(decomposedRes)
  dev.off()
  
  stlRes <- stl(dta_ts, s.window = "periodic")
  
  acfRes <- acf(dta_ts)
  pacfRes <- pacf(dta_ts)
#}


dta$month <- format(as.Date(dta$full_date, format="%d-%m-%Y"),"%m")

seasonality.func = function(df){
  lm.1 = lm(nox ~ factor(month), data = dta)
  summary(lm.1)
  p.vals = summary(lm.1)$coefficients[,4]
  p.vals.lt.01 = as.numeric(sum(p.vals<.01)>0)
  return(p.vals.lt.01)
}

seasonality.func(df)

## all of the data
## collapse over site and month
sub_agg_month <- aggregate(sub$nox, by=list(sub$month), FUN=mean, na.rm=T)
names(sub_agg_month) <- c('month', 'avg_nox')
sub_agg_month$month_name <- c('Jan', 'Feb', 'Mar', 'April', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')
pdf("figures/avg_nox_concentration_monthly.pdf", height=4, width=8)
plot(sub_agg_month$month, sub_agg_month$avg_nox, pch=16, col=rgb(0,0,0,0.8),
     main="Monthly Average NOx Concentration" ,xlab="", ylab="Average NOx Concentration", xaxt='n',
     ylim=c(0,max(sub_agg_month$avg_nox)+20))
axis(side=1, at=sub_agg_month$month, label=sub_agg_month$month_name)
dev.off()

#lm.2 = lm(nox ~ factor(month), data=sub)
#summary(lm.2)


# compute the Fourier Transform
#p = periodogram(dta_site$nox)
p = periodogram(sub$nox)
dd = data.frame(freq=p$freq, spec=p$spec)
order = dd[order(-dd$spec),]
top2 = head(order, 2)
top2
# convert frequency to time periods
time = 1/top2$f
time
