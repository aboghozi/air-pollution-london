## clear workspace
rm(list=ls())

## setwd
setwd("~/Documents/MIT/Spring2019/6.862/project/air-pollution-london/")

library(openair)
library(TSA)
library(forecast)
library(ggplot2)

data_m2 <- c(0.293738121763502, 0.45728805929725613, 0.9980547846555808,
            0.14154699613486757, 0.1458077967763497, 0.3122744859095631,
            0.1588594575454301, 0.19108466579058989, NA, 0.2656939430133758,
            0.3799254185984513, NA)

data_m1 <- c(0.8456863378328031, 0.9498865601879187, 0.8622644623424907,
             0.938434528957358, 0.8891895921499919, 0.9099270364887586,
             0.9036238395028933, 0.9413148346599836, 1.0006019083124036,
             0.9437203826510108, 0.881269262930133, 0.865659369444384)

data_m3 <- c(0.28575088342736316, 1.6687370138120152, 0.20234739028743054,
             0.12071800470906409, 0.11270673779635344, 0.2643669979671157,
             0.10461483930076629, 0.1414502763958607, 0.11038204107636336,
             0.21791853695775756, 0.3481292191320313, 0.25517511767401957)

months=c('Jan', 'Feb', 'March', 'April',
         'May', 'June', 'July', 'Aug', 'Sept',
         'Oct', 'Nov', 'Dec')

df <- data.frame(month=rep(months,3))
df$x <- rep(seq(1,12,1),3)
df$mse <- c(data_m1, data_m2, data_m3)
df$model <- c(rep("Model 1", 12), rep("Model 2", 12), rep("Model 3", 12))

## plot all lines
pdf("figures/mse_months_all_models.pdf", height=5, width=10)
ggplot(data=df, aes(x=x, y=mse, group=model, colour=model)) +
  geom_line() + labs(x="month",y="mse") +
  scale_x_continuous(breaks = seq(1,12,1), labels=c(months))
dev.off()


