##*********************************
##
## Model
##
##*********************************

## clear workspace
rm(list=ls())

library(kernlab)
library(gptk)
library(caret)
library(robustHD)
library(doParallel)
rCluster <- makePSOCKcluster(4)
registerDoParallel(rCluster)
on.exit(stopCluster(rCluster))

## setwd
setwd("~/Documents/MIT/Spring2019/6.862/project/air-pollution-london/")

##===============================================================================
## read in munged data

## collapsed over time
dta_winter_col <- readRDS("data/kcl_london_model_data_winter_collapsed.rds")
dta_nowinter_col <- readRDS("data/kcl_london_model_data_nowinter_collapsed.rds")

## aggregated over time
dta_winter_agg <- readRDS("data/kcl_london_model_data_winter_agg_time.rds")
dta_nowinter_agg <- readRDS("data/kcl_london_model_data_nowinter_agg_time.rds")

dta_monthly_agg <- readRDS("data/kcl_london_model_data_monthly.rds")

##===============================================================================
## split into train and test

## SHOULD I INCLUDE OTHER PARAMETERS?
## standardize nox (mean 0, variance 1)
## is it making mistakes on very high and very low or generically?
## make bandwidth small
## try radial basis kernel (RBF)
## nugget? noise parameter point-wise error -- allows for noise at the same point (quality of the point-wise)
## generate estimation for sigma values
## could add in site_type (could start as a regression) (make a new kernel, product of the kernel)


##===============================================================================
## define model
gp <- list(type = "Regression",
              library = "kernlab")

prm <- data.frame(parameter = c("var", "sigma"),
                  class = rep("numeric", 2),
                  label = c("Noise", "Sigma"))
gp$parameters <- prm

gpGrid <- function(x, y, len = NULL, search = "grid", sigma_hold = TRUE, var_hold = FALSE) {
  ## This produces low, middle and high values for sigma 
  ## (i.e. a vector with 3 elements). 
  srange <- kernlab::sigest(nox~., data=training, frac=0.5, scaled=TRUE, na.action=na.omit)
  sigma_vals <- seq(srange[1], srange[3], by=(srange[3]-srange[1])/len)
  ## hold sigma to be mean if indicated
  if (sigma_hold == TRUE){
    #sigmas <- mean(as.vector(sigma_vals[-2]))
    sigmas <-  2.371746
  } else{
    sigmas <- sigma_vals
  }
  if (var_hold == TRUE){
    var_value <- 0.001
  } else{
    var_max <- 0.5
    var_min <- 0.001
    var_value <- seq(from=var_min, to=var_max, by=(var_max-var_min)/len)
  }
  
  ## To use grid search:
  if(search == "grid") {
    out <- expand.grid(sigma = sigmas,
                       var = var_value)
  }
  out
}
gp$grid <- gpGrid

gpFit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) { 
  kernlab::gausspr(
    x = as.matrix(x), y = y,
    kernel = "rbfdot",
    kpar = list(sigma = param$sigma),
    var = param$var,
    prob.model = classProbs,
    ...
  )
}
gp$fit <- gpFit

gpPred <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  kernlab::predict(modelFit, newdata)
gp$predict <- gpPred

gpProb <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  kernlab::predict(modelFit, newdata, type = "probabilities")
gp$prob <- gpProb

##===============================================================================

x <- c('longitude', 'latitude')
y <- c('nox')
## monthly data with time
dta <- dta_monthly_agg[which(dta_monthly_agg$month == '01'),]
## winter data with time
#dta <- dta_winter_agg
dta <- dta[,c(x,y)]

set.seed(998)
inTraining <- createDataPartition(dta[,y], p=.75, list=FALSE)
training <- dta[inTraining,]
## standardize separately to not corrupt testing and training
training$nox <- standardize(training$nox)
testing  <- dta[-inTraining,]
testing$nox <- standardize(testing$nox)

fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  ## repeated ten times
  repeats = 10)

set.seed(825)
system.time(
gpTrainRBF <- train(nox ~ ., data = training, 
                   method = gp, 
                   #preProc = c("center", "scale"),
                   tuneLength = 10,
                   trControl = fitControl)
)
gpTrainRBF

## plot RMSE
pdf("figures/gpTrainRBF_noise.pdf", height=4, width=6)
trellis.par.set(caretTheme())
plot(gpTrainRBF, metric = "RMSE")
dev.off()

## testing data
gp_rbf_pred <- predict(gpTrainRBF, newdata=testing)
postResample(pred = gp_rbf_pred, obs = testing$nox)

## Jan best sigma: 2.371746


## create customized kernel
k <- function(x,y,z,sigma=0.05){
  return(exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y))))
}
class(k) <- "kernel"



set.seed(777)
gp2 <- gausspr(nox ~ ., data=training, kernel="rbfdot", kpar=list(sigma=1), var=0.001)
gp2



##===============================================================================

## can add doMC package for faster processing time
## http://philipmgoddard.com/R/custom_models_with_caret

set.seed(667)
system.time({gpFitRBF_test <- train(nox ~ ., data = training, 
                         method = "gaussprRadial", 
                         trControl = fitControl,
                         preProcess="scale",
                         tuneGrid=gpGrid)
  gpFitRBF_test
})

set.seed(667)
gpFitRBF_test <- train(nox ~ ., data = training, 
                                    method = "gaussprRadial", 
                                    trControl = fitControl,
                                    preProcess="scale",
                                    #tuneGrid=gpGrid,
                                    var = 0.001,
                                    sigma=3)
gpFitRBF_test

gpFitRBF <- train(training[,x], training[,y],
                  method='gaussprRadial', preProcess="scale",
                  trControl=fitControl, tuneGrid=sigmaValues)
gpFitRBF

## next: keep sigma and re-run with noise parameter -- done

## lat, long, add in time
## 
## option #1) stay with lat, long -- could stick with the same kernel -- done

## option #2) add time component -- keep kernel parameters from above (first cut)
## try break up by month, winter/non-winter
## A) multiply RBF in space * RBF in time
## break up by month, winter/non-winter
## B) linear/ polynomial (to capture linear trend in time) * squared exponential in space


## periodic kernel in time (every month) * linear/polynomial (linear trend) * squared exponential in space


#predict and variance
x = c(-4, -3, -2, -1,  0, 0.5, 1, 2)
y = c(-2,  0,  -0.5,1,  2, 1, 0, -1)
plot(x,y)
foo2 <- gausspr(x, y, variance.model = TRUE)
xtest <- seq(-4,2,0.2)
lines(xtest, predict(foo2, xtest))
lines(xtest,
      predict(foo2, xtest)+2*predict(foo2,xtest, type="sdeviation"),
      col="red")
lines(xtest,
      predict(foo2, xtest)-2*predict(foo2,xtest, type="sdeviation"),
      col="red")




##-----------------------------------------------------------
## experimental design (new data)
