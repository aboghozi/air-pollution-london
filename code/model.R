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

gpGrid <- function(x, y, len = NULL, search = "grid", sigma_hold = TRUE, var_hold = TRUE) {
  ## This produces low, middle and high values for sigma 
  ## (i.e. a vector with 3 elements). 
  srange <- kernlab::sigest(nox~., data=training, frac=0.5, scaled=TRUE, na.action=na.omit)
  sigma_vals <- seq(srange[1], srange[3], by=(srange[3]-srange[1])/len)
  ## hold sigma to be mean if indicated
  if (sigma_hold == TRUE){
    #sigmas <- mean(as.vector(sigma_vals[-2]))
    sigmas <-  2.371746 ## Jan (RMSE=0.7746, R^2=0.4310823)
    #sigmas <- 2.854153 ## winter_agg (RMSE=0.6409, R^2=0.5837)
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
    #kernel = "rbfdot",
    kernel = k_rbf_test,
    kpar = list(sigma = param$sigma),
    var = param$var,
    variance_model=TRUE,
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
fake_dta <- 0
x <- c('longitude', 'latitude', 'year')
y <- c('nox')
## monthly data with time
dta <- dta_monthly_agg[which(dta_monthly_agg$month == '01'),]
## winter data with time
#dta <- dta_winter_agg
## winter data collapsed over time
#dta <- dta_winter_col
dta <- dta[,c(x,y)]
## rescale year
if ('year' %in% colnames(dta)){
  dta$year <- as.numeric(dta$year) - 2000
}

set.seed(998)
inTraining <- createDataPartition(dta[,y], p=.75, list=FALSE)
training <- dta[inTraining,]
## standardize separately to not corrupt testing and training
training$nox <- standardize(training$nox)
testing  <- dta[-inTraining,]
testing$nox <- standardize(testing$nox)

if (fake_dta == 1){
  ## fake data
  x1 <- seq(-0.4929032,0.2054607,(0.2054607+0.4929032)/400)
  x1 <- x1[sample(1:length(x1))]
  x2 <- seq(51.35866,51.66864,(51.66864-51.35864)/length(x1))
  x2 <- x2[sample(1:length(x2))]
  x <- seq(-20,20,.1)
  y <- sin(x)/x + rnorm(401,sd=0.03)
  y <- y[sample(1:length(y))]
  dta <- data.frame(longitude=x1, latitude=x2, nox=y)
  inTraining <- createDataPartition(dta[,"nox"], p=.75, list=FALSE)
  training <- dta[inTraining,]
  ## standardize separately to not corrupt testing and training
  training$nox <- standardize(training$nox)
  testing  <- dta[-inTraining,]
  testing$nox <- standardize(testing$nox)
}

#if (plotting == 1){
  ## read in london boundaries
  ldnb <- readOGR(dsn = "data/statistical-gis-boundaries-london/ESRI/",
                  layer="London_Borough_Excluding_MHW") %>% spTransform(CRS("+proj=longlat +datum=WGS84"))
  ## plot both periods
  ## winter
  ## transform to spatial data object
  coords <- cbind(longitude=training$longitude, latitude=training$latitude)
  sub_training_sp <- SpatialPointsDataFrame(coords, data=training, proj4string = CRS("+proj=longlat +datum=WGS84"))
  ## assign color based on nox value
  val <- training$nox
  valcol <- (val + abs(min(training$nox)))/max(training$nox + abs(min(training$nox)))
  #pdf("figures/training_data.pdf")
  plot(ldnb)
  points(sub_training_sp, pch=16, col=rgb(1,0,0,valcol))
  title(main = list("Training Data", cex=0.8))
  title(sub = list(paste("Number of Sites: ", nrow(training), sep=""), cex=0.8))
  #dev.off()
#}


fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  ## repeated ten times
  repeats = 10)

## customized kernel with RBF(space)*RBF(time)
k_rbf_space_time <- function(x,y,sigma=0.05){
  rbf <- function(x, y, sigma) {
    exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
  }
  x_a <- x[c(1,2)]
  y_a <- y[c(1,2)]
  x_b <- x[3]
  y_b <- y[3]
  rbf_space <- rbf(x_a, y_a, sigma)
  rbf_time <- rbf(x_b, y_b, sigma)
  return(rbf_space*rbf_time)
}
class(k_rbf_space_time) <- "kernel"





set.seed(825)
system.time(
#gpTrainRBF <- train(nox ~ ., data = training, 
gpTrainRBF <- train(nox ~ ., data = training,
                   method = gp, 
                   #preProc = c("center", "scale"),
                   tuneLength = 10,
                   trControl = fitControl)
)
gpTrainRBF

## plot RMSE
pdf("figures/gpTrainRBF_jan_noise&sigma.pdf", height=4, width=6)
trellis.par.set(caretTheme())
plot(gpTrainRBF, metric = "RMSE")
dev.off()

## testing data
gp_rbf_pred <- predict(gpTrainRBF, newdata=testing)
postResample(pred = gp_rbf_pred, obs = testing$nox)

## Jan best sigma: 2.371746, var=0.01


## create customized kernel (https://stackoverflow.com/questions/17056080/how-to-customize-a-kernel-function-in-ksvm-of-kernlab-package?rq=1)
## customized kernel with RBF(space)*RBF(time)
k_rbf_space_time <- function(x,y,sigma=0.05){
  rbf <- function(x, y, sigma) {
    exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
  }
  x_a <- x[c(1,2)]
  y_a <- y[c(1,2)]
  x_b <- x[3]
  y_b <- y[3]
  rbf_space <- rbf(x_a, y_a, sigma)
  rbf_time <- rbf(x_b, y_b, sigma)
  return(rbf_space*rbf_time)
}
class(k_rbf_space_time) <- "kernel"


x <- c(-0.3457146, 51.4243, 0)
y <- c(-0.012381, 51.4725, 10)
k(x, y)


k_rbf_test <- function(x,y,sigma=1){
#  rbf <- function(x, y, sigma) {
#    exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
#  }
  x_a <- x[c(1,2)]
  y_a <- y[c(1,2)]
  exp(sigma*(2*crossprod(x_a,y_a) - crossprod(x_a) - crossprod(y_a)))
  #rbf_space <- rbf(x, y, sigma)
  #return(rbf_space)
}
class(k_rbf_test) <- "kernel"


#exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
#exp(-0.5*norm((as.matrix(x)-as.matrix(y)),"f")^2) 



set.seed(777)
training <- training[,c('longitude', 'latitude', 'nox')]
gp3 <- gausspr(nox ~ ., data=training, kernel="rbfdot", kpar=list(sigma=1), var=0.001)
gp2 <- gausspr(nox ~ ., data=training, kernel=k_rbf_test, kpar=list(sigma=1), var=0.001)
gp2
gp3



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

## option #3) all data
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
