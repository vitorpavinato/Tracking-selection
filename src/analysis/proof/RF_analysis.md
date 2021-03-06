Tracking Selection with ABC-RF
================

ABC Random Forests Analysis
---------------------------

#### Install required packages

``` r
#install.packages(c("DataExplorer", 
#                   "dplyr",
#                   "abcrf",
#                   "weights",
#                   "grDevices"),
#                 dependencies = T)
```

#### Load packages and source file

``` r
library(DataExplorer)
library(dplyr)
library(abcrf)
library(weights)
library(grDevices)

source("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_fun.R")
```

#### Upload the reference table files

``` r
# Super-batch 1
load("/home/pavinato/My_repositories/Tracking-selection/results/pipeline_v4_PopGroup_results_backup/Tracking-selection-1.0/results/pooled_reftable_1.RData")
pooled_raw_reftable_1 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

``` r
# Super-batch 2
load("/home/pavinato/My_repositories/Tracking-selection/results/pipeline_v4_PopGroup_results_backup/Tracking-selection-1.1/results/pooled_reftable_2.RData")
pooled_raw_reftable_2 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

``` r
# Super-batch 3
load("/home/pavinato/My_repositories/Tracking-selection/results/pipeline_v4_PopGroup_results_backup/Tracking-selection-1.2/results/pooled_reftable_3.RData")
pooled_raw_reftable_3 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

``` r
# Super-batch 4
load("/home/pavinato/My_repositories/Tracking-selection/results/pipeline_v4_PopGroup_results_backup/Tracking-selection-1.3/results/pooled_reftable_4.RData")
pooled_raw_reftable_4 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

#### Pool together the reftables from super-batches

``` r
pooled_raw_reftable_model_1 <- rbind(pooled_raw_reftable_1,pooled_raw_reftable_2,pooled_raw_reftable_3,pooled_raw_reftable_4)

rm(pooled_raw_reftable_1)
rm(pooled_raw_reftable_2)
rm(pooled_raw_reftable_3)
rm(pooled_raw_reftable_4)

save(pooled_raw_reftable_model_1, file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/pooled_raw_reftable_model_1",".RData"))
#load(file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/pooled_raw_reftable_model_1",".RData"))
```

#### Prepare the reftable

``` r
# Check missing data
plot_missing(pooled_raw_reftable_model_1) 
missing_count <- sapply(pooled_raw_reftable_model_1, function(x) sum(is.na(x)))

# produce a data.frame with only summary statistics that contain missing values higher than a threshold
pooled_reftable_model_1 <- pooled_raw_reftable_model_1[,missing_count < nrow(pooled_raw_reftable_model_1)/10]

# check it!
plot_missing(pooled_reftable_model_1) 

# remove remaining rows with missing data
pooled_reftable_model_1 <- pooled_reftable_model_1[complete.cases(pooled_reftable_model_1), ]
save(pooled_reftable_model_1, file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/pooled_reftable_model_1",".RData"))
#load(file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/pooled_reftable_model_1",".RData"))

# check it again
plot_missing(pooled_reftable_model_1)
#missing_count <- sapply(pooled_reftable_model_1, function(x) sum(is.na(x)))

# Check != numeric values
#infinite_count <- sapply(pooled_reftable_model_1, function(x) sum(is.infinite(x)))

# Sample 1000 random simulations to use as PODs
random_pods <- sample(1:dim(pooled_reftable_model_1)[1], 1000, replace = F)
save(random_pods, file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/random_pods",".RData"))
#load(file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/random_pods",".RData"))

# Prepare the global summary statistics table by removing columns that do not contain summary stats
global_sumstats <- pooled_reftable_model_1[-random_pods, -c(1:169)]

save(global_sumstats, file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/global_sumstats",".RData"))
load(file=paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/global_sumstats",".RData"))
```

### ABC-RF Regression for Selection: Growing Trees

#### Genetic Load

``` r
averageGenLoad <- pooled_reftable_model_1[-random_pods, "averageGeneticLoad"]
#hist(averageGenLoad, freq = TRUE, xlab = "averaged genetic load", main = "", col = "#bdbdbd")

logaverageGenload <- -log10(1-averageGenLoad)
#hist(logaverageGenload, freq = TRUE, xlab = expression(-log[10]("averaged genetic load")), main = "", col = "#bdbdbd")

# abf-rf regression
reg_averageGenLoad <- regAbcrf(formula = logaverageGenload~.,
                               data    = data.frame(logaverageGenload, global_sumstats),
                               ntree   = 1000,
                               paral   = T,
                               ncores  = 28)

save(reg_averageGenLoad, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_averageGenLoad",".RData"))

# variable Importance plot
plot(x = reg_averageGenLoad, n.var = 20)

#head(sort(reg_averageGenLoad$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_averageGenLoad$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_averageGenLoad,
#             training = data.frame(logaverageGenload, global_sumstats),
#             paral    = T,
#             ncores   = 28)

# oob predictions vs true value
plot(x    = logaverageGenload,
     y    = reg_averageGenLoad$model.rf$predictions,
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,3),
     ylim = c(0,3),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### FDR 5% of the FST hypothesis of Ns&gt;1

``` r
fstfdrNs05 <- pooled_reftable_model_1[-random_pods, "FSTfdrNS05"]
#hist(fstfdrNs05, freq = TRUE, xlab = "5% FST FDR Ns>1", main = "", col = "#bdbdbd")

# abf-rf regression
reg_fstfdrNs05 <- regAbcrf(formula = fstfdrNs05~.,
                           data    = data.frame(fstfdrNs05, global_sumstats),
                           ntree   = 1000,
                           paral   = T,
                           ncores  = 28)

save(reg_fstfdrNs05, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_fstfdrNs05",".RData"))

# variable Importance plot
plot(x = reg_fstfdrNs05, n.var = 20)

#head(sort(reg_fstfdrNs05$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_fstfdrNs05$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_fstfdrNs05,
#             training = data.frame(fstfdrNs05, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = fstfdrNs05,
     y    = reg_fstfdrNs05$model.rf$predictions,
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,1),
     ylim = c(0,1),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Proportion of selected regions

``` r
## proportion of selected regions
prRSel <- pooled_reftable_model_1$PrGWSel[-random_pods] * pooled_reftable_model_1$PrMSel[-random_pods]

# PrRsel = prior (PrGWSel * PrMSel)
#hist(prRSel, freq = TRUE, xlab = "proportion of selected regions", main = "", col = "#bdbdbd")

# correlation between the compound parameter PrRSel and the calculated PrPOPMSel
cor(prRSel, pooled_reftable_model_1$PrPOPMSel[-random_pods])

plot(prRSel, pooled_reftable_model_1$PrPOPMSel[-random_pods], 
     xlab = "proportion of selected regions",
     ylab = "PrPOPMSel")
abline(lm(pooled_reftable_model_1$PrPOPMSel[-random_pods] ~ prRSel), col="#cb181d")

# correlation between the compound parameter PrRSel and the calculated PrSAMMSel
cor(prRSel, pooled_reftable_model_1$PrSAMMSel[-random_pods])

plot(prRSel, pooled_reftable_model_1$PrSAMMSel[-random_pods], 
     xlab = "proportion of selected regions",
     ylab = "PrSAMMSel")
abline(lm(pooled_reftable_model_1$PrSAMMSel[-random_pods] ~ prRSel), col="#cb181d")

# correlation between the the calculated PrPOPMSel and PrSAMMSel
cor(pooled_reftable_model_1$PrPOPMSel[-random_pods], pooled_reftable_model_1$PrSAMMSel[-random_pods])

plot(pooled_reftable_model_1$PrPOPMSel[-random_pods], pooled_reftable_model_1$PrSAMMSel[-random_pods], 
     xlab = "PrPOPMSel",
     ylab = "PrSAMMSel")
abline(lm(pooled_reftable_model_1$PrSAMMSel[-random_pods] ~ pooled_reftable_model_1$PrPOPMSel[-random_pods]), col="#cb181d")

# abf-rf regression
reg_prRSel <- regAbcrf(formula = prRSel~.,
                       data    = data.frame(prRSel, global_sumstats),
                       ntree   = 1000,
                       paral   = T,
                       ncores  = 28)

save(reg_prRSel, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_prRSel",".RData"))

# variable Importance plot
plot(x = reg_prRSel, n.var = 20)

#head(sort(reg_prRSel$model.rf$variable.importance, decreasing = T), n=10)

# prediction error
reg_prRSel$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_prRSel,
#             training = data.frame(prRSel, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = prRSel,
     y    = reg_prRSel$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,1),
     ylim = c(0,1),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## log prRSel
logrsel <- log10(prRSel)

# PrRsel = log prior (PrGWSel * PrMSel)
#hist(logrsel, freq = TRUE, xlab = expression(log[10]("proportion of selected regions")), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logrsel <- regAbcrf(formula = logrsel~.,
                        data    = data.frame(logrsel, global_sumstats),
                        ntree   = 1000,
                        paral   = T,
                        ncores  = 28)

save(reg_logrsel, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logrsel",".RData"))

# variable Importance plot
plot(x = reg_logrsel, n.var = 20)

#head(sort(reg_logrsel$model.rf$variable.importance, decreasing = T), n=10)

# prediction error
reg_logrsel$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logrsel,
#             training = data.frame(logrsel, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logrsel,
     y    = reg_logrsel$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-5,0),
     ylim = c(-5,0),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## PrPOPMSel
prPOPMSel <- pooled_reftable_model_1[-random_pods, "PrPOPMSel"]

#hist(prPOPMSel, freq = TRUE, xlab = "PrPOPMSel", main = "", col = "#bdbdbd")

# abf-rf regression
reg_prPOPMSel <- regAbcrf(formula = prPOPMSel~.,
                          data    = data.frame(prPOPMSel, global_sumstats),
                          ntree   = 1000,
                          paral   = T,
                          ncores  = 28)

save(reg_prPOPMSel, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_prPOPMSel",".RData"))

# variable Importance plot
plot(x = reg_prPOPMSel, n.var = 20)

#head(sort(reg_prPOPMSel$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_prPOPMSel$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_prPOPMSel,
#             training = data.frame(prPOPMSel, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = prPOPMSel,
     y    = reg_prPOPMSel$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,1),
     ylim = c(0,1),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## log PrPOPMSel
logpopmsel <- log10(prPOPMSel) 

#hist(logpopmsel, freq = TRUE, xlab = expression(log[10]("PrPOPMSel")), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logpopmsel <- regAbcrf(formula = logpopmsel~.,
                           data    = data.frame(logpopmsel, global_sumstats),
                           ntree   = 1000,
                           paral   = T,
                           ncores  = 28)

save(reg_logpopmsel, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logpopmsel",".RData"))

# variable Importance plot
plot(x = reg_logpopmsel, n.var = 20)

#head(sort(reg_logpopmsel$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logpopmsel$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logpopmsel,
#             training = data.frame(logpopmsel, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logpopmsel,
     y    = reg_logpopmsel$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-5,0),
     ylim = c(-5,0),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Proportion of strongly selected regions

``` r
# correlation between the the calculated PrPOPStrongMSel and PrSAMStrongMSel
cor(pooled_reftable_model_1$PrPOPStrongMSel[-random_pods], pooled_reftable_model_1$PrSAMStrongMSel[-random_pods] )

plot(x = pooled_reftable_model_1$PrPOPStrongMSel[-random_pods], 
     y = pooled_reftable_model_1$PrSAMStrongMSel[-random_pods], 
     xlab = "PrPOPStrongMSel",
     ylab = "PrSAMStrongMSel")
abline(lm(pooled_reftable_model_1$PrSAMStrongMSel[-random_pods] ~ pooled_reftable_model_1$PrPOPStrongMSel[-random_pods]), col="#cb181d")

# PrPOPStrongMSel
prPOPStrongMSel <- pooled_reftable_model_1[-random_pods, "PrPOPStrongMSel"]
#hist(prPOPStrongMSel, freq = TRUE, xlab = "PrPOPStrongMSel", main = "", col = "#bdbdbd")

# abf-rf regression
reg_prPOPStrongMSel <- regAbcrf(formula = prPOPStrongMSel~.,
                                data    = data.frame(prPOPStrongMSel, global_sumstats),
                                ntree   = 1000,
                                paral   = T,
                                ncores  = 28)

save(reg_prPOPStrongMSel, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_prPOPStrongMSel",".RData"))

# variable Importance plot
plot(x = reg_prPOPStrongMSel, n.var = 20)

#head(sort(reg_prPOPStrongMSel$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_prPOPStrongMSel$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_prPOPStrongMSel,
#             training = data.frame(prPOPStrongMSel, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = prPOPStrongMSel,
     y    = reg_prPOPStrongMSel$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,1),
     ylim = c(0,1),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## log10 PrPOPStrongMSel
logpopstrongmsel <- log10(pooled_reftable_model_1[-random_pods, "PrPOPStrongMSel"])
#logpopstrongmsel[which(logpopstrongmsel == -Inf)] <- -6

infindex <- which(logpopstrongmsel == -Inf)
logpopstrongmsel <- logpopstrongmsel[-infindex]

global_sumstats_logpopstrongmsel <- global_sumstats[-infindex,]

#hist(logpopstrongmsel, freq = TRUE, xlab = expression(log[10]("PrPOPStrongMSel")), main = "", col = "#bdbdbd")

reg_logpopstrongmsel <- regAbcrf(formula = logpopstrongmsel~.,
                                 data    = data.frame(logpopstrongmsel, global_sumstats_logpopstrongmsel),
                                 ntree   = 1000,
                                 paral   = T,
                                 ncores  = 28)

save(reg_logpopstrongmsel, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logpopstrongmsel",".RData"))

# variable Importance plot
plot(x = reg_logpopstrongmsel, n.var = 20)

#head(sort(reg_logpopstrongmsel$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logpopstrongmsel$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logpopstrongmsel,
#             training = data.frame(logpopstrongmsel, global_sumstats_logpopstrongmsel),
#             paral    = T,
#            ncores   = 28)

plot(x    = logpopstrongmsel,
     y    = reg_logpopstrongmsel$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-6,0),
     ylim = c(-6,0),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### DFE Gamma mean

``` r
gammamean <- pooled_reftable_model_1[-random_pods, "gammaMean"]

#hist(gammamean, freq = TRUE, xlab = expression(paste("gamma ", kappa, theta)), main = "", col = "#bdbdbd")

# abf-rf regression
reg_gammamean <- regAbcrf(formula = gammamean~.,
                          data    = data.frame(gammamean, global_sumstats),
                          ntree   = 1000,
                          paral   = T,
                          ncores  = 28)

save(reg_gammamean, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_gammamean",".RData"))

# variable Importance plot
plot(x = reg_gammamean, n.var = 20)

#head(sort(reg_gammamean$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_gammamean$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_gammamean,
#             training = data.frame(gammamean, global_sumstats),
#            paral    = T,
#             ncores   = 28)

plot(x    = gammamean,
     y    = reg_gammamean$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,1),
     ylim = c(0,1),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")


## log gammamean
loggammamean <- log10(gammamean)

#hist(loggammamean, freq = TRUE, xlab = expression(log[10](paste("gamma ", kappa, theta))), main = "", col = "#bdbdbd")

# abf-rf regression
reg_loggammamean <- regAbcrf(formula = loggammamean~.,
                             data    = data.frame(loggammamean, global_sumstats),
                             ntree   = 1000,
                             paral   = T,
                             ncores  = 28)

save(reg_loggammamean, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_loggammamean",".RData"))

# variable Importance plot
plot(x = reg_loggammamean, n.var = 20)

#head(sort(reg_loggammamean$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_loggammamean$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_loggammamean,
#             training = data.frame(loggammamean, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = loggammamean,
     y    = reg_loggammamean$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-3,0),
     ylim = c(-3,0),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Ns

``` r
logNs <- log10(pooled_reftable_model_1$PedigreeNetotal[-random_pods] * pooled_reftable_model_1[-random_pods, "gammaMean"])

#hist(logNs, freq = TRUE, xlab = expression(log[10](N[s])), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logNs <- regAbcrf(formula = logNs~.,
                      data    = data.frame(logNs, global_sumstats),
                      ntree   = 1000,
                      paral   = T,
                      ncores  = 28)

save(reg_logNs, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logNs",".RData"))

# variable Importance plot
plot(x = reg_logNs, n.var = 20)

#head(sort(reg_logNs$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logNs$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logNs,
#             training = data.frame(logNs, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logNs,
     y    = reg_logNs$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-3,3),
     ylim = c(-3,3),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### ThetaS

``` r
genomeS = 100e+5

logthetaS <- log10(4 * pooled_reftable_model_1$PedigreeNetotal[-random_pods] * (pooled_reftable_model_1$mu[-random_pods]*prRSel) * genomeS)

hist(logthetaS, freq = TRUE, xlab = expression(log[10](theta[S])), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logthetaS <- regAbcrf(formula = logthetaS~.,
                          data    = data.frame(logthetaS, global_sumstats),
                          ntree   = 1000,
                          paral   = T,
                          ncores  = 28)

save(reg_logthetaS, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logthetaS",".RData"))

# variable Importance plot
plot(x = reg_logthetaS, n.var = 20)

#head(sort(reg_logthetaS$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logthetaS$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logthetaS,
#             training = data.frame(logthetaS, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logthetaS,
     y    = reg_logthetaS$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-3,5),
     ylim = c(-3,5),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

### ABC-RF Regression for Demography: Growing Trees

#### Per site mutation rate mu

``` r
logmu <- log10(pooled_reftable_model_1$mu[-random_pods])

#hist(logmu, freq = TRUE, xlab = expression(log[10](mu)), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logmu <- regAbcrf(formula = logmu~.,
                      data    = data.frame(logmu, global_sumstats),
                      ntree   = 1000,
                      paral   = T,
                      ncores  = 28)

save(reg_logmu, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logmu",".RData"))

# variable Importance plot
plot(x = reg_logmu, n.var = 20)

#head(sort(reg_logmu$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logmu$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logmu,
#             training = data.frame(logmu, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logmu,
     y    = reg_logmu$model.rf$predictions,
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-8,-4),
     ylim = c(-8,-4),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Theta

``` r
logtheta <- log10(4 * pooled_reftable_model_1$PedigreeNetotal[-random_pods] * (pooled_reftable_model_1$mu[-random_pods]) * genomeS)

#hist(logtheta, freq = TRUE, xlab = expression(log[10](theta)), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logtheta <- regAbcrf(formula = logtheta~.,
                         data    = data.frame(logtheta, global_sumstats),
                         ntree   = 1000,
                         paral   = T,
                         ncores  = 28)

save(reg_logtheta, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logtheta",".RData"))

# variable Importance plot
plot(x = reg_logtheta, n.var = 20)

#head(sort(reg_logtheta$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logtheta$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logtheta,
#             training = data.frame(logtheta, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logtheta,
     y    = reg_logtheta$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-1,6),
     ylim = c(-1,6),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Per site recombination rate

``` r
logrr <- log10(pooled_reftable_model_1$rr[-random_pods])

#hist(logrr, freq = TRUE, xlab = expression(log[10](r)), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logrr <- regAbcrf(formula = logrr~.,
                      data    = data.frame(logrr, global_sumstats),
                      ntree   = 1000,
                      paral   = T,
                      ncores  = 28)

save(reg_logrr, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logrr",".RData"))

# variable Importance plot
plot(x = reg_logrr, n.var = 20)

#head(sort(reg_logmu$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logrr$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logrr,
#             training = data.frame(logrr, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logrr,
     y    = reg_logrr$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-7,-4),
     ylim = c(-7,-4),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Rho = Nr

``` r
logrho <- log10(pooled_reftable_model_1$rr[-random_pods] * pooled_reftable_model_1$PedigreeNetotal[-random_pods])

#hist(logrho, freq = TRUE, xlab = expression(log[10](rho)), main = "", col = "#bdbdbd")

# abf-rf regression
reg_logrho <- regAbcrf(formula = logrho~.,
                       data    = data.frame(logrho, global_sumstats),
                       ntree   = 1000,
                       paral   = T,
                       ncores  = 28)

save(reg_logrho, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logrho",".RData"))

# variable Importance plot
plot(x = reg_logrho, n.var = 20)

#head(sort(reg_logrho$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logrho$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logrho,
#             training = data.frame(logrho, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logrho,
     y    = reg_logrho$model.rf$predictions, 
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-9,1),
     ylim = c(-9,1),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Pedigree Effective Population Size - PedigreeNetotal

``` r
pedigreeNetotal <- pooled_reftable_model_1[-random_pods, "PedigreeNetotal"]

logpedigreeNetotal <- log10(pedigreeNetotal)

#hist(logpedigreeNetotal, freq = TRUE, xlab = expression(log[10]("pedigree N"[e])), main = "", col = "#bdbdbd")

reg_logpedigreeNetotal <- regAbcrf(formula = logpedigreeNetotal~.,
                                   data    = data.frame(logpedigreeNetotal, global_sumstats),
                                   ntree   = 1000,
                                   paral   = T,
                                   ncores  = 28)

save(reg_logpedigreeNetotal, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logpedigreeNetotal",".RData"))

# Variable Importance plot
plot(x = reg_logpedigreeNetotal, n.var = 20)

#head(sort(reg_logpedigreeNetotal$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logpedigreeNetotal$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logpedigreeNetotal,
#             training = data.frame(logpedigreeNetotal, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logpedigreeNetotal,
     y    = reg_logpedigreeNetotal$model.rf$predictions,
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(-2,4),
     ylim = c(-2,4),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

#### Population Census Size for the Sampling Phase - N

``` r
n <- pooled_reftable_model_1[-random_pods, "N"]

logn <- log10(n)

#hist(logn, freq = TRUE, xlab = expression(log[10]("N"[c])), main = "", col = "#bdbdbd")

reg_logn <- regAbcrf(formula = logn~.,
                        data = data.frame(logn, global_sumstats),
                       ntree = 1000,
                       paral = T,
                      ncores = 28)

save(reg_logn, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/reg_logn",".RData"))

# Variable Importance plot
plot(x = reg_logn, n.var = 20)

#head(sort(reg_logn$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logn$model.rf$prediction.error

# error rate plot
#err.regAbcrf(object   = reg_logn,
#             training = data.frame(logn, global_sumstats),
#             paral    = T,
#             ncores   = 28)

plot(x    = logn,
     y    = reg_logn$model.rf$predictions,
     xlab = "True value",
     ylab = "OOB predictions",
     xlim = c(0,3),
     ylim = c(0,3),
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")
```

### ABC-RF Regression for Selection: Prediction

#### RANDOM PODS

``` r
global_sumstats_randompods <- pooled_reftable_model_1[random_pods, -c(1:169)]

save(global_sumstats_randompods, file = paste0("/home/pavinato/My_repositories/Tracking-selection/src/analysis/RF_analysis_files/global_sumstats_randompods",".RData"))
```

#### Load PODs

``` r
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/PODs/model_Neutral/reference_table_neutral.RData")

neutral_pods_reftable <- raw_reftable[, names(raw_reftable) %in% names(pooled_reftable_model_1), drop=FALSE]
  
global_sumstats_neutralpods <- neutral_pods_reftable[, -c(1:169)]
```

``` r
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/PODs/model_DN/LSelection/reference_table_DN_LSelection.RData")

dnlsel_pods_reftable <- raw_reftable[, names(raw_reftable) %in% names(pooled_reftable_model_1), drop=FALSE]
  
global_sumstats_dnlselpods <- dnlsel_pods_reftable[, -c(1:169)]
```

``` r
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/PODs/model_DN/HSelection/reference_table_DN_HSelection.RData")

dnhsel_pods_reftable <- raw_reftable[, names(raw_reftable) %in% names(pooled_reftable_model_1), drop=FALSE]
  
global_sumstats_dnhselpods <- dnhsel_pods_reftable[, -c(1:169)]
```

#### Genetic Load

``` r
## Random PODs
log_radnompods_averageGenLoad <- -log10(1-pooled_reftable_model_1[random_pods, "averageGeneticLoad"])

posterior_radnompods_averageGenLoad <- predict(object     = reg_averageGenLoad,
                                               obs        = global_sumstats_randompods,
                                               training   = data.frame(logaverageGenload, global_sumstats),
                                               quantiles  = c(0.025,0.5,0.975),
                                               paral      = T,
                                               ncores     = 28,
                                               rf.weights = T)
(posterior_radnompods_averageGenLoad)

# Expected True vs estimated
plot(x    = log_radnompods_averageGenLoad, 
     y    = posterior_radnompods_averageGenLoad$expectation,
     xlim = c(0,3),
     ylim = c(0,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_averageGenLoad, 
     y = posterior_radnompods_averageGenLoad$med,
     xlim = c(0,3),
     ylim = c(0,3),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## Fixed PODs 

# Neutral PODs
log_neutralpods_averageGenLoad <- -log10(1-neutral_pods_reftable[, "averageGeneticLoad"])

posterior_neutralpods_averageGenLoad <- predict(object     = reg_averageGenLoad,
                                                obs        = global_sumstats_neutralpods,
                                                training   = data.frame(logaverageGenload, global_sumstats),
                                                quantiles  = c(0.025,0.5,0.975),
                                                paral      = T,
                                                ncores     = 28,
                                                rf.weights = T)
(posterior_neutralpods_averageGenLoad)

# Low Selection PODs
log_dnlselpods_averageGenLoad <- -log10(1-dnlsel_pods_reftable[, "averageGeneticLoad"])

posterior_dnlselpods_averageGenLoad <- predict(object     = reg_averageGenLoad,
                                               obs        = global_sumstats_dnlselpods,
                                               training   = data.frame(logaverageGenload, global_sumstats),
                                               quantiles  = c(0.025,0.5,0.975),
                                               paral      = T,
                                               ncores     = 28,
                                               rf.weights = T)
(posterior_dnlselpods_averageGenLoad)

# High Selection PODs
log_dnhselpods_averageGenLoad <- -log10(1-dnhsel_pods_reftable[, "averageGeneticLoad"])

posterior_dnhselpods_averageGenLoad <- predict(object     = reg_averageGenLoad,
                                               obs        = global_sumstats_dnhselpods,
                                               training   = data.frame(logaverageGenload, global_sumstats),
                                               quantiles  = c(0.025,0.5,0.975),
                                               paral      = T,
                                               ncores     = 28,
                                               rf.weights = T)
(posterior_dnhselpods_averageGenLoad)

# Expected True vs estimated
plot(x    = c(log_neutralpods_averageGenLoad,log_dnlselpods_averageGenLoad, log_dnhselpods_averageGenLoad), 
     y    = c(posterior_neutralpods_averageGenLoad$expectation, posterior_dnlselpods_averageGenLoad$expectation, posterior_dnhselpods_averageGenLoad$expectation),
     xlim = c(0,2),
     ylim = c(0,2),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x    = c(log_neutralpods_averageGenLoad,log_dnlselpods_averageGenLoad, log_dnhselpods_averageGenLoad), 
     y    = c(posterior_neutralpods_averageGenLoad$med, posterior_dnlselpods_averageGenLoad$med, posterior_dnhselpods_averageGenLoad$med),
     xlim = c(0,2),
     ylim = c(0,2),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

# Boxplot

## boxplot
pred_averageGenLoad <- data.frame(model  = c(rep("Neutral", 300), rep("Low_selection", 300), rep("High_selection", 300)), 
                                  type   = c(rep("True",100), rep("Expected",100), rep("Median",100), 
                                             rep("True",100), rep("Expected",100), rep("Median",100),
                                             rep("True",100), rep("Expected",100), rep("Median",100)),
                                  values = c(log_neutralpods_averageGenLoad, posterior_neutralpods_averageGenLoad$expectation, posterior_neutralpods_averageGenLoad$med,
                                             log_dnlselpods_averageGenLoad, posterior_dnlselpods_averageGenLoad$expectation, posterior_dnlselpods_averageGenLoad$med,
                                             log_dnhselpods_averageGenLoad, posterior_dnhselpods_averageGenLoad$expectation, posterior_dnhselpods_averageGenLoad$med))

boxplot(pred_averageGenLoad$values ~ pred_averageGenLoad$type * pred_averageGenLoad$model,
        lwd = 2,
        ylim = c(0,2),
        ylab = expression(-log[10]*(1 - "average genetic load")),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3),rep("#969696",3)))
legend("topright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:9, labels = FALSE)
text(x = 1:9, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

#### FDR 5% of the FST hypothesis of Ns&gt;1

``` r
## Random PODs
radnompods_fstfdrNs05 <- pooled_reftable_model_1[random_pods, "FSTfdrNS05"]

posterior_radnompods_fstfdrNs05 <- predict(object     = reg_fstfdrNs05,
                                           obs        = global_sumstats_randompods,
                                           training   = data.frame(fstfdrNs05, global_sumstats),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 28,
                                           rf.weights = T)
(posterior_radnompods_fstfdrNs05)

# Expected True vs estimated
plot(x    = radnompods_fstfdrNs05, 
     y    = posterior_radnompods_fstfdrNs05$expectation,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = radnompods_fstfdrNs05, 
     y = posterior_radnompods_fstfdrNs05$med,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## Fixed PODs 

# Neutral PODs
neutralpods_fstfdrNs05 <- neutral_pods_reftable[, "FSTfdrNS05"]

posterior_neutralpods_fstfdrNs05 <- predict(object     = reg_fstfdrNs05,
                                            obs        = global_sumstats_neutralpods,
                                            training   = data.frame(fstfdrNs05, global_sumstats),
                                            quantiles  = c(0.025,0.5,0.975),
                                            paral      = T,
                                            ncores     = 28,
                                            rf.weights = T)
(posterior_neutralpods_fstfdrNs05)

# Low Selection PODs
dnlselpods_fstfdrNs05 <- dnlsel_pods_reftable[, "FSTfdrNS05"]

posterior_dnlselpods_fstfdrNs05 <- predict(object     = reg_fstfdrNs05,
                                           obs        = global_sumstats_dnlselpods,
                                           training   = data.frame(fstfdrNs05, global_sumstats),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 28,
                                           rf.weights = T)
(posterior_dnlselpods_fstfdrNs05)

# High Selection PODs
dnhselpods_fstfdrNs05 <- dnhsel_pods_reftable[, "FSTfdrNS05"]

posterior_dnhselpods_fstfdrNs05 <- predict(object     = reg_fstfdrNs05,
                                           obs        = global_sumstats_dnhselpods,
                                           training   = data.frame(fstfdrNs05, global_sumstats),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 28,
                                           rf.weights = T)
(posterior_dnhselpods_fstfdrNs05)

# Expected True vs estimated
plot(x    = c(neutralpods_fstfdrNs05,dnlselpods_fstfdrNs05, dnhselpods_fstfdrNs05), 
     y    = c(posterior_neutralpods_fstfdrNs05$expectation, posterior_dnlselpods_fstfdrNs05$expectation, posterior_dnhselpods_fstfdrNs05$expectation),
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x    = c(neutralpods_fstfdrNs05,dnlselpods_fstfdrNs05, dnhselpods_fstfdrNs05), 
     y    = c(posterior_neutralpods_fstfdrNs05$med, posterior_dnlselpods_fstfdrNs05$med, posterior_dnhselpods_fstfdrNs05$med),
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

## boxplot
pred_fstfdrNs05 <- data.frame(model  = c(rep("Neutral", 300), rep("Low_selection", 300), rep("High_selection", 300)), 
                              type   = c(rep("True",100), rep("Expected",100), rep("Median",100), 
                                         rep("True",100), rep("Expected",100), rep("Median",100),
                                         rep("True",100), rep("Expected",100), rep("Median",100)),
                              values = c(neutralpods_fstfdrNs05, posterior_neutralpods_fstfdrNs05$expectation, posterior_neutralpods_fstfdrNs05$med,
                                         dnlselpods_fstfdrNs05, posterior_dnlselpods_fstfdrNs05$expectation, posterior_dnlselpods_fstfdrNs05$med,
                                         dnhselpods_fstfdrNs05, posterior_dnhselpods_fstfdrNs05$expectation, posterior_dnhselpods_fstfdrNs05$med))

boxplot(pred_fstfdrNs05$values ~ pred_fstfdrNs05$type * pred_fstfdrNs05$model,
        lwd = 2,
        ylim = c(0,1),
        ylab = "FDR 5% of the FST hypothesis of Ns>1",
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3),rep("#969696",3)))
legend("topright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:9, labels = FALSE)
text(x = 1:9, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

#### Proportion of selected mutations

``` r
## Random PODs
log_radnompods_popmsel <- log10(pooled_reftable_model_1[random_pods, "PrPOPMSel"])

posterior_radnompods_logpopmsel <- predict(object     = reg_logpopmsel,
                                           obs        = global_sumstats_randompods,
                                           training   = data.frame(logpopmsel, global_sumstats),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 28,
                                           rf.weights = T)
(posterior_radnompods_logpopmsel)

# Expected True vs estimated
plot(x    = log_radnompods_popmsel, 
     y    = posterior_radnompods_logpopmsel$expectation,
     xlim = c(-5,0),
     ylim = c(-5,0),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_popmsel, 
     y = posterior_radnompods_logpopmsel$med,
     xlim = c(-5,0),
     ylim = c(-5,0),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## Fixed PODs 

# Low Selection PODs
log_dnlselpods_popmsel <- log10(dnlsel_pods_reftable[, "PrPOPMSel"])

posterior_dnlselpods_logpopmsel <- predict(object     = reg_logpopmsel,
                                           obs        = global_sumstats_dnlselpods,
                                           training   = data.frame(logpopmsel, global_sumstats),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 28,
                                           rf.weights = T)
(posterior_dnlselpods_logpopmsel)

# High Selection PODs
log_dnhselpods_popmsel <- log10(dnhsel_pods_reftable[, "PrPOPMSel"])

posterior_dnhselpods_logpopmsel <- predict(object     = reg_logpopmsel,
                                           obs        = global_sumstats_dnhselpods,
                                           training   = data.frame(logpopmsel, global_sumstats),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 28,
                                           rf.weights = T)
(posterior_dnhselpods_logpopmsel)

# Expected True vs estimated
plot(x    = c(log_dnlselpods_popmsel, log_dnhselpods_popmsel), 
     y    = c(posterior_dnlselpods_logpopmsel$expectation, posterior_dnhselpods_logpopmsel$expectation),
     xlim = c(-5,0),
     ylim = c(-5,0),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")


# Median True vs estimated
plot(x    = c(log_dnlselpods_popmsel, log_dnhselpods_popmsel), 
     y    = c(posterior_dnlselpods_logpopmsel$med, posterior_dnhselpods_logpopmsel$med),
     xlim = c(-5,0),
     ylim = c(-5,0),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

## boxplot
pred_logpopmsel <- data.frame(model  = c(rep("Low_selection", 300), rep("High_selection", 300)), 
                              type   = c(rep("True",100), rep("Expected",100), rep("Median",100),
                                         rep("True",100), rep("Expected",100), rep("Median",100)),
                              values = c(log_dnlselpods_popmsel, posterior_dnlselpods_logpopmsel$expectation, posterior_dnlselpods_logpopmsel$med,
                                         log_dnhselpods_popmsel, posterior_dnhselpods_logpopmsel$expectation, posterior_dnhselpods_logpopmsel$med))

boxplot(pred_logpopmsel$values ~ pred_logpopmsel$type * pred_logpopmsel$model,
        lwd = 2,
        ylab = expression(-log[10]*("Pr Selected mutations")),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3)))
legend("topright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:6, labels = FALSE)
text(x = 1:6, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

#### Proportion of strongly selected mutations

``` r
## Random PODs
log_radnompods_popstrongmsel <- log10(pooled_reftable_model_1[random_pods, "PrPOPStrongMSel"])

infindex_randompods <- which(log_radnompods_popstrongmsel == -Inf)
log_radnompods_popstrongmsel <- log_radnompods_popstrongmsel[-infindex_randompods]

global_sumstats_randompods_logpopstrongmsel <- global_sumstats_randompods[-infindex_randompods,]


posterior_radnompods_logpopstrongmsel <- predict(object     = reg_logpopstrongmsel,
                                                 obs        = global_sumstats_randompods_logpopstrongmsel,
                                                 training   = data.frame(logpopstrongmsel, global_sumstats_logpopstrongmsel),
                                                 quantiles  = c(0.025,0.5,0.975),
                                                 paral      = T,
                                                 ncores     = 28,
                                                 rf.weights = T)
(posterior_radnompods_logpopstrongmsel)

# Expected True vs estimated
plot(x    = log_radnompods_popstrongmsel, 
     y    = posterior_radnompods_logpopstrongmsel$expectation,
     xlim = c(-6,1),
     ylim = c(-6,1),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_popstrongmsel, 
     y = posterior_radnompods_logpopstrongmsel$med,
     xlim = c(-6,1),
     ylim = c(-6,1),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## Fixed PODs 

# Low Selection PODs
log_dnlselpods_popstrongmsel <- log10(dnlsel_pods_reftable[, "PrPOPStrongMSel"])

posterior_dnlselpods_logpopstrongmsel <- predict(object     = reg_logpopstrongmsel,
                                                 obs        = global_sumstats_dnlselpods,
                                                 training   = data.frame(logpopstrongmsel, global_sumstats_logpopstrongmsel),
                                                 quantiles  = c(0.025,0.5,0.975),
                                                 paral      = T,
                                                 ncores     = 28,
                                                 rf.weights = T)
(posterior_dnlselpods_logpopstrongmsel)

# High Selection PODs
log_dnhselpods_popstrongmsel <- log10(dnhsel_pods_reftable[, "PrPOPStrongMSel"])

posterior_dnhselpods_logpopstrongmsel <- predict(object     = reg_logpopstrongmsel,
                                                 obs        = global_sumstats_dnhselpods,
                                                 training   = data.frame(logpopstrongmsel, global_sumstats_logpopstrongmsel),
                                                 quantiles  = c(0.025,0.5,0.975),
                                                 paral      = T,
                                                 ncores     = 28,
                                                 rf.weights = T)
(posterior_dnhselpods_logpopstrongmsel)

# Expected True vs estimated
plot(x    = c(log_dnlselpods_popstrongmsel, log_dnhselpods_popstrongmsel), 
     y    = c(posterior_dnlselpods_logpopstrongmsel$expectation, posterior_dnhselpods_logpopstrongmsel$expectation),
     xlim = c(-6,1),
     ylim = c(-6,1),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")


# Median True vs estimated
plot(x    = c(log_dnlselpods_popstrongmsel, log_dnhselpods_popstrongmsel), 
     y    = c(posterior_dnlselpods_logpopstrongmsel$med, posterior_dnhselpods_logpopstrongmsel$med),
     xlim = c(-6,1),
     ylim = c(-6,1),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

# Boxplot

## boxplot
pred_logpopstrongmsel <- data.frame(model  = c(rep("Low_selection", 300), rep("High_selection", 300)), 
                                    type   = c(rep("True",100), rep("Expected",100), rep("Median",100),
                                               rep("True",100), rep("Expected",100), rep("Median",100)),
                                    values = c(log_dnlselpods_popstrongmsel, posterior_dnlselpods_logpopstrongmsel$expectation, posterior_dnlselpods_logpopstrongmsel$med,
                                               log_dnhselpods_popstrongmsel, posterior_dnhselpods_logpopstrongmsel$expectation, posterior_dnhselpods_logpopstrongmsel$med))

boxplot(pred_logpopstrongmsel$values ~ pred_logpopstrongmsel$type * pred_logpopstrongmsel$model,
        lwd = 2,
        ylim = c(-4,0),
        ylab = expression(log[10]*("Pr Strongly selected mutations")),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3)))
legend("topright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:6, labels = FALSE)
text(x = 1:6, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

#### Ns

``` r
## Random PODs
log_radnompods_Ns <- log10(pooled_reftable_model_1$PedigreeNetotal[random_pods] * pooled_reftable_model_1$gammaMean[random_pods])


posterior_radnompods_logNs <- predict(object     = reg_logNs,
                                      obs        = global_sumstats_randompods,
                                      training   = data.frame(logNs, global_sumstats),
                                      quantiles  = c(0.025,0.5,0.975),
                                      paral      = T,
                                      ncores     = 28,
                                      rf.weights = T)

(posterior_radnompods_logNs)

# Expected True vs estimated
plot(x    = log_radnompods_Ns, 
     y    = posterior_radnompods_logNs$expectation,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_Ns, 
     y = posterior_radnompods_logNs$med,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")


## Fixed PODs 

# Low Selection PODs
log_dnlselpods_Ns <- log10(dnlsel_pods_reftable$PedigreeNetotal * dnlsel_pods_reftable$gammaMean)

posterior_dnlselpods_logNs <- predict(object     = reg_logNs,
                                      obs        = global_sumstats_dnlselpods,
                                      training   = data.frame(logNs, global_sumstats),
                                      quantiles  = c(0.025,0.5,0.975),
                                      paral      = T,
                                      ncores     = 28,
                                      rf.weights = T)
(posterior_dnlselpods_logNs)

# High Selection PODs
log_dnhselpods_Ns <- log10(dnhsel_pods_reftable$PedigreeNetotal * dnhsel_pods_reftable$gammaMean)

posterior_dnhselpods_logNs <- predict(object     = reg_logNs,
                                      obs        = global_sumstats_dnhselpods,
                                      training   = data.frame(logNs, global_sumstats),
                                      quantiles  = c(0.025,0.5,0.975),
                                      paral      = T,
                                      ncores     = 28,
                                      rf.weights = T)

(posterior_dnhselpods_logNs)

# Expected True vs estimated
plot(x    = c(log_dnlselpods_Ns, log_dnhselpods_Ns), 
     y    = c(posterior_dnlselpods_logNs$expectation, posterior_dnhselpods_logNs$expectation),
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")


# Median True vs estimated
plot(x    = c(log_dnlselpods_Ns, log_dnhselpods_Ns), 
     y    = c(posterior_dnlselpods_logNs$med, posterior_dnhselpods_logNs$med),
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

## boxplot
pred_logNs <- data.frame(model  = c(rep("Low_selection", 300), rep("High_selection", 300)), 
                         type   = c(rep("True",100), rep("Expected",100), rep("Median",100),
                                    rep("True",100), rep("Expected",100), rep("Median",100)),
                                    values = c(log_dnlselpods_Ns, posterior_dnlselpods_logNs$expectation, posterior_dnlselpods_logNs$med,
                                               log_dnhselpods_Ns, posterior_dnhselpods_logNs$expectation, posterior_dnhselpods_logNs$med))

boxplot(pred_logNs$values ~ pred_logNs$type * pred_logNs$model,
        lwd = 2,
        ylim = c(0,3),
        ylab = expression(log[10]*(Ns)),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3)))
legend("topright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:6, labels = FALSE)
text(x = 1:6, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

#### ThetaS

``` r
## Random PODs
prRSel_randompods <- pooled_reftable_model_1$PrGWSel[random_pods] * pooled_reftable_model_1$PrMSel[random_pods]

log_radnompods_thetaS <- log10(4 * pooled_reftable_model_1$PedigreeNetotal[random_pods] * (pooled_reftable_model_1$mu[random_pods] * prRSel_randompods) * genomeS)


posterior_radnompods_logthetaS <- predict(object     = reg_logthetaS,
                                          obs        = global_sumstats_randompods,
                                          training   = data.frame(logthetaS, global_sumstats),
                                          quantiles  = c(0.025,0.5,0.975),
                                          paral      = T,
                                          ncores     = 28,
                                          rf.weights = T)

(posterior_radnompods_logthetaS)

# Expected True vs estimated
plot(x    = log_radnompods_thetaS, 
     y    = posterior_radnompods_logthetaS$expectation,
     xlim = c(-5,5),
     ylim = c(-5,5),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_thetaS, 
     y = posterior_radnompods_logthetaS$med,
     xlim = c(-5,5),
     ylim = c(-5,5),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")


## Fixed PODs 

# Low Selection PODs
log_dnlselpods_thetaS <- log10(4 * dnlsel_pods_reftable$PedigreeNetotal * (dnlsel_pods_reftable$mu * dnlsel_pods_reftable$PrGWSel * dnlsel_pods_reftable$PrMSel) * genomeS)


posterior_dnlselpods_logthetaS <- predict(object     = reg_logthetaS,
                                          obs        = global_sumstats_dnlselpods,
                                          training   = data.frame(logthetaS, global_sumstats),
                                          quantiles  = c(0.025,0.5,0.975),
                                          paral      = T,
                                          ncores     = 28,
                                          rf.weights = T)
(posterior_dnlselpods_logthetaS)

# High Selection PODs
log_dnhselpods_thetaS <- log10(4 * dnhsel_pods_reftable$PedigreeNetotal * (dnhsel_pods_reftable$mu * dnhsel_pods_reftable$PrGWSel * dnhsel_pods_reftable$PrMSel) * genomeS)

posterior_dnhselpods_logthetaS <- predict(object     = reg_logthetaS,
                                          obs        = global_sumstats_dnhselpods,
                                          training   = data.frame(logthetaS, global_sumstats),
                                          quantiles  = c(0.025,0.5,0.975),
                                          paral      = T,
                                          ncores     = 28,
                                          rf.weights = T)

(posterior_dnhselpods_logthetaS)

# Expected True vs estimated
plot(x    = c(log_dnlselpods_thetaS, log_dnhselpods_thetaS), 
     y    = c(posterior_dnlselpods_logthetaS$expectation, posterior_dnhselpods_logthetaS$expectation),
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")


# Median True vs estimated
plot(x    = c(log_dnlselpods_thetaS, log_dnhselpods_thetaS), 
     y    = c(posterior_dnlselpods_logthetaS$med, posterior_dnhselpods_logthetaS$med),
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

## boxplot
pred_logthetaS <- data.frame(model  = c(rep("Low_selection", 300), rep("High_selection", 300)), 
                             type   = c(rep("True",100), rep("Expected",100), rep("Median",100),
                                        rep("True",100), rep("Expected",100), rep("Median",100)),
                             values = c(log_dnlselpods_thetaS, posterior_dnlselpods_logthetaS$expectation, posterior_dnlselpods_logthetaS$med,
                                        log_dnhselpods_thetaS, posterior_dnhselpods_logthetaS$expectation, posterior_dnhselpods_logthetaS$med))

boxplot(pred_logthetaS$values ~ pred_logthetaS$type * pred_logthetaS$model,
        lwd = 2,
        ylim = c(0,3),
        ylab = expression(log[10]*(theta[s])),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3)))
legend("bottomright", legend = c("1.0% Selection", "2.5% Selection"), col = c("#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:6, labels = FALSE)
text(x = 1:6, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

### ABC-RF Regression for Demography: Prediction

#### Pedigree Effective Population Size - PedigreeNetotal

``` r
## Random PODs
log_radnompods_pedigreeNetotal <- log10(pooled_reftable_model_1[random_pods, "PedigreeNetotal"])

posterior_radnompods_logpedigreeNetotal <- predict(object     = reg_logpedigreeNetotal,
                                                   obs        = global_sumstats_randompods,
                                                   training   = data.frame(logpedigreeNetotal, global_sumstats),
                                                   quantiles  = c(0.025,0.5,0.975),
                                                   paral      = T,
                                                   ncores     = 28,
                                                   rf.weights = T)
(posterior_radnompods_logpedigreeNetotal)

# Expected True vs estimated
plot(x    = log_radnompods_pedigreeNetotal, 
     y    = posterior_radnompods_logpedigreeNetotal$expectation,
     xlim = c(-2,3),
     ylim = c(-2,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_pedigreeNetotal, 
     y = posterior_radnompods_logpedigreeNetotal$med,
     xlim = c(-2,3),
     ylim = c(-2,3),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## Fixed PODs 

# Neutral PODs
log_neutralpods_pedigreeNetotal <- log10(neutral_pods_reftable[, "PedigreeNetotal"])

posterior_neutralpods_logpedigreeNetotal <- predict(object     = reg_logpedigreeNetotal,
                                                    obs        = global_sumstats_neutralpods,
                                                    training   = data.frame(logpedigreeNetotal, global_sumstats),
                                                    quantiles  = c(0.025,0.5,0.975),
                                                    paral      = T,
                                                    ncores     = 28,
                                                    rf.weights = T)
(posterior_neutralpods_logpedigreeNetotal)

# Low Selection PODs
log_dnlselpods_pedigreeNetotal <- log10(dnlsel_pods_reftable[, "PedigreeNetotal"])

posterior_dnlselpods_logpedigreeNetotal <- predict(object     = reg_logpedigreeNetotal,
                                                   obs        = global_sumstats_dnlselpods,
                                                   training   = data.frame(logpedigreeNetotal, global_sumstats),
                                                   quantiles  = c(0.025,0.5,0.975),
                                                   paral      = T,
                                                   ncores     = 28,
                                                   rf.weights = T)
(posterior_dnlselpods_logpedigreeNetotal)

# High Selection PODs
log_dnhselpods_pedigreeNetotal <- log10(dnhsel_pods_reftable[, "PedigreeNetotal"])

posterior_dnhselpods_logpedigreeNetotal <- predict(object     = reg_logpedigreeNetotal,
                                                   obs        = global_sumstats_dnhselpods,
                                                   training   = data.frame(logpedigreeNetotal, global_sumstats),
                                                   quantiles  = c(0.025,0.5,0.975),
                                                   paral      = T,
                                                   ncores     = 28,
                                                   rf.weights = T)
(posterior_dnhselpods_logpedigreeNetotal)

# Expected True vs estimated
plot(x    = c(log_neutralpods_pedigreeNetotal,log_dnlselpods_pedigreeNetotal, log_dnhselpods_pedigreeNetotal), 
     y    = c(posterior_neutralpods_logpedigreeNetotal$expectation, posterior_dnlselpods_logpedigreeNetotal$expectation, posterior_dnhselpods_logpedigreeNetotal$expectation),
     xlim = c(1,3),
     ylim = c(1,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x    = c(log_neutralpods_pedigreeNetotal,log_dnlselpods_pedigreeNetotal, log_dnhselpods_pedigreeNetotal), 
     y    = c(posterior_neutralpods_logpedigreeNetotal$med, posterior_dnlselpods_logpedigreeNetotal$med, posterior_dnhselpods_logpedigreeNetotal$med),
     xlim = c(1,3),
     ylim = c(1,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

## boxplot
pred_logpedigreeNetotal <- data.frame(model  = c(rep("Neutral", 300), rep("Low_selection", 300), rep("High_selection", 300)), 
                                  type   = c(rep("True",100), rep("Expected",100), rep("Median",100), 
                                             rep("True",100), rep("Expected",100), rep("Median",100),
                                             rep("True",100), rep("Expected",100), rep("Median",100)),
                                  values = c(log_neutralpods_pedigreeNetotal, posterior_neutralpods_logpedigreeNetotal$expectation, posterior_neutralpods_logpedigreeNetotal$med,
                                             log_dnlselpods_pedigreeNetotal, posterior_dnlselpods_logpedigreeNetotal$expectation, posterior_dnlselpods_logpedigreeNetotal$med,
                                             log_dnhselpods_pedigreeNetotal, posterior_dnhselpods_logpedigreeNetotal$expectation, posterior_dnhselpods_logpedigreeNetotal$med))

boxplot(pred_logpedigreeNetotal$values ~ pred_logpedigreeNetotal$type * pred_logpedigreeNetotal$model,
        lwd = 2,
        ylim = c(1,3),
        ylab = expression(log[10]*("N"[e])),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3),rep("#969696",3)))
abline(h=log10(500), col = "green", lty = 3)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:9, labels = FALSE)
text(x = 1:9, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)
```

#### Population Census Size for the Sampling Phase - N

``` r
log_radnompods_n <- log10(pooled_reftable_model_1[random_pods, "N"])

posterior_radnompods_logn <- predict(object     = reg_logn,
                                     obs        = global_sumstats_randompods,
                                     training   = data.frame(logn, global_sumstats),
                                     quantiles  = c(0.025,0.5,0.975),
                                     paral      = T,
                                     ncores     = 28,
                                     rf.weights = T)
(posterior_radnompods_logn)

# Expected True vs estimated
plot(x    = log_radnompods_n, 
     y    = posterior_radnompods_logn$expectation,
     xlim = c(0,3),
     ylim = c(0,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Random PODs",
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x = log_radnompods_n, 
     y = posterior_radnompods_logn$med,
     xlim = c(0,3),
     ylim = c(0,3),
     xlab="True values", 
     ylab="Estimated values",
     main="Median Random PODs",
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
abline(a=0, b=1, col = "#cb181d")

## Fixed PODs 

# Neutral PODs
log_neutralpods_n <- log10(neutral_pods_reftable[, "N"])

posterior_neutralpods_logn <- predict(object     = reg_logn,
                                      obs        = global_sumstats_neutralpods,
                                      training   = data.frame(logn, global_sumstats),
                                      quantiles  = c(0.025,0.5,0.975),
                                      paral      = T,
                                      ncores     = 28,
                                      rf.weights = T)
(posterior_neutralpods_logn)

# Low Selection PODs
log_dnlselpods_n <- log10(dnlsel_pods_reftable[, "N"])

posterior_dnlselpods_logn <- predict(object     = reg_logn,
                                     obs        = global_sumstats_dnlselpods,
                                     training   = data.frame(logn, global_sumstats),
                                     quantiles  = c(0.025,0.5,0.975),
                                     paral      = T,
                                     ncores     = 28,
                                     rf.weights = T)
(posterior_dnlselpods_logn)

# High Selection PODs
log_dnhselpods_n <- log10(dnhsel_pods_reftable[, "N"])

posterior_dnhselpods_logn <- predict(object     = reg_logn,
                                     obs        = global_sumstats_dnhselpods,
                                     training   = data.frame(logn, global_sumstats),
                                     quantiles  = c(0.025,0.5,0.975),
                                     paral      = T,
                                     ncores     = 28,
                                     rf.weights = T)
(posterior_dnhselpods_logpn)

# Expected True vs estimated
plot(x    = c(log_neutralpods_n,log_dnlselpods_n, log_dnhselpods_n), 
     y    = c(posterior_neutralpods_logn$expectation, posterior_dnlselpods_logn$expectation, posterior_dnhselpods_logn$expectation),
     xlim = c(0,3),
     ylim = c(0,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Expected Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch=1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

# Median True vs estimated
plot(x    = c(log_neutralpods_n,log_dnlselpods_n, log_dnhselpods_n), 
     y    = c(posterior_neutralpods_logn$med, posterior_dnlselpods_logn$med, posterior_dnhselpods_logn$med),
     xlim = c(0,3),
     ylim = c(0,3),
     xlab = "True values", 
     ylab = "Estimated values",
     main = "Median Fixed PODs",
     col  = c(rep("#969696",100), rep("#08519c",100), rep("#a50f15",100)),
     pch  = 1,
     cex.axis = 1.2,
     cex.lab  = 1.5)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 1, bty = "n", cex = 1.2)
abline(a=0, b=1, col = "#cb181d")

## boxplot
pred_logn <- data.frame(model  = c(rep("Neutral", 300), rep("Low_selection", 300), rep("High_selection", 300)), 
                                  type   = c(rep("True",100), rep("Expected",100), rep("Median",100), 
                                             rep("True",100), rep("Expected",100), rep("Median",100),
                                             rep("True",100), rep("Expected",100), rep("Median",100)),
                                  values = c(log_neutralpods_n, posterior_neutralpods_logn$expectation, posterior_neutralpods_logn$med,
                                             log_dnlselpods_n, posterior_dnlselpods_logn$expectation, posterior_dnlselpods_logn$med,
                                             log_dnhselpods_n, posterior_dnhselpods_logn$expectation, posterior_dnhselpods_logn$med))

boxplot(pred_logn$values ~ pred_logn$type * pred_logn$model,
        lwd = 2,
        ylim = c(1,3),
        ylab = expression(log[10]*("N"[c])),
        xlab = "",
        xaxt = "n",
        cex.axis = 1.2,
        cex.lab  = 1.2,
        col = c(rep("#a50f15",3),rep("#08519c",3),rep("#969696",3)))
abline(h=log10(500), col = "green", lty = 3)
legend("bottomright", legend = c("Neutral", "1.0% Selection", "2.5% Selection"), col = c("#969696", "#08519c", "#a50f15"), pch = 15, bty = "n", cex = 1.2)
axis(side = 1, at = 1:9, labels = FALSE)
text(x = 1:9, 
     y = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), 
     srt=90, 
     adj=1, 
     xpd=TRUE, 
     labels = paste(rep(c("Expected", "Median", "True"),3)),
     cex = 1.2)


## Multiple histogram-like plots
#par(mfrow=c(1,3))
# Neutral
#hist(x      = logn,
#     breaks = seq(0,3,0.05),
#     col    = "grey",
#     freq   = FALSE,
#     ylim   = c(0,6),
#     main   = "Neutral",
#     xlab   = expression(log[10](N[c])),
#     ylab   = "Probability density")
#wtd.hist(x      = logn,
#         breaks = seq(0,3,0.05),
#         col    = adjustcolor( "red", alpha.f = 0.2), freq=FALSE, add=TRUE,
#         weight = posterior_neutralpod_n$weights)
#
# 1.0% Selection
#hist(x      = logn,
#     breaks = seq(0,3,0.05),
#     col    = "grey",
#     freq   = FALSE,
#     ylim   = c(0,6),
#     main   = "1.0% of Selection",
#     xlab   = expression(log[10](N[c])),
#     ylab   = "Probability density")
#wtd.hist(x      = logn,
#         breaks = seq(0,3,0.05),
#         col    = adjustcolor( "red", alpha.f = 0.2), freq=FALSE, add=TRUE,
#         weight = posterior_dnlselpod_n$weights)
#
# 2.5% Selection
#hist(x      = logn,
#     breaks = seq(0,3,0.05),
#    col    = "grey",
#     freq   = FALSE,
#     ylim   = c(0,6),
#     main   = "2.5% of Selection",
#     xlab   = expression(log[10](N[c])),
#     ylab   = "Probability density")
#wtd.hist(x      = logn,
#         breaks = seq(0,3,0.05),
#         col    = adjustcolor( "red", alpha.f = 0.2), freq=FALSE, add=TRUE,
#         weight = posterior_dnhselpod_n$weights)
#
#
# Multiple density plot
#par(mfrow=c(1,3))
## Neutral
#densityPlot(object    = reg_logn,
#            obs       = neutralpod_sumstats,
#            training  = data.frame(logn, global_sumstats),
#            main      = "Neutral",
#            xlab      = expression(log[10](N[c])),
#            ylab      = "Probability density", 
#            xlim      = c(0,3),
#            ylim      = c(0,3),
#            paral     = T, 
#            ncores    = 28)
#abline(v=c(logneutralpod_n,
#           posterior_neutralpod_n$expectation,
#           posterior_neutralpod_n$quantiles[1], 
#           posterior_neutralpod_n$quantiles[3]),
#       col="red",
#       lty=c(1, rep(2,3)))
#
## 1.0% Selection
#densityPlot(object    = reg_logn,
#           obs       = dnlselpod_sumstats,
#           training  = data.frame(logn, global_sumstats),
#           main      = "1.0% of Selection",
#           xlab      = expression(log[10](N[c])),
#           ylab      = "Probability density", 
#            xlim      = c(0,3),
#            ylim      = c(0,3),
#            paral     = T, 
#            ncores    = 28)
#abline(v=c(logdnlselpod_n,
#           posterior_dnlselpod_n$expectation,
#           posterior_dnlselpod_n$quantiles[1], 
#           posterior_dnlselpod_n$quantiles[3]),
#       col="red",
#       lty=c(1, rep(2,3)))
#
## 2.5% Selection
#densityPlot(object    = reg_logn,
#            obs       = dnhselpod_sumstats,
#            training  = data.frame(logn, global_sumstats),
#            main      = "2.5% of Selection",
#            xlab      = expression(log[10](N[c])),
#           ylab      = "Probability density", 
#            xlim      = c(0,3),
#            ylim      = c(0,3),
#            paral     = T, 
#            ncores    = 28)
#abline(v=c(logdnhselpod_n,
#           posterior_dnhselpod_n$expectation,
#           posterior_dnhselpod_n$quantiles[1], 
#           posterior_dnhselpod_n$quantiles[3]),
#       col="red",
#       lty=c(1, rep(2,3)))
```
