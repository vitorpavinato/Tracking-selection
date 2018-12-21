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

#### Load packages

``` r
library(DataExplorer)
library(dplyr)
library(abcrf)
library(weights)
library(grDevices)
```

#### Upload the reference table files

``` r
# Super-batch 1
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/Tracking-selection-1.0/results/pooled_reftable_1.RData")
pooled_raw_reftable_1 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

``` r
# Super-batch 2
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/Tracking-selection-1.1/results/pooled_reftable_2.RData")
pooled_raw_reftable_2 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

``` r
# Super-batch 3
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/Tracking-selection-1.2/results/pooled_reftable_3.RData")
pooled_raw_reftable_3 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

``` r
# Super-batch 4
load("/home/pavinato/My_repositories/Tracking-selection/results/lastPipelineModels/Tracking-selection-1.3/results/pooled_reftable_4.RData")
pooled_raw_reftable_4 <- pooled_raw_reftable
rm(pooled_raw_reftable)
```

#### Pool together the reftables from super-batches

``` r
pooled_raw_reftable_model_1 <- rbind(pooled_raw_reftable_1,pooled_raw_reftable_2,pooled_raw_reftable_3,pooled_raw_reftable_4)
```

#### Prepare the reftable

``` r
# Check missing data
plot_missing(pooled_raw_reftable_model_1) 
missing_count <- sapply(pooled_raw_reftable_model_1, function(x) sum(is.na(x)))

# produce a data.frame with only summary statistics that contain missing values higher than a threshold
ss2remove <- pooled_raw_reftable_model_1[, which(missing_count > 2000)]

# use this data.frame to remove the summary statistics in the reftable 
pooled_reftable_model_1 <- pooled_raw_reftable_model_1[, !pooled_raw_reftable_model_1 %in% ss2remove, drop=FALSE]

# check it!
plot_missing(pooled_reftable_model_1) 

# remove remaining rows with missing data
pooled_reftable_model_1 <- pooled_reftable_model_1[complete.cases(pooled_reftable_model_1), ]

# check it again
plot_missing(pooled_reftable_model_1)
missing_count <- sapply(pooled_reftable_model_1, function(x) sum(is.na(x)))

# Check != numeric values
infinite_count <- sapply(pooled_reftable_model_1, function(x) sum(is.infinite(x)))

# Prepare the global summary statistics table by removing columns that do not contain summary stats
global_sumstats <- pooled_reftable_model_1[, -c(1:166)]
```

### ABC-RF Regression for Selection: Growing Trees

#### Genetic Load

``` r
averageGenLoad <- pooled_reftable_model_1[, "averageGeneticLoad"]
hist(averageGenLoad, freq = TRUE, xlab = "Averaged Genetic Load", main = "", col = "#bdbdbd")

logaverageGenload <- -log10(1-averageGenLoad)
hist(logaverageGenload, freq = TRUE, xlab = "Averaged Genetic Load", main = "", col = "#bdbdbd")

# abf-rf regression
reg_averageGenLoad <- regAbcrf(formula = logaverageGenload~.,
                               data = data.frame(logaverageGenload, global_sumstats),
                               ntree = 1000,
                               paral = T,
                               ncores = 28)

# variable Importance plot
plot(x = reg_averageGenLoad, n.var = ncol(global_sumstats))

head(sort(reg_averageGenLoad$model.rf$variable.importance, decreasing = T), n=10)

# prediction error
reg_averageGenLoad$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_averageGenLoad,
             training = data.frame(logaverageGenload, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### FDR 5% of the FST hypothesis of Ns&gt;1

``` r
fstfdrNs05 <- pooled_reftable_model_1[, "FSTfdrNS05"]
hist(fstfdrNs05, freq = TRUE, xlab = "5% FST FDR Ns>1", main = "", col = "#bdbdbd")

# abf-rf regression
reg_fstfdrNs05 <- regAbcrf(formula = fstfdrNs05~.,
                              data = data.frame(fstfdrNs05, global_sumstats),
                             ntree = 1000,
                             paral = T,
                            ncores = 28)

# variable Importance plot
plot(x = reg_fstfdrNs05, n.var = ncol(global_sumstats))

head(sort(reg_fstfdrNs05$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_fstfdrNs05$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_fstfdrNs05,
             training = data.frame(fstfdrNs05, global_sumstats),
             paral    = T,
             ncores   = 28)
```

### ABC-RF Regression for Demography: Growing Trees

#### Genome-wide Inbreeding Effective Population Size - IBDNeGWtotal

``` r
ibdNeGWtotal <- pooled_reftable_model_1[, "IBDNeGWtotal"]

# replace negative values with 10*max value
ibdNeGWtotal <- ifelse(ibdNeGWtotal < 0, 10*max(ibdNeGWtotal), ibdNeGWtotal)
ibdNeGWtotal <- ifelse(ibdNeGWtotal == 0, 1, ibdNeGWtotal)

logibdNeGWtotal <- log10(ibdNeGWtotal)

hist(logibdNeGWtotal, freq = TRUE, xlab = "Genome-wide IBD Ne", main = "", col = "#bdbdbd")

reg_logibdNeGWtotal <- regAbcrf(formula = logibdNeGWtotal~.,
                           data = data.frame(logibdNeGWtotal, global_sumstats),
                           ntree = 1000,
                           paral = T,
                           ncores = 28)

# variable Importance plot
plot(x = reg_logibdNeGWtotal, n.var = ncol(global_sumstats))

head(sort(reg_logibdNeGWtotal$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logibdNeGWtotal$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logibdNeGWtotal,
             training = data.frame(logibdNeGWtotal, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### EXTRA Chromosome Inbreeding Effective Population Size - IBDNeChrtotal

``` r
ibdNeChrtotal <- pooled_reftable_model_1[, "IBDNeChrtotal"]

# replace negative values with 10*max value
ibdNeChrtotal <- ifelse(ibdNeChrtotal < 0, 10*max(ibdNeChrtotal), ibdNeChrtotal)
ibdNeChrtotal <- ifelse(ibdNeChrtotal == 0, 1, ibdNeChrtotal)

logibdNeChrtotal <- log10(ibdNeChrtotal)

hist(logibdNeChrtotal, freq = TRUE, xlab = "EXTRA Chromosome IBD Ne", main = "", col = "#bdbdbd")

reg_logibdNeChrtotal <- regAbcrf(formula = logibdNeChrtotal~.,
                                data = data.frame(logibdNeChrtotal, global_sumstats),
                                ntree = 1000,
                                paral = T,
                                ncores = 28)

# variable Importance plot
plot(x = reg_logibdNeChrtotal, n.var = ncol(global_sumstats))

head(sort(reg_logibdNeChrtotal$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logibdNeChrtotal$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logibdNeChrtotal,
             training = data.frame(logibdNeChrtotal, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### Genome-wide Variance Effective Population Size - VARNeGWtotal

``` r
varNeGWtotal <- pooled_reftable_model_1[, "VARNeGWtotal"]

# replace negative values with 10*max value
varNeGWtotal <- ifelse(varNeGWtotal < 0, 10*max(varNeGWtotal), varNeGWtotal)
varNeGWtotal <- ifelse(varNeGWtotal == 0, 1, varNeGWtotal)
logvarNeGWtotal <- log10(varNeGWtotal)

hist(logvarNeGWtotal, freq = TRUE, xlab = "Genome-wide variance Ne", main = "", col = "#bdbdbd")

reg_logvarNeGWtotal <- regAbcrf(formula = logvarNeGWtotal~.,
                                data = data.frame(logvarNeGWtotal, global_sumstats),
                                ntree = 1000,
                                paral = T,
                                ncores = 28)

# variable Importance plot
plot(x  = reg_logvarNeGWtotal, n.var = ncol(global_sumstats))

head(sort(reg_logvarNeGWtotal$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logvarNeGWtotal$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logvarNeGWtotal,
             training = data.frame(logvarNeGWtotal, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### EXTRA Chromosome Variance Effective Population Size - VARNeChrtotal

``` r
varNeChrtotal <- pooled_reftable_model_1[, "VARNeChrtotal"]

# replace negative values with 10*max value
varNeChrtotal <- ifelse(varNeChrtotal < 0, 10*max(varNeChrtotal), varNeChrtotal)
varNeChrtotal <- ifelse(varNeChrtotal == 0, 1, varNeChrtotal)

logvarNeChrtotal <- log10(varNeChrtotal)

hist(logvarNeChrtotal, freq = TRUE, xlab = "EXTRA Chromosome variance Ne", main = "", col = "#bdbdbd")

reg_logvarNeChrtotal <- regAbcrf(formula = logvarNeChrtotal~.,
                                 data = data.frame(logvarNeChrtotal, global_sumstats),
                                 ntree = 1000,
                                 paral = T,
                                 ncores = 28)

# Variable Importance plot
plot(x = reg_logvarNeChrtotal, n.var = ncol(global_sumstats))

head(sort(reg_logvarNeChrtotal$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logvarNeChrtotal$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logvarNeChrtotal,
             training = data.frame(logvarNeChrtotal, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### Pedigree Effective Population Size - PedigreeNetotal

``` r
pedigreeNetotal <- pooled_reftable_model_1[, "PedigreeNetotal"]

logpedigreeNetotal <- log10(pedigreeNetotal)

hist(logpedigreeNetotal, freq = TRUE, xlab = "Pedigree Effective Population Size", main = "", col = "#bdbdbd")

reg_logpedigreeNetotal <- regAbcrf(formula = logpedigreeNetotal~.,
                                 data = data.frame(logpedigreeNetotal, global_sumstats),
                                 ntree = 1000,
                                 paral = T,
                                 ncores = 28)

# Variable Importance plot
plot(x = reg_logpedigreeNetotal, n.var = ncol(global_sumstats))

head(sort(reg_logpedigreeNetotal$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logpedigreeNetotal$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logpedigreeNetotal,
             training = data.frame(logpedigreeNetotal, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### Population Census Size for the Sampling Phase - N

``` r
n <- pooled_reftable_model_1[, "N"]

logn <- log10(n)

hist(logn, freq = TRUE, xlab = "Population Census Size for the Sampling Phase", main = "", col = "#bdbdbd")

reg_logn <- regAbcrf(formula = logn~.,
                        data = data.frame(logn, global_sumstats),
                       ntree = 1000,
                       paral = T,
                      ncores = 28)

# Variable Importance plot
plot(x = reg_logn, n.var = ncol(global_sumstats))

head(sort(reg_logn$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logn$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logn,
             training = data.frame(logn, global_sumstats),
             paral    = T,
             ncores   = 28)
```

#### Population Census Size for the Equilibrium Phase - Neq

``` r
neq <- pooled_reftable_model_1[, "Neq"]

logneq <- log10(neq)

hist(logneq, freq = TRUE, xlab = "Population Census Size for the Equilibrium Phase", main = "", col = "#bdbdbd")

reg_logneq <- regAbcrf(formula = logneq~.,
                        data = data.frame(logneq, global_sumstats),
                       ntree = 1000,
                       paral = T,
                      ncores = 28)

# Variable Importance plot
plot(x = reg_logneq, n.var = ncol(global_sumstats))

head(sort(reg_logneq$model.rf$variable.importance, decreasing = T), n=20)

# prediction error
reg_logneq$model.rf$prediction.error

# error rate plot
err.regAbcrf(object   = reg_logneq,
             training = data.frame(logneq, global_sumstats),
             paral    = T,
             ncores   = 28)
```

### ABC-RF Regression for Selection: Prediction

#### Load PODs

``` r
load("results/lastPipelineModels/PODs/Selection/reference_table_selection.RData")

selectedpod_sumstats <- raw_reftable[, -c(1:217)]
selectedpod_sumstats <- selectedpod_sumstats[!sapply(selectedpod_sumstats, function(x) all(is.na(x)))]
```

#### Genetic Load

``` r
selectedpod_averageGenLoad <- raw_reftable[, "averageGeneticLoad"]

# log transformation
logselectedpodaverageGenload <- -log10(1-selectedpod_averageGenLoad)

posterior_averageGenLoad <- predict(   object  = reg_averageGenLoad,
                                       obs     = selectedpod_sumstats,
                                     training  = data.frame(logaverageGenload, global_sumstats),
                                     quantiles = c(0.025,0.5,0.975),
                                     paral     = T,
                                     ncores    = 28,
                                    rf.weights = T)
(posterior_averageGenLoad)

# histrogram-like plot
hist(x      = logaverageGenload,
     breaks = seq(0,3,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,6),
     main   = "",
     xlab   = "Genetic load",
     ylab   = "Probability density")
wtd.hist(x      = logaverageGenload,
         breaks = seq(0,3,0.05),
         col    = adjustcolor( "red", alpha.f = 0.2), freq=FALSE, add=TRUE,
         weight = posterior_averageGenLoad$weights)

# density plot
densityPlot(object    = reg_averageGenLoad,
            obs       = selectedpod_sumstats,
            training  = data.frame(logaverageGenload, global_sumstats),
            main      = "Genetic load",
            xlab      = expression(-log[10]*(1 - "genetic load")),
            ylab      = "Probability density", 
            paral     = T, 
            ncores    = 28)
abline(v=c(logselectedpodaverageGenload,
           posterior_averageGenLoad$expectation,
           posterior_averageGenLoad$quantiles[1], 
           posterior_averageGenLoad$quantiles[3]),
       col="red",
       lty=c(1, rep(2,3)))
```

### ABC-RF Regression for Demography: Prediction

#### Pedigree Effective Population Size - PedigreeNetotal

``` r
selectedpod_pedigreeNetotal <- raw_reftable[, "PedigreeNetotal"]

logselectedpodpedigreeNetotal <- log10(selectedpod_pedigreeNetotal)

posterior_pedigreeNetotal <- predict(  object  = reg_logpedigreeNetotal,
                                       obs     = selectedpod_sumstats,
                                     training  = data.frame(logpedigreeNetotal, global_sumstats),
                                     quantiles = c(0.025,0.5,0.975),
                                     paral     = T,
                                     ncores    = 28,
                                    rf.weights = T)
(posterior_pedigreeNetotal)

# histrogram-like plot
hist(x      = logpedigreeNetotal,
     breaks = seq(-2,4,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = "Effective population size",
     ylab   = "Probability density")
wtd.hist(x      = logpedigreeNetotal,
         breaks = seq(-2,4,0.05),
         col    = adjustcolor( "red", alpha.f = 0.2), freq=FALSE, add=TRUE,
         weight = posterior_pedigreeNetotal$weights)

# density plot
densityPlot(object    = reg_logpedigreeNetotal,
            obs       = selectedpod_sumstats,
            training  = data.frame(logpedigreeNetotal, global_sumstats),
            main      = "Effective Population Size",
            xlab      = expression(log[10]*"Ne"),
            ylab      = "Probability density", 
            xlim =c(0,3),
            paral     = T, 
            ncores    = 28)
abline(v=c(logselectedpodpedigreeNetotal,
           posterior_pedigreeNetotal$expectation,
           posterior_pedigreeNetotal$quantiles[1], 
           posterior_pedigreeNetotal$quantiles[3]),
       col="red",
       lty=c(1, rep(2,3)))
```

#### Population Census Size for the Sampling Phase - N

``` r
selectedpod_n <- raw_reftable[, "N"]

logselectedpodn <- log10(selectedpod_global_n)

posterior_logn <- predict(    object = reg_logn,
                             obs     = selectedpod_sumstats,
                           training  = data.frame(logn, global_sumstats),
                           quantiles = c(0.025,0.5,0.975),
                           paral     = T,
                           ncores    = 28,
                          rf.weights = T)
(posterior_logn)

# histrogram-like plot
hist(x      = logn,
     breaks = seq(0,3,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = "Census size N",
     ylab   = "Probability density")
wtd.hist(x      = logn,
         breaks = seq(0,3,0.05),
         col    = adjustcolor( "red", alpha.f = 0.2), freq=FALSE, add=TRUE,
         weight = posterior_logn$weights)

# density plot
densityPlot(object    = reg_logn,
            obs       = selectedpod_global_sumstats,
            training  = data.frame(logn, global_sumstats),
            main      = "Population Census Size",
            xlab      = expression(log[10]*"N"),
            ylab      = "Probability density", 
            paral     = T, 
            ncores    = 28)
abline(v=c(logselectedpodn,
           posterior_logn$expectation,
           posterior_logn$quantiles[1], 
           posterior_logn$quantiles[3]),
       col="red",
       lty=c(1, rep(2,3)))
```
