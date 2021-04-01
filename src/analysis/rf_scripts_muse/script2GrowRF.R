### Script to grow trees

# clean memory
rm(list=ls())
ls()

# load packages
library(abcrf)

# load reftable
load(file=paste0("/nfs/work/select/rf_analysis/pooled_reftable_total",".RData"))

# load global summary stats
load(file=paste0("/nfs/work/select/rf_analysis/global_sumstats",".RData"))

### RF for mean selection coefficient

## POPSelMean

## Preparation
##--------------------------------------------------

# Load POPPrMSel vector
#prPOPPrMSel <- pooled_reftable_total[, "POPPrMSel"]

# Mark simulations with POPPrMSel == 0
#zero_prmsel <- which(prPOPPrMSel == 0)

# Prepare global sumstats with same dim
#global_sumstats_popmsel <- global_sumstats[-zero_prmsel,]

# Load the vector of target parameter
#meanPOPSel <- pooled_reftable_total[, "POPSelMean"]

# Remove simulations with POPPrMSel == 0
#meanPOPSel_2 <- meanPOPSel[-zero_prmsel]
#logmeanpopsel_1 <- log10(meanPOPSel_2)

## abf-rf regression
##--------------------------------------------------

#reg_logmeanpopsel_1 <- regAbcrf(formula = logmeanpopsel_1~.,
#                              data    = data.frame(logmeanpopsel_1, global_sumstats_popmsel),
#                              ntree   = 1000,
#                              paral   = T,
#                              ncores  = 20)

#save(reg_logmeanpopsel_1, 
#     file = paste0("/nfs/work/select/rf_analysis/reg_logmeanpopsel_1",".RData"))
     
## POPStrongSelMean

## Preparation
##--------------------------------------------------

# Load PrPOPStrongMSel vector
#prPOPStrongMSel <- pooled_reftable_total[, "POPPrStrongMSel"]

# Mark simulations with PrPOPStrongMSel == 0
#neutralsims <- which(prPOPStrongMSel == 0)

# Prepare global sumstats with same dim
#global_sumstats_popstrongmsel <- global_sumstats[-neutralsims,]

# Load the vector of target parameter
#meanstrongPOPSel <- pooled_reftable_total[, "POPStrongSelMean"]

# Remove simulations with PrPOPStrongMSel == 0
#meanstrongPOPSel_2 <- meanstrongPOPSel[-neutralsims]
#logmeanpopsel_2 <- log10(meanstrongPOPSel_2)     
     
## abf-rf regression
##--------------------------------------------------

#reg_logmeanpopsel_2 <- regAbcrf(formula = logmeanpopsel_2~.,
#                              data    = data.frame(logmeanpopsel_2, global_sumstats_popstrongmsel),
#                              ntree   = 1000,
#                              paral   = T,
#                              ncores  = 20)

#save(reg_logmeanpopsel_2, 
#     file = paste0("/nfs/work/select/rf_analysis/reg_logmeanpopsel_2",".RData"))     
 

## Ns

## preparation
##--------------------------------------------------
#logNs <- log10(pooled_reftable_total[, "meanNe2"] * pooled_reftable_total[, "gammaMean"])

# abf-rf regression
##--------------------------------------------------
#reg_logNs <- regAbcrf(formula = logNs~.,
#                      data    = data.frame(logNs, global_sumstats),
#                      ntree   = 1000,
#                      paral   = T,
#                      ncores  = 28)

#save(reg_logNs, 
#     file = paste0("/nfs/work/select/rf_analysis/reg_logNs",".RData"))


## Theta 2

## preparation
##--------------------------------------------------
#genomeS = 100e+6
#logtheta2 <- log10(4 * pooled_reftable_total[, "meanNe2"] * (pooled_reftable_total[, "mu"]) * genomeS)


## abc-rf regression
##--------------------------------------------------
#reg_logtheta2 <- regAbcrf(formula = logtheta2~.,
#                          data    = data.frame(logtheta2, global_sumstats),
#                          ntree   = 1000,
#                          paral   = T,
#                          ncores  = 28)

#save(reg_logtheta2, 
#     file = paste0("/nfs/work/select/rf_analysis/reg_logtheta2",".RData")) 

gc()






