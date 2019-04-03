#### Script to grow trees

# clean memory
rm(list=ls())
ls()

# load packages
library(abcrf)

### Reference Table
# load reftable
load(file=paste0("/nfs/work/select/rf_analysis/pooled_reftable_total",".RData"))

# load global summary stats
load(file=paste0("/nfs/work/select/rf_analysis/global_sumstats",".RData"))

### Random PODs
# load reftable
load(file=paste0("/nfs/work/select/rf_analysis/pooled_reftable_pods",".RData"))

# load global summary stats
load(file=paste0("/nfs/work/select/rf_analysis/global_sumstats_pods",".RData"))

### Fixed PODs
# load reftable
load(file=paste0("/nfs/work/select/rf_analysis/pooled_reftable_fixedpods",".RData"))

# load global summary stats
load(file=paste0("/nfs/work/select/rf_analysis/global_sumstats_fixedpods",".RData"))

### RF Prediction

## ThetaS_2 - Proportion of selecte mutations POPPrMSel
load(file=paste0("/nfs/work/select/rf_analysis/reg_logthetaS_2",".RData"))

prPOPPrMSel <- pooled_reftable_total[, "POPPrMSel"]
zero_prmsel <- which(prPOPPrMSel == 0)
prPOPPrMSel_2 <- prPOPPrMSel[-zero_prmsel]
global_sumstats_popmsel <- global_sumstats[-zero_prmsel,]
logthetaS_2 <- log10(4 * pooled_reftable_total[-zero_prmsel, "meanNe2"] * (pooled_reftable_total[-zero_prmsel, "mu"] * prPOPPrMSel_2) * genomeS)

# RANDOM PODs
posterior_radnompods_logthetaS_2 <- predict(object     = reg_logthetaS_2,
                                            obs        = global_sumstats_pods,
                                            training   = data.frame(logthetaS_2, global_sumstats_popmsel),
                                            quantiles  = c(0.025,0.5,0.975),
                                            paral      = T,
                                            ncores     = 20,
                                            rf.weights = T)
save(posterior_radnompods_logthetaS_2, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_logthetaS_2",".RData"))

# FIXED PODs
posterior_fixedpods_logthetaS_2 <- predict(object     = reg_logthetaS_2,
                                           obs        = global_sumstats_fixedpods,
                                           training   = data.frame(logthetaS_2, global_sumstats_popmsel),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 20,
                                           rf.weights = T)

save(posterior_fixedpods_logthetaS_2, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logthetaS_2",".RData"))

## ThetaS_3 - Proportion of strongly selecte mutations POPPrStrongMSel
load(file=paste0("/nfs/work/select/rf_analysis/reg_logthetaS_3",".RData"))

prPOPStrongMSel <- pooled_reftable_total[, "POPPrStrongMSel"]
neutralsims <- which(prPOPStrongMSel == 0)
prPOPStrongMSel_2 <- prPOPStrongMSel[-neutralsims]
global_sumstats_popstrongmsel <- global_sumstats[-neutralsims,]
logthetaS_3 <- log10(4 * pooled_reftable_total[-neutralsims, "meanNe2"] * (pooled_reftable_total[-neutralsims, "mu"] * prPOPStrongMSel_2) * genomeS)

# RANDOM PODs
posterior_radnompods_logthetaS_3 <- predict(object     = reg_logthetaS_3,
                                            obs        = global_sumstats_pods,
                                            training   = data.frame(logthetaS_3, global_sumstats_popstrongmsel),
                                            quantiles  = c(0.025,0.5,0.975),
                                            paral      = T,
                                            ncores     = 20,
                                            rf.weights = T)
save(posterior_radnompods_logthetaS_3, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_logthetaS_3",".RData"))

# FIXED PODs
posterior_fixedpods_logthetaS_3 <- predict(object     = reg_logthetaS_3,
                                           obs        = global_sumstats_fixedpods,
                                           training   = data.frame(logthetaS_3, global_sumstats_popstrongmsel),
                                           quantiles  = c(0.025,0.5,0.975),
                                           paral      = T,
                                           ncores     = 20,
                                           rf.weights = T)

save(posterior_fixedpods_logthetaS_3, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logthetaS_3",".RData"))

## logpopstrongmsel
load(file=paste0("/nfs/work/select/rf_analysis/reg_logpopstrongmsel",".RData"))

logpopstrongmsel <- log10(prPOPStrongMSel_2)

# RANDOM PODs
posterior_radnompods_logpopstrongmsel <- predict(object     = reg_logpopstrongmsel,
                                                 obs        = global_sumstats_pods,
                                                 training   = data.frame(logpopstrongmsel, global_sumstats_popstrongmsel),
                                                 quantiles  = c(0.025,0.5,0.975),
                                                 paral      = T,
                                                 ncores     = 20,
                                                 rf.weights = T)

save(posterior_radnompods_logpopstrongmsel, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_logpopstrongmsel",".RData"))

# FIXED PODs
posterior_fixedpods_logpopstrongmsel <- predict(object     = reg_logpopstrongmsel,
                                                obs        = global_sumstats_fixedpods,
                                                training   = data.frame(logpopstrongmsel, global_sumstats_popstrongmsel),
                                                quantiles  = c(0.025,0.5,0.975),
                                                paral      = T,
                                                ncores     = 20,
                                                rf.weights = T)

save(posterior_fixedpods_logpopstrongmsel, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logpopstrongmsel",".RData"))

## logitpopstrongmsel
load(file=paste0("/nfs/work/select/rf_analysis/reg_logitpopstrongmsel",".RData"))

logitpopstrongmsel <- log(prPOPStrongMSel_2/(1-prPOPStrongMSel_2))

# RANDOM PODs
posterior_radnompods_logitpopstrongmsel <- predict(object     = reg_logitpopstrongmsel,
                                                   obs        = global_sumstats_pods,
                                                   training   = data.frame(logitpopstrongmsel,global_sumstats_popstrongmsel),
                                                   quantiles  = c(0.025,0.5,0.975),
                                                   paral      = T,
                                                   ncores     = 28,
                                                   rf.weights = T)

save(posterior_radnompods_logitpopstrongmsel, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_logitpopstrongmsel",".RData"))

# FIXED PODs
posterior_fixedpods_logitpopstrongmsel <- predict(object     = reg_logitpopstrongmsel,
                                                   obs        = global_sumstats_fixedpods,
                                                   training   = data.frame(logitpopstrongmsel,global_sumstats_popstrongmsel),
                                                   quantiles  = c(0.025,0.5,0.975),
                                                   paral      = T,
                                                   ncores     = 28,
                                                   rf.weights = T)

save(posterior_fixedpods_logitpopstrongmsel, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logitpopstrongmsel",".RData"))


## DFE Gamma mean
load(file=paste0("/nfs/work/select/rf_analysis/reg_loggammamean",".RData"))

loggammamean <- log10(pooled_reftable_total[, "gammaMean"])

# RANDOM PODs
posterior_radnompods_loggammamean <- predict(object     = reg_loggammamean,
                                             obs        = global_sumstats_pods,
                                             training   = data.frame(loggammamean,global_sumstats),
                                             quantiles  = c(0.025,0.5,0.975),
                                             paral      = T,
                                             ncores     = 28,
                                             rf.weights = T)

save(posterior_radnompods_loggammamean, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_loggammamean",".RData"))

# FIXED PODs
posterior_fixedpods_loggammamean <- predict(object     = reg_loggammamean,
                                            obs        = global_sumstats_fixedpods,
                                            training   = data.frame(loggammamean,global_sumstats),
                                            quantiles  = c(0.025,0.5,0.975),
                                            paral      = T,
                                            ncores     = 28,
                                            rf.weights = T)

save(posterior_fixedpods_loggammamean, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_loggammamean",".RData"))

## DFE Gamma - POPSelMean
load(file=paste0("/nfs/work/select/rf_analysis/reg_meanPOPSel_1",".RData"))
meanPOPSel <- pooled_reftable_total[, "POPSelMean"]
meanPOPSel_2 <- meanPOPSel[-zero_prmsel]

# RANDOM PODs
posterior_radnompods_meanPOPSel_1 <- predict(object     = reg_meanPOPSel_1,
                                             obs        = global_sumstats_pods,
                                             training   = data.frame(meanPOPSel_2,global_sumstats_popmsel),
                                             quantiles  = c(0.025,0.5,0.975),
                                             paral      = T,
                                             ncores     = 28,
                                             rf.weights = T)

save(posterior_radnompods_meanPOPSel_1, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_meanPOPSel_1",".RData"))

# FIXED PODs
posterior_fixedpods_meanPOPSel_1 <- predict(object     = reg_meanPOPSel_1,
                                            obs        = global_sumstats_fixedpods,
                                            training   = data.frame(meanPOPSel_2,global_sumstats_popmsel),
                                            quantiles  = c(0.025,0.5,0.975),
                                            paral      = T,
                                            ncores     = 28,
                                            rf.weights = T)

save(posterior_fixedpods_meanPOPSel_1, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_meanPOPSel_1",".RData"))

## DFE Gamma - POPStrongSelMean
load(file=paste0("/nfs/work/select/rf_analysis/reg_meanPOPSel_2",".RData"))

meanstrongPOPSel <- pooled_reftable_total[, "POPStrongSelMean"]
meanstrongPOPSel_2 <- meanstrongPOPSel[-neutralsims]

# RANDOM PODs
posterior_radnompods_meanPOPSel_2 <- predict(object     = reg_meanPOPSel_2,
                                             obs        = global_sumstats_pods,
                                             training   = data.frame(meanstrongPOPSel_2,global_sumstats_popstrongmsel),
                                             quantiles  = c(0.025,0.5,0.975),
                                             paral      = T,
                                             ncores     = 28,
                                             rf.weights = T)

save(posterior_radnompods_meanPOPSel_2, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_meanPOPSel_2",".RData"))

# FIXED PODs
posterior_fixedpods_meanPOPSel_2 <- predict(object     = reg_meanPOPSel_2,
                                            obs        = global_sumstats_fixedpods,
                                            training   = data.frame(meanstrongPOPSel_2,global_sumstats_popstrongmsel),
                                            quantiles  = c(0.025,0.5,0.975),
                                            paral      = T,
                                            ncores     = 28,
                                            rf.weights = T)

save(posterior_fixedpods_meanPOPSel_2, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_meanPOPSel_2",".RData"))

## DFE Gamma - logmeanpopsel_1
load(file=paste0("/nfs/work/select/rf_analysis/reg_logmeanpopsel_1",".RData"))

logmeanpopsel_1 <- log10(meanPOPSel_2)

# RANDOM PODs
posterior_radnompods_logmeanpopsel_1 <- predict(object     = reg_logmeanpopsel_1,
                                                obs        = global_sumstats_pods,
                                                training   = data.frame(logmeanpopsel_1,global_sumstats_popmsel),
                                                quantiles  = c(0.025,0.5,0.975),
                                                paral      = T,
                                                ncores     = 28,
                                                rf.weights = T)

save(posterior_radnompods_logmeanpopsel_1, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_logmeanpopsel_1",".RData"))

# FIXED PODs
posterior_fixedpods_logmeanpopsel_1 <- predict(object     = reg_logmeanpopsel_1,
                                               obs        = global_sumstats_fixedpods,
                                               training   = data.frame(logmeanpopsel_1,global_sumstats_popmsel),
                                               quantiles  = c(0.025,0.5,0.975),
                                               paral      = T,
                                               ncores     = 28,
                                               rf.weights = T)

save(posterior_fixedpods_logmeanpopsel_1, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logmeanpopsel_1",".RData"))

## DFE Gamma - logmeanpopsel_2
load(file=paste0("/nfs/work/select/rf_analysis/reg_logmeanpopsel_2",".RData"))

logmeanpopsel_2 <- log10(meanstrongPOPSel_2)

# RANDOM PODs
posterior_radnompods_logmeanpopsel_2 <- predict(object     = reg_logmeanpopsel_2,
                                                obs        = global_sumstats_pods,
                                                training   = data.frame(logmeanpopsel_2,global_sumstats_popstrongmsel),
                                                quantiles  = c(0.025,0.5,0.975),
                                                paral      = T,
                                                ncores     = 28,
                                                rf.weights = T)

save(posterior_radnompods_logmeanpopsel_2, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_radnompods_logmeanpopsel_2",".RData"))

# FIXED PODs
posterior_fixedpods_logmeanpopsel_2 <- predict(object     = reg_logmeanpopsel_2,
                                               obs        = global_sumstats_fixedpods,
                                               training   = data.frame(logmeanpopsel_2,global_sumstats_popstrongmsel),
                                               quantiles  = c(0.025,0.5,0.975),
                                               paral      = T,
                                               ncores     = 28,
                                               rf.weights = T)

save(posterior_fixedpods_logmeanpopsel_2, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logmeanpopsel_2",".RData"))

## Census size Ncs
load(file=paste0("/nfs/work/select/rf_analysis/reg_logncs",".RData"))

logncs <- log10(pooled_reftable_total[, "Ncs"])

# FIXED PODs
posterior_fixedpods_logncs <- predict(object     = reg_logncs,
                                      obs        = global_sumstats_fixedpods,
                                      training   = data.frame(logncs, global_sumstats),
                                      quantiles  = c(0.025,0.5,0.975),
                                      paral      = T,
                                      ncores     = 28,
                                      rf.weights = T)

save(posterior_fixedpods_logncs, 
     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logncs",".RData"))

## Theta1
#load(file=paste0("/nfs/work/select/rf_analysis/reg_logtheta1",".RData"))
#
#logtheta1 <- log10(4 * pooled_reftable_total[, "meanNe1"] * (pooled_reftable_total[, "mu"]) * genomeS)
#
# RANDOM PODs
#posterior_randompods_logtheta1 <- predict(object     = reg_logtheta1,
#                                          obs        = global_sumstats_pods,
#                                          training   = data.frame(logtheta1, global_sumstats),
#                                          quantiles  = c(0.025,0.5,0.975),
#                                          paral      = T,
#                                          ncores     = 28,
#                                          rf.weights = T)
#
#save(posterior_randompods_logtheta1, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_randompods_logtheta1",".RData"))
#
# FIXED PODs
#posterior_fixedpods_logtheta1 <- predict(object     = reg_logtheta1,
#                                         obs        = global_sumstats_fixedpods,
#                                         training   = data.frame(logtheta1, global_sumstats),
#                                         quantiles  = c(0.025,0.5,0.975),
#                                         paral      = T,
#                                         ncores     = 28,
#                                         rf.weights = T)
#
#save(posterior_fixedpods_logtheta1, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logtheta1",".RData"))
#
## Per site mutation rate mu
#load(file=paste0("/nfs/work/select/rf_analysis/reg_logmu",".RData"))
#
#logmu <- log10(pooled_reftable_total[, "mu"])
#
# RANDOM PODs
#posterior_randompods_logmu <- predict(object     = reg_logmu,
#                                      obs        = global_sumstats_pods,
#                                      training   = data.frame(logmu, global_sumstats),
#                                      quantiles  = c(0.025,0.5,0.975),
#                                      paral      = T,
#                                      ncores     = 28,
#                                      rf.weights = T)
#
#save(posterior_randompods_logmu, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_randompods_logmu",".RData"))
#
# FIXED PODs
#posterior_fixedpods_logmu <- predict(object     = reg_logmu,
#                                     obs        = global_sumstats_fixedpods,
#                                     training   = data.frame(logmu, global_sumstats),
#                                     quantiles  = c(0.025,0.5,0.975),
#                                     paral      = T,
#                                     ncores     = 28,
#                                     rf.weights = T)
#
#save(posterior_fixedpods_logmu, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logmu",".RData"))
#
## Ratio Effective Size / Census Size for the Sampling Phase - MeanNe2/Ncs
#load(file=paste0("/nfs/work/select/rf_analysis/reg_logmeanNe2ncs",".RData"))
#
#logmeanNe2ncs <- log10(pooled_reftable_total[, "meanNe2"]/pooled_reftable_total[, "Ncs"])
#
# RANDOM PODs
#posterior_randompods_logmeanNe2ncs <- predict(object     = reg_logmeanNe2ncs,
#                                              obs        = global_sumstats_pods,
#                                              training   = data.frame(logmu, global_sumstats),
#                                              quantiles  = c(0.025,0.5,0.975),
#                                              paral      = T,
#                                              ncores     = 28,
#                                              rf.weights = T)
#
#save(posterior_randompods_logmeanNe2ncs, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_randompods_logmeanNe2ncs",".RData"))
#
# FIXED PODs
#posterior_fixedpods_logmeanNe2ncs <- predict(object     = reg_logmeanNe2ncs,
#                                             obs        = global_sumstats_fixedpods,
#                                             training   = data.frame(logmeanNe2ncs, global_sumstats),
#                                             quantiles  = c(0.025,0.5,0.975),
#                                             paral      = T,
#                                             ncores     = 28,
#                                             rf.weights = T)
#
#save(posterior_fixedpods_logmeanNe2ncs, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logmeanNe2ncs",".RData"))
#
#
## Rho = Nr
#load(file=paste0("/nfs/work/select/rf_analysis/reg_logrho",".RData"))
#
#logrho <- log10(pooled_reftable_total[, "meanNe2"] * pooled_reftable_total[, "rr"])
#
# RANDOM PODs
#posterior_randompods_logrho <- predict(object     = reg_logrho,
#                                       obs        = global_sumstats_pods,
#                                       training   = data.frame(logrho, global_sumstats),
#                                       quantiles  = c(0.025,0.5,0.975),
#                                       paral      = T,
#                                       ncores     = 28,
#                                       rf.weights = T)
#
#save(posterior_randompods_logrho, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_randompods_logrho",".RData"))
#
# FIXED PODs
#posterior_fixedpods_logrho <- predict(object     = reg_logrho,
#                                      obs        = global_sumstats_fixedpods,
#                                      training   = data.frame(logrho, global_sumstats),
#                                      quantiles  = c(0.025,0.5,0.975),
#                                      paral      = T,
#                                      ncores     = 28,
#                                      rf.weights = T)
#
#save(posterior_fixedpods_logrho, 
#     file = paste0("/nfs/work/select/rf_analysis/posterior_fixedpods_logrho",".RData"))    
#     
gc()






