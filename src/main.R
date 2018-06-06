#!/usr/bin/env Rscript

##########################################################################################################
##                            Tracking Selection in Time-Series Data with                               ## 
##                                      Forward-time Simulations                                        ##
##########################################################################################################

## Script to perform an ABC analysis 
## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

### Remove last features
rm(list=ls())
ls()

################
## R PACKAGES ##
################
#install.packages(c("KScorrect", "moments", "foreach", "parallel", "doParallel"), dependencies = T)

###########################################
##           GLOBAL SETTINGS             ##
###########################################

nsim               <- 3                              # number of simulations 1000
model              <- "src/models/model_2_v2.1.slim"
slim_output_folder <- "results/slim_output/"
egglib_input       <- "results/egglib_input/"
egglib_output      <- "results/egglib_output/"
reftable_file      <- "results/reference_table"       # reference table file name
seed               <- 1234
set.seed(seed,"Mersenne-Twister")
parallel_sims    <- TRUE
num_of_threads   <- 3
remove_files     <- FALSE

##########################################
##       EGGLIB SUMMSTAT SETTINGS       ##
##########################################

python_path     <- "/home/pavinato/py-egglib-3.0.0b19/bin/python"
#python_path     <- "/home/pavinato/py-egglib-3.0.0b21/bin/python"

egglib_summstat <- "bin/summstats_0.0.py" # this version works with egglib-3.0.0b19
#egglib_summstat <- "bin/summstats_1.0.py" # this version works with egglib-3.0.0b21

wss_wspan = 100
sfs_bins = 5

############################################
##             SAMPLED VALUES             ##
## DEFINE PRIORS (upper and lower limits) ##
############################################

# THETA:
# We are going to sample thetas from a log uniform distribution; however SLiM is set up on population size Ne;
# Given the mutation rate, theta is calculated as THETA = Ne*mutation_rate
# For now, only sampling theta and calculating Ne with known mutation rate mu

mu_rate = 0 # 0 = "FIXED"; 1 = "RANDOM" sample from prior
mu_min  = 1e-8
mu_max  = 1e-5

ne0_min = 1
ne0_max = 1000

ne1_min = 1
ne1_max = 1000

pge2_min = 0
pge2_max = 0.5

mpb_min = 0
mpb_max = 1

gamma_min = 0.00001 # this defines a lower and an upper limits of a uniform distribution where gamma MEAN and SHAPE (K) values;
gamma_max = 0.1

rr_rate = 0 # 0 = "FIXED"; 1 = "RANDOM" sample from prior
rr_min  = 4.2 * 1e-8
rr_max  = 4.2 * 1e-5

############################################
##              FIXED VALUES              ##
############################################

SS1 = 80
SS2 = 115

ts2 = 8

###########################################
##    LOAD REQUIRED FUNCTIONS/PACKAGES   ##
###########################################

# load required libraries
library(KScorrect) # for the log uniform distribution
#library(plyr)     # for sample markers in chromossome 2 (sample_n)
#library(dplyr)
#library(mon)
library(moments)
if(parallel_sims){
  library(foreach)    #
  library(parallel)   # -> fo    r simulating in parallel
  library(doParallel) #
}

# load other functions (from file distributed together with present file)
source("src/fun.R")

###########################################
##              SIMULATION               ##
###########################################

# Parallelization with doParallel and foreach
system.time(if(parallel_sims){
  cl <- makeCluster(num_of_threads)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(KScorrect))
  ref_table <- foreach(sim=seq_len(nsim),.combine=rbind) %dopar% {
                                                                        #library(plyr)
                                                                        #library(dplyr)
                                                                        library(moments)
                                                                        do_sim(sim, nsim, model,
                                                                               mu_rate, mu_min, mu_max, 
                                                                               ne0_min, ne0_max, ne1_min, ne1_max, pge2_min, pge2_max, mpb_min, mpb_max, 
                                                                               gamma_min, gamma_max, rr_rate, rr_min, rr_max, SS1, SS2, ts2,
                                                                               slim_output_folder, egglib_input, save_extra_data, extra_output,  
                                                                               python_path, egglib_summstat, wss_wspan, sfs_bins,
                                                                               egglib_output, 
                                                                               remove_files
                                                                               )
    }

  stopCluster(cl)
  
}else{
  ref_table <- vector("list", nsim)
  for(sim in 1:nsim){
    ref_table[[sim]] <- do_sim(sim, nsim, model,
                               mu_rate, mu_min, mu_max, 
                               ne0_min, ne0_max, ne1_min, ne1_max, pge2_min, pge2_max, mpb_min, mpb_max, 
                               gamma_min, gamma_max, rr_rate, rr_min, rr_max, SS1, SS2, ts2,
                               slim_output_folder, egglib_input, save_extra_data, extra_output,  
                               python_path, egglib_summstat, wss_wspan, sfs_bins,
                               egglib_output, 
                               remove_files
                               )
  }
  ref_table <- do.call(rbind, ref_table)
})
gc()

write.table(ref_table,
            file      = paste0(reftable_file,".txt"),
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

ref_table <- read.table(paste0(reftable_file,".txt"),header=T)
#dim(ref_table)
#head(ref_table)
#summary(ref_table)
save(ref_table,file=paste0(reftable_file,".RData"))

cat("\n Simulations finished\n\n")
