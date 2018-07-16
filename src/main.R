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

###########################################
##              R PACKAGES               ##
###########################################

#install.packages(c("moments", "foreach", 
#                   "doParallel", "parallel"),dependencies = T)

###########################################
##            EXTERNAL TOOLS             ##
###########################################

# bgzip & tabix - https://github.com/TabakoffLab/General/wiki/Install-vcftools#Install_tabix
# bcftools - https://github.com/samtools/bcftools

###########################################
##           GLOBAL SETTINGS             ##
###########################################

nsim                    <- 5
slim_model              <- paste0("src/models/model", ".slim")
path_to_slim            <- "/usr/local/bin/slim"
slim_output_folder      <- "results/slim_output/"
path_to_bgzip           <- "/usr/local/bin/bgzip"
path_to_tabix           <- "/usr/local/bin/tabix"
path_to_bcftools        <- "/usr/local/bin/bcftools"
egglib_input_folder     <- "results/egglib_input/"
egglib_output_folder    <- "results/egglib_output/"
egglib_input_data      <- "data/egglib_input"
path_to_python          <- "/home/pavinato/py-egglib-3.0.0b21/bin/python"
path_to_egglib_summstat <- "bin/summstats_1.0.py" # this version works with egglib-3.0.0b21
reftable_file_folder    <- "results/reference_table"
arg <- commandArgs(TRUE)
#seed                    <- arg
seed                    <- 1234
set.seed(seed,"Mersenne-Twister")
parallel_sims           <- TRUE
num_of_threads          <- 25
remove_files            <- TRUE

############################################
##            SLiM SIMULATION             ##
##             SAMPLED VALUES             ##
## DEFINE PRIORS (upper and lower limits) ##
############################################

# THETA:
# We are going to sample thetas from a log uniform distribution; however SLiM is set up on population size Ne;
# Given the mutation rate, theta is calculated as THETA = Ne*mutation_rate
# For now, only sampling theta and calculating Ne with known mutation rate mu

# Mutation Rate
mu_rate = 1 # 0 = "FIXED"; 1 = "RANDOM" sample from prior
mu_min  = 1e-8
mu_max  = 1e-5

# Theta's effective population size 0 - Ne0
ne0_min = 1
ne0_max = 1000

# EFFECTIVE POPULATION SIZE 1 - Ne1
ne1_min = 1
ne1_max = 1000

# GENOME-WIDE DFE FOR BENEFICIAL MUTATIONS 
gammaM_gammak = TRUE # if TRUE, rate=1, then only gammaM will be sample from prior

gammaM_min = 0.001
gammaM_max = 1

gammak_min = 0.001 # this defines a lower and an upper limits of a uniform distribution where gamma MEAN and SHAPE (K) values;
gammak_max = 1

# PROPORTION OF THE GENOME THAT HOLDS BENEFICIAL MUTATIONS
PrGWSel_min = 0.00001 
PrGWSel_max = 1

# PROPORTION OF BENEFICIAL MUTATIONS
prbe_min = 0.00001 
prbe_max = 1

# DOMINANCE FOR GENOME-WIDE MUTATIONS
# Neutral mutations - m1 and m2
domN = 0 # 0 = "FIXED"; 1 = "RANDOM" sample from prior
domN_min = 0.5
domN_max = 1

# Beneficial mutations - m3
domB = 0 # 0 = "FIXED"; 1 = "RANDOM" sample from prior
domB_min = 0.5
dom_max = 1

# GENOME-WIDE RECOMBINATION RATE
rr_rate = 1 # 0 = "FIXED"; 1 = "RANDOM" sample from prior
rr_min  = 4.2 * 1e-8
rr_max  = 4.2 * 1e-5

############################################
##            SLiM SIMULATION             ##
##              FIXED VALUES              ##
############################################

SS1 = 80
SS2 = 115
ts2 = 8

genomeS = 135e+6
fragS   = 4.5e+6 
chrN    = 4

###########################################
##      EGGLIB SUMMSTAT SETTINGS         ##
###########################################

wss_wspan = 1000
sfs_bins = 10

###########################################
##    LOAD REQUIRED FUNCTIONS/PACKAGES   ##
###########################################

# load required libraries
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
if(parallel_sims){
  cl <- makeCluster(num_of_threads)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(moments))
  raw_reftable <- foreach(sim=seq_len(nsim),.combine=rbind) %dopar% {
                                                                        do_sim(sim, nsim, slim_model, path_to_slim, slim_output_folder,
                                                                               path_to_bgzip, path_to_tabix, path_to_bcftools,
                                                                               egglib_input_folder, egglib_output_folder, egglib_input_data,
                                                                               path_to_python, path_to_egglib_summstat,remove_files,  
                                                                               mu_rate, mu_min, mu_max, ne0_min, ne0_max, ne1_min, ne1_max,  
                                                                               gammaM_gammak, gammaM_min, gammaM_max, gammak_min, gammak_max, 
                                                                               PrGWSel_min, PrGWSel_max, prbe_min, prbe_max, 
                                                                               domN, domN_min, domN_max, domB, domB_min, domB_max, 
                                                                               rr_rate, rr_min, rr_max,
                                                                               SS1, SS2, ts2, genomeS, fragS, chrN,
                                                                               wss_wspan, sfs_bins
                                                                               )
    }

  stopCluster(cl)
  
}else{
  raw_reftable <- vector("list", nsim)
  for(sim in 1:nsim){
    raw_reftable[[sim]] <- do_sim(sim, nsim, slim_model, path_to_slim, slim_output_folder,
                               path_to_bgzip, path_to_tabix, path_to_bcftools,
                               egglib_input_folder, egglib_output_folder, egglib_input_data,
                               path_to_python, path_to_egglib_summstat,remove_files,  
                               mu_rate, mu_min, mu_max, ne0_min, ne0_max, ne1_min, ne1_max,  
                               gammaM_gammak, gammaM_min, gammaM_max, gammak_min, gammak_max, 
                               PrGWSel_min, PrGWSel_max, prbe_min, prbe_max, 
                               domN, domN_min, domN_max, domB, domB_min, domB_max, 
                               rr_rate, rr_min, rr_max,
                               SS1, SS2, ts2, genomeS, fragS, chrN,
                               wss_wspan, sfs_bins
                               )
  }
  raw_reftable <- do.call(rbind, raw_reftable)
}
gc()

#write.table(ref_table,
#            file      = paste0(reftable_file,".txt"),
#            sep       = "\t",
#            quote     = FALSE,
#            col.names = TRUE,
#            row.names = FALSE,
#            append    = FALSE)

#ref_table <- read.table(paste0(reftable_file,".txt"),header=T)
#dim(ref_table)
#head(ref_table)
#summary(ref_table)
save(raw_reftable,file=paste0(reftable_file,".RData"))

cat("\n Simulations finished\n\n")
