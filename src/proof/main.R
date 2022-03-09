#!/usr/bin/env Rscript

##########################################################################################################
##                            Tracking Selection in Time-Series Data with                               ## 
##                                      Forward-time Simulations                                        ##
##                                      PIPELINE FOR SIMULATIONS                                        ##
##########################################################################################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

### Remove last features
rm(list=ls())
ls()

###########################################
##              R PACKAGES               ##
###########################################

#install.packages(c("moments", "ROCR"
#                   "foreach", "doParallel", "parallel"),dependencies = T)

###########################################
##            EXTERNAL TOOLS             ##
###########################################

# bcftools - https://github.com/samtools/bcftools
# bgzip & tabix - https://github.com/TabakoffLab/General/wiki/Install-vcftools#Install_tabix
# (usually bgzip and tabis come with bcftools installation)
# egglib-3.0.0b22

###########################################
##           GLOBAL SETTINGS             ##
###########################################

nsim                    <- 10
path_to_slim_model      <- "src/models/"
slim_model_prefix       <- "model"
path_to_slim            <- "bin/slim" #SLiM v 3.1
slim_output_folder      <- "results/slim_output/"
path_to_bgzip           <- "bin/bgzip"    
path_to_tabix           <- "bin/tabix"    
path_to_bcftools        <- "bin/bcftools" 
egglib_input_folder     <- "results/egglib_input/"
egglib_output_folder    <- "results/egglib_output/"
path_to_python          <- "bin/python"
path_to_egglib_summstat <- "bin/summstats.py" 
reftable_file           <- "results/reference_table"
arg <- commandArgs(TRUE)
seed                    <- arg
#seed                    <- 12347319 # this is here only for small tests and for fixed pods (proof model)
set.seed(seed,"Mersenne-Twister")
parallel_sims           <- FALSE
num_of_threads          <- 6
remove_files            <- FALSE # it is not necessary to grow RF
debug_sim               <- TRUE
debug_output_folder     <- "results/debug_output/"
wfabc_input_file        <- FALSE # it is not necessary to grow RF
wfabc_input_folder      <- "results/wfabc_input/"
egglib_input_selcoeff   <- FALSE # it is not necessary to grow RF

############################################
##            SLiM SIMULATION             ##
##             FIXED VALUES               ##
############################################

# MODEL SELECTION
model_type = 1             # 1 = de novo beneficial mutations ("DN"); 
                           # 2 = background selection ("BS"); 
                           # 3 = selection on standing variation ("SV");

# GENOME SPECIFICATION
genomeS = 100e+6           # genomeS => Genome Size;
fragS   = 5e+4             # fragS   => Fragment size to define g1 or g2 elements;
chrN    = 1                # chrN    => Chromosome number to define independent genome blocks;

chrTAG = FALSE             # chrTAG  => if TRUE, SNPs are tagged with chromosome ID based on its position;

# DATASET SPECIFICATION 
data_type = 1              # 1 = Whole genome sequencing (WGS);
                           # 2 = RADseq;
radseq_readL = 100         # radseq_readL => the read length produced by the RAD library (bp);
radseq_cov = 0.20          # radseq_cov   => the proportion of the genome covered by the RADseq reads;

missing_data = 0.0         # missing_data => specify the proportion of missing genotypes per locus;

haplotype = FALSE          # haplotype => define how the homozygotes genotypes in the sample will be processed;

SS1 = 100                  # SS1 => Sample Size for T=1;
SS2 = 100                  # SS2 => Sample Size for T=2; Total sample size (SS1+SS2) should be at least > 3;
tau = 10                   # tau => Time between samples;

############################################
##            SLiM SIMULATION             ##
##             SAMPLED VALUES             ##
## DEFINE PRIORS (upper and lower limits) ##
############################################

# MUTATION RATE
mu_random = TRUE
mu_rate <- 1e-7
mu_min  = 1e-10
mu_max  = 1e-6

# POPULATION SIZE EQUILIBRIUM PHASE - Neq
neq_random = TRUE
neq_value <- 500
neq_min = 1
neq_max = 2000

# POPULATION CENSUS SIZE - Ncs
ncs_random = TRUE
ncs_value <- 500
ncs_min = 1
ncs_max = 2000

# GENOME-WIDE DFE FOR BENEFICIAL MUTATIONS 
# gamma mean
gammaM_random = TRUE
gammaM_value <- 0.1
gammaM_min = 0.001
gammaM_max = 1              

# gamma shape
gammak_random = TRUE
gammak_value <- 0.1         # gamma shape k must be positive;
gammak_min = 0.001          # this defines a lower and an upper limits of a uniform distribution where gamma MEAN and SHAPE (K) values;
gammak_max = 1

gammaM_gammak = TRUE        # if TRUE, rate=1, then only gammaM will be sample from prior;

# PROPORTION OF THE GENOME THAT CONTAINS BENEFICIAL MUTATIONS - G2 ELEMENTS
PrGWSel_random = TRUE
PrGWSel_value <- 0
PrGWSel_min = 0 
PrGWSel_max = 1

# PROPORTION OF BENEFICIAL MUTATION IN G2 ELEMENTS
prbe_random = TRUE
prbe_value <- 0.1
prbe_min = 1e-05
prbe_max = 1

# DOMINANCE FOR GENOME-WIDE MUTATIONS
# Neutral mutations - m1 and m2
domN_random = FALSE
domN <- 0.5
domN_min = 0.5
domN_max = 1

# Beneficial mutations - m3
domB_random = FALSE 
domB <- 0.5
domB_min = 0.5
domB_max = 1

# PER BASE RECOMBINATION RATE
rr_random = TRUE
rr_rate <- 5 * 1e-7
rr_min  = 5 * 1e-10
rr_max  = 5 * 1e-7

# SELFING RATE
selfing_random = FALSE
selfing_rate = 0.0         
selfing_min = 0.80
selfing_max = 1.00

# TIME OF SELECTION CHANGE
tc_random = FALSE
tc_value = 0             # from = 0 to = tau     

###########################################
##      EGGLIB SUMMSTAT SETTINGS         ##
###########################################

wss_wspan_run = 500 #300
sfs_bins_run = 10

add_WSSwspan_SFSbins_1 = TRUE
add_wss_wspan_1 = 5000 #500
add_sfs_bins_1 = 15

add_WSSwspan_SFSbins_2 = TRUE
add_wss_wspan_2 = 10000 #1000
add_sfs_bins_2 = 20

###########################################
##    LOAD REQUIRED FUNCTIONS/PACKAGES   ##
###########################################

# load required libraries
library(moments)
library(ROCR)
if(parallel_sims){
  library(foreach)    
  library(parallel)   
  library(doParallel)
}

# load other functions (from file distributed together with present file)
source("Tracking-selection/src/proof/fun.R")
#source("src/fun.R")

###########################################
##              SIMULATION               ##
###########################################

# Parallelization with doParallel and foreach
if(parallel_sims){
  cl <- makeCluster(num_of_threads)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(moments))
  raw_reftable <- foreach(sim=seq_len(nsim),.combine=rbind) %dopar% {   library(ROCR)
                                                                        do_sim(sim, nsim, 
                                                                               path_to_slim_model, slim_model_prefix,
                                                                               path_to_slim, slim_output_folder,
                                                                               path_to_bgzip, path_to_tabix, path_to_bcftools,
                                                                               egglib_input_folder, egglib_output_folder,
                                                                               path_to_python, path_to_egglib_summstat,
                                                                               remove_files, debug_sim, debug_output_folder,
                                                                               wfabc_input_file, wfabc_input_folder, egglib_input_selcoeff,  # add egglib_input_selcoeff 25/Oct/19
                                                                               model_type, model_title, genomeS, fragS, chrN, chrS,
                                                                               chrTAG, chromtagging, data_type, 
                                                                               radseq_readL, radseq_readN, radseq_cov,
                                                                               missing_data, haplotype, SS1, SS2, tau,
                                                                               mu_rate, mu_random, mu_min, mu_max, 
                                                                               neq_value, neq_random, neq_min, neq_max,
                                                                               ncs_value, ncs_random, ncs_min, ncs_max,  
                                                                               gammaM_value, gammak_value, gammaM_gammak, 
                                                                               gammaM_random, gammaM_min, gammaM_max, 
                                                                               gammak_random, gammak_min, gammak_max, 
                                                                               PrGWSel_value, PrGWSel_random, PrGWSel_min, PrGWSel_max, 
                                                                               prbe_value, prbe_random, prbe_min, prbe_max, 
                                                                               domN, domN_random, domN_min, domN_max, 
                                                                               domB, domB_random, domB_min, domB_max, 
                                                                               rr_rate, rr_random, rr_min, rr_max, rr_limits,
                                                                               selfing_rate, selfing_random, selfing_min, selfing_max,
                                                                               tc_random, tc_value, wss_wspan_run, sfs_bins_run,
                                                                               add_WSSwspan_SFSbins_1, add_wss_wspan_1, add_sfs_bins_1,
                                                                               add_WSSwspan_SFSbins_2, add_wss_wspan_2, add_sfs_bins_2
                                                                               )
    }

  stopCluster(cl)
  
}else{
  raw_reftable <- vector("list", nsim)
  for(sim in 1:nsim){
    raw_reftable[[sim]] <- do_sim(sim, nsim, 
                                  path_to_slim_model, slim_model_prefix,
                                  path_to_slim, slim_output_folder,
                                  path_to_bgzip, path_to_tabix, path_to_bcftools,
                                  egglib_input_folder, egglib_output_folder,
                                  path_to_python, path_to_egglib_summstat,
                                  remove_files, debug_sim, debug_output_folder,
                                  wfabc_input_file, wfabc_input_folder, egglib_input_selcoeff,  # add egglib_input_selcoeff 25/Oct/19
                                  model_type, model_title, genomeS, fragS, chrN, chrS,
                                  chrTAG, chromtagging, data_type, 
                                  radseq_readL, radseq_readN, radseq_cov,
                                  missing_data, haplotype, SS1, SS2, tau,
                                  mu_rate, mu_random, mu_min, mu_max, 
                                  neq_value, neq_random, neq_min, neq_max,
                                  ncs_value, ncs_random, ncs_min, ncs_max,  
                                  gammaM_value, gammak_value, gammaM_gammak, 
                                  gammaM_random, gammaM_min, gammaM_max, 
                                  gammak_random, gammak_min, gammak_max, 
                                  PrGWSel_value, PrGWSel_random, PrGWSel_min, PrGWSel_max, 
                                  prbe_value, prbe_random, prbe_min, prbe_max, 
                                  domN, domN_random, domN_min, domN_max, 
                                  domB, domB_random, domB_min, domB_max, 
                                  rr_rate, rr_random, rr_min, rr_max, rr_limits,
                                  selfing_rate, selfing_random, selfing_min, selfing_max,
                                  tc_random, tc_value, wss_wspan_run, sfs_bins_run,
                                  add_WSSwspan_SFSbins_1, add_wss_wspan_1, add_sfs_bins_1,
                                  add_WSSwspan_SFSbins_2, add_wss_wspan_2, add_sfs_bins_2
                                  )
  }
  raw_reftable <- do.call(rbind, raw_reftable)
}

raw_reftable_header <- c("model","seed","sim","mu","rr","selfing","Neq","Ncs","gammaMean","gammaShape","Tc","PrGWSel","PrMSel","meanNe1","meanNe2",paste0("timesNe2_",seq(from=0,to=(tau))),"averageGeneticLoad","lastGeneticLoad",
                         "POPPrMSel", "POPSelMean","POPSelSd","POPPrStrongMSel","POPStrongSelMean","POPStrongSelSd","SAMPrMSel","SAMSelMean","SAMSelSd","SAMPrStrongMSel","SAMStrongSelMean","SAMStrongSelSd","FSTfdrNS005","FSTfdrNS01","FSTfdrNS02","FSTfdrNS05","FSTfdrNS10",
                         "ID","MID","MT","S","DOM","GO","SAAF1","SAAF2","LSS_He","LSS_Dj","LSS_WCst","LSSp_He1","LSSp_He2","LSSp_Dj1","LSSp_Dj2",paste0("WSS",wss_wspan_run,"_He"),paste0("WSS",wss_wspan_run,"_Dj"),paste0("WSS",wss_wspan_run,"_WCst"),
                         paste0("WSS",wss_wspan_run,"_S"),paste0("WSS",wss_wspan_run,"_thetaW"),paste0("WSS",wss_wspan_run,"_Pi"),paste0("WSS",wss_wspan_run,"_D"),paste0("WSS",wss_wspan_run,"_Da"),paste0("WSS",wss_wspan_run,"_ZZ"),
                         paste0("WSS",wss_wspan_run,"_ZnS"),paste0("WSSp",wss_wspan_run,"_He1"),paste0("WSSp",wss_wspan_run,"_He2"),paste0("WSSp",wss_wspan_run,"_Dj1"),paste0("WSSp",wss_wspan_run,"_Dj2"),paste0("WSSp",wss_wspan_run,"_S1"),
                         paste0("WSSp",wss_wspan_run,"_S2"),paste0("WSSp",wss_wspan_run,"_thetaW1"),paste0("WSSp",wss_wspan_run,"_thetaW2"),paste0("WSSp",wss_wspan_run,"_Pi1"),paste0("WSSp",wss_wspan_run,"_Pi2"),paste0("WSSp",wss_wspan_run,"_D1"),
                         paste0("WSSp",wss_wspan_run,"_D2"),paste0("WSSp",wss_wspan_run,"_ZZ1"),paste0("WSSp",wss_wspan_run,"_ZZ2"),paste0("WSSp",wss_wspan_run,"_ZnS1"),paste0("WSSp",wss_wspan_run,"_ZnS2"),
                         "GSS_He","GSS_Dj","GSS_WCst","GSS_S","GSS_thetaW","GSS_Pi","GSS_D","GSS_Da","GSSp_He1","GSSp_He2","GSSp_S1","GSSp_S2","GSSp_thetaW1","GSSp_thetaW2","GSSp_Pi1","GSSp_Pi2","GSSp_D1","GSSp_D2","GSSp_Da1","GSSp_Da2",
                         paste0("SFSbin",sfs_bins_run,"_",seq(from=1,to=sfs_bins_run)),"MEAN_LSS_He","MEAN_LSS_Dj","MEAN_LSS_WCst","MEAN_LSSp_He1","MEAN_LSSp_He2","MEAN_LSSp_Dj1","MEAN_LSSp_Dj2",
                         paste0("MEAN_WSS",wss_wspan_run,"_He"),paste0("MEAN_WSS",wss_wspan_run,"_Dj"),paste0("MEAN_WSS",wss_wspan_run,"_WCst"),
                         paste0("MEAN_WSS",wss_wspan_run,"_S"),paste0("MEAN_WSS",wss_wspan_run,"_thetaW"),paste0("MEAN_WSS",wss_wspan_run,"_Pi"),
                         paste0("MEAN_WSS",wss_wspan_run,"_D"),paste0("MEAN_WSS",wss_wspan_run,"_Da"),paste0("MEAN_WSS",wss_wspan_run,"_ZZ"),
                         paste0("MEAN_WSS",wss_wspan_run,"_ZnS"),paste0("MEAN_WSSp",wss_wspan_run,"_He1"),paste0("MEAN_WSSp",wss_wspan_run,"_He2"),
                         paste0("MEAN_WSSp",wss_wspan_run,"_Dj1"),paste0("MEAN_WSSp",wss_wspan_run,"_Dj2"),paste0("MEAN_WSSp",wss_wspan_run,"_S1"),
                         paste0("MEAN_WSSp",wss_wspan_run,"_S2"),paste0("MEAN_WSSp",wss_wspan_run,"_thetaW1"),paste0("MEAN_WSSp",wss_wspan_run,"_thetaW2"),
                         paste0("MEAN_WSSp",wss_wspan_run,"_Pi1"),paste0("MEAN_WSSp",wss_wspan_run,"_Pi2"),paste0("MEAN_WSSp",wss_wspan_run,"_D1"),
                         paste0("MEAN_WSSp",wss_wspan_run,"_D2"),paste0("MEAN_WSSp",wss_wspan_run,"_ZZ1"),paste0("MEAN_WSSp",wss_wspan_run,"_ZZ2"),
                         paste0("MEAN_WSSp",wss_wspan_run,"_ZnS1"),paste0("MEAN_WSSp",wss_wspan_run,"_ZnS2"),
                         "VAR_LSS_He","VAR_LSS_Dj","VAR_LSS_WCst","VAR_LSSp_He1","VAR_LSSp_He2","VAR_LSSp_Dj1","VAR_LSSp_Dj2",
                         paste0("VAR_WSS",wss_wspan_run,"_He"),paste0("VAR_WSS",wss_wspan_run,"_Dj"),paste0("VAR_WSS",wss_wspan_run,"_WCst"),
                         paste0("VAR_WSS",wss_wspan_run,"_S"),paste0("VAR_WSS",wss_wspan_run,"_thetaW"),paste0("VAR_WSS",wss_wspan_run,"_Pi"),
                         paste0("VAR_WSS",wss_wspan_run,"_D"),paste0("VAR_WSS",wss_wspan_run,"_Da"),paste0("VAR_WSS",wss_wspan_run,"_ZZ"),
                         paste0("VAR_WSS",wss_wspan_run,"_ZnS"),paste0("VAR_WSSp",wss_wspan_run,"_He1"),paste0("VAR_WSSp",wss_wspan_run,"_He2"),
                         paste0("VAR_WSSp",wss_wspan_run,"_Dj1"),paste0("VAR_WSSp",wss_wspan_run,"_Dj2"),paste0("VAR_WSSp",wss_wspan_run,"_S1"),
                         paste0("VAR_WSSp",wss_wspan_run,"_S2"),paste0("VAR_WSSp",wss_wspan_run,"_thetaW1"),paste0("VAR_WSSp",wss_wspan_run,"_thetaW2"),
                         paste0("VAR_WSSp",wss_wspan_run,"_Pi1"),paste0("VAR_WSSp",wss_wspan_run,"_Pi2"),paste0("VAR_WSSp",wss_wspan_run,"_D1"),
                         paste0("VAR_WSSp",wss_wspan_run,"_D2"),paste0("VAR_WSSp",wss_wspan_run,"_ZZ1"),paste0("VAR_WSSp",wss_wspan_run,"_ZZ2"),
                         paste0("VAR_WSSp",wss_wspan_run,"_ZnS1"),paste0("VAR_WSSp",wss_wspan_run,"_ZnS2"),
                         "KURT_LSS_He","KURT_LSS_Dj","KURT_LSS_WCst","KURT_LSSp_He1","KURT_LSSp_He2","KURT_LSSp_Dj1","KURT_LSSp_Dj2",
                         paste0("KURT_WSS",wss_wspan_run,"_He"),paste0("KURT_WSS",wss_wspan_run,"_Dj"),paste0("KURT_WSS",wss_wspan_run,"_WCst"),
                         paste0("KURT_WSS",wss_wspan_run,"_S"),paste0("KURT_WSS",wss_wspan_run,"_thetaW"),paste0("KURT_WSS",wss_wspan_run,"_Pi"),
                         paste0("KURT_WSS",wss_wspan_run,"_D"),paste0("KURT_WSS",wss_wspan_run,"_Da"),paste0("KURT_WSS",wss_wspan_run,"_ZZ"),
                         paste0("KURT_WSS",wss_wspan_run,"_ZnS"),paste0("KURT_WSSp",wss_wspan_run,"_He1"),paste0("KURT_WSSp",wss_wspan_run,"_He2"),
                         paste0("KURT_WSSp",wss_wspan_run,"_Dj1"),paste0("KURT_WSSp",wss_wspan_run,"_Dj2"),paste0("KURT_WSSp",wss_wspan_run,"_S1"),
                         paste0("KURT_WSSp",wss_wspan_run,"_S2"),paste0("KURT_WSSp",wss_wspan_run,"_thetaW1"),paste0("KURT_WSSp",wss_wspan_run,"_thetaW2"),
                         paste0("KURT_WSSp",wss_wspan_run,"_Pi1"),paste0("KURT_WSSp",wss_wspan_run,"_Pi2"),paste0("KURT_WSSp",wss_wspan_run,"_D1"),
                         paste0("KURT_WSSp",wss_wspan_run,"_D2"),paste0("KURT_WSSp",wss_wspan_run,"_ZZ1"),paste0("KURT_WSSp",wss_wspan_run,"_ZZ2"),
                         paste0("KURT_WSSp",wss_wspan_run,"_ZnS1"),paste0("KURT_WSSp",wss_wspan_run,"_ZnS2"),
                         "SKEW_LSS_He","SKEW_LSS_Dj","SKEW_LSS_WCst","SKEW_LSSp_He1","SKEW_LSSp_He2","SKEW_LSSp_Dj1","SKEW_LSSp_Dj2",
                         paste0("SKEW_WSS",wss_wspan_run,"_He"),paste0("SKEW_WSS",wss_wspan_run,"_Dj"),paste0("SKEW_WSS",wss_wspan_run,"_WCst"),
                         paste0("SKEW_WSS",wss_wspan_run,"_S"),paste0("SKEW_WSS",wss_wspan_run,"_thetaW"),paste0("SKEW_WSS",wss_wspan_run,"_Pi"),
                         paste0("SKEW_WSS",wss_wspan_run,"_D"),paste0("SKEW_WSS",wss_wspan_run,"_Da"),paste0("SKEW_WSS",wss_wspan_run,"_ZZ"),
                         paste0("SKEW_WSS",wss_wspan_run,"_ZnS"),paste0("SKEW_WSSp",wss_wspan_run,"_He1"),paste0("SKEW_WSSp",wss_wspan_run,"_He2"),
                         paste0("SKEW_WSSp",wss_wspan_run,"_Dj1"),paste0("SKEW_WSSp",wss_wspan_run,"_Dj2"),paste0("SKEW_WSSp",wss_wspan_run,"_S1"),
                         paste0("SKEW_WSSp",wss_wspan_run,"_S2"),paste0("SKEW_WSSp",wss_wspan_run,"_thetaW1"),paste0("SKEW_WSSp",wss_wspan_run,"_thetaW2"),
                         paste0("SKEW_WSSp",wss_wspan_run,"_Pi1"),paste0("SKEW_WSSp",wss_wspan_run,"_Pi2"),paste0("SKEW_WSSp",wss_wspan_run,"_D1"),
                         paste0("SKEW_WSSp",wss_wspan_run,"_D2"),paste0("SKEW_WSSp",wss_wspan_run,"_ZZ1"),paste0("SKEW_WSSp",wss_wspan_run,"_ZZ2"),
                         paste0("SKEW_WSSp",wss_wspan_run,"_ZnS1"),paste0("SKEW_WSSp",wss_wspan_run,"_ZnS2"),
                         "Q05_LSS_He","Q05_LSS_Dj","Q05_LSS_WCst","Q05_LSSp_He1","Q05_LSSp_He2","Q05_LSSp_Dj1","Q05_LSSp_Dj2",
                         paste0("Q05_WSS",wss_wspan_run,"_He"),paste0("Q05_WSS",wss_wspan_run,"_Dj"),paste0("Q05_WSS",wss_wspan_run,"_WCst"),
                         paste0("Q05_WSS",wss_wspan_run,"_S"),paste0("Q05_WSS",wss_wspan_run,"_thetaW"),paste0("Q05_WSS",wss_wspan_run,"_Pi"),
                         paste0("Q05_WSS",wss_wspan_run,"_D"),paste0("Q05_WSS",wss_wspan_run,"_Da"),paste0("Q05_WSS",wss_wspan_run,"_ZZ"),
                         paste0("Q05_WSS",wss_wspan_run,"_ZnS"),paste0("Q05_WSSp",wss_wspan_run,"_He1"),paste0("Q05_WSSp",wss_wspan_run,"_He2"),
                         paste0("Q05_WSSp",wss_wspan_run,"_Dj1"),paste0("Q05_WSSp",wss_wspan_run,"_Dj2"),paste0("Q05_WSSp",wss_wspan_run,"_S1"),
                         paste0("Q05_WSSp",wss_wspan_run,"_S2"),paste0("Q05_WSSp",wss_wspan_run,"_thetaW1"),paste0("Q05_WSSp",wss_wspan_run,"_thetaW2"),
                         paste0("Q05_WSSp",wss_wspan_run,"_Pi1"),paste0("Q05_WSSp",wss_wspan_run,"_Pi2"),paste0("Q05_WSSp",wss_wspan_run,"_D1"),
                         paste0("Q05_WSSp",wss_wspan_run,"_D2"),paste0("Q05_WSSp",wss_wspan_run,"_ZZ1"),paste0("Q05_WSSp",wss_wspan_run,"_ZZ2"),
                         paste0("Q05_WSSp",wss_wspan_run,"_ZnS1"),paste0("Q05_WSSp",wss_wspan_run,"_ZnS2"),
                         "Q95_LSS_He","Q95_LSS_Dj","Q95_LSS_WCst","Q95_LSSp_He1","Q95_LSSp_He2","Q95_LSSp_Dj1","Q95_LSSp_Dj2",
                         paste0("Q95_WSS",wss_wspan_run,"_He"),paste0("Q95_WSS",wss_wspan_run,"_Dj"),paste0("Q95_WSS",wss_wspan_run,"_WCst"),
                         paste0("Q95_WSS",wss_wspan_run,"_S"),paste0("Q95_WSS",wss_wspan_run,"_thetaW"),paste0("Q95_WSS",wss_wspan_run,"_Pi"),
                         paste0("Q95_WSS",wss_wspan_run,"_D"),paste0("Q95_WSS",wss_wspan_run,"_Da"),paste0("Q95_WSS",wss_wspan_run,"_ZZ"),
                         paste0("Q95_WSS",wss_wspan_run,"_ZnS"),paste0("Q95_WSSp",wss_wspan_run,"_He1"),paste0("Q95_WSSp",wss_wspan_run,"_He2"),
                         paste0("Q95_WSSp",wss_wspan_run,"_Dj1"),paste0("Q95_WSSp",wss_wspan_run,"_Dj2"),paste0("Q95_WSSp",wss_wspan_run,"_S1"),
                         paste0("Q95_WSSp",wss_wspan_run,"_S2"),paste0("Q95_WSSp",wss_wspan_run,"_thetaW1"),paste0("Q95_WSSp",wss_wspan_run,"_thetaW2"),
                         paste0("Q95_WSSp",wss_wspan_run,"_Pi1"),paste0("Q95_WSSp",wss_wspan_run,"_Pi2"),paste0("Q95_WSSp",wss_wspan_run,"_D1"),
                         paste0("Q95_WSSp",wss_wspan_run,"_D2"),paste0("Q95_WSSp",wss_wspan_run,"_ZZ1"),paste0("Q95_WSSp",wss_wspan_run,"_ZZ2"),
                         paste0("Q95_WSSp",wss_wspan_run,"_ZnS1"),paste0("Q95_WSSp",wss_wspan_run,"_ZnS2")
                         )

add_WSSwspan_header_1 <- c(paste0("WSS",add_wss_wspan_1,"_He"),paste0("WSS",add_wss_wspan_1,"_Dj"),paste0("WSS",add_wss_wspan_1,"_WCst"),
                         paste0("WSS",add_wss_wspan_1,"_S"),paste0("WSS",add_wss_wspan_1,"_thetaW"),paste0("WSS",add_wss_wspan_1,"_Pi"),
                         paste0("WSS",add_wss_wspan_1,"_D"),paste0("WSS",add_wss_wspan_1,"_Da"),paste0("WSS",add_wss_wspan_1,"_ZZ"),
                         paste0("WSS",add_wss_wspan_1,"_ZnS"),paste0("WSSp",add_wss_wspan_1,"_He1"),paste0("WSSp",add_wss_wspan_1,"_He2"),
                         paste0("WSSp",add_wss_wspan_1,"_Dj1"),paste0("WSSp",add_wss_wspan_1,"_Dj2"),paste0("WSSp",add_wss_wspan_1,"_S1"),
                         paste0("WSSp",add_wss_wspan_1,"_S2"),paste0("WSSp",add_wss_wspan_1,"_thetaW1"),paste0("WSSp",add_wss_wspan_1,"_thetaW2"),
                         paste0("WSSp",add_wss_wspan_1,"_Pi1"),paste0("WSSp",add_wss_wspan_1,"_Pi2"),paste0("WSSp",add_wss_wspan_1,"_D1"),
                         paste0("WSSp",add_wss_wspan_1,"_D2"),paste0("WSSp",add_wss_wspan_1,"_ZZ1"),paste0("WSSp",add_wss_wspan_1,"_ZZ2"),
                         paste0("WSSp",add_wss_wspan_1,"_ZnS1"),paste0("WSSp",add_wss_wspan_1,"_ZnS2")
                         )

add_WSSwspan_header_2 <- c(paste0("WSS",add_wss_wspan_2,"_He"),paste0("WSS",add_wss_wspan_2,"_Dj"),paste0("WSS",add_wss_wspan_2,"_WCst"),
                         paste0("WSS",add_wss_wspan_2,"_S"),paste0("WSS",add_wss_wspan_2,"_thetaW"),paste0("WSS",add_wss_wspan_2,"_Pi"),
                         paste0("WSS",add_wss_wspan_2,"_D"),paste0("WSS",add_wss_wspan_2,"_Da"),paste0("WSS",add_wss_wspan_2,"_ZZ"),
                         paste0("WSS",add_wss_wspan_2,"_ZnS"),paste0("WSSp",add_wss_wspan_2,"_He1"),paste0("WSSp",add_wss_wspan_2,"_He2"),
                         paste0("WSSp",add_wss_wspan_2,"_Dj1"),paste0("WSSp",add_wss_wspan_2,"_Dj2"),paste0("WSSp",add_wss_wspan_2,"_S1"),
                         paste0("WSSp",add_wss_wspan_2,"_S2"),paste0("WSSp",add_wss_wspan_2,"_thetaW1"),paste0("WSSp",add_wss_wspan_2,"_thetaW2"),
                         paste0("WSSp",add_wss_wspan_2,"_Pi1"),paste0("WSSp",add_wss_wspan_2,"_Pi2"),paste0("WSSp",add_wss_wspan_2,"_D1"),
                         paste0("WSSp",add_wss_wspan_2,"_D2"),paste0("WSSp",add_wss_wspan_2,"_ZZ1"),paste0("WSSp",add_wss_wspan_2,"_ZZ2"),
                         paste0("WSSp",add_wss_wspan_2,"_ZnS1"),paste0("WSSp",add_wss_wspan_2,"_ZnS2")
                         )

add_SFSbins_header_1 <- c(paste0("SFSbin", add_sfs_bins_1, "_", seq(from=1, to=add_sfs_bins_1)),
                          paste0("MEAN_WSS",add_wss_wspan_1,"_He"),paste0("MEAN_WSS",add_wss_wspan_1,"_Dj"),paste0("MEAN_WSS",add_wss_wspan_1,"_WCst"),
                          paste0("MEAN_WSS",add_wss_wspan_1,"_S"),paste0("MEAN_WSS",add_wss_wspan_1,"_thetaW"),paste0("MEAN_WSS",add_wss_wspan_1,"_Pi"),
                          paste0("MEAN_WSS",add_wss_wspan_1,"_D"),paste0("MEAN_WSS",add_wss_wspan_1,"_Da"),paste0("MEAN_WSS",add_wss_wspan_1,"_ZZ"),
                          paste0("MEAN_WSS",add_wss_wspan_1,"_ZnS"),paste0("MEAN_WSSp",add_wss_wspan_1,"_He1"),paste0("MEAN_WSSp",add_wss_wspan_1,"_He2"),
                          paste0("MEAN_WSSp",add_wss_wspan_1,"_Dj1"),paste0("MEAN_WSSp",add_wss_wspan_1,"_Dj2"),paste0("MEAN_WSSp",add_wss_wspan_1,"_S1"),
                          paste0("MEAN_WSSp",add_wss_wspan_1,"_S2"),paste0("MEAN_WSSp",add_wss_wspan_1,"_thetaW1"),paste0("MEAN_WSSp",add_wss_wspan_1,"_thetaW2"),
                          paste0("MEAN_WSSp",add_wss_wspan_1,"_Pi1"),paste0("MEAN_WSSp",add_wss_wspan_1,"_Pi2"),paste0("MEAN_WSSp",add_wss_wspan_1,"_D1"),
                          paste0("MEAN_WSSp",add_wss_wspan_1,"_D2"),paste0("MEAN_WSSp",add_wss_wspan_1,"_ZZ1"),paste0("MEAN_WSSp",add_wss_wspan_1,"_ZZ2"),
                          paste0("MEAN_WSSp",add_wss_wspan_1,"_ZnS1"),paste0("MEAN_WSSp",add_wss_wspan_1,"_ZnS2"),
                          paste0("VAR_WSS",add_wss_wspan_1,"_He"),paste0("VAR_WSS",add_wss_wspan_1,"_Dj"),paste0("VAR_WSS",add_wss_wspan_1,"_WCst"),
                          paste0("VAR_WSS",add_wss_wspan_1,"_S"),paste0("VAR_WSS",add_wss_wspan_1,"_thetaW"),paste0("VAR_WSS",add_wss_wspan_1,"_Pi"),
                          paste0("VAR_WSS",add_wss_wspan_1,"_D"),paste0("VAR_WSS",add_wss_wspan_1,"_Da"),paste0("VAR_WSS",add_wss_wspan_1,"_ZZ"),
                          paste0("VAR_WSS",add_wss_wspan_1,"_ZnS"),paste0("VAR_WSSp",add_wss_wspan_1,"_He1"),paste0("VAR_WSSp",add_wss_wspan_1,"_He2"),
                          paste0("VAR_WSSp",add_wss_wspan_1,"_Dj1"),paste0("VAR_WSSp",add_wss_wspan_1,"_Dj2"),paste0("VAR_WSSp",add_wss_wspan_1,"_S1"),
                          paste0("VAR_WSSp",add_wss_wspan_1,"_S2"),paste0("VAR_WSSp",add_wss_wspan_1,"_thetaW1"),paste0("VAR_WSSp",add_wss_wspan_1,"_thetaW2"),
                          paste0("VAR_WSSp",add_wss_wspan_1,"_Pi1"),paste0("VAR_WSSp",add_wss_wspan_1,"_Pi2"),paste0("VAR_WSSp",add_wss_wspan_1,"_D1"),
                          paste0("VAR_WSSp",add_wss_wspan_1,"_D2"),paste0("VAR_WSSp",add_wss_wspan_1,"_ZZ1"),paste0("VAR_WSSp",add_wss_wspan_1,"_ZZ2"),
                          paste0("VAR_WSSp",add_wss_wspan_1,"_ZnS1"),paste0("VAR_WSSp",add_wss_wspan_1,"_ZnS2"),
                          paste0("KURT_WSS",add_wss_wspan_1,"_He"),paste0("KURT_WSS",add_wss_wspan_1,"_Dj"),paste0("KURT_WSS",add_wss_wspan_1,"_WCst"),
                          paste0("KURT_WSS",add_wss_wspan_1,"_S"),paste0("KURT_WSS",add_wss_wspan_1,"_thetaW"),paste0("KURT_WSS",add_wss_wspan_1,"_Pi"),
                          paste0("KURT_WSS",add_wss_wspan_1,"_D"),paste0("KURT_WSS",add_wss_wspan_1,"_Da"),paste0("KURT_WSS",add_wss_wspan_1,"_ZZ"),
                          paste0("KURT_WSS",add_wss_wspan_1,"_ZnS"),paste0("KURT_WSSp",add_wss_wspan_1,"_He1"),paste0("KURT_WSSp",add_wss_wspan_1,"_He2"),
                          paste0("KURT_WSSp",add_wss_wspan_1,"_Dj1"),paste0("KURT_WSSp",add_wss_wspan_1,"_Dj2"),paste0("KURT_WSSp",add_wss_wspan_1,"_S1"),
                          paste0("KURT_WSSp",add_wss_wspan_1,"_S2"),paste0("KURT_WSSp",add_wss_wspan_1,"_thetaW1"),paste0("KURT_WSSp",add_wss_wspan_1,"_thetaW2"),
                          paste0("KURT_WSSp",add_wss_wspan_1,"_Pi1"),paste0("KURT_WSSp",add_wss_wspan_1,"_Pi2"),paste0("KURT_WSSp",add_wss_wspan_1,"_D1"),
                          paste0("KURT_WSSp",add_wss_wspan_1,"_D2"),paste0("KURT_WSSp",add_wss_wspan_1,"_ZZ1"),paste0("KURT_WSSp",add_wss_wspan_1,"_ZZ2"),
                          paste0("KURT_WSSp",add_wss_wspan_1,"_ZnS1"),paste0("KURT_WSSp",add_wss_wspan_1,"_ZnS2"),
                          paste0("SKEW_WSS",add_wss_wspan_1,"_He"),paste0("SKEW_WSS",add_wss_wspan_1,"_Dj"),paste0("SKEW_WSS",add_wss_wspan_1,"_WCst"),
                          paste0("SKEW_WSS",add_wss_wspan_1,"_S"),paste0("SKEW_WSS",add_wss_wspan_1,"_thetaW"),paste0("SKEW_WSS",add_wss_wspan_1,"_Pi"),
                          paste0("SKEW_WSS",add_wss_wspan_1,"_D"),paste0("SKEW_WSS",add_wss_wspan_1,"_Da"),paste0("SKEW_WSS",add_wss_wspan_1,"_ZZ"),
                          paste0("SKEW_WSS",add_wss_wspan_1,"_ZnS"),paste0("SKEW_WSSp",add_wss_wspan_1,"_He1"),paste0("SKEW_WSSp",add_wss_wspan_1,"_He2"),
                          paste0("SKEW_WSSp",add_wss_wspan_1,"_Dj1"),paste0("SKEW_WSSp",add_wss_wspan_1,"_Dj2"),paste0("SKEW_WSSp",add_wss_wspan_1,"_S1"),
                          paste0("SKEW_WSSp",add_wss_wspan_1,"_S2"),paste0("SKEW_WSSp",add_wss_wspan_1,"_thetaW1"),paste0("SKEW_WSSp",add_wss_wspan_1,"_thetaW2"),
                          paste0("SKEW_WSSp",add_wss_wspan_1,"_Pi1"),paste0("SKEW_WSSp",add_wss_wspan_1,"_Pi2"),paste0("SKEW_WSSp",add_wss_wspan_1,"_D1"),
                          paste0("SKEW_WSSp",add_wss_wspan_1,"_D2"),paste0("SKEW_WSSp",add_wss_wspan_1,"_ZZ1"),paste0("SKEW_WSSp",add_wss_wspan_1,"_ZZ2"),
                          paste0("SKEW_WSSp",add_wss_wspan_1,"_ZnS1"),paste0("SKEW_WSSp",add_wss_wspan_1,"_ZnS2"),
                          paste0("Q05_WSS",add_wss_wspan_1,"_He"),paste0("Q05_WSS",add_wss_wspan_1,"_Dj"),paste0("Q05_WSS",add_wss_wspan_1,"_WCst"),
                          paste0("Q05_WSS",add_wss_wspan_1,"_S"),paste0("Q05_WSS",add_wss_wspan_1,"_thetaW"),paste0("Q05_WSS",add_wss_wspan_1,"_Pi"),
                          paste0("Q05_WSS",add_wss_wspan_1,"_D"),paste0("Q05_WSS",add_wss_wspan_1,"_Da"),paste0("Q05_WSS",add_wss_wspan_1,"_ZZ"),
                          paste0("Q05_WSS",add_wss_wspan_1,"_ZnS"),paste0("Q05_WSSp",add_wss_wspan_1,"_He1"),paste0("Q05_WSSp",add_wss_wspan_1,"_He2"),
                          paste0("Q05_WSSp",add_wss_wspan_1,"_Dj1"),paste0("Q05_WSSp",add_wss_wspan_1,"_Dj2"),paste0("Q05_WSSp",add_wss_wspan_1,"_S1"),
                          paste0("Q05_WSSp",add_wss_wspan_1,"_S2"),paste0("Q05_WSSp",add_wss_wspan_1,"_thetaW1"),paste0("Q05_WSSp",add_wss_wspan_1,"_thetaW2"),
                          paste0("Q05_WSSp",add_wss_wspan_1,"_Pi1"),paste0("Q05_WSSp",add_wss_wspan_1,"_Pi2"),paste0("Q05_WSSp",add_wss_wspan_1,"_D1"),
                          paste0("Q05_WSSp",add_wss_wspan_1,"_D2"),paste0("Q05_WSSp",add_wss_wspan_1,"_ZZ1"),paste0("Q05_WSSp",add_wss_wspan_1,"_ZZ2"),
                          paste0("Q05_WSSp",add_wss_wspan_1,"_ZnS1"),paste0("Q05_WSSp",add_wss_wspan_1,"_ZnS2"),
                          paste0("Q95_WSS",add_wss_wspan_1,"_He"),paste0("Q95_WSS",add_wss_wspan_1,"_Dj"),paste0("Q95_WSS",add_wss_wspan_1,"_WCst"),
                          paste0("Q95_WSS",add_wss_wspan_1,"_S"),paste0("Q95_WSS",add_wss_wspan_1,"_thetaW"),paste0("Q95_WSS",add_wss_wspan_1,"_Pi"),
                          paste0("Q95_WSS",add_wss_wspan_1,"_D"),paste0("Q95_WSS",add_wss_wspan_1,"_Da"),paste0("Q95_WSS",add_wss_wspan_1,"_ZZ"),
                          paste0("Q95_WSS",add_wss_wspan_1,"_ZnS"),paste0("Q95_WSSp",add_wss_wspan_1,"_He1"),paste0("Q95_WSSp",add_wss_wspan_1,"_He2"),
                          paste0("Q95_WSSp",add_wss_wspan_1,"_Dj1"),paste0("Q95_WSSp",add_wss_wspan_1,"_Dj2"),paste0("Q95_WSSp",add_wss_wspan_1,"_S1"),
                          paste0("Q95_WSSp",add_wss_wspan_1,"_S2"),paste0("Q95_WSSp",add_wss_wspan_1,"_thetaW1"),paste0("Q95_WSSp",add_wss_wspan_1,"_thetaW2"),
                          paste0("Q95_WSSp",add_wss_wspan_1,"_Pi1"),paste0("Q95_WSSp",add_wss_wspan_1,"_Pi2"),paste0("Q95_WSSp",add_wss_wspan_1,"_D1"),
                          paste0("Q95_WSSp",add_wss_wspan_1,"_D2"),paste0("Q95_WSSp",add_wss_wspan_1,"_ZZ1"),paste0("Q95_WSSp",add_wss_wspan_1,"_ZZ2"),
                          paste0("Q95_WSSp",add_wss_wspan_1,"_ZnS1"),paste0("Q95_WSSp",add_wss_wspan_1,"_ZnS2")
                          )

add_SFSbins_header_2 <- c(paste0("SFSbin", add_sfs_bins_2, "_", seq(from=1, to=add_sfs_bins_2)),
                          paste0("MEAN_WSS",add_wss_wspan_2,"_He"),paste0("MEAN_WSS",add_wss_wspan_2,"_Dj"),paste0("MEAN_WSS",add_wss_wspan_2,"_WCst"),
                          paste0("MEAN_WSS",add_wss_wspan_2,"_S"),paste0("MEAN_WSS",add_wss_wspan_2,"_thetaW"),paste0("MEAN_WSS",add_wss_wspan_2,"_Pi"),
                          paste0("MEAN_WSS",add_wss_wspan_2,"_D"),paste0("MEAN_WSS",add_wss_wspan_2,"_Da"),paste0("MEAN_WSS",add_wss_wspan_2,"_ZZ"),
                          paste0("MEAN_WSS",add_wss_wspan_2,"_ZnS"),paste0("MEAN_WSSp",add_wss_wspan_2,"_He1"),paste0("MEAN_WSSp",add_wss_wspan_2,"_He2"),
                          paste0("MEAN_WSSp",add_wss_wspan_2,"_Dj1"),paste0("MEAN_WSSp",add_wss_wspan_2,"_Dj2"),paste0("MEAN_WSSp",add_wss_wspan_2,"_S1"),
                          paste0("MEAN_WSSp",add_wss_wspan_2,"_S2"),paste0("MEAN_WSSp",add_wss_wspan_2,"_thetaW1"),paste0("MEAN_WSSp",add_wss_wspan_2,"_thetaW2"),
                          paste0("MEAN_WSSp",add_wss_wspan_2,"_Pi1"),paste0("MEAN_WSSp",add_wss_wspan_2,"_Pi2"),paste0("MEAN_WSSp",add_wss_wspan_2,"_D1"),
                          paste0("MEAN_WSSp",add_wss_wspan_2,"_D2"),paste0("MEAN_WSSp",add_wss_wspan_2,"_ZZ1"),paste0("MEAN_WSSp",add_wss_wspan_2,"_ZZ2"),
                          paste0("MEAN_WSSp",add_wss_wspan_2,"_ZnS1"),paste0("MEAN_WSSp",add_wss_wspan_2,"_ZnS2"),
                          paste0("VAR_WSS",add_wss_wspan_2,"_He"),paste0("VAR_WSS",add_wss_wspan_2,"_Dj"),paste0("VAR_WSS",add_wss_wspan_2,"_WCst"),
                          paste0("VAR_WSS",add_wss_wspan_2,"_S"),paste0("VAR_WSS",add_wss_wspan_2,"_thetaW"),paste0("VAR_WSS",add_wss_wspan_2,"_Pi"),
                          paste0("VAR_WSS",add_wss_wspan_2,"_D"),paste0("VAR_WSS",add_wss_wspan_2,"_Da"),paste0("VAR_WSS",add_wss_wspan_2,"_ZZ"),
                          paste0("VAR_WSS",add_wss_wspan_2,"_ZnS"),paste0("VAR_WSSp",add_wss_wspan_2,"_He1"),paste0("VAR_WSSp",add_wss_wspan_2,"_He2"),
                          paste0("VAR_WSSp",add_wss_wspan_2,"_Dj1"),paste0("VAR_WSSp",add_wss_wspan_2,"_Dj2"),paste0("VAR_WSSp",add_wss_wspan_2,"_S1"),
                          paste0("VAR_WSSp",add_wss_wspan_2,"_S2"),paste0("VAR_WSSp",add_wss_wspan_2,"_thetaW1"),paste0("VAR_WSSp",add_wss_wspan_2,"_thetaW2"),
                          paste0("VAR_WSSp",add_wss_wspan_2,"_Pi1"),paste0("VAR_WSSp",add_wss_wspan_2,"_Pi2"),paste0("VAR_WSSp",add_wss_wspan_2,"_D1"),
                          paste0("VAR_WSSp",add_wss_wspan_2,"_D2"),paste0("VAR_WSSp",add_wss_wspan_2,"_ZZ1"),paste0("VAR_WSSp",add_wss_wspan_2,"_ZZ2"),
                          paste0("VAR_WSSp",add_wss_wspan_2,"_ZnS1"),paste0("VAR_WSSp",add_wss_wspan_2,"_ZnS2"),
                          paste0("KURT_WSS",add_wss_wspan_2,"_He"),paste0("KURT_WSS",add_wss_wspan_2,"_Dj"),paste0("KURT_WSS",add_wss_wspan_2,"_WCst"),
                          paste0("KURT_WSS",add_wss_wspan_2,"_S"),paste0("KURT_WSS",add_wss_wspan_2,"_thetaW"),paste0("KURT_WSS",add_wss_wspan_2,"_Pi"),
                          paste0("KURT_WSS",add_wss_wspan_2,"_D"),paste0("KURT_WSS",add_wss_wspan_2,"_Da"),paste0("KURT_WSS",add_wss_wspan_2,"_ZZ"),
                          paste0("KURT_WSS",add_wss_wspan_2,"_ZnS"),paste0("KURT_WSSp",add_wss_wspan_2,"_He1"),paste0("KURT_WSSp",add_wss_wspan_2,"_He2"),
                          paste0("KURT_WSSp",add_wss_wspan_2,"_Dj1"),paste0("KURT_WSSp",add_wss_wspan_2,"_Dj2"),paste0("KURT_WSSp",add_wss_wspan_2,"_S1"),
                          paste0("KURT_WSSp",add_wss_wspan_2,"_S2"),paste0("KURT_WSSp",add_wss_wspan_2,"_thetaW1"),paste0("KURT_WSSp",add_wss_wspan_2,"_thetaW2"),
                          paste0("KURT_WSSp",add_wss_wspan_2,"_Pi1"),paste0("KURT_WSSp",add_wss_wspan_2,"_Pi2"),paste0("KURT_WSSp",add_wss_wspan_2,"_D1"),
                          paste0("KURT_WSSp",add_wss_wspan_2,"_D2"),paste0("KURT_WSSp",add_wss_wspan_2,"_ZZ1"),paste0("KURT_WSSp",add_wss_wspan_2,"_ZZ2"),
                          paste0("KURT_WSSp",add_wss_wspan_2,"_ZnS1"),paste0("KURT_WSSp",add_wss_wspan_2,"_ZnS2"),
                          paste0("SKEW_WSS",add_wss_wspan_2,"_He"),paste0("SKEW_WSS",add_wss_wspan_2,"_Dj"),paste0("SKEW_WSS",add_wss_wspan_2,"_WCst"),
                          paste0("SKEW_WSS",add_wss_wspan_2,"_S"),paste0("SKEW_WSS",add_wss_wspan_2,"_thetaW"),paste0("SKEW_WSS",add_wss_wspan_2,"_Pi"),
                          paste0("SKEW_WSS",add_wss_wspan_2,"_D"),paste0("SKEW_WSS",add_wss_wspan_2,"_Da"),paste0("SKEW_WSS",add_wss_wspan_2,"_ZZ"),
                          paste0("SKEW_WSS",add_wss_wspan_2,"_ZnS"),paste0("SKEW_WSSp",add_wss_wspan_2,"_He1"),paste0("SKEW_WSSp",add_wss_wspan_2,"_He2"),
                          paste0("SKEW_WSSp",add_wss_wspan_2,"_Dj1"),paste0("SKEW_WSSp",add_wss_wspan_2,"_Dj2"),paste0("SKEW_WSSp",add_wss_wspan_2,"_S1"),
                          paste0("SKEW_WSSp",add_wss_wspan_2,"_S2"),paste0("SKEW_WSSp",add_wss_wspan_2,"_thetaW1"),paste0("SKEW_WSSp",add_wss_wspan_2,"_thetaW2"),
                          paste0("SKEW_WSSp",add_wss_wspan_2,"_Pi1"),paste0("SKEW_WSSp",add_wss_wspan_2,"_Pi2"),paste0("SKEW_WSSp",add_wss_wspan_2,"_D1"),
                          paste0("SKEW_WSSp",add_wss_wspan_2,"_D2"),paste0("SKEW_WSSp",add_wss_wspan_2,"_ZZ1"),paste0("SKEW_WSSp",add_wss_wspan_2,"_ZZ2"),
                          paste0("SKEW_WSSp",add_wss_wspan_2,"_ZnS1"),paste0("SKEW_WSSp",add_wss_wspan_2,"_ZnS2"),
                          paste0("Q05_WSS",add_wss_wspan_2,"_He"),paste0("Q05_WSS",add_wss_wspan_2,"_Dj"),paste0("Q05_WSS",add_wss_wspan_2,"_WCst"),
                          paste0("Q05_WSS",add_wss_wspan_2,"_S"),paste0("Q05_WSS",add_wss_wspan_2,"_thetaW"),paste0("Q05_WSS",add_wss_wspan_2,"_Pi"),
                          paste0("Q05_WSS",add_wss_wspan_2,"_D"),paste0("Q05_WSS",add_wss_wspan_2,"_Da"),paste0("Q05_WSS",add_wss_wspan_2,"_ZZ"),
                          paste0("Q05_WSS",add_wss_wspan_2,"_ZnS"),paste0("Q05_WSSp",add_wss_wspan_2,"_He1"),paste0("Q05_WSSp",add_wss_wspan_2,"_He2"),
                          paste0("Q05_WSSp",add_wss_wspan_2,"_Dj1"),paste0("Q05_WSSp",add_wss_wspan_2,"_Dj2"),paste0("Q05_WSSp",add_wss_wspan_2,"_S1"),
                          paste0("Q05_WSSp",add_wss_wspan_2,"_S2"),paste0("Q05_WSSp",add_wss_wspan_2,"_thetaW1"),paste0("Q05_WSSp",add_wss_wspan_2,"_thetaW2"),
                          paste0("Q05_WSSp",add_wss_wspan_2,"_Pi1"),paste0("Q05_WSSp",add_wss_wspan_2,"_Pi2"),paste0("Q05_WSSp",add_wss_wspan_2,"_D1"),
                          paste0("Q05_WSSp",add_wss_wspan_2,"_D2"),paste0("Q05_WSSp",add_wss_wspan_2,"_ZZ1"),paste0("Q05_WSSp",add_wss_wspan_2,"_ZZ2"),
                          paste0("Q05_WSSp",add_wss_wspan_2,"_ZnS1"),paste0("Q05_WSSp",add_wss_wspan_2,"_ZnS2"),
                          paste0("Q95_WSS",add_wss_wspan_2,"_He"),paste0("Q95_WSS",add_wss_wspan_2,"_Dj"),paste0("Q95_WSS",add_wss_wspan_2,"_WCst"),
                          paste0("Q95_WSS",add_wss_wspan_2,"_S"),paste0("Q95_WSS",add_wss_wspan_2,"_thetaW"),paste0("Q95_WSS",add_wss_wspan_2,"_Pi"),
                          paste0("Q95_WSS",add_wss_wspan_2,"_D"),paste0("Q95_WSS",add_wss_wspan_2,"_Da"),paste0("Q95_WSS",add_wss_wspan_2,"_ZZ"),
                          paste0("Q95_WSS",add_wss_wspan_2,"_ZnS"),paste0("Q95_WSSp",add_wss_wspan_2,"_He1"),paste0("Q95_WSSp",add_wss_wspan_2,"_He2"),
                          paste0("Q95_WSSp",add_wss_wspan_2,"_Dj1"),paste0("Q95_WSSp",add_wss_wspan_2,"_Dj2"),paste0("Q95_WSSp",add_wss_wspan_2,"_S1"),
                          paste0("Q95_WSSp",add_wss_wspan_2,"_S2"),paste0("Q95_WSSp",add_wss_wspan_2,"_thetaW1"),paste0("Q95_WSSp",add_wss_wspan_2,"_thetaW2"),
                          paste0("Q95_WSSp",add_wss_wspan_2,"_Pi1"),paste0("Q95_WSSp",add_wss_wspan_2,"_Pi2"),paste0("Q95_WSSp",add_wss_wspan_2,"_D1"),
                          paste0("Q95_WSSp",add_wss_wspan_2,"_D2"),paste0("Q95_WSSp",add_wss_wspan_2,"_ZZ1"),paste0("Q95_WSSp",add_wss_wspan_2,"_ZZ2"),
                          paste0("Q95_WSSp",add_wss_wspan_2,"_ZnS1"),paste0("Q95_WSSp",add_wss_wspan_2,"_ZnS2")
                          )

if (add_WSSwspan_SFSbins_1 & add_WSSwspan_SFSbins_2){
  
  raw_reftable_header <- c(raw_reftable_header[1:(15+(tau+1)+60)],add_WSSwspan_header_1,add_WSSwspan_header_2,
                           raw_reftable_header[((15+(tau+1)+60)+1):length(raw_reftable_header)],add_SFSbins_header_1,add_SFSbins_header_2)

} else if (add_WSSwspan_SFSbins_1){
  
  raw_reftable_header <- c(raw_reftable_header[1:(15+(tau+1)+60)],add_WSSwspan_header_1,
                           raw_reftable_header[((15+(tau+1)+60)+1):length(raw_reftable_header)],add_SFSbins_header_1)
  
} else if (add_WSSwspan_SFSbins_2){
  
  raw_reftable_header <- c(raw_reftable_header[1:(15+(tau+1)+60)],add_WSSwspan_header_2,
                           raw_reftable_header[((15+(tau+1)+60)+1):length(raw_reftable_header)],add_SFSbins_header_2)
}

colnames(raw_reftable) <- raw_reftable_header

save(raw_reftable,file=paste0(reftable_file,".RData"))

gc()

cat("\n Simulations finished\n\n")
