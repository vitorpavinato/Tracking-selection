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
# (usually bgzip and tabix come with bcftools installation)
# egglib-3.0.0b22

###########################################
##           GLOBAL SETTINGS             ##
###########################################

nsim                    <- 5
path_to_slim_model      <- "src/models/bees/"
slim_model_prefix       <- "model"
path_to_slim            <- "/usr/local/bin/slim"   
slim_output_folder      <- "results/slim_output/"
path_to_bgzip           <- "/usr/local/bin/bgzip" 
path_to_tabix           <- "/usr/local/bin/tabix" 
path_to_bcftools        <- "/usr/local/bin/bcftools" 
egglib_input_folder     <- "results/egglib_input/"
egglib_output_folder    <- "results/egglib_output/"
path_to_python          <- "/Users/vitorpavinato/myPython2venvs/python-egglib3.0.0b22/bin/python" 
path_to_egglib_summstat <- "bin/summstats.py" 
reftable_file           <- "results/reference_table"
arg <- commandArgs(TRUE)
seed                    <- arg
#seed                    <- 1234 # this is here only for small tests
set.seed(seed,"Mersenne-Twister")
parallel_sims           <- FALSE
num_of_threads          <- 4
remove_files            <- TRUE
debug_sim               <- TRUE
debug_output_folder     <- "results/debug_output/"
wfabc_input_file        <- FALSE
wfabc_input_folder      <- "results/wfabc_input/"
egglib_input_selcoeff   <- FALSE 

############################################
##            SLiM SIMULATION             ##
##             FIXED VALUES               ##
############################################

# MODEL SELECTION
model_type = 1             # 1 = de novo beneficial mutations ("DN"); 
                           # 2 = background selection ("BS"); 
                           # 3 = selection on standing variation ("SV");

# GENOME SPECIFICATION
genomeS = 250e+6           # genomeS => Genome Size;
fragS   = 5e+4             # fragS   => Fragment size to define g1 or g2 elements;
chrN    = 1                # chrN    => Chromosome number to define independent genome blocks;

chrTAG = FALSE             # chrTAG  => if TRUE, SNPs are tagged with chromosome ID based on its position;

# DATASET SPECIFICATION 
data_type = 1              # 1 = Whole genome sequencing (WGS);
                           # 2 = RADseq;
radseq_readL   = 100       # radseq_readL => the read length produced by the RAD library (bp);
radseq_cov     = 1.5e-4    # radseq_cov   => the proportion of the genome covered by the RADseq reads;
one_snp_radseq = TRUE      # one_snp_radseq => if true, it allows to take one SNP per RADseq loci;

missing_data = 0.00        # missing_data => specify the proportion of missing genotypes per locus;

haplotype = FALSE          # haplotype => define how the homozygotes genotypes in the sample will be processed;

SSs = c(2,6,2,2,5,5,8)     # SSs => a string with the sample size for each TS;
tau = c(55,57,66,86,89,104) # tau => a string with the Time between samples;

############################################
##            SLiM SIMULATION             ##
##             SAMPLED VALUES             ##
## DEFINE PRIORS (upper and lower limits) ##
############################################

# MUTATION RATE
mu_random = TRUE
mu_rate <- 1e-7
mu_prior <- 2              # 1 = log-uniform; 2 = log-normal; 
mu_min  = 1e-10            # min for the log-uniform, otherwise NULL;
mu_max  = 1e-4             # max for the log-uniform, otherwise NULL;
mu_mean = 3.4e-9           # mean for the log-normal, otherwise NULL;
mu_sd   = 0.5              # sd for the log-normal, otherwise NULL;

# POPULATION SIZE EQUILIBRIUM PHASE - Neq
neq_random = TRUE
neq_value <- 50
neq_min = 2
neq_max = 10000

# POPULATION CENSUS SIZE - Ncs
ncs_random = TRUE
ncs_value <- 50
ncs_min = 2
ncs_max = 10000

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
PrGWSel_value <- 0.25
PrGWSel_min = 0 
PrGWSel_max = 1

# PROPORTION OF BENEFICIAL MUTATION IN G2 ELEMENTS
prbe_random = TRUE
prbe_value <- 0.25
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
rr_rate <- 1e-7
rr_prior <- 1             # 1 = log-uniform; 2 = log-normal; 
rr_min  = 1e-8            # min for the log-uniform, otherwise NULL; 
rr_max  = 1e-4            # max for the log-uniform, otherwise NULL;
rr_mean = 2e-7            # mean for the log-normal, otherwise NULL;
rr_sd   = 0.5             # sd for the log-normal, otherwise NULL;

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

wss_wspan_run = 10e+3
sfs_bins_run = 10

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
source("src/bees/fun_bees.R")

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
                                                                               wfabc_input_file, wfabc_input_folder, egglib_input_selcoeff,               # add egglib_input_selcoeff
                                                                               model_type, model_title, genomeS, fragS, 
                                                                               chrN, chrTAG, chromtagging, chrs_lowers,                                   # removeer chrS
                                                                               data_type, radseq_readL, radseq_interval, radseq_sampled, one_snp_radseq,  # add chrs_lowers  #radseq_readN to radseq_lociN ## remove radseq_readN ## remove radseq_cov
                                                                               missing_data, haplotype, SSs, tau,                                         # add radseq_sampled and radseq_interval and one_snp_radseq; SS1 and SS2 to SSs
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
                                                                               tc_random, tc_value, wss_wspan_run, sfs_bins_run
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
                                  wfabc_input_file, wfabc_input_folder, egglib_input_selcoeff,               # add egglib_input_selcoeff
                                  model_type, model_title, genomeS, fragS, 
                                  chrN, chrTAG, chromtagging, chrs_lowers,                                   # removeer chrS
                                  data_type, radseq_readL, radseq_interval, radseq_sampled, one_snp_radseq,  # add chrs_lowers  #radseq_readN to radseq_lociN ## remove radseq_readN ## remove radseq_cov
                                  missing_data, haplotype, SSs, tau,                                         # add radseq_sampled and radseq_interval and one_snp_radseq; SS1 and SS2 to SSs
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
                                  tc_random, tc_value, wss_wspan_run, sfs_bins_run
                                  )
  }
  raw_reftable <- do.call(rbind, raw_reftable)
}

reftable_header_incommon <- c("meanNe2","averageGenLoad","POPPrMSel", "POPSelMean","POPSelSd","POPPrStrongMSel","POPStrongSelMean","POPStrongSelSd",
                              "SAMPrMSel","SAMSelMean","SAMSelSd","SAMPrStrongMSel","SAMStrongSelMean","SAMStrongSelSd",
                              "FSTfdrNS005","FSTfdrNS01","FSTfdrNS02","FSTfdrNS05","FSTfdrNS10",
                              "ID","MID","MT","S","DOM","GO","SAAF1","SAAF2","LSS_He","LSS_Dj","LSS_WCst","LSSp_He1","LSSp_He2","LSSp_Dj1","LSSp_Dj2",
                              paste0("WSS",wss_wspan_run,"_He"),paste0("WSS",wss_wspan_run,"_Dj"),paste0("WSS",wss_wspan_run,"_WCst"),
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
                              paste0("Q95_WSSp",wss_wspan_run,"_ZnS1"),paste0("Q95_WSSp",wss_wspan_run,"_ZnS2"))

raw_reftable_header_1 <- c("model","seed","sim","mu","rr","selfing","Neq","Ncs","gammaMean","gammaShape","Tc","PrGWSel","PrMSel",
                           "meanNe1",paste0("timesNe2_",seq(from=0,to=(tau[6]))))

raw_reftable_header_2 <- c(paste0(reftable_header_incommon, "_Ava"),
                           paste0(reftable_header_incommon, "_Hum"),
                           paste0(reftable_header_incommon, "_Dav"),
                           paste0(reftable_header_incommon, "_Sta"),
                           paste0(reftable_header_incommon, "_Ste"),
                           paste0(reftable_header_incommon, "_Riv"),
                           paste0(reftable_header_incommon, "_Pla"))

colnames(raw_reftable) <- c(raw_reftable_header_1, raw_reftable_header_2)

save(raw_reftable,file=paste0(reftable_file,".RData"))

warnings()

gc()

cat("\n Simulations finished\n\n")
