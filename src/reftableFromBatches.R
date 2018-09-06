#!/usr/bin/env Rscript

##########################################################################################################
##                            Tracking Selection in Time-Series Data with                               ## 
##                                      Forward-time Simulations                                        ##
##                                      SCRIPT TO POOL REFTABLE                                         ##
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

#install.packages(c("foreach", "doParallel", "parallel"),dependencies = T)

###########################################
##           GLOBAL SETTINGS             ##
###########################################

super_batch_number     <- 1
number_of_batches      <- 1000 
batch_size             <- 20
working_dir            <- "home/pavinato/Desktop/Tracking-selection-02.1/"
reftable_file_folder   <- "batch"
reftable_file          <- "reference_table"
pooled_reftable_output <- "results/" 
pooled_reftable_file   <- "pooled_reftable"
parallel_pool          <- TRUE
num_of_threads         <- 28

###########################################
##         LOAD REQUIRED PACKAGES        ##
###########################################

# load required libraries
if(parallel_pool){
  library(foreach)    #
  library(parallel)   # -> fo    r simulating in parallel
  library(doParallel) #
}


###########################################
##           POOLING REFTABLES           ##
###########################################

## FUNCTION FOR PARALLEL PROCESSING
#-------------------------------------------
poolFromBatches <- function(i){
    if (file.exists(paste0("/", working_dir, reftable_file_folder, ".", i, "/", reftable_file, ".RData"))){
      load(paste0("/", working_dir, reftable_file_folder, ".", i, "/", reftable_file, ".RData"))
    
      sim <- raw_reftable$sim + (i - 1) * batch_size
      raw_reftable <- data.frame(sim, raw_reftable[, c(2:11, 95:113, 12:94)])
      #raw_reftable <- data.frame(sim, raw_reftable[, -c(1)])
      
      return(raw_reftable)
    }
}


if(parallel_pool){
  cl <- makeCluster(num_of_threads)
  registerDoParallel(cl)
  #clusterEvalQ(cl)
  
  ## PARALLEL PROCESSING
  #------------------------------------------
  pooled_raw_reftable <- foreach(i=1:number_of_batches,.combine=rbind, 
                                                       .errorhandling="pass", .verbose = TRUE) %dopar% {
                                                                                                          poolFromBatches(i=i)
    }
  stopCluster(cl)
    
} else {
  ## SERIAL PROCESSING
  #------------------------------------------
  pooled_raw_reftable <- NULL
  for (i in 1:number_of_batches){
    if (file.exists(paste0("/", working_dir, reftable_file_folder, ".", i, "/", reftable_file, ".RData"))){
      load(paste0("/", working_dir, reftable_file_folder, ".", i, "/", reftable_file, ".RData"))
    
      sim <- raw_reftable$sim + (i - 1) * batch_size
      raw_reftable <- data.frame(sim, raw_reftable[, c(2:11, 95:113, 12:94)])
      #raw_reftable <- data.frame(sim, raw_reftable[, -c(1)])
      pooled_raw_reftable <- rbind(pooled_raw_reftable, raw_reftable)
    }
  }
}

if (!file_test("-d", paste0("/", working_dir, pooled_reftable_output))){
  dir.create(file.path(paste0("/", working_dir, pooled_reftable_output)))
}

save(pooled_raw_reftable,
     file=paste0("/", working_dir, pooled_reftable_output, pooled_reftable_file, "_", super_batch_number, ".RData"))

gc()

cat("\n Pooling batches from Super-batch #", super_batch_number,  "finished !!!\n\n")
