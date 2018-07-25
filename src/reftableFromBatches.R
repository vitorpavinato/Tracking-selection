#!/usr/bin/env Rscript

rm(list=ls())
ls()

super_batch_number     <- 1
number_of_batches      <- 1000 
working_dir            <- "home/pavinato/project/Tracking-selection/"
reftable_file_folder   <- "batch"
reftable_file          <- "reference_table"
pooled_reftable_folder <- "pooled_reference_table/" 
pooled_global_reftable <- "pooled_global_reftable"
pooled_locus_reftable  <- "pooled_locus_reftable"

pooled_raw_global_reftable <- NULL
pooled_raw_locus_reftable <- NULL
for (i in 1:number_of_batches){
  if (file.exists(paste0("/", working_dir, reftable_file_folder, ".", i, "/", reftable_file, ".RData"))){
    load(paste0("/", working_dir, reftable_file_folder, ".", i, "/", reftable_file, ".RData"))
    
    raw_global_reftable <- raw_reftable[, 1:94]
    raw_global_reftable <- raw_global_reftable[1, ]
    pooled_raw_global_reftable <- rbind(pooled_raw_global_reftable, raw_global_reftable)
    
    raw_locus_reftable <- raw_reftable[, c(1, 95:length(names(raw_reftable)))]
    pooled_raw_locus_reftable <- rbind(pooled_raw_locus_reftable, raw_locus_reftable)
  }
}

pooled_reftable <- rbind(pooled_reftable,ref_table)

if (!file_test("-d", paste0("/", working_dir, pooled_reftable_folder))){
  dir.create(file.path(paste0("/", working_dir, pooled_reftable_folder)))
}

save(pooled_raw_global_reftable,
     file=paste0("/", working_dir, pooled_reftable_folder, pooled_global_reftable, "_", super_batch_number, ".RData"))

save(pooled_raw_locus_reftable,
     file=paste0("/", working_dir, pooled_reftable_folder, pooled_loucus_reftable, "_", super_batch_number, ".RData"))

gc()

cat("\n Pooling batches from Super-batch #", super_batch_number,  "finished !!!\n\n")
