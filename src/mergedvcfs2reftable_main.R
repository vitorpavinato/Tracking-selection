#!/usr/bin/env Rscript

##########################################################################################################
##                            Tracking Selection in Time-Series Data with                               ## 
##                                      Forward-time Simulations                                        ##
##########################################################################################################

## Script to:
# 1 - convert merged vcfs data (table) to egglib input file;
# 2 - calculate the summary statistics with egglib summstats.py;
# 3 - prepare reference table to be used in ABC-RF analysis

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

### Remove last features
rm(list=ls())
ls()

###########################################
##           GLOBAL SETTINGS             ##
###########################################

dataset_path <- "data/ApisMellifera/merged_file/"
dataset_name <- "placerita_merged.txt"
pop_prefix <- "Pla"
snpsTable2egglib <- TRUE
snptTable2wfabc <- TRUE
wfacbInput_path <- "data/ApisMellifera/wfabc/"
wfabcInput_file <- "placerita_wfabc_input.txt"
globalStats_path <- "data/ApisMellifera/globalStats_reftable/"
globalStats_file <- "placerita_globalStats_reftable"

##########################################
##       EGGLIB SUMMSTAT SETTINGS       ##
##########################################

python_path     <- "/Users/vitorpavinato/myPython2venvs/python-egglib3.0.0b22/bin/python"
egglib_summstat <- "bin/summstats.py" # this version works with egglib-3.0.0b22

wss_wspan = 10e+3
sfs_bins = 10

snps2egglib_path <- "data/ApisMellifera/egglib/"
snps2egglib_file <- "placerita_egglib_input.txt"
egglib2summstats_path <- "data/ApisMellifera/egglib/"
egglib2summstats_file <- "placerita_summstats.txt"

############################################
##              SAMPLE SIZE               ##
############################################

SS1 = 5
SS2 = 6
tau = 15

###########################################
##    LOAD REQUIRED FUNCTIONS/PACKAGES   ##
###########################################

# load required libraries
#list.of.packages <- c("moments")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
library(moments)

source("src/mergedvcfs2reftable_fun.R")

###########################################
##        DATASET REFERENCE TABLE        ##
###########################################

data2globalStats <- do_data(dataset_path, 
                            dataset_name,
                            pop_prefix,
                            snpsTable2egglib,
                            python_path, 
                            egglib_summstat, 
                            wss_wspan, 
                            sfs_bins,
                            SS1,
                            SS2,
                            snps2egglib_path,
                            snps2egglib_file,
                            egglib2summstats_path,
                            egglib2summstats_file,
                            snptTable2wfabc,
                            wfacbInput_path,
                            wfabcInput_file,
                            tau
                          )

save(data2globalStats,file=paste0(globalStats_path,globalStats_file,".RData"))

