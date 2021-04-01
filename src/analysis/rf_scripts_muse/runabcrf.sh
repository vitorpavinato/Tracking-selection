#!/bin/bash

#SBATCH -J trackselABCRf
#SBATCH -o trackselABCRf.o%A
#SBATCH -e trackselABCRf.e%A
#SBATCH --partition=imag

module purge
module load cv-standard R

#Rscript ./scriptGrowRF.R
Rscript ./scriptPredictRF.R
