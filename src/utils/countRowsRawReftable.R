setwd("/fs/scratch/PAS1715/TrackSel_analysis/Tracking-selection")
dir()

N=2000

l <- matrix(NA, nrow=N, ncol = 1) 
for (i in 1296:N)
{
  load(paste0("batch.",i,"/reference_table.RData")); l[i] <- dim(raw_reftable)[1]
  #rm(raw_reftable)
}
table(l)
l[1:10]

