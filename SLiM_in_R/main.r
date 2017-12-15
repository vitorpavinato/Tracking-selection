#################################################
## Tracking Selection in Time-Series Data with ## 
##          Forward-time Simulations           ##
#################################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#install.packages(c("abc","abcrf"));
#library(abc);
#library(abcrf);

#sessionInfo()$otherPkgs$abc$Version;
#sessionInfo()$otherPkgs$abcrf$Version;

recordSessionInfo <- sessionInfo();
setwd("/Users/vitorpavinato/Documents/My_repositories/Tracking-selection/SLiM_in_R");
save(recordSessionInfo,file="results/sessionInfo.RData");
save.image(file = "results/workspaceFile.RData");

#################################################
##       Runing Forward-time Simulations       ##
#################################################

## Creating SLiM infiles

source("src/create_slim_infile.r");

# Simulating N from prior;
N_size <- as.integer(runif(3, 100, 300));

for (n in 1:length(N_size)){
cat(paste("You simulated sample",n,"of size N =", N_size[n],"\n"));
};

if (file.exists(dir('data/', pattern = '.slim$', full.names = TRUE))){
  file_removed <- file.remove(dir('data/', pattern = '.slim$', full.names = TRUE))
}

model1(ne = N_size,
       filename = 'infile_N_size', 
       folder = 'data/');  ## ate aqui esta ok;

## Running SLiM simulations 

source("src/slim.r");

# Random seeds;
random_seeds <- as.integer(runif(length(N_size), 100000, 900000));

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters and the output in the same file;
# This function is one I worked last time;
slim(    parm = N_size, 
         seed = random_seeds,
     filename = 'outfile_N_size',
       folder = 'data/',
       infile = 'infile_N_size');

# SLiM2 runs v2:
# It runs SLiM2 and saves only the output;
run_slim_infile2(parm = N_size, 
                seed = random_seeds,
                filename = 'outfile_N_size',
                folder = 'data/',
                infile = 'infile_N_size');

tmp <- system2(command = '/usr/local/bin/slim', 
               args = paste0('-s', ' ', 755552, ' ', 'data/', 'infile_N_size', '_', 1293, '.slim'), 
               stdout = TRUE);

tmp[17:270]
tmp[271:length(tmp)]
