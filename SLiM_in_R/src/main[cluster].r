#################################################
## Tracking Selection in Time-Series Data with ## 
##          Forward-time Simulations           ##
#################################################

## Vitor Pavinato & Miguel Navascues
## vitor.pavinato@supagro.fr
## INRA

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#install.packages(c("abc","abcrf","plyr"));
#library(abc);
#library(abcrf);
#library(plyr); # for the convertEggLib function & join function
#library(dplyr)

#sessionInfo()$otherPkgs$abc$Version;
#sessionInfo()$otherPkgs$abcrf$Version;
#sessionInfo()$otherPkgs$plyr$Version;

recordSessionInfo <- sessionInfo();
working_dir <- "/home/pavinato/trackingsel_1"

if (is.na(match(working_dir,getwd()))){
  setwd(working_dir)
} else {
  print("Working directory have already been set")
};

save(recordSessionInfo,file="results/sessionInfo.RData");
save.image(file = "results/workspaceFile.RData");

#################################################
##       Runing Forward-time Simulations       ##
#################################################

#load(file = "results/workspaceFile.RData");

## Creating SLiM infiles

# Set the number of simulations

nsim <- 10000; 

# Set values for the prior - only for Ne

sim_ne <- as.integer(runif(nsim, 100, 10000));
sim_seed <- as.integer(runif(nsim, 100000000, 900000000));

priors <- as.data.frame(cbind(sim=1:nsim, seed=sim_seed,ne=sim_ne));
#hist(priors$ne)

# Create SLiM infiles using values taken from the prior

source("src/create_slim_infile.r");

ifelse(file.exists(dir("results/slim-infiles/", pattern = ".slim$", full.names = TRUE)), 
       file.remove(dir("results/slim-infiles/", pattern = ".slim$", full.names = TRUE)), NULL);
      
createslim <- model1(      ne = sim_ne,
                     filename = "infile_slim", 
                       folder = "results/slim-infiles/",
                       output = "output_slim",
                    folderout = "results/slim-outputs");

priors <- cbind(priors, createslim); 
priors <- priors[,-4]

write.table(priors, file = "results/ref-tables/sim_priors.txt", quote = FALSE, row.names = FALSE)

## Running SLiM simulations 

source("src/slim.r");

ifelse(file.exists(dir('results/slim-outputs/', pattern = "output_slim_", full.names = TRUE)), 
       file.remove(dir('results/slim-outputs/', pattern = "output_slim_", full.names = TRUE)), NULL);

ifelse(file.exists(dir('results/slim-outputs/', pattern = "outfile_slim_", full.names = TRUE)), 
       file.remove(dir('results/slim-outputs/', pattern = "outfile_slim_", full.names = TRUE)), NULL);

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters of each simulation in an outfile;
# Thanks to Eidos code each simulation now generate 2 files for each sampling;

slim(     seed = sim_seed,
       outfile = 'outfile_slim',
     folderout = 'results/slim-outputs/',
        infile = 'infile_slim',
      folderin = 'results/slim-infiles/');

# SLiM2 runs v2:
# It runs SLiM2 and saves only the outputs;

# slimclean(    seed = sim_seed, 
#             infile = 'infile_slim',      
#           folderin = 'infiles/');

save.image(file = "results/workspaceFile.RData");

#################################################
##     Calculating the Summary Statistics      ##
#################################################

#load(file = "results/workspaceFile.RData");

## Converting SLiM outputs to EggLib inputs

#source("src/create_ref_table.r");

#ifelse(file.exists(dir('results/egglib-inputs/', pattern = "input_egglib_", full.names = TRUE)), 
#       file.remove(dir('results/egglib-inputs/', pattern = "input_egglib_", full.names = TRUE)), NULL)

#dataconvertion <- convertEggLib(       nsim = 10,
#                                select_freq = 0.75,
#                                     output = "input_egglib" , 
#                                  folderout = "results/egglib-inputs/", 
#                                      input = "output_slim", 
#                                   folderin = "results/slim-outputs/",
#                                     remove = FALSE);
#
## Running EggLib

#ifelse(file.exists(dir('results/egglig-outputs/', pattern = "output_egglib_", full.names = TRUE)), 
#       file.remove(dir('results/egglig-outputs/', pattern = "output_egglib_", full.names = TRUE)), NULL)

#runEggLib(    output = "output_egglib", 
#           folderout = "results/egglig-outputs/", 
#               input = "input_egglib",
#            folderin = "results/egglib-inputs/",
#          pythonpath = '/Users/vitorpavinato/anaconda/bin/python',
#          sourcepath = 'src/summstats.py',
#                 lss = c("He", "Dj", "WCst"),
#                 wss = c("S","thetaW","D", "Da", "ZZ"),
#                 gss = c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"),
#               wspan = 100, select = "all");

#save.image(file = "results/workspaceFile.RData");

#################################################
##                    ABC-RF                   ##
#################################################

#load(file = "results/workspaceFile.RData");

## Creating the Reference Tables

# Global Reference Table - GRT

#globalrf <- createGRT(sim_priors = priors, 
#                          output = "output_egglib", 
#                       outfolder = "results/egglig-outputs/",
#                        gss_list = c("GSS.He", "GSS.Dj","GSS.WCst", "GSS.S", "GSS.thetaW","GSS.D", "GSS.Da","GSS.ZZ"));

#write.table(globalrf, file = "results/ref-tables/globalrf.txt", quote = FALSE, row.names = FALSE)

# Locus-specific Reference Table - LRF

#locusrf <- createLRT(sim_priors = priors,
#                           data = dataconvertion,
#                          input = "input_egglib",
#                       infolder = "results/egglib-inputs/",
#                         output = "output_egglib", 
                      outfolder = "results/egglig-outputs/",
                         remove = FALSE)

#write.table(locusrf, file = "results/ref-tables/locusrf.txt", quote = FALSE, row.names = FALSE)


#colnames(priors) <- c("Sim", "NE", "sampleS", "murate", "neuDom", "neuFitness", 
#                      "benDom", "benFitness", "sampleT1","Time1", "INT12", "sampleT2", "Time2");

#save.image(file = "results/workspaceFile.RData");
