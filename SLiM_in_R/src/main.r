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
library(plyr); # for the convertEggLib function 

#sessionInfo()$otherPkgs$abc$Version;
#sessionInfo()$otherPkgs$abcrf$Version;
sessionInfo()$otherPkgs$plyr$Version;

recordSessionInfo <- sessionInfo();
working_dir <- "/Users/vitorpavinato/Documents/My_repositories/Tracking-selection/SLiM_in_R"

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

## Creating SLiM infiles

# Set the number of simulations

nsim <- 10; 

# Set values for the prior - only for Ne
sim_ne <- as.integer(runif(nsim, 100, 1000));
sim_seed <- as.integer(runif(nsim, 100000000, 900000000));

priors <- as.data.frame(cbind(sim=1:nsim, seed=sim_seed,ne=sim_ne));

# Create SLiM infiles using values taken from the prior

source("src/create_slim_infile.r");

ifelse(file.exists(dir("results/slim-infiles/", pattern = ".slim$", full.names = TRUE)), 
       file.remove(dir("results/slim-infiles/", pattern = ".slim$", full.names = TRUE)), NULL);
      
createslim <- model1(      ne = sim_ne,
                     filename = "infile_slim", 
                       folder = "results/slim-infiles/",
                       output = "output_slim",
                    folderout = "results/slim-outputs");

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
# It runs SLiM2 and saves only the output;

# slimclean(    seed = sim_seed, 
#             infile = 'infile_slim',      
#           folderin = 'infiles/');

## Converting SLiM outputs to EggLib inputs

ifelse(file.exists(dir('results/egglib-inputs/', pattern = "input_egglib_", full.names = TRUE)), 
       file.remove(dir('results/egglib-inputs/', pattern = "input_egglib_", full.names = TRUE)), NULL)

system.time(converteddata <- convertEggLib(     nsim = 10, 
                                              output = "input_egglib" , 
                                           folderout = "results/egglib-inputs/", 
                                               input = "output_slim", 
                                            folderin = "results/slim-outputs/"));

## Running EggLib

ifelse(file.exists(dir('results/egglig-outputs/', pattern = "output_egglib_", full.names = TRUE)), 
       file.remove(dir('results/egglig-outputs/', pattern = "output_egglib_", full.names = TRUE)), NULL)

system.time(runEggLib(      nsim = 10, 
              output = "output_egglib", 
           folderout = "results/egglig-outputs/", 
               input = "input_egglib", 
            folderin = "results/egglib-inputs/"));

#  user  system elapsed 
# 5.255   0.497   5.794

## Creating the Reference Tables

# Global Simulations Parameters

simparameters <- createslim

colnames(priors) <- c("Sim", "NE", "sampleS", "murate", "neuDom", "neuFitness", 
                      "benDom", "benFitness", "sampleT1","Time1", "INT12", "sampleT2", "Time2");

# Global Summary Statistics

system.time(globalss <- createGSS(    nsim = 1, 
                                   inputss = "output_egglib", 
                                  folderin = "data/"));

# user  system elapsed 
# 0.069   0.001   0.070 

# Locus-specific Summary Statistics

system.time(locusss <- createLSS(    nsim = 1,
                                    input = "input_egglib",
                                  inputss = "output_egglib", 
                                 folderin = "data/"));

# user  system elapsed 
# 0.419   0.003   0.423 