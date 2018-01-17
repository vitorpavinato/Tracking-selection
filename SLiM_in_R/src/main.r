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

#install.packages(c("abc","abcrf","plyr"));
#library(abc);
#library(abcrf);
library(plyr);

#sessionInfo()$otherPkgs$abc$Version;
#sessionInfo()$otherPkgs$abcrf$Version;

recordSessionInfo <- sessionInfo();
working_dir <- "/Users/vitorpavinato/Documents/My_repositories/Tracking-selection/SLiM_in_R"

if (is.na(match(working_dir,getwd()))){
  setwd(working_dir)
} else {
  print("Working directory already set up")
}

save(recordSessionInfo,file="results/sessionInfo.RData");
save.image(file = "results/workspaceFile.RData");

#################################################
##       Runing Forward-time Simulations       ##
#################################################

## Creating SLiM infiles

source("src/create_slim_infile.r");

ifelse(file.exists(dir("infiles/", pattern = ".slim$", full.names = TRUE)), 
       file.remove(dir("infiles/", pattern = ".slim$", full.names = TRUE)), NULL)
      
createslim <- model1(      nsim = 10,
                         ne_min = 100, ne_max = 500,
                         filename = "infile_slim", 
                         folder = "infiles/"); 
print(createslim);

## Running SLiM simulations 

source("src/slim.r");

ifelse(file.exists(dir('data/', pattern = "output_", full.names = TRUE)), 
       file.remove(dir('data/', pattern = "output_", full.names = TRUE)), NULL)

ifelse(file.exists(dir('data/', pattern = "outfile_", full.names = TRUE)), 
       file.remove(dir('data/', pattern = "outfile_", full.names = TRUE)), NULL)

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters of each simulation in an outfile;
# Thanks to Eidos code each simulation now generate 2 files for each sampling;

runslim1 <- slim(     nsim = 10,
                   outfile = 'outfile_slim',
                 folderout = 'data/',
                    infile = 'infile_slim',
                  folderin = 'infiles/');
print(runslim1);

# system time: 10 simulation Ne < 1000
# user  system elapsed 
# 17.502   0.146  17.690

# SLiM2 runs v2:
# It runs SLiM2 and saves only the output;

runslim2 <- slimclean(    nsim = 10, 
                        infile = 'infile_slim',      
                      folderin = 'infiles/');
print(runslim2)

# system time: 10 simulation Ne < 1000
# user  system elapsed 
# 18.895   0.155  19.081

## Converting SLiM outputs to EggLib inputs

convertEggLib(nsim = 10, 
              output = "input_egglib" , 
              folderout = "data/", 
              input = "output_slim", 
              folderin = "data/")

## Running EggLib

runEggLib (nsim = 10, 
           output = "output_egglib", 
           folderout = "data/", 
           input = "input_egglib", 
           folderin = "data/")
