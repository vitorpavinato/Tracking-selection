##########################################################################################
# RUNNING SLiM                                                                           #
# This function runs each SLiM infile and saves each outfile to a .txt file.             #
#                                                                                        #
# slim(parm = parm, seed = random_seed, filename = 'outfile_slim', folderout = 'data/',  #
#                                       infile = 'infile_slim', folderin = 'infiles');   #
#                                                                                        #
#                                                                                        #
# parm = a vector or a  matrix containing the simulated parameters;                      #
# seed = seed generated randomly (e.g using a uniform distribution);                     #   
# filename = prefix of the outfiles;                                                     #
# folderout = path where to save the generated SLiM2 outfiles;                           #
# infile = prefix of the infiles;                                                        #
# folderin = path to the infiles.                                                        #            
##########################################################################################

slim <- function(seed = 845295086, 
                 outfile, folderout, 
                 infile, folderin)
{
    parm <- as.data.frame(seed);
    
    for (i in 1:dim(parm)[1]){
      system2(command = '/usr/local/bin/slim', 
                 args = paste0('-s', ' ', parm$seed[i], ' ', folderin, infile, '_', i, '.slim'), 
               stdout = paste0(folderout, outfile, '_', i, '.txt')); 
    };
};

slimclean <- function(seed = 845295086, 
                      infile, folderin)
{
    parm <- as.data.frame(seed);
  
    for (i in 1:dim(parm)[1]){
      system2(command = '/usr/local/bin/slim', 
                 args = paste0('-s', ' ', parm$seed[i], ' ', folderin, infile, '_', i, '.slim'), 
               stdout = FALSE);
    };
};