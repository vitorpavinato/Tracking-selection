##########################################################################################
# slim(parm = parm, seed = random_seed, filename = 'outfile_slim', folderout = 'data/',  #
#                                       infile = 'infile_slim', folderin = 'infiles');   #
#                                                                                        #
# This function runs each SLiM infile and saves each outfile to a .txt file.             #
#                                                                                        #
# parm = a matrix containing the simulated parameters;                                   #
# seed = seed generated randomly (e.g using a uniform distribution);                     #   
# filename = prefix of the outfiles;                                                     #
# folderout = path where to save the generated SLiM2 outfiles;                           #
# infile = prefix of the infiles;                                                        #
# folderin = path to the infiles.                                                        #            
##########################################################################################

slim <- function(parm, seed, filename, folderout, infile, folderin){
  for (i in 1:length(parm)){
    system2(command = '/usr/local/bin/slim', 
               args = paste0('-s', ' ', seed[i], ' ', folderin, infile, '_N_', parm[i], '.slim'), 
             stdout = paste0(folderout, filename, '_N_', parm[i], '.txt')); 
  };
};

slimclean <- function(parm, seed, infile, folderin){
  for (i in 1:length(parm)){
    system2(command = '/usr/local/bin/slim', 
               args = paste0('-s', ' ', seed[i], ' ', folderin, infile, '_N_', parm[i], '.slim'), 
             stdout = FALSE);
  };
};