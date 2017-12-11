##########################################################################################
# slim(parm = N_size, seed = random_seed, filename = 'outfile_N_size',                   #
#                         folder = 'data/', infile = 'infile_N_size');                   #
#                                                                                        #
# This function runs each SLiM infile and saves its output to a .txt file.               #
#                                                                                        #
# parm = a vector containing population size N to be simulated;                          #
# seed = seed generated randomly (e.g using a uniform distribution);                     #   
# filename = prefix of the outfiles;                                                     #
# folder = path where to save the generated SLiM2 infiles;                               #
# infile = prefix of the infiles.                                                        #
##########################################################################################

slim <- function(parm, seed, filename, folder, infile){
  for (i in 1:length(parm)){
    system2(command = '/usr/local/bin/slim', 
               args = paste0('-s', ' ', seed[i], ' ', folder, infile, '_', parm[i], '.slim'), 
             stdout = paste0(folder, filename, '_', parm[i], '.txt')); 
  };
};

# This run slim in the terminal and send the stdout to a file specified;
# The output consist of the commands for the simulation and the output specified in late() call;


#run_slim_infile2 <- function(parm, seed, filename, folder, infile){
#  for (i in 1:length(parm)){
#      tmp <- system2(command = '/usr/local/bin/slim', 
#                        args = paste0('-s', ' ', seed[i], ' ', folder, infile, '_', parm[i], '.slim'), 
#                      stdout = TRUE);
#      
#      write(x = as.numeric(tmp[length(tmp)]), 
#            file = paste0(folder, filename, '_', parm[i], '_', 'output', '.txt'));
#  };
#};