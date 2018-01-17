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

slim <- function(nsim, 
                 seed = 845295086, 
                 outfile, folderout, 
                 infile, folderin)
{
     	    simnumb <- seq(from = 1, to = nsim, by = 1);
     	     if (nsim > 1){
     	      seed <- as.integer(runif(nsim, 100000000, 900000000));
     	    } else {
     	      seed = seed;
     	    }
     	    
     	    parm <- cbind(simnumb, seed);
     	    
     	      for (i in 1:nsim){
     	            system2(command = '/usr/local/bin/slim', 
     	                       args = paste0('-s', ' ', seed[i], ' ', folderin, infile, '_', simnumb[i], '.slim'), 
     	                     stdout = paste0(folderout, outfile, '_', simnumb[i], '.txt')); 
     	      };
     	    return(parm);
};

slimclean <- function(nsim,
                      seed = 845295086, 
                      infile, folderin)
{
            simnumb <- seq(1,nsim,1);
             if (nsim > 1){
              seed <- as.integer(runif(nsim, 100000000, 900000000));
            } else {
              seed = seed;
            }
  
            parm <- cbind(simnumb, seed);
  
              for (i in 1:nsim){
                    system2(command = '/usr/local/bin/slim', 
                               args = paste0('-s', ' ', seed[i], ' ', folderin, infile, '_', simnumb[i], '.slim'), 
                             stdout = FALSE);
              };
            return(parm);
};

##########################################################################################
# FILE CONVERSION                                                                        #
# This function runs each SLiM infile and saves each outfile to a .txt file.             #
#                                                                                        #
#                                                                                        #            
##########################################################################################

convertEggLib <- function(nsim, 
                          output, folderout, 
                          input, folderin)
{
       		simnumb <- seq(1,nsim,1);
       		 for(i in 1:nsim){
       		 	  f1 <- paste0(input,"_",simnumb[i],"_","mutInfo_t1.txt");
       		 	    inft1 <- read.table(paste0(folderin,f1), header = TRUE);
       		 	  
       		 	  f2 <- paste0(input,"_",simnumb[i],"_","genotypes_t1.txt");
       		 	    gnt1 <- read.table(paste0(folderin,f2), header = FALSE);
       		 	    gnt1 <- t(gnt1);
       		 	  
       		 	    snpindex1 = seq(from=1, to=dim(gnt1)[1], by=1); rownames(gnt1) <- paste(snpindex1);
       		 	    iindex1 = seq(from=1, to=dim(gnt1)[2], by=1); colnames(gnt1) <- paste0("indiv",iindex1,"@pop1", "")
       		 	  
       		 	        sp1 <- cbind(inft1,gnt1)
       		 	  
       		 	  f3 <- paste0(input,"_",simnumb[i],"_","mutInfo_t2.txt")
       		 	    inft2 <- read.table(paste0(folderin,f3), header = TRUE)
       		 	  
       		 	  f4 <- paste0(input,"_",simnumb[i],"_","genotypes_t2.txt")
       		 	    gnt2 <- read.table(paste0(folderin,f4), header = FALSE)
       		 	    gnt2 <- t(gnt2)
       		 	  
       		 	    snpindex2 = seq(from=1, to=dim(gnt2)[1], by=1); rownames(gnt2) <- paste(snpindex2)
       		 	    iindex2 = seq(from=1, to=dim(gnt2)[2], by=1); colnames(gnt2) <- paste0("indiv",iindex2,"@pop2", "")
       		 	  
       		 	        sp2 <- cbind(inft2,gnt2)
       		 	  
       		 	            mrgdata <- merge(sp1, sp2, by.x= c(1,2,3,4,5), by.y = c(1,2,3,4,5), all=TRUE)
       		 	            mrgdata[is.na(mrgdata)] <- "11" # it might fix
       		 	  
       		 	                colnames(mrgdata)[3] <- "status"
       		 	                mrgdata[,3]<- rep("NS", times=dim(mrgdata)[1])
       		 	                
       		 	                colnames(mrgdata)[5] <- "alleles"
       		 	                mrgdata[,5]<- rep("A,T", times=dim(mrgdata)[1])
       		 	                
       		 	                out <- paste0(output,"_",simnumb[i],".txt");
       		 	                write.table(mrgdata, file = paste0(folderout,out), quote=FALSE, sep="\t", row.names = FALSE)
       		 };
};

##########################################################################################
# RUNNING EggLib                                                                         #
# This function runs each dataset in Egglib and saves each output to a .txt file.        #
#                                                                                        #
#                                                                                        #            
##########################################################################################

runEggLib <- function(nsim, 
                      output, folderout, 
                      input, folderin,
                      workdir = paste0(working_dir),
                      pythonpath = '/Users/vitorpavinato/anaconda/bin/python',
                      sourcepath = 'src/summstats.py')
{
       		simnumb <- seq(1,nsim,1);
       		  for(i in 1:nsim){
       		    
       		      args <- c(paste0(workdir,'/', sourcepath),
       		                paste0('input-file=', workdir, '/', folderin, input, '_', simnumb[i], '.txt'),
       		                paste0('output-file=', workdir, '/', folderout, output, '_', simnumb[i], '.txt'),
       		                paste0('LSS=','He,Dj,WCst'),
       		                paste0('WSS=','S,thetaW,D,Da,ZZ'),
       		                paste0('GSS=','S,thetaW,D,Da,ZZ'),
       		                paste0('wspan=','10'),
       		                paste0('select=','all'))
       		      
       		      system2(command = pythonpath, args = args, stdout = TRUE, wait = FALSE);
       		    
       		   } 
  
}