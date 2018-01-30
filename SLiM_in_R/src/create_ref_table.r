##########################################################################################
# FILE CONVERSION                                                                        #
# This function runs each SLiM infile and saves each outfile to a .txt file.             #
#                                                                                        #
# Working here 25/01                                                                                        #            
##########################################################################################

convertEggLib <- function(nsim=1,
                          select_freq = 1, 
                          output, folderout, 
                          input, folderin)
{
  simnumb <- seq(1, nsim, 1)
  out <- vector("list", nsim)
  for(i in 1:nsim){
    f1 <- paste0(input,"_", simnumb[i],"_","mutInfo_t1.txt");
    inft1 <- read.table(paste0(folderin,f1), header = TRUE);
    
    f2 <- paste0(input,"_", simnumb[i],"_","genotypes_t1.txt");
    gnt1 <- read.table(paste0(folderin,f2), header = FALSE);
    gnt1 <- t(gnt1);
    
    snpindex1 = seq(from=1, to=dim(gnt1)[1], by=1); rownames(gnt1) <- paste(snpindex1);
    iindex1 = seq(from=1, to=dim(gnt1)[2], by=1); colnames(gnt1) <- paste0("indiv",iindex1,"@pop1", "")
    
    sp1 <- cbind(inft1,gnt1); sp1 <- sp1[,-4];
    
    f3 <- paste0(input,"_", simnumb[i],"_","mutInfo_t2.txt")
    inft2 <- read.table(paste0(folderin,f3), header = TRUE)
    
    f4 <- paste0(input,"_", simnumb[i],"_","genotypes_t2.txt")
    gnt2 <- read.table(paste0(folderin,f4), header = FALSE)
    gnt2 <- t(gnt2)
    
    snpindex2 = seq(from=1, to=dim(gnt2)[1], by=1); rownames(gnt2) <- paste(snpindex2)
    iindex2 = seq(from=1, to=dim(gnt2)[2], by=1); colnames(gnt2) <- paste0("indiv",iindex2,"@pop2", "")
    
    sp2 <- cbind(inft2,gnt2); sp2 <- sp2[,-4];
    
    rawmrgdata <- merge(sp1, sp2, by.x= c(1,2,3,4,5,6), by.y = c(1,2,3,4,5,6), all=TRUE); 
    mrgdata <- rawmrgdata[,-3]
    colnames(mrgdata)[4] <- "selection"
    mrgdata[is.na(mrgdata)] <- "11"
    
    mn <- matrix(data = NA, nrow = dim(mrgdata)[1], ncol = 1)
    tm <- mrgdata[,-c(1,2,3,4,5)]
    for(j in 1:length(mn)){
      s <- sum(tm[j,] == 11)
      mn[j] <- s
    }
    numcol = dim(tm)[2]
    rm <- ifelse(mn < numcol & mn !=0, TRUE, FALSE)
    newdata <- mrgdata[rm, ]
    
    ss <- newdata$selection
    ss[which(ss != 0 )] <- "Y"
    ss[sample(which(ss ==0), round(select_freq*length(which(ss==0))))] <- "Y"
    ss[which(ss == 0)] <- "N"
    newdata$selection <- ss
    
    out[[i]] <- rawmrgdata[rm, 1:6];
    
    outfile <- paste0(output,"_",simnumb[i],".txt");
    write.table(newdata, file = paste0(folderout,outfile), quote=FALSE, sep="\t", row.names = FALSE)
  };
  return(out);
};

##########################################################################################
# RUNNING EggLib                                                                         #
# This function runs each dataset in Egglib and saves each output to a .txt file.        #
#                                                                                        #
#                                                                                        #            
##########################################################################################

runEggLib <- function(output = "out", folderout = ".", 
                      input = NULL, folderin = NULL,
                      workdir = ".",
                      pythonpath = NULL,
                      sourcepath = NULL,
                      lss = "He", wss = "S", gss = "WCst",
                      wspan = 1, SFS_bins = NULL,
                      select = "all", select_num = NULL, select_freq = NULL
)
{
  
  if (is.null(input))
    stop("You should specify the prefix of the input files")
  if (is.null(folderin))
    stop("You should specify the path where the inputs are")
  if (is.null(pythonpath))
    stop("You should specify the path to python")
  if (is.null(sourcepath))
    stop("You should specify the path to EggLib summstats.py")
  if (is.na(pmatch(sourcepath, "src/summstats.py")))
    stop("Wrong path to EggLib summstats.py! (should be src/summstats.py")
  
  files <- length(list.files(path= paste0(folderin), pattern=glob2rx(paste0(input, "_", "*", ".txt"))))        
  
  for(i in 1:files){
    args <- c(paste0(workdir,'/', sourcepath),
              paste0('input-file=', workdir, '/', folderin, input, '_', i, '.txt'),
              paste0('output-file=', workdir, '/', folderout, output, '_', i, '.txt'),
              paste0('LSS=', paste0(lss, collapse = ",")),
              paste0('WSS=', paste0(wss, collapse = ",")),
              paste0('GSS=', paste0(gss, collapse = ",")),
              paste0('wspan=', wspan),
              paste0('select=', select));
    
    system2(command = pythonpath, args = args, stdout = TRUE, wait = FALSE);
    
  }; 
  
};

createGRT <- function(sim_priors = NULL,
                      output, outfolder,
                      gss_list = "GSS.He")
{
  if (is.null(sim_priors))
    stop("You should specify the data table that contain the prior information")
  
  gss_matrix <- matrix(ncol = length(gss_list), nrow = dim(sim_priors)[1]);
  for(i in 1:dim(sim_priors)[1]){
    f <- paste0(output,"_", i,".txt");
    fdata <- read.table(file = paste0(outfolder,f), header = TRUE, na.strings = "NA");
    gss_matrix[i,] <- as.matrix(fdata[1, gss_list]); 
  };
  colnames(gss_matrix) <- gss_list
  grt <- cbind(sim_priors, gss_matrix);
  return(grt);
};

createLRT <- function(sim_priors = NULL,
                      data = NULL,
                      input, infolder,
                      output, outfolder)
{
  if (is.null(sim_priors))
    stop("You should specify the data table that contain the prior information")
  if (is.null(data))
    stop("You should specify the object that contain the information data converted")
  
  out <- vector("list",dim(sim_priors)[1])
  for(i in 1:dim(sim_priors)[1]){
    bii <- data[[i]][which(data[[i]]$selCoeff != 0.0), ]; 
    if (is.data.frame(bii) && nrow(bii)==0){
      bii <- NA; bai <- NA
    } else {
      bii <- bii; bai <- paste0(bii$chrom, ":", bii$position);
    }  
    
    f1 <- paste0(input,"_",i, ".txt");
    fdata1 <- read.table(paste0(infolder,f1), header = TRUE, na.strings = "NA");
    
    f2 <- paste0(output,"_",i, ".txt");
    fdata2 <- read.table(paste0(outfolder,f2), header = TRUE, na.strings = "NA");
    
    nii <- fdata1[which(fdata1$chrom == "chr2" & fdata1$selection == "Y"), ];
      nis <- sample_n(nii, size = 1); nai <- paste0(nis$chrom, ":", nis$position)
    
    if (is.na(bai)){
        neu <- cbind(sim_priors$sim[i], fdata2[which(fdata2$ID == nai),])
        bAn <- neu;
      } else {
        ben <- cbind(sim_priors$sim[i], fdata2[which(fdata2$ID == bai),])
        neu <- cbind(sim_priors$sim[i], fdata2[which(fdata2$ID == nai),])
        bAn <- rbind(ben,neu);
      }
    rmv <- bAn[grepl("[0-10]", names(bAn))]; drop <- names(rmv)
    bAn <- bAn[,!(names(bAn) %in% drop)]
    
    out[[i]] <- bAn
  };
  do.call(rbind, out)
};