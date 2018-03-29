do_sim <- function(sim, nsim, model,
                   mu_rate, mu_min, mu_max, 
                   ne_min, ne_max,
                   slim_output, egglib_input, save_extra_data, extra_output,  
                   python_path, egglib_summstat, egglib_output, remove_files){
  
  # write progress of simulation on screen 
  if (sim==1 | sim%%10==0 | sim==nsim){
    cat(paste(sim,"of",nsim))  
  }
  
  ## SAMPLE FROM PRIOR
  ##-------------------
  
  # Mutation Rate: mu
  if (mu_rate == 0){
    mu <- 1e-7
  } else {
    mu <- rlunif(1, mu_min=1e-8, mu_max=1e-5, base=exp(10))
  }
  
  # Theta 
  theta_min = 4*ne_min*mu
  theta_max = 4*ne_max*mu
  theta <- rlunif(1, min = theta_min, max =theta_max, base = exp(10))
  ne = (theta/(4*mu))
  
  # Simulation Seed - SLiM seed
  sim_seed  <- as.integer(runif(1, 100000000, 900000000));
  
  ## RUM SLiM2
  ##-------------------
  
  # check if the folder exists
  if (!file_test("-d", slim_output)){
    dir.create(file.path(slim_output))
  }
  
  # generate text with slim command
  slim_run <- paste( "/usr/local/bin/slim",
                     "-s", sim_seed,         # seed number
                     paste0("-d simID=", sim),
                     paste0("-d theta=", theta),    # theta 
                     paste0("-d mu=", mu),            # mutation rate mu
                     paste0(model, ".slim")) # modelo
  
  # rum slim on system
  system(slim_run)
  
  ## READ SIMULATED DATA AND CONVERT TO EGGLIB FORMAT
  ##-----------------------------------------------------
  
  ts_1 <- read.csv(file = paste0(slim_output,"slim_output_t1_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  ts_2 <- read.csv(file = paste0(slim_output,"slim_output_t2_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  
  rwm <- merge.data.frame(ts_1, ts_2, by.x = c(1,2,3,4,5,106,107,108,109), 
                          by.y = c(1,2,3,4,5,106,107,108,109), 
                          all = T, sort = T); 
  
  prm <- rwm[, -c(6,7,8,9,110,211)]
  prm[is.na(prm)] <- "11"
  
  tempm <- prm[,-c(1,2,3,4,5)]
  filter1 <- apply(tempm, 1, function(x){sum(x == 11)})
  rm <- filter1 < dim(tempm)[2] & filter1 !=0
  m <- prm[rm, ]
  m <- m[!duplicated(m[ ,1:2]), ]
  
  if (!file_test("-d", egglib_input)){
    dir.create(file.path(egglib_input))
  }
  
  conv_file <- paste0("egglib_input", "_", sim, ".txt");
  write.table(m, file = paste0(egglib_input, conv_file), quote=FALSE, sep="\t", row.names = FALSE)
  
  pre_extradata <- rwm[, c(1,2,6,7,8,9,110,211)]
  extradata <- pre_extradata[rm, ]
  extradata <- extradata[!duplicated(extradata[, 1:2]), ]
  
  if (save_extra_data) {
    
    if (!file_test("-d", extra_output)){
      dir.create(file.path(extra_output))
    }
    
    extra_file <- paste0("extra_output", "_", sim, ".txt");
    write.table(extradata, file = paste0(extra_output, extra_file), quote=FALSE, sep="\t", row.names = FALSE)
  }
  
  ## READ EGGLIB INPUT AND RUN EGGLIB SUMSTAT
  ##----------------------------------------------------- 
  
  # check if the folder exists
  if (!file_test("-d", egglib_output)){
    dir.create(file.path(egglib_output))
  }
  
  # generate text with slim command  
  egglib_run <- paste(python_path,
                      paste0(getwd(), "/", egglib_summstat),
                      paste0("input-file=", egglib_input, conv_file),
                      paste0("output-file=", egglib_output, "egglib_output", "_", sim, ".txt"),
                      paste0("LSS=", paste0(c("He", "Dj", "WCst"), collapse = ",")),
                      paste0("WSS=", paste0(c("S","thetaW","D", "Da", "ZZ"), collapse = ",")),
                      paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"), collapse = ",")),
                      paste0("wspan=", 100),
                      paste0("select=", "all"));
  
  # rum egglib summstat on system
  if(.Platform$OS.type == "unix") {
    slim_run <- paste( ".", egglib_run, sep="")
    system(egglib_run)
  }
  
  ## READ EGGLIB OUTPUT AND CREATE THE REFERENCE TABLE
  ##-----------------------------------------------------
  
  summstats <- read.csv(file = paste0(egglib_output,"egglib_output_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  summstats <- summstats[, unique(names(summstats))]
  colnames(summstats) <- gsub(":", "_", names(summstats))
  
  extradata <- cbind(ID=paste0(extradata$chrom, ":", extradata$position), extradata[, -c(1,2)])
  
  # global reference table
  glb_summstats <- summstats[1 , grepl( "GSS" , unique(names(summstats)) ) ]
  glb_reftable  <- cbind(sim=sim, seed=sim_seed, theta=theta, mu=mu, Ne=ne, glb_summstats)
  
  # locus specific reference table
  bai <- extradata[which(extradata$muSel != 0.0), ] 
  if (is.data.frame(bai) && nrow(bai)==0){
    bai <- NA
  } else {
    bai <- bai$ID;
  } 
  
  nii <- m[which(m$chrom == "chr2"), c("chrom", "position")];
  nis <- sample_n(nii, size = 1); nai <- paste0(nis$chrom, ":", nis$position);
  
  if (is.na(bai)){
    neutral <- cbind(lc_classif="N", 
                     extradata[which(extradata$ID == nai), -1],
                     summstats[which(summstats$ID == nai), -1])
    loc_reftable <- neutral;
  } else {
    benefit <- cbind(lc_classif="A",
                     extradata[which(extradata$ID == bai), -1],
                     summstats[which(summstats$ID == bai), -1])
    neutral <- cbind(lc_classif="N", 
                     extradata[which(extradata$ID == nai), -1],
                     summstats[which(summstats$ID == nai), -1])
    loc_reftable <- rbind(benefit,neutral);
  }
  
  if (remove_files){
    file.remove(paste0(slim_output,"slim_output_t1_", sim, ".txt"))
    file.remove(paste0(slim_output,"slim_output_t2_", sim, ".txt"))
    file.remove(paste0(egglib_input, conv_file))
    file.remove(paste0(egglib_output,"egglib_output_", sim, ".txt"))
  }
  
  #res <- vector(mode = "list", length = 2)
  #res[[1]] <- glb_reftable 
  #res[[2]] <- loc_reftable
  
  #return(res)
  
  res <- suppressWarnings(cbind(glb_reftable, loc_reftable))
  return(res)
  
}
