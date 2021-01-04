do_data <- function(dataset_path, 
                    dataset_name,
                    pop_prefix,
                    snpsTable2egglib,
                    python_path, 
                    egglib_summstat, 
                    wss_wspan, 
                    sfs_bins,
                    SS1,
                    SS2,
                    snps2egglib_path,
                    snps2egglib_file,
                    egglib2summstats_path,
                    egglib2summstats_file,
                    snptTable2wfabc,
                    wfacbInput_path,
                    wfabcInput_file,
                    tau){
  
  ## HANDLE SLiM2 OUTPUT 3 - GENETIC DATA
  ##-----------------------------------------------------
  
  # file name
  input_data        <- paste0(dataset_path, dataset_name)
  
  # assembly the header for DATA
  header_1       <- c("chrom", "position", "alleles", "selection")
  sample_names_1 <- paste0("indiv", seq(from=1, to=SS1, by=1), "@pop1", "")
  sample_names_2 <- paste0("indiv", seq(from=1, to=SS2, by=1), "@pop2", "")
  
  full_header <- c(header_1, sample_names_1, sample_names_2)
  
  # imported DATA
  input_data_raw <- read.table(file = input_data, header = F, col.names = full_header, check.names = F, na.strings = c("./."))
  
  # remove non bi-allelic SNPs
  #---------------------------------------
  
  # split the allele column
  splited <- strsplit(as.character(input_data_raw$alleles), split= ",", fixed = TRUE)
  
  # count the number of allele of each SNP
  number_alleles <- lengths(splited)
  
  # keep only biallelic markers
  keep_biallelic <- number_alleles == 2
  input_data_biallelic <- input_data_raw[keep_biallelic, ]
  splited_biallelic <- splited[keep_biallelic]
  
  # remove "complex" SNPs
  splited2table <- do.call(rbind, splited_biallelic)
  colnames(splited2table) <- c("REF", "ALT")
  
  ref.simple <- nchar(splited2table[,1]) == 1
  alt.simple <- nchar(splited2table[,2]) == 1
  
  keep_simple <- ref.simple & alt.simple
  input_data_biallelic_simple <- input_data_biallelic[keep_simple, ]
  splited2table_simple <- splited2table[keep_simple,]
  
  # split the data
  input_data_gen <- input_data_biallelic_simple[, 5:(SS1+SS2+4)]
  input_data_snp <- input_data_biallelic_simple[, 1:4]
  
  # change the genotypes annotation
  input_data_gen <- as.matrix(input_data_gen)
  input_data_gen[is.na(input_data_gen)]   <- "00"
  input_data_gen[input_data_gen == "0/0"] <- "11"
  input_data_gen[input_data_gen =="1/1"]  <- "22"
  input_data_gen[input_data_gen == "0/1" | input_data_gen=="1/0"] <- "12"
  
  input_data_gen <- as.data.frame(input_data_gen)
  
  # mark monomophormic mutations (all "11" + "00" or "22")
  count_ref_gen    <- apply(input_data_gen, 1, function(x){all(x == "11" | x == "00")})
  count_alt_gen    <- apply(input_data_gen, 1, function(x){all(x == "22" | x == "00")})
  keep_snps         <- !(count_ref_gen | count_alt_gen) # MARK MONOMORPHIC MUTATIONS
  
  input_data <- cbind(input_data_snp, input_data_gen)
  input_data <- input_data[keep_snps, ]
  
  # make WFABC input file
  if (snptTable2wfabc){
    
    # declare a function that make WFABC input
    countgen4wfabc <- function(input, t_points = 2){
      
      d = input
      g = vector("list", length = t_points)
      n_chrom = NULL
      n_Aalleles = NULL
      for (t in seq(t_points)){
        g[[t]] <- d[grepl(paste0("@pop", t), names(d))]
        n_chromT <- 2 * sum(g[[t]] != "00")
        n_AallelesT <- sum(g[[t]] == "12" | g[[t]] == "21") + 2 * sum(g[[t]] == "22")
        
        n_chrom <- cbind(n_chrom, n_chromT)
        n_Aalleles <- cbind(n_Aalleles, n_AallelesT)
      }
      return(rbind(n_chrom, n_Aalleles))
    }
    
    data2wfabc <- input_data[, -c(1:length(header_1))]
    data2wfabc <- as.data.frame(t(data2wfabc))
    
    wfabc_data  <- do.call(rbind, sapply(data2wfabc, countgen4wfabc, t_points=2, simplify = F))
    
    if (!file_test("-d", wfacbInput_path)){
      dir.create(file.path(wfacbInput_path))
    }
    
    write(paste(dim(data2wfabc)[2], 2),file=paste0(wfacbInput_path, wfabcInput_file)) 
    write(paste(0, (tau), sep = ","),file=paste0(wfacbInput_path, wfabcInput_file), append = TRUE) 
    write.table(wfabc_data, file=paste0(wfacbInput_path, wfabcInput_file),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE) 
  }
   
  # prepare egglib input data
  data_to_egglib <- data.frame(chrom=input_data$chrom, 
                               position=input_data$position, 
                               status="NS",
                               selection=input_data$selection, 
                               alleles=input_data$alleles)
  
  # assembly final egglib input
  data_to_egglib <- cbind(data_to_egglib, input_data[, (length(header_1)+1):ncol(input_data)])
  
  if (snpsTable2egglib){
    
    if (!file_test("-d", snps2egglib_path)){
      dir.create(file.path(snps2egglib_path))
    }
    
    write.table(data_to_egglib, file = paste0(snps2egglib_path, snps2egglib_file), quote=FALSE, sep="\t", row.names = FALSE)
  }
    
  ## READ EGGLIB INPUT AND RUN EGGLIB SUMSTAT
  ##----------------------------------------------------- 
  
  # generate text with slim command  
  egglib_run <- paste(python_path,
                       paste0(getwd(), "/", egglib_summstat),
                       paste0("input-file=", snps2egglib_path, snps2egglib_file),
                       paste0("output-file=", egglib2summstats_path, egglib2summstats_file),
                       paste0("LSS=", paste0(c("He", "Dj", "WCst"), collapse = ",")),
                       paste0("LSSp=", paste0(c("He", "Dj"), collapse = ",")),
                       paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "Pi", "D", "Da", "ZZ", "ZnS"), collapse = ",")),
                       paste0("WSSp=", paste0(c("He", "Dj", "S", "thetaW", "Pi", "D", "ZZ", "ZnS"), collapse = ",")),
                       paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW" , "Pi", "D", "Da", "SFS"), collapse = ",")),
                       paste0("GSSp=",paste0(c("He", "S", "thetaW", "Pi", "D", "Da"), collapse = ",")),
                       paste0("wspan=", wss_wspan),
                       paste0("SFS-bins=", sfs_bins),
                       paste0("select=", "all"));
  
  # rum egglib summstat on system
  if(.Platform$OS.type == "unix") {
    slim_run <- paste( ".", egglib_run, sep="")
    system(egglib_run)
  }
  
  ## READ EGGLIB OUTPUT AND CREATE THE REFERENCE TABLE
  ##-----------------------------------------------------
  
  # import egglib output
  
  summstats_data <- read.csv(file = paste0(egglib2summstats_path, egglib2summstats_file), 
                            header = T, sep = "\t", check.names = F)
  
  
  # remove redundant summary statistics
  summstats_data_noredundant <- summstats_data[, unique(names(summstats_data))]
  
  # rename the summary statistics
  colnames(summstats_data_noredundant) <- gsub(":", "_", names(summstats_data_noredundant))
  
  # egglib calculated GLOBAL statistics
  global_stats <- summstats_data_noredundant[1 , grepl("^GSS" , unique(names(summstats_data_noredundant)))]
  
  global_SFS   <- summstats_data_noredundant[1 , grepl("^SFS" , unique(names(summstats_data_noredundant)))]
  
  # calculate additional GLOBAL summary statistics
  mean_locus_stats <- apply(summstats_data_noredundant[,-c(1, which(grepl("^GSS" , unique(names(summstats_data_noredundant))))[1]:length(summstats_data_noredundant))], 2, function(x){mean(x, na.rm=T)})
  
  var_locus_stats <- apply(summstats_data_noredundant[,-c(1, which(grepl("^GSS" , unique(names(summstats_data_noredundant))))[1]:length(summstats_data_noredundant))], 2, function(x){var(x, na.rm=T)})
  
  kurt_locus_stats <- apply(summstats_data_noredundant[,-c(1, which(grepl("^GSS" , unique(names(summstats_data_noredundant))))[1]:length(summstats_data_noredundant))], 2, function(x){kurtosis(x, na.rm=T)})
  
  skew_locus_stats <- apply(summstats_data_noredundant[,-c(1, which(grepl("^GSS" , unique(names(summstats_data_noredundant))))[1]:length(summstats_data_noredundant))], 2, function(x){skewness(x, na.rm=T)})
  
  q05_locus_stats <- apply(summstats_data_noredundant[,-c(1, which(grepl("^GSS" , unique(names(summstats_data_noredundant))))[1]:length(summstats_data_noredundant))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
  
  q95_locus_stats <- apply(summstats_data_noredundant[,-c(1, which(grepl("^GSS" , unique(names(summstats_data_noredundant))))[1]:length(summstats_data_noredundant))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
  
  # assemble additional GLOBAL summary statistics
  add_global_stats <-cbind(as.data.frame(t(mean_locus_stats)), 
                           as.data.frame(t(var_locus_stats )), 
                           as.data.frame(t(kurt_locus_stats)), 
                           as.data.frame(t(skew_locus_stats)), 
                           as.data.frame(t(q05_locus_stats )), 
                           as.data.frame(t(q95_locus_stats ))
                          )
  
  # ASSEMBLY default GLOBAL summary statistics
  global_summstat_reftable <- cbind(global_stats, global_SFS, add_global_stats)
  
  ## HEADERS
  ##-----------------------------------------------------
  
  header_global_stats <- c("GSS_He","GSS_Dj","GSS_WCst","GSS_S","GSS_thetaW","GSS_Pi","GSS_D","GSS_Da","GSSp_He1","GSSp_He2",
                           "GSSp_S1","GSSp_S2","GSSp_thetaW1","GSSp_thetaW2","GSSp_Pi1","GSSp_Pi2","GSSp_D1","GSSp_D2","GSSp_Da1","GSSp_Da2")
  
  header_global_SFS <- c(paste0("SFSbin",sfs_bins,"_",seq(from=1,to=sfs_bins)))
  
  header_add_global_stats <- c("MEAN_LSS_He","MEAN_LSS_Dj","MEAN_LSS_WCst","MEAN_LSSp_He1","MEAN_LSSp_He2","MEAN_LSSp_Dj1","MEAN_LSSp_Dj2",
                               paste0("MEAN_WSS",wss_wspan,"_He"),paste0("MEAN_WSS",wss_wspan,"_Dj"),paste0("MEAN_WSS",wss_wspan,"_WCst"),
                               paste0("MEAN_WSS",wss_wspan,"_S"),paste0("MEAN_WSS",wss_wspan,"_thetaW"),paste0("MEAN_WSS",wss_wspan,"_Pi"),
                               paste0("MEAN_WSS",wss_wspan,"_D"),paste0("MEAN_WSS",wss_wspan,"_Da"),paste0("MEAN_WSS",wss_wspan,"_ZZ"),
                               paste0("MEAN_WSS",wss_wspan,"_ZnS"),paste0("MEAN_WSSp",wss_wspan,"_He1"),paste0("MEAN_WSSp",wss_wspan,"_He2"),
                               paste0("MEAN_WSSp",wss_wspan,"_Dj1"),paste0("MEAN_WSSp",wss_wspan,"_Dj2"),paste0("MEAN_WSSp",wss_wspan,"_S1"),
                               paste0("MEAN_WSSp",wss_wspan,"_S2"),paste0("MEAN_WSSp",wss_wspan,"_thetaW1"),paste0("MEAN_WSSp",wss_wspan,"_thetaW2"),
                               paste0("MEAN_WSSp",wss_wspan,"_Pi1"),paste0("MEAN_WSSp",wss_wspan,"_Pi2"),paste0("MEAN_WSSp",wss_wspan,"_D1"),
                               paste0("MEAN_WSSp",wss_wspan,"_D2"),paste0("MEAN_WSSp",wss_wspan,"_ZZ1"),paste0("MEAN_WSSp",wss_wspan,"_ZZ2"),
                               paste0("MEAN_WSSp",wss_wspan,"_ZnS1"),paste0("MEAN_WSSp",wss_wspan,"_ZnS2"),
                               "VAR_LSS_He","VAR_LSS_Dj","VAR_LSS_WCst","VAR_LSSp_He1","VAR_LSSp_He2","VAR_LSSp_Dj1","VAR_LSSp_Dj2",
                               paste0("VAR_WSS",wss_wspan,"_He"),paste0("VAR_WSS",wss_wspan,"_Dj"),paste0("VAR_WSS",wss_wspan,"_WCst"),
                               paste0("VAR_WSS",wss_wspan,"_S"),paste0("VAR_WSS",wss_wspan,"_thetaW"),paste0("VAR_WSS",wss_wspan,"_Pi"),
                               paste0("VAR_WSS",wss_wspan,"_D"),paste0("VAR_WSS",wss_wspan,"_Da"),paste0("VAR_WSS",wss_wspan,"_ZZ"),
                               paste0("VAR_WSS",wss_wspan,"_ZnS"),paste0("VAR_WSSp",wss_wspan,"_He1"),paste0("VAR_WSSp",wss_wspan,"_He2"),
                               paste0("VAR_WSSp",wss_wspan,"_Dj1"),paste0("VAR_WSSp",wss_wspan,"_Dj2"),paste0("VAR_WSSp",wss_wspan,"_S1"),
                               paste0("VAR_WSSp",wss_wspan,"_S2"),paste0("VAR_WSSp",wss_wspan,"_thetaW1"),paste0("VAR_WSSp",wss_wspan,"_thetaW2"),
                               paste0("VAR_WSSp",wss_wspan,"_Pi1"),paste0("VAR_WSSp",wss_wspan,"_Pi2"),paste0("VAR_WSSp",wss_wspan,"_D1"),
                               paste0("VAR_WSSp",wss_wspan,"_D2"),paste0("VAR_WSSp",wss_wspan,"_ZZ1"),paste0("VAR_WSSp",wss_wspan,"_ZZ2"),
                               paste0("VAR_WSSp",wss_wspan,"_ZnS1"),paste0("VAR_WSSp",wss_wspan,"_ZnS2"),
                               "KURT_LSS_He","KURT_LSS_Dj","KURT_LSS_WCst","KURT_LSSp_He1","KURT_LSSp_He2","KURT_LSSp_Dj1","KURT_LSSp_Dj2",
                               paste0("KURT_WSS",wss_wspan,"_He"),paste0("KURT_WSS",wss_wspan,"_Dj"),paste0("KURT_WSS",wss_wspan,"_WCst"),
                               paste0("KURT_WSS",wss_wspan,"_S"),paste0("KURT_WSS",wss_wspan,"_thetaW"),paste0("KURT_WSS",wss_wspan,"_Pi"),
                               paste0("KURT_WSS",wss_wspan,"_D"),paste0("KURT_WSS",wss_wspan,"_Da"),paste0("KURT_WSS",wss_wspan,"_ZZ"),
                               paste0("KURT_WSS",wss_wspan,"_ZnS"),paste0("KURT_WSSp",wss_wspan,"_He1"),paste0("KURT_WSSp",wss_wspan,"_He2"),
                               paste0("KURT_WSSp",wss_wspan,"_Dj1"),paste0("KURT_WSSp",wss_wspan,"_Dj2"),paste0("KURT_WSSp",wss_wspan,"_S1"),
                               paste0("KURT_WSSp",wss_wspan,"_S2"),paste0("KURT_WSSp",wss_wspan,"_thetaW1"),paste0("KURT_WSSp",wss_wspan,"_thetaW2"),
                               paste0("KURT_WSSp",wss_wspan,"_Pi1"),paste0("KURT_WSSp",wss_wspan,"_Pi2"),paste0("KURT_WSSp",wss_wspan,"_D1"),
                               paste0("KURT_WSSp",wss_wspan,"_D2"),paste0("KURT_WSSp",wss_wspan,"_ZZ1"),paste0("KURT_WSSp",wss_wspan,"_ZZ2"),
                               paste0("KURT_WSSp",wss_wspan,"_ZnS1"),paste0("KURT_WSSp",wss_wspan,"_ZnS2"),
                               "SKEW_LSS_He","SKEW_LSS_Dj","SKEW_LSS_WCst","SKEW_LSSp_He1","SKEW_LSSp_He2","SKEW_LSSp_Dj1","SKEW_LSSp_Dj2",
                               paste0("SKEW_WSS",wss_wspan,"_He"),paste0("SKEW_WSS",wss_wspan,"_Dj"),paste0("SKEW_WSS",wss_wspan,"_WCst"),
                               paste0("SKEW_WSS",wss_wspan,"_S"),paste0("SKEW_WSS",wss_wspan,"_thetaW"),paste0("SKEW_WSS",wss_wspan,"_Pi"),
                               paste0("SKEW_WSS",wss_wspan,"_D"),paste0("SKEW_WSS",wss_wspan,"_Da"),paste0("SKEW_WSS",wss_wspan,"_ZZ"),
                               paste0("SKEW_WSS",wss_wspan,"_ZnS"),paste0("SKEW_WSSp",wss_wspan,"_He1"),paste0("SKEW_WSSp",wss_wspan,"_He2"),
                               paste0("SKEW_WSSp",wss_wspan,"_Dj1"),paste0("SKEW_WSSp",wss_wspan,"_Dj2"),paste0("SKEW_WSSp",wss_wspan,"_S1"),
                               paste0("SKEW_WSSp",wss_wspan,"_S2"),paste0("SKEW_WSSp",wss_wspan,"_thetaW1"),paste0("SKEW_WSSp",wss_wspan,"_thetaW2"),
                               paste0("SKEW_WSSp",wss_wspan,"_Pi1"),paste0("SKEW_WSSp",wss_wspan,"_Pi2"),paste0("SKEW_WSSp",wss_wspan,"_D1"),
                               paste0("SKEW_WSSp",wss_wspan,"_D2"),paste0("SKEW_WSSp",wss_wspan,"_ZZ1"),paste0("SKEW_WSSp",wss_wspan,"_ZZ2"),
                               paste0("SKEW_WSSp",wss_wspan,"_ZnS1"),paste0("SKEW_WSSp",wss_wspan,"_ZnS2"),
                               "Q05_LSS_He","Q05_LSS_Dj","Q05_LSS_WCst","Q05_LSSp_He1","Q05_LSSp_He2","Q05_LSSp_Dj1","Q05_LSSp_Dj2",
                               paste0("Q05_WSS",wss_wspan,"_He"),paste0("Q05_WSS",wss_wspan,"_Dj"),paste0("Q05_WSS",wss_wspan,"_WCst"),
                               paste0("Q05_WSS",wss_wspan,"_S"),paste0("Q05_WSS",wss_wspan,"_thetaW"),paste0("Q05_WSS",wss_wspan,"_Pi"),
                               paste0("Q05_WSS",wss_wspan,"_D"),paste0("Q05_WSS",wss_wspan,"_Da"),paste0("Q05_WSS",wss_wspan,"_ZZ"),
                               paste0("Q05_WSS",wss_wspan,"_ZnS"),paste0("Q05_WSSp",wss_wspan,"_He1"),paste0("Q05_WSSp",wss_wspan,"_He2"),
                               paste0("Q05_WSSp",wss_wspan,"_Dj1"),paste0("Q05_WSSp",wss_wspan,"_Dj2"),paste0("Q05_WSSp",wss_wspan,"_S1"),
                               paste0("Q05_WSSp",wss_wspan,"_S2"),paste0("Q05_WSSp",wss_wspan,"_thetaW1"),paste0("Q05_WSSp",wss_wspan,"_thetaW2"),
                               paste0("Q05_WSSp",wss_wspan,"_Pi1"),paste0("Q05_WSSp",wss_wspan,"_Pi2"),paste0("Q05_WSSp",wss_wspan,"_D1"),
                               paste0("Q05_WSSp",wss_wspan,"_D2"),paste0("Q05_WSSp",wss_wspan,"_ZZ1"),paste0("Q05_WSSp",wss_wspan,"_ZZ2"),
                               paste0("Q05_WSSp",wss_wspan,"_ZnS1"),paste0("Q05_WSSp",wss_wspan,"_ZnS2"),
                               "Q95_LSS_He","Q95_LSS_Dj","Q95_LSS_WCst","Q95_LSSp_He1","Q95_LSSp_He2","Q95_LSSp_Dj1","Q95_LSSp_Dj2",
                               paste0("Q95_WSS",wss_wspan,"_He"),paste0("Q95_WSS",wss_wspan,"_Dj"),paste0("Q95_WSS",wss_wspan,"_WCst"),
                               paste0("Q95_WSS",wss_wspan,"_S"),paste0("Q95_WSS",wss_wspan,"_thetaW"),paste0("Q95_WSS",wss_wspan,"_Pi"),
                               paste0("Q95_WSS",wss_wspan,"_D"),paste0("Q95_WSS",wss_wspan,"_Da"),paste0("Q95_WSS",wss_wspan,"_ZZ"),
                               paste0("Q95_WSS",wss_wspan,"_ZnS"),paste0("Q95_WSSp",wss_wspan,"_He1"),paste0("Q95_WSSp",wss_wspan,"_He2"),
                               paste0("Q95_WSSp",wss_wspan,"_Dj1"),paste0("Q95_WSSp",wss_wspan,"_Dj2"),paste0("Q95_WSSp",wss_wspan,"_S1"),
                               paste0("Q95_WSSp",wss_wspan,"_S2"),paste0("Q95_WSSp",wss_wspan,"_thetaW1"),paste0("Q95_WSSp",wss_wspan,"_thetaW2"),
                               paste0("Q95_WSSp",wss_wspan,"_Pi1"),paste0("Q95_WSSp",wss_wspan,"_Pi2"),paste0("Q95_WSSp",wss_wspan,"_D1"),
                               paste0("Q95_WSSp",wss_wspan,"_D2"),paste0("Q95_WSSp",wss_wspan,"_ZZ1"),paste0("Q95_WSSp",wss_wspan,"_ZZ2"),
                               paste0("Q95_WSSp",wss_wspan,"_ZnS1"),paste0("Q95_WSSp",wss_wspan,"_ZnS2"))
  
  header_stats_final <- c(header_global_stats,
                          header_global_SFS,
                          header_add_global_stats)
  
  
  if (is.null(pop_prefix)) {
    
    colnames(global_summstat_reftable) <- header_stats_final
   
  } else {
    
    colnames(global_summstat_reftable) <- paste0(header_stats_final, "_", pop_prefix)
    
  }
  
  ## OUTPUT FINAL DATA REFTABLE
  return(global_summstat_reftable)
  
}
