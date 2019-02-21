# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

fdrthreshold <- function(x, tau, output_snps = TRUE){
  
  sim = x["sim"]
  ne2 = x["meanNe2"]
  fstfdr = x["True_fdr05"]
  threshold = x["Infe_fdr05"]
  
  
  # load selection coefficients files
  filenames_genome <- paste0("results/pipeline_v5/Cluster_data/pooled_reftable_pods/batch",".",sim,"/slim_output/", 
                             "slim_output_pmuts_t", seq(from=0, to=tau, by=1), "_", sim, ".txt")
  
  datalist_genome <- lapply(filenames_genome, function(x){read.table(file= x, header=T, na.strings = "NA")})
  
  merged_genome <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome)
  
  merged_genome <- data.frame(ID = paste0("chr1",":", merged_genome$POS), merged_genome)
  #merged_genome[which(duplicated(merged_genome$ID)), ]
  
  merged_genome <- merged_genome[!duplicated(merged_genome$ID), ]
  
  # load egglib summary statistics
  filename_egglib_output <- paste0("results/pipeline_v5/Cluster_data/pooled_reftable_pods/batch", ".", sim, "/egglib_output/", 
                                   "egglib_output_sample", "_", sim, ".txt")
  
  data_sumstats <- read.table(file= filename_egglib_output, header=T, na.strings = "NA")
  data_sumstats <- data_sumstats[, c("ID", "LSS.WCst")]
  
  # combined data
  snps_selcoeff <- merged_genome[merged_genome$ID %in% data_sumstats$ID, c("ID", "S")]
  
  snp_fst_table <- merge(snps_selcoeff, data_sumstats, by = 1, all = TRUE)
  
  snp_fst_table <- snp_fst_table[order(-snp_fst_table$LSS.WCst), ]
  
  n_snps = nrow(snp_fst_table)
  
  tp_snps = snp_fst_table[which(snp_fst_table[, "LSS.WCst"] > threshold & ne2*snp_fst_table[, "S"] > 1), ]
  
  n_tp = nrow(tp_snps)
  
  n_fp = nrow(snp_fst_table[which(snp_fst_table[, "LSS.WCst"] > threshold & ne2*snp_fst_table[, "S"] <= 1), ])
  
  n_sel = nrow(snp_fst_table[which(ne2*snp_fst_table[, "S"] > 1), ])
  
  results = data.frame(sim=sim, meanNe2=ne2, total_snps=n_snps, total_sel=n_sel,
                       True_fdr05=fstfdr, Infe_fdr05=threshold,
                       ppv=n_tp/(n_tp+n_fp), fdr= 1-(n_tp/(n_tp+n_fp)))
  
  if (output_snps){
    
    output_folder = paste0("results/pipeline_v5/Cluster_data/pooled_reftable_pods/batch",".",sim,"/fdrthreshold_output/")
    
    if (!file_test("-d", output_folder)){
      dir.create(file.path(output_folder))
    }
    
    write.table(cbind("snp", "ID", "S", "LSS_WCst"), 
                file=paste0(output_folder, "true_positive_snps_",sim,".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    write.table(tp_snps, 
                file=paste0(output_folder, "true_positive_snps_",sim,".txt"), row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  }
  
  return(results)
}

selstrongthreshold <- function(x, tau){
  
  #x=popstrongmsel_applic_data[1,]
  
  sim = as.numeric(x["sim"])
  ne2 = as.numeric(x["meanNe2"])
  strongmsel = as.numeric(x["True_popstrongmsel"])
  threshold = as.numeric(x["Infe_popstrongmsel"])
  
  
  # load selection coefficients files
  filenames_genome <- paste0("results/pipeline_v5/Cluster_data/pooled_reftable_pods/batch",".",sim,"/slim_output/", 
                             "slim_output_pmuts_t", seq(from=0, to=tau, by=1), "_", sim, ".txt")
  
  datalist_genome <- lapply(filenames_genome, function(x){read.table(file= x, header=T, na.strings = "NA")})
  
  merged_genome <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome)
  
  merged_genome <- data.frame(ID = paste0("chr1",":", merged_genome$POS), merged_genome)
  #merged_genome[which(duplicated(merged_genome$ID)), ]
  
  merged_genome <- merged_genome[!duplicated(merged_genome$ID), ]
  
  # load egglib summary statistics
  filename_egglib_output <- paste0("results/pipeline_v5/Cluster_data/pooled_reftable_pods/batch", ".", sim, "/egglib_output/", 
                                   "egglib_output_sample", "_", sim, ".txt")
  
  data_sumstats <- read.table(file= filename_egglib_output, header=T, na.strings = "NA")
  data_sumstats <- data_sumstats[, c("ID", "LSS.WCst")]
  
  # combined data
  snps_selcoeff <- merged_genome[merged_genome$ID %in% data_sumstats$ID, c("ID", "S")]
  
  snp_fst_table <- merge(snps_selcoeff, data_sumstats, by = 1, all = TRUE)
  
  snp_fst_table <- snp_fst_table[order(-snp_fst_table$LSS.WCst), ]
  
  n_snps = nrow(snp_fst_table)
  
  snp_id_threshold <- snp_fst_table[round(threshold*n_snps)+1, "ID"]
  snp_fst_threshold <- snp_fst_table[round(threshold*n_snps)+1, "LSS.WCst"]
  
  snp_fst_table$Ns_test <- ifelse(ne2* snp_fst_table[,"S"] > 1, 1, 0)
  
  if (any(snp_fst_table$Ns_test != 0)){
    if(!all(snp_fst_table$Ns_test == 1)){
      
      pred_fst <- prediction(predictions = snp_fst_table[,"LSS.WCst"], labels = snp_fst_table[, "Ns_test"])
      perf_fst_1 <- performance(pred_fst, "ppv", "fpr")
      perf_fst_2 <- performance(pred_fst, "tnr", "fnr")
      perf_fst_3 <- performance(pred_fst, "auc")
      
      perf_temp_table <- data.frame(fdr=1-perf_fst_1@y.values[[1]], fpr=perf_fst_1@x.values[[1]],
                                    tnr=perf_fst_2@y.values[[1]], fnr=perf_fst_2@x.values[[1]])
      
      perf_snp_fst_table <- cbind(snp_fst_table[1:dim(perf_temp_table)[1], ], perf_temp_table)
      
      if (round(threshold*n_snps)+1 <= nrow(perf_snp_fst_table)){
        perf_snp_fst_out <- perf_snp_fst_table[which(perf_snp_fst_table$ID == snp_id_threshold), ]
      } else {
        perf_snp_fst_out <- perf_snp_fst_table[nrow(perf_snp_fst_table), ]
      }
      
      
      fdr = perf_snp_fst_out$fdr
      fpr = perf_snp_fst_out$fpr
      tnr = perf_snp_fst_out$tnr
      fnr = perf_snp_fst_out$fnr
      auc = perf_fst_3@y.values[[1]]
      
    } else {
      fdr <- as.numeric(NA)
      fpr <- as.numeric(NA)
      tnr <- as.numeric(NA)
      fnr <- as.numeric(NA)
      auc <- as.numeric(NA)
    }
  } else {
    fdr <- as.numeric(NA)
    fpr <- as.numeric(NA)
    tnr <- as.numeric(NA)
    fnr <- as.numeric(NA)
    auc <- as.numeric(NA)
  }
    
  results = data.frame(sim=sim, meanNe2=ne2, total_snps=n_snps,
                       True_popstrongmsel=strongmsel, Infe_popstrongmsel=threshold,
                       fst_threshold=snp_fst_threshold,
                       fdr = fdr, fpr=fpr,
                       tnr = tnr, fnr=fnr,
                       auc = auc)
  
  return(results)
}

#loopPrediction <- function(pods_reftable, reg_ABC_obj, vector_par, global_reftable){
#  
#  pods_prediction = NULL
#  for (p in pods_reftable$sim){
#    
#    pod_sumstats <- pods_reftable[p, -c(1:217)] 
#    pod_sumstats <- pod_sumstats[!sapply(pod_sumstats, function(x) all(is.na(x)))]
#    
#   pods_posterior <- predict(object    = reg_ABC_obj,
#                             obs       = pod_sumstats,
#                              training  = data.frame(vector_par, global_reftable),
#                              quantiles = c(0.025,0.5,0.975),
#                              paral     = T,
#                              ncores    = 28)
#    
#    pods_prediction <- rbind(pods_prediction, cbind(pods_posterior$expectation, pods_posterior$med, 
#                                                    pods_posterior$quantiles[1], pods_posterior$quantiles[3]))
# }
#  colnames(pods_prediction) <- c("expectation" ,"median", "q1", "q3")
#  pods_prediction <- as.data.frame(pods_prediction)
#  return(pods_prediction)
#}