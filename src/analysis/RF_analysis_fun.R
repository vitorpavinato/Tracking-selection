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
  threshold = x["threshold"]
  fstfdr = x[3]
  
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
  
  results = data.frame(sim=sim, n_snps=n_snps, fstfdr, threshold=threshold,
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