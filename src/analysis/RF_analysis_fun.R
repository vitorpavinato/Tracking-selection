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


loopPrediction <- function(pods_reftable, reg_ABC_obj, vector_par, global_reftable){
  
  pods_prediction = NULL
  for (p in pods_reftable$sim){
    
    pod_sumstats <- pods_reftable[p, -c(1:217)] 
    pod_sumstats <- pod_sumstats[!sapply(pod_sumstats, function(x) all(is.na(x)))]
    
    pods_posterior <- predict(object    = reg_ABC_obj,
                              obs       = pod_sumstats,
                              training  = data.frame(vector_par, global_reftable),
                              quantiles = c(0.025,0.5,0.975),
                              paral     = T,
                              ncores    = 28)
    
    pods_prediction <- rbind(pods_prediction, cbind(pods_posterior$expectation, pods_posterior$med, 
                                                    pods_posterior$quantiles[1], pods_posterior$quantiles[3]))
  }
  colnames(pods_prediction) <- c("expectation" ,"median", "q1", "q3")
  pods_prediction <- as.data.frame(pods_prediction)
  return(pods_prediction)
}