
datatable <- read.table("thetaS_Pstrong_class_table.txt")

head(datatable)

error <- (datatable$obs!=datatable$infer)

summary(datatable)

tol<-0.2
local_error <- array(NA,nrow(datatable))
for (i in seq_len(nrow(datatable))){
  
    distance <- abs(datatable$logthetaS_1-datatable$logthetaS_1[i])
  
  # calculate weigths from epachnikov kernel
  nacc <- ceiling(length(distance) * tol)
  ds   <- sort(distance)[nacc]
  weights <- 1 - (distance/ds)^2
  weights[which(weights<0)]<-0

  # calculated weighthed proportion of error
  local_error[i]<-sum(error*weights)/sum(weights)
}
datatable<-data.frame(thetaS=datatable$logthetaS_1,error=local_error)

datatable <- datatable[order(datatable$thetaS),]


plot(datatable$thetaS,
     datatable$error,
     xlab=expression(log[10]*theta[S]),
     ylab="error rate",
     type="l",lwd=2,col="red")
abline(v=0)

