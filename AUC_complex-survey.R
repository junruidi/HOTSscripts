##########################################################################################
# This script runs logistic models using replicate-weights survey design  
#
# and calculates AUC using 20 seleted cutpoints 
#
# Logistic model needs to be specified to adapt to user's purpose 
##########################################################################################
#
#  Parameters :
#
#    mydata           A data.frame with ID in column named 'id', PSU in column named 'PSU',  
#                     strata in column named 'strata', dependent variable in column 'depvar' 
#                     and independent variable in column 'indvar'
#
################################################################################


library(survey)
options(survey.lonely.psu = "adjust")




# Select cutpoints
cutpts <- function(x, N) {
  k <- length(x)
  b <- k/(N-1)
  xcut <- rep(NA, length=N+1)
  xcut[1] <- 0
  xcut[2] <- x[1]
  for (i in 1: (N-2)) {
    j <- round(i*b)
    xcut[i+2] <- x[j]
  }
  xcut[N+1] <- x[k]
  xcut <- sort(xcut, decreasing=TRUE)
  return(xcut)
}
    

auc <- function(weight, data=mydata) {
    
  sdesign <- svydesign(id=~PSU, weights=weight, strata=~strata, data=mydata, nest=TRUE)

  logistic <- svyglm(depvar ~ indvar, data=mydata, design=sdesign, family=quasibinomial)

  pred <- predict(logistic, newdata=subset(mydata, is.na(mydata$depvar)==FALSE), type="response")
  mydata <- merge(mydata, pred, by=0, all=TRUE) 
  mydata <- mydata[order(mydata$id),]
  
  sdesign <- svydesign(id=~PSU, weights=weight, strata=~strata, data=mydata, nest=TRUE)
  
  cutoff <- unique(sort(mydata[,c("response")], decreasing=FALSE))
  
  
  # Calculate AUC
  auccalc <- function(x, N) {
    
    xcut <- cutpts(x, N)
    
    x1 <- NULL
    x2 <- NULL
    for (j in xcut) {
      a <- svymean(~I(response>j), subset(sdesign, depvar==TRUE), data=mydata, na.rm=TRUE)
      b <- svymean(~I(response>j), subset(sdesign, depvar==FALSE), data=mydata, na.rm=TRUE)
      x1 <- c(x1, a[2])
      x2 <- c(x2, b[2])
    }
    auc <-  sum((x1[1:N-1] + x1[2:(N)])/2 * (x2[2:(N)] - x2[1:N-1]))
    print(auc)
    return(auc)
  }
  auccalc(x=cutoff, N=20)
  
}


#bsweight is the bootstrap replicate-weights
#withReplicates(bsweight, auc)














                                                            