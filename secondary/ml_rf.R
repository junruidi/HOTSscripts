rm(list = ls())
setwd("ml0612")
load("data/survivaforall.rda")

library(caret)

sc2 = c(paste0("sc2_",c(1:12)))
hosvdsc3 = c(paste0("hosvd_sc3_",c(1:12)))
hosvdsc4 = c(paste0("hosvd_sc4_",c(1:12)))
robsc3 = c(paste0("rob_sc3_",c(1:12)))
robsc4 = c(paste0("rob_sc4_",c(1:12)))
survival_all$yr5_mort = ifelse(survival_all$yr5_mort == 1, "deceased","alive")

outcomes = c("Cancer","Diabetes","yr5_mort","CHF","CHD")


control= trainControl(method="repeatedcv", number=10, search="random",
                      repeats = 5, 
                      summaryFunction=twoClassSummary,classProbs = TRUE)

vImp = data.frame()
model.tune = data.frame()
for(i in 1:5){
  x1 = survival_all[,c(outcomes[i],sc2)]
  x2 = survival_all[,c(outcomes[i],sc2,hosvdsc3,hosvdsc4)]
  x3 = survival_all[,c(outcomes[i],sc2,robsc3,robsc4)]
  names(x1)[1] = names(x2)[1] = names(x3)[1] = "class"
  
  
  rf_random_1 = train(class~., data = x1, metric = "ROC",method="rf", tuneLength=15, trControl=control)
  rf_random_2 = train(class~., data = x2, metric = "ROC",method="rf", tuneLength=15, trControl=control)
  rf_random_3 = train(class~., data = x3, metric = "ROC",method="rf", tuneLength=15, trControl=control)
  
  dat.i = data.frame(outcome = rep(outcomes[i],3),model = c(1,2,3), n.var = c(NA,NA,NA), roc = c(NA,NA,NA))
   for(j in 1:3){
     rf = get(paste0("rf_random_",j))
     dat.i$n.var[j] = rf$bestTune
     dat.i$roc[j] = min(rf$results$ROC)
     varimpj = varImp(rf)$importance
     varimpj$var = row.names(varimpj)
     row.names(varimpj) = NULL
     varimpj$model = j
     varimpj$outcome = outcomes[i]
     vImp = rbind(vImp, varimpj)
   }
  model.tune = rbind(model.tune,dat.i)
}

save(vImp, model.tune)