rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/survivaforall.rda")

library(caret) 
library(pROC)

sc2 = c(paste0("sc2_",c(1:12)))
hosvdsc3 = c(paste0("hosvd_sc3_",c(1:12)))
hosvdsc4 = c(paste0("hosvd_sc4_",c(1:12)))
robsc3 = c(paste0("rob_sc3_",c(1:12)))
robsc4 = c(paste0("rob_sc4_",c(1:12)))
survival_all$yr5_mort = ifelse(survival_all$yr5_mort == 1, "Yes","No")
outcomes = c("Cancer","Diabetes","yr5_mort","CHF","CHD")

trctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                      classProbs = TRUE,summaryFunction = twoClassSummary)

model.tune = data.frame()
vimp_test = data.frame()
vimp_insamp = data.frame()

for(i in 1:5){
  print(paste0("begin outcome ",i))
  x1 = scale(survival_all[,c(sc2)])
  x2 = scale(survival_all[,c(sc2,hosvdsc3,hosvdsc4)])
  x3 = scale(survival_all[,c(sc2,robsc3,robsc4)])
  
  x1 = data.frame(class = survival_all[,outcomes[i]],x1)
  x2 = data.frame(class = survival_all[,outcomes[i]],x2)
  x3 = data.frame(class = survival_all[,outcomes[i]],x3)
  
  
  grid_radial = expand.grid(C=2^(-2:7), sigma = seq(from=0.5,to=0.1,by=-0.1) )
  
  
  m.i = data.frame(m = c(1,2,3), roc_test = c(0,0,0), roc_insample = c(0,0,0))
  for(m in 1:3){
    print(paste0("begin features ",m))
    x = get(paste0("x",m))
    set.seed(3033)
    intrain = createDataPartition(y = x$class, p = 0.7,list = F)
    training = x[intrain,]
    testing = x[-intrain,]
    svm_Radial_Grid = train(class ~., data = training, method = "svmRadial",
                            trControl=trctrl,
                            tuneGrid = grid_radial,
                            metric = "ROC")
    test_pred_Radial_Grid = predict(svm_Radial_Grid, newdata = testing,type = "prob")
    
    gc_pROC = roc(response = testing$class, predictor = test_pred_Radial_Grid$No)$auc
    
    svm_Radial_Grid_insample = train(class ~., data = x, method = "svmRadial",
                                     trControl=trctrl,
                                     tuneGrid = grid_radial,
                                     metric = "ROC")
    m.i$roc_test[m] = gc_pROC
    m.i$roc_insample[m] = max(svm_Radial_Grid_insample$results$ROC)
    
    vpt = varImp(svm_Radial_Grid)$importance
    vpt$var = row.names(vpt)
    vpi = varImp(svm_Radial_Grid_insample)$importance
    vpi$var = row.names(vpi)
    vpt$m = vpi$m = m
    vpt$outcome = vpi$outcome = outcomes[i]
    vimp_test = rbind(vimp_test,vpt)
    vimp_insamp = rbind(vimp_insamp,vpi)
  }
  model.tune = rbind(model.tune,m.i)
}

save(model.tune,vimp_test,vimp_insamp, file = "result.rda")

