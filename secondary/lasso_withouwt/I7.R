rm(list = ls())
setwd("lasso0509_2")
load("survivaforall.rda")
library(glmnet)
library(nhanesdata)
library(pROC)
library(caret)

i = 5 #3:7
outcome = names(survival_all)[i]
survival_all$yr5_mort = ifelse(survival_all$yr5_mort == 1, "Yes", "No")


repeated.cv.glmnet = function(x = x, y = y, n = 50){
  aucs = NULL
  for (i in 1:n){
    cv = cv.glmnet(x = x, y = y, type.measure = "auc",standardize = T,
                   family = "binomial", alpha = 1)
    
    la.i = data.frame(lambda = cv$lambda, auc = cv$cvm)
    aucs = rbind(aucs,la.i)
  }
  mean_aucs = aggregate(aucs$auc, by = list(aucs$lambda), mean)
  best.auc = max(mean_aucs$x)
  best.lambda = mean_aucs$Group.1[which.max(mean_aucs$x)]
  return(list("best.auc" = best.auc,"best.lambda" = best.lambda))
}


set.seed(1234)
intrain = createDataPartition(y = survival_all[,i], p = 0.7,list = F)
train = survival_all[intrain,]
test = survival_all[-intrain,]

insamp1 = model.matrix(~.,survival_all[,paste0("sc2_",c(1:12))])[,-1]
insamp2 = model.matrix(~.,survival_all[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
insamp3 = model.matrix(~.,survival_all[,c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))])[,-1]


train1 = model.matrix(~.,train[,paste0("sc2_",c(1:12))])[,-1]
train2 = model.matrix(~.,train[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
train3 = model.matrix(~.,train[,c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))])[,-1]

test1 = model.matrix(~.,test[,paste0("sc2_",c(1:12))])[,-1]
test2 = model.matrix(~.,test[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
test3 = model.matrix(~.,test[,c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))])[,-1]

y_train = train[,outcome]
y_test = test[,outcome]
y_all = survival_all[,outcome]

fit_yi_x1_cv = repeated.cv.glmnet(x = train1, y = y_train, n = 50)
fit_y1_x1_final = glmnet(x = train1, y = y_train, standardize = T, lambda = fit_yi_x1_cv$best.lambda,
                         family = "binomial", alpha = 1)
predict_x1 = predict.glmnet(fit_y1_x1_final,newx = test1,s = fit_yi_x1_cv$best.lambda,type = "response")[,1]
auc_1 = roc(y_test,predict_x1)$auc[1]

coef_x1 = as.data.frame(as.matrix(coef(fit_y1_x1_final)))
coef_x1$var = row.names(coef_x1)
coef_x1$auc = auc_1
coef_x1 = coef_x1[-c(1,which(coef_x1$s0 == 0)),]
coef_x1$outcome = outcome
coef_x1$lambda = fit_yi_x1_cv$best.lambda

fit_yi_x1_cv_is = repeated.cv.glmnet(x = insamp1, y = y_all, n = 50)
fit_y1_x1_final_is = glmnet(x = insamp1, y = y_all, standardize = T, lambda = fit_yi_x1_cv_is$best.lambda,
                            family = "binomial", alpha = 1)
predict_x1_is = predict.glmnet(fit_y1_x1_final_is,newx = insamp1,s = fit_yi_x1_cv_is$best.lambda,type = "response")[,1]
auc_1_is = roc(y_all,predict_x1_is)$auc[1]

coef_x1$insampauc = auc_1_is


fit_yi_x2_cv = repeated.cv.glmnet(x = train2, y = y_train, n = 50)
fit_y1_x2_final = glmnet(x = train2, y = y_train, standardize = T, lambda = fit_yi_x2_cv$best.lambda,
                         family = "binomial", alpha = 1)
predict_x2 = predict.glmnet(fit_y1_x2_final,newx = test2,s = fit_yi_x2_cv$best.lambda,type = "response")[,1]
auc_2 = roc(y_test,predict_x2)$auc[1]

coef_x2 = as.data.frame(as.matrix(coef(fit_y1_x2_final)))
coef_x2$var = row.names(coef_x2)
coef_x2$auc = auc_2
coef_x2 = coef_x2[-c(1,which(coef_x2$s0 == 0)),]
coef_x2$outcome = outcome
coef_x2$lambda = fit_yi_x2_cv$best.lambda

fit_yi_x2_cv_is = repeated.cv.glmnet(x = insamp2, y = y_all, n = 50)
fit_y1_x2_final_is = glmnet(x = insamp2, y = y_all, standardize = T, lambda = fit_yi_x2_cv_is$best.lambda,
                            family = "binomial", alpha = 1)
predict_x2_is = predict.glmnet(fit_y1_x2_final_is,newx = insamp2,s = fit_yi_x2_cv_is$best.lambda,type = "response")[,1]
auc_2_is = roc(y_all,predict_x2_is)$auc[1]

coef_x2$insampauc = auc_2_is


fit_yi_x3_cv = repeated.cv.glmnet(x = train3, y = y_train, n = 50)
fit_y1_x3_final = glmnet(x = train3, y = y_train, standardize = T, lambda = fit_yi_x3_cv$best.lambda,
                         family = "binomial", alpha = 1)
predict_x3 = predict.glmnet(fit_y1_x3_final,newx = test3,s = fit_yi_x3_cv$best.lambda,type = "response")[,1]
auc_3 = roc(y_test,predict_x3)$auc[1]

coef_x3 = as.data.frame(as.matrix(coef(fit_y1_x3_final)))
coef_x3$var = row.names(coef_x3)
coef_x3$auc = auc_3
coef_x3 = coef_x3[-c(1,which(coef_x3$s0 == 0)),]
coef_x3$outcome = outcome
coef_x3$lambda = fit_yi_x3_cv$best.lambda

fit_yi_x3_cv_is = repeated.cv.glmnet(x = insamp3, y = y_all, n = 50)
fit_y1_x3_final_is = glmnet(x = insamp3, y = y_all, standardize = T, lambda = fit_yi_x3_cv_is$best.lambda,
                            family = "binomial", alpha = 1)
predict_x3_is = predict.glmnet(fit_y1_x3_final_is,newx = insamp3,s = fit_yi_x3_cv_is$best.lambda,type = "response")[,1]
auc_3_is = roc(y_all,predict_x3_is)$auc[1]
coef_x3$insampauc = auc_3_is


I7_2 = coef_x1
I7_h234 = coef_x2
I7_r234 = coef_x3
save(I7_2, I7_h234, I7_r234, file = "I7.rda")