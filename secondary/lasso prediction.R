####################################################
##           LASSO Prediction Modeling            ##
##               J Di 05/03/2018                  ##
####################################################

rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/survivaforall.rda")
library(glmnet)
library(pROC)
library(nhanesdata)
library(plotmo)

# 1. In-sample: trajecotry of variables -----------------------------------
survival_all = reweight_accel(survival_all)
scores = survival_all[,8:295]
scores = as.data.frame(scale(scores, center = F,scale = T))

x1 = model.matrix(~.,scores[,paste0("sc2_",c(1:12))])[,-1]
x2 = model.matrix(~.,scores[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
x3 = model.matrix(~.,scores[,c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))])[,-1]
wt = survival_all$NormWts



png(file = "results/exploration/lasso_trajectory2.png",units = "in",res = 300, height = 15,width = 12)
par(mfrow = c(5,3))
for(i in 3:7){
  y = survival_all[,i]
  y.name = names(survival_all)[i]
  
  fit_yi_x1 = glmnet(x = x1, y = y, standardize = F, 
                     weights = wt, family = "binomial", alpha = 1)
  fit_yi_x2 = glmnet(x = x2, y = y, standardize = F, 
                     weights = wt, family = "binomial", alpha = 1)
  fit_yi_x3 = glmnet(x = x3, y = y, standardize = F, 
                     weights = wt, family = "binomial", alpha = 1)
  
  # vn1 = c(paste0("sc2_",c(1:12)))
  # vn2 = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))
  # vnat1=coef(fit_yi_x1)
  # vnat1=vnat1[-1,ncol(vnat1)] 
  # vnat2=coef(fit_yi_x2)
  # vnat2=vnat2[-1,ncol(vnat2)] 
  
  plot_glmnet(fit_yi_x1,xvar = "norm",main = y.name, label = 12)
  # axis(2, at=vnat1,line=-.5,label=vn1,las=1,tick=FALSE, cex.axis=0.5) 
  
  plot_glmnet(fit_yi_x2,xvar = "norm",main = y.name,label = 12)
  # axis(2, at=vnat2,line=-.5,label=vn2,las=1,tick=FALSE, cex.axis=0.5) 
  plot_glmnet(fit_yi_x3,xvar = "norm",main = y.name,label = 12)
}

dev.off()


# 2. In-sample: CV --------------------------------------------------------
survival_all = reweight_accel(survival_all)
scores = survival_all[,8:295]
scores = as.data.frame(scale(scores, center = T, scale = T))
x1 = model.matrix(~.,scores[,paste0("sc2_",c(1:12))])[,-1]
x2 = model.matrix(~.,scores[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
wt = survival_all$NormWts

repeated.cv.glmnet = function(x = x, y = y, n = 50){
  aucs = NULL
  for (i in 1:n){
    cv = cv.glmnet(x = x, y = y, type.measure = "auc",standardize = F,
                   weights = wt, family = "binomial", alpha = 1)
    
    la.i = data.frame(lambda = cv$lambda, auc = cv$cvm)
    aucs = rbind(aucs,la.i)
  }
  mean_aucs = aggregate(aucs$auc, by = list(aucs$lambda), mean)
  best.auc = max(mean_aucs$x)
  best.lambda = mean_aucs$Group.1[which.max(mean_aucs$x)]
  return(list("best.auc" = best.auc,"best.lambda" = best.lambda))
}

PC2 = NULL
PC234 = NULL
for(i in 3:7){
  y = survival_all[,i]
  y.name = names(survival_all)[i]

  fit_yi_x1_cv = repeated.cv.glmnet(x = x1, y = y, n = 50)
  fit_y1_x1_final = glmnet(x = x1, y = y, standardize = F, lambda = fit_yi_x1_cv$best.lambda,
                           weights = wt, family = "binomial", alpha = 1)
  
  coef_x1 = as.data.frame(as.matrix(coef(fit_y1_x1_final)))
  coef_x1$var = row.names(coef_x1)
  coef_x1$auc = fit_yi_x1_cv$best.auc
  coef_x1 = coef_x1[-c(1,which(coef_x1$s0 == 0)),]
  coef_x1$outcome = y.name
  coef_x1$lambda = fit_yi_x1_cv$best.lambda
  PC2 = rbind(PC2,coef_x1)
  
  fit_yi_x2_cv = repeated.cv.glmnet(x = x2, y = y, n = 50)
  fit_y1_x2_final = glmnet(x = x2, y = y, standardize = F, lambda = fit_yi_x2_cv$best.lambda,
                           weights = wt, family = "binomial", alpha = 1)
  
  coef_x2 = as.data.frame(as.matrix(coef(fit_y1_x2_final)))
  coef_x2$var = row.names(coef_x2)
  coef_x2$auc = fit_yi_x2_cv$best.auc
  coef_x2 = coef_x2[-c(1,which(coef_x2$s0 == 0)),]
  coef_x2$outcome = y.name
  coef_x2$lambda = fit_yi_x2_cv$best.lambda
  PC234 = rbind(PC234,coef_x2)
}



# 3. Outsample - auc ------------------------------------------------------
i = 3 #3:7
outcome = names(survival_all)[i]
scores = survival_all[,8:295]
scores = as.data.frame(scale(scores, center = F,scale = T))
survival_all[,8:295] = scores

split.data = function(data, p.train, by){
  ot = data[,by]
  split.data = split(data, f = ot)
  dat_boot_control = split.data[[1]]
  dat_boot_case = split.data[[2]]
  N_control = nrow(dat_boot_control)
  N_case = nrow(dat_boot_case)
  
  train_ind_control = sample(1:N_control, replace=FALSE, size=floor(N_control*p.train))
  train_ind_case = sample(1:N_case, replace=FALSE, size=floor(N_case*p.train))
  
  train_control = dat_boot_control[train_ind_control,]
  test_control = dat_boot_control[-train_ind_control,]
  
  train_case = dat_boot_case[train_ind_case,]
  test_case = dat_boot_case[-train_ind_case,]
  
  train_dat = rbind(train_control,train_case)
  test_dat = rbind(test_control,test_case)
  
  return(list("train" = train_dat, "test" = test_dat))
}

repeated.cv.glmnet = function(x = x, y = y, n = 50){
  aucs = NULL
  for (i in 1:n){
    cv = cv.glmnet(x = x, y = y, type.measure = "auc",standardize = F,
                   weights = wt, family = "binomial", alpha = 1)
    coef(cv,s = "lambda.min")
    la.i = data.frame(lambda = cv$lambda, auc = cv$cvm)
    aucs = rbind(aucs,la.i)
  }
  mean_aucs = aggregate(aucs$auc, by = list(aucs$lambda), mean)
  best.auc = max(mean_aucs$x)
  best.lambda = mean_aucs$Group.1[which.max(mean_aucs$x)]
  return(list("best.auc" = best.auc,"best.lambda" = best.lambda))
}

sp = split.data(data = survival_all,p.train = 0.9, by = outcome)

train = sp$train
test = sp$test

train = reweight_accel(train)
test = reweight_accel(test)

train1 = model.matrix(~.,train[,paste0("sc2_",c(1:12))])[,-1]
train2 = model.matrix(~.,train[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
train3 = model.matrix(~.,train[,c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))])[,-1]

test1 = model.matrix(~.,test[,paste0("sc2_",c(1:12))])[,-1]
test2 = model.matrix(~.,test[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
test3 = model.matrix(~.,test[,c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))])[,-1]


y_train = train[,outcome]
y_test = test[,outcome]
wt = train$NormWts

fit_yi_x1_cv = repeated.cv.glmnet(x = train1, y = y_train, n = 50)
fit_y1_x1_final = glmnet(x = train1, y = y_train, standardize = F, lambda = fit_yi_x1_cv$best.lambda,
                         weights = wt, family = "binomial", alpha = 1)
predict_x1 = predict.glmnet(fit_y1_x1_final,newx = test1,s = fit_yi_x1_cv$best.lambda,type = "response")[,1]
auc_1 = roc(y_test,predict_x1)$auc[1]

coef_x1 = as.data.frame(as.matrix(coef(fit_y1_x1_final)))
coef_x1$var = row.names(coef_x1)
coef_x1$auc = auc_1
coef_x1 = coef_x1[-c(1,which(coef_x1$s0 == 0)),]
coef_x1$outcome = outcome
coef_x1$lambda = fit_yi_x1_cv$best.lambda


fit_yi_x2_cv = repeated.cv.glmnet(x = train2, y = y_train, n = 50)
fit_y1_x2_final = glmnet(x = train2, y = y_train, standardize = F, lambda = fit_yi_x2_cv$best.lambda,
                         weights = wt, family = "binomial", alpha = 1)
predict_x2 = predict.glmnet(fit_y1_x2_final,newx = test2,s = fit_yi_x2_cv$best.lambda,type = "response")[,1]
auc_2 = roc(y_test,predict_x2)$auc[1]

coef_x2 = as.data.frame(as.matrix(coef(fit_y1_x2_final)))
coef_x2$var = row.names(coef_x2)
coef_x2$auc = auc_2
coef_x2 = coef_x2[-c(1,which(coef_x2$s0 == 0)),]
coef_x2$outcome = outcome
coef_x2$lambda = fit_yi_x2_cv$best.lambda

fit_yi_x3_cv = repeated.cv.glmnet(x = train3, y = y_train, n = 50)
fit_y1_x3_final = glmnet(x = train3, y = y_train, standardize = F, lambda = fit_yi_x3_cv$best.lambda,
                         weights = wt, family = "binomial", alpha = 1)
predict_x3 = predict.glmnet(fit_y1_x3_final,newx = test3,s = fit_yi_x3_cv$best.lambda,type = "response")[,1]
auc_3 = roc(y_test,predict_x3)$auc[1]

coef_x3 = as.data.frame(as.matrix(coef(fit_y1_x2_final)))
coef_x3$var = row.names(coef_x2)
coef_x3$auc = auc_3
coef_x3 = coef_x3[-c(1,which(coef_x3$s0 == 0)),]
coef_x3$outcome = outcome
coef_x3$lambda = predict_x3$best.lambda