#########################################################################
##        Prediction using nonlinear machine learning techniques       ##
##                            05/01/0218                               ##
#########################################################################

rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/survivaforall.rda")
library(caret)
library(pROC)

# create training and test data set
N = nrow(survival_all)
inx_train = sample(1:N, replace=F, size=floor(N*0.8), prob=survival_all$wt4yr_norm)

training = survival_all[inx_train,]
testing = survival_all[-inx_train,]

# 1. predict 5yr mortality
y5mort_train = ifelse(training$yr5_mort == 1, "deceased", "alive")
y5mort_test = ifelse(testing$yr5_mort == 1, "deceased", "alive")
ctrl_svm = trainControl(method = "repeatedcv", number = 10, repeats = 5)

# 1.1 svm
X_2ndpc_train = training[,paste0("sc2_",c(1:12))]
X_2ndpc_test = testing[,paste0("sc2_",c(1:12))]
svm_model1 = train(x = X_2ndpc_train, y = y5mort_train, trControl = ctrl_svm,
                   method = "svmRadial",preProcess = c("scale"),weights = training$wt4yr_norm)

svm_model1_pred = predict(svm_model1,newdata = X_2ndpc_test)
confusionMatrix(y5mort_test,svm_model1_pred)

X_hosvd_train = training[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))]
X_hosvd_test = testing[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))]
svm_model2 = train(x = X_hosvd_train, y = y5mort_train, trControl = ctrl_svm,
                   method = "svmRadial",preProcess = c("scale"))

svm_model2_pred = predict(svm_model2,newdata = X_hosvd_test)
confusionMatrix(y5mort_test,svm_model2_pred)


