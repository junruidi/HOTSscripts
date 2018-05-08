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

survival_all = survival_all
scores = survival_all[,8:295]

x1 = model.matrix(~.,scores[,paste0("sc2_",c(1:12))])[,-1]
x2 = model.matrix(~.,scores[,c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))])[,-1]
wt = survival_all$NormWts
# 1. 5yr mortality --------------------------------------------------------
y = survival_all$yr5_mort

fit_y1_x1 = glmnet(x = x1, y = y, standardize = F, 
                         weights = wt, family = "binomial", alpha = 1)
vn = paste0("sc2_",c(1:12))
par(mar=c(4.5,4.5,1,4))
plot(fit_y1_x1,xvar = "lambda")
vnat=coef(fit_y1_x1)
vnat=vnat[-1,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path
axis(2, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.5) 

fit_y1_x1_cv = cv.glmnet(x = x1, y = y, type.measure = "auc",standardize = F,
                         weights = wt, family = "binomial", alpha = 1)
aucs = NULL
for (i in 1:50){
  cv = cv.glmnet(x = x1, y = y, type.measure = "auc",standardize = F,
                  weights = wt, family = "binomial", alpha = 1)
  la.i = data.frame(lambda = cv$lambda, auc = cv$cvm)
  aucs = rbind(aucs,la.i)
}

mean_aucs = aggregate(aucs$auc, by = list(aucs$lambda), mean)
best.lambda = mean_aucs$Group.1[which.max(mean_aucs$x)]
fit_y1_x1_final = glmnet(x = x1, y = y, standardize = F, lambda = best.lambda,
                         weights = wt, family = "binomial", alpha = 1)
coef(fit_y1_x1_final)


fit_y1_x2 = glmnet(x = x2, y = y, standardize = F, 
                   weights = wt, family = "binomial", alpha = 1)
vn = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))
par(mar=c(4.5,4.5,1,4))
plot(fit_y1_x2,xvar = "lambda")
vnat=coef(fit_y1_x2)
vnat=vnat[-1,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path
axis(2, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.5) 
