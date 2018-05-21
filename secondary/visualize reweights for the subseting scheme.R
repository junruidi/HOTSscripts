# check the weights in the case group
rm(list = ls())
library(nhanesdata)
load("data/survivaforall.rda")

par(mfrow = c(1,3))
for(outcome in c("yr5_mort","diabetes","cancer","CHF","CHD")){
  dat_boot = survival_all[,c("ID","wave","age",outcome)]
  dat_boot = reweight_accel(dat_boot)
  ot = dat_boot[,outcome]
  split.data = split(dat_boot, f = ot)
  dat_boot_control = split.data[[1]]
  dat_boot_case = split.data[[2]]
  
  N_control = nrow(dat_boot_control)
  N_case = nrow(dat_boot_case)
  
  train_ind_control = sample(1:N_control, replace=FALSE, size=floor(N_control*0.9))
  train_ind_case = sample(1:N_case, replace=FALSE, size=floor(N_case*0.9))
  
  train_control = dat_boot_control[train_ind_control,]
  test_control = dat_boot_control[-train_ind_control,]
  
  train_case = dat_boot_case[train_ind_case,]
  test_case = dat_boot_case[-train_ind_case,]
  
  train_dat = rbind(train_control,train_case)
  test_dat = rbind(test_control,test_case)
  
  train_dat = reweight_accel(train_dat)
  test_dat = reweight_accel(test_dat)
  
  plot(dat_boot$age,dat_boot$NormWts, xlab = "age", ylab = "weight",main = paste0(outcome,": oringinal sample, N = ", nrow(dat_boot)))
  plot(train_dat$age, train_dat$NormWts,xlab = "age", ylab = "weight",main = paste0(outcome,": train sample, N = ", nrow(train_dat)))
  plot(test_dat$age,test_dat$NormWts, xlab = "age", ylab = "weight",main = paste0(outcome,": test sample, N = ", nrow(test_dat)))
  
}
