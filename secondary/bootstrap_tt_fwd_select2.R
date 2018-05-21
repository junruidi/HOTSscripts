library(pROC)
library(nhanesaccel)

bootstrap_tt_fwd_select2 = function(data, ind_vars, outcome, nboot, seed.start){
  inc_vars = c()
  auc_mat_full = data.frame("Variable" = rep(NA_character_,length(ind_vars)),
                            "AUC" = rep(NA_real_,length(ind_vars)),
                            "AIC" = rep(NA_real_,length(ind_vars)),
                            "EPIC" = rep(NA_real_,length(ind_vars)),
                            stringsAsFactors = FALSE)
  dat_boot = data[,c("ID","wave",outcome,ind_vars)]
  ot = dat_boot[,outcome]
  split.data = split(dat_boot, f = ot)
  dat_boot_control = split.data[[1]]
  dat_boot_case = split.data[[2]]
  N_control = nrow(dat_boot_control)
  N_case = nrow(dat_boot_case)
  if(N_control <= N_case){
    print("need to check the categorization")
  }
  
  
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
  
  rm(list = c("train_ind_control","train_ind_case","train_control",
              "test_control","train_case","test_case"))
  
  N_train = nrow(train_dat)
  N_test = nrow(test_dat)
  
  for(i in 1:length(ind_vars)){
    exc_vars = setdiff(ind_vars, inc_vars)
    ret    = cbind.data.frame(exc_vars, NA_real_, NA_real_, NA_real_,stringsAsFactors=FALSE)
    auc_ij = aic_ij = epic_ij = matrix(NA,nrow=nboot,ncol=length(exc_vars))
    for(j in 1:nboot){
      set.seed(j + seed.start)
      
      ## use same bootstrapped data set for each iteration of the forward selection procedure
      inx      = sample(1:N_train, replace=TRUE, size=N_train, prob=train_dat$NormWts)
      inx_pred = sample(1:N_test, replace=TRUE, size=N_test, prob=test_dat$NormWts)
      
      dat_tmp = train_dat[inx,]
      dat_tmp_pred = test_dat[inx_pred,]
      
      for(k in 1:length(exc_vars)){
        var_tmp = exc_vars[k]
        form = paste(c(inc_vars, var_tmp), collapse=" +")
        fit_tmp = glm(as.formula(paste(outcome, form, sep = "~")), data=dat_tmp,family=binomial())
        
        nparam = length(fit_tmp$coefficients)#number of parameters in the model
        LL = logLik(fit_tmp)[1]
        
        aic_ij[j,k]  = round(-2*LL + 2*nparam, 4)
        epic_ij[j,k] = round(-2*LL + 4*nparam, 4)
        
        auc_ij[j,k] = roc(dat_tmp_pred[,outcome],predict(fit_tmp, newdata=dat_tmp_pred,type='response'))$auc[1]
        
        rm(list=c("fit_tmp","form","nparam","var_tmp"))
      }
      if(j %% 200 == 0) print(j)
    }
    auc_j = data.frame( exc_vars, colMeans(auc_ij), colMeans(auc_ij), colMeans(epic_ij), stringsAsFactors = FALSE)
    auc_mat_full[i,] = auc_j[which.max(auc_j[,2]),]
    if(i == 1) {
      auc_mat_1 = auc_j
      auc_mat_1$Outcome = outcome
    }
    
    inc_vars <- c(inc_vars, auc_mat_full[i,1])
    print(paste(i, " Independent variable finished"))
  }
  auc_mat_full$Outcome = outcome
  result = list("auc_mat_full" = auc_mat_full, "auc_mat_1" = auc_mat_1)
  return(result)
}





