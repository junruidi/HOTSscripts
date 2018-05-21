library(pROC)

bootstrap_tt_fwd_select = function(data, ind_vars, outcome, nboot, seed.start){
  inc_vars = c()
  auc_mat_full = data.frame("Variable" = rep(NA_character_,length(ind_vars)),
                            "AUC" = rep(NA_real_,length(ind_vars)),
                            "AIC" = rep(NA_real_,length(ind_vars)),
                            "EPIC" = rep(NA_real_,length(ind_vars)),
                            stringsAsFactors = FALSE)
  dat_boot = data[,c("ID","wave",outcome,ind_vars)]
  ot = dat_boot[,outcome]
  split.data = split(dat_boot, f = ot)
  dat_boot_control = reweight_accel(split.data[[1]])
  dat_boot_case = reweight_accel(split.data[[2]])
  N_control = nrow(dat_boot_control)
  N_case = nrow(dat_boot_case)
  if(N_contro00l <= N_case){
    print("need to check the categorization")
  }
  
  for(i in 1:length(ind_vars)){
    exc_vars = setdiff(ind_vars, inc_vars)
    ret    = cbind.data.frame(exc_vars, NA_real_, NA_real_, NA_real_,stringsAsFactors=FALSE)
    auc_ij = aic_ij = epic_ij = matrix(NA,nrow=nboot,ncol=length(exc_vars))
    for(j in 1:nboot){
      set.seed(j + seed.start)
      
      inx_control = sample(1:N_control, replace=TRUE, size=N_control, prob=dat_boot_control$NormWts)
      inx_case = sample(1:N_case, replace=TRUE, size=N_case, prob=dat_boot_case$NormWts)
      
      
      dat_tmp_control = dat_boot_control[inx_control,]
      dat_tmp_case = dat_boot_case[inx_case,]
      
      train_ind_control = sample(1:N_control, replace=FALSE, size=floor(N_control*0.9))
      train_ind_case = sample(1:N_case, replace=FALSE, size=floor(N_case*0.9))
      
      
      dat_tmp_control_train = dat_tmp_control[train_ind_control,]
      dat_tmp_control_test = dat_tmp_control[-train_ind_control,]
      dat_tmp_case_train = dat_tmp_case[train_ind_case,]
      dat_tmp_case_test = dat_tmp_case[-train_ind_case,]
      
      dat_tmp = rbind(dat_tmp_control_train,dat_tmp_case_train)
      dat_tmp$ID = dat_tmp$wave = NULL
      dat_tmp_pred = rbind(dat_tmp_control_test,dat_tmp_case_test)
      dat_tmp_pred$ID = dat_tmp_pred$wave = NULL
      
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





