fwd_select = function(data, ind_vars, outcome, nboot, seed.start){
  library("pROC")
  inc_vars = c()
  auc_mat_full = data.frame("Variable" = rep(NA_character_,length(ind_vars)),
                             "AUC" = rep(NA_real_,length(ind_vars)),
                             "AIC" = rep(NA_real_,length(ind_vars)),
                             "EPIC" = rep(NA_real_,length(ind_vars)),
                             stringsAsFactors = FALSE)
  N = nrow(data)
  dat_boot = data[,c(outcome,"wt4yr_norm",ind_vars)]
  for(i in 1:length(ind_vars)){
    exc_vars = setdiff(ind_vars, inc_vars)
    ret    = cbind.data.frame(exc_vars, NA_real_, NA_real_, NA_real_,stringsAsFactors=FALSE)
    auc_ij = aic_ij = epic_ij = matrix(NA,nrow=nboot,ncol=length(exc_vars))
    for(j in 1:nboot){
      set.seed(j + seed.start)
      
      ## use same bootstrapped data set for each iteration of the forward selection procedure
      inx      = sample(1:N, replace=TRUE, size=N, prob=dat_boot$wt4yr_norm)
      inx_pred = sample(1:N, replace=TRUE, size=N, prob=dat_boot$wt4yr_norm)
      
      dat_tmp = dat_boot[inx,]
      dat_tmp_pred = dat_boot[inx_pred,]
      
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
    auc_j = data.frame( exc_vars, colMeans(auc_ij), colMeans(aic_ij), colMeans(epic_ij), stringsAsFactors = FALSE)
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





