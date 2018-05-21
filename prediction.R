###################################################################
##                  Prediction Mortality                         ##
##                        04/05/2018                             ##
###################################################################

rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/survivaforall.rda")
source("HOTSscripts/fwd_select.R")

pred = paste0("sc2_",c(1:12))
pred = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))
pred = c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))


result = fwd_select(data = survival_all,ind_vars = hosvd_sc, outcome = "y5r-mort",nboot = 1000,seed.start = 1234,outcome = "CHD")
warnings()




###################################################################
##                      Model Selection                          ##
##                        04/05/2018                             ##
###################################################################

rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")


rm(list = ls())
for(i in 1:5){
  for(j in 1:5){
    nm.i = paste0("model",i,j,".rda")
    load(nm.i)
  }
}
rm(list = c("i","j","nm.i"))
models = ls()

cut_to_model = function(data){
  ind.model = which.min(diff(data$AUC)>0)
  return(data[1:ind.model,])
}
fwd.selection = list()
auc_uinvariate = list()

for(i in 1:length(models)){
  x = get(models[i])
  uni = x$auc_mat_1
  multi = x$auc_mat_full
  
  uni_order = uni[order(uni$colMeans.auc_ij.,decreasing=T),]
  multi.final = cut_to_model(multi)
  
  fwd.selection[[i]] = multi.final
  auc_uinvariate[[i]] = uni_order
}

save(fwd.selection,auc_uinvariate, file = "finalmodel.rda")

for(j in 1:25){
  print(fwd.selection[[j]][nrow(fwd.selection[[j]]),])
}

############
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/finalmodel0518.rda")
for(i in 1:9){
  fwd.selection[[i]]$Order = c(1:nrow(fwd.selection[[i]]))
}
for(i in c(1,4,7)){
  x1 = fwd.selection[[i]]
  x2 = fwd.selection[[i+1]]
  x3 = fwd.selection[[i+2]]
  ori = rbind(x1,x2,x3)
  ori = subset(ori, select = c("Variable","AUC","Outcome","Order"))
  ori$AUC = as.character(round(ori$AUC,4))
  ot.i = unique(ori$Outcome)
  nm.ot.i = paste0(ot.i,"_","ori")
  write.csv(ori, file = paste0("results/final_models/",nm.ot.i,".csv"),row.names = F)
  
}

# accumulate the data model selection
AUC_fwd = data.frame()
for(i in 1:25){
  x = fwd.selection[[i]][,c(1,2,5)]
  x = rbind(x,c(NA,NA,NA))
  AUC_fwd = rbind(AUC_fwd,x)
}
ref = data.frame(var = unique(na.omit(AUC_fwd$Outcome)), labels = c("5yr Mortality","Diabetes","Cancer","CHF","CHD"))
