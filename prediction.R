###################################################################
##                  Prediction Mortality                         ##
##                        04/05/2018                             ##
###################################################################

rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/survivaforall.rda")
source("HOTSscripts/fwd_select.R")


hosvd_sc = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))
hosvd_sc_ct = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc_y_convert3_",c(1:12)),paste0("hosvd_sc_y_convert4_",c(1:12)))
rob_sc = c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))
rob_sc_ct = c(paste0("sc2_",c(1:12)),paste0("rob_sc_y_convert3_",c(1:12)),paste0("rob_sc_y_convert4_",c(1:12)))


result_hosvd_sc = fwd_select(data = survival_all,ind_vars = hosvd_sc, nboot = 1000,seed.start = 1234,outcome = "CHD")
warnings()
save(result_hosvd_sc,file = "result_hosvd_sc.rda")

result_hosvd_sc_ct = fwd_select(data = survival_all,ind_vars = hosvd_sc_ct, nboot = 1000,seed.start = 1234)
warnings()
save(result_hosvd_sc_ct,file = "result_hosvd_sc_ct.rda")

result_rob_sc = fwd_select(data = survival_all,ind_vars = rob_sc, nboot = 1000,seed.start = 1234)
warnings()
save(result_rob_sc,file = "result_rob_sc.rda")

result_rob_sc_ct = fwd_select(data = survival_all,ind_vars = rob_sc_ct, nboot = 1000,seed.start = 1234)
warnings()
save(result_rob_sc_ct,file = "result_rob_sc_ct.rda")



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
load("data/finalmodel.rda")

# accumulate the data model selection
AUC_fwd = data.frame()
for(i in 1:25){
  x = fwd.selection[[i]][,c(1,2,5)]
  x = rbind(x,c(NA,NA,NA))
  AUC_fwd = rbind(AUC_fwd,x)
}
ref = data.frame(var = unique(na.omit(AUC_fwd$Outcome)), labels = c("5yr Mortality","Diabetes","Cancer","CHF","CHD"))

pdf(file = "results/exploration/model_selection_visulization.pdf",width = 16,height = 20)
par(mfrow = c(5,2))
for(i in c(1,6,11,16,21)){
  outcome.i = unique(fwd.selection[[i]]$Outcome)
  main.i = as.character(ref$labels[which(ref$var == outcome.i)])
  
  x1 = fwd.selection[[i]]$AUC
  x2 = fwd.selection[[i+1]]$AUC
  x3 = fwd.selection[[i+3]]$AUC
  
  lowlim = floor(min(c(x1,x2,x3)) / 0.005) * 0.005
  uplim = ceiling(max(c(x1,x2,x3)) / 0.005) * 0.005
  
  plot(x1,xlim = c(1,max(length(x1),length(x2),length(x3))),ylim = c(lowlim,uplim),type="b",pch=NA,xaxt = "n",
       ylab = "AUC", xlab = "Number of PCs", main = paste(main.i,"Oringal Scores", sep = ": "))
  axis(1, at = c(1:max(length(x1),length(x2),length(x3))))
  text(c(1:length(x1)), x1,labels=c("Model1"),cex=0.75,col="black")
  lines(x2, type = "b",col = "red", pch = NA)
  text(c(1:length(x2)), x2,labels=c("Model2"),cex=0.75,col="red")
  lines(x3, type = "b", col = "blue", pch = NA)
  text(c(1:length(x3)), x3,labels=c("Model3"),cex=0.75,col="blue")
  
  
  x1 = fwd.selection[[i]]$AUC
  x2 = fwd.selection[[i+2]]$AUC
  x3 = fwd.selection[[i+4]]$AUC
  
  lowlim = floor(min(c(x1,x2,x3)) / 0.005) * 0.005
  uplim = ceiling(max(c(x1,x2,x3)) / 0.005) * 0.005
  
  plot(x1,xlim = c(1,max(length(x1),length(x2),length(x3))),ylim = c(lowlim,uplim),type="b",pch=NA,xaxt = "n",
       ylab = "AUC", xlab = "Number of PCs", main = paste(main.i,"Converted Scores", sep = ": "))
  axis(1, at = c(1:max(length(x1),length(x2),length(x3))))
  text(c(1:length(x1)), x1,labels=c("Model1"),cex=0.75,col="black")
  lines(x2, type = "b",col = "red", pch = NA)
  text(c(1:length(x2)), x2,labels=c("Model2"),cex=0.75,col="red")
  lines(x3, type = "b", col = "blue", pch = NA)
  text(c(1:length(x3)), x3,labels=c("Model3"),cex=0.75,col="blue")
  
}
dev.off()
