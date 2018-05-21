###################################################################
#       Data processing for Activity and Covariance Data         ##
##                             04/05/2018                        ##
###################################################################

# 1.create half hourly level activity data
rm(list = ls())
load("D:/Dropbox/NHANES Data/Processed/wt10/act.rda")
act10 = na.omit(act10)
wear10 = na.omit(wear10)

x = act10[,-c(1:4)]
x = act10[,c(421:1380)]
d = c(1:960)
d_s = split(d, ceiling(seq_along(d)/30))

hr_sum = data.frame(matrix(NA,nrow = nrow(act10),ncol = 32))
names(hr_sum) = paste0("HR",c(1:32))
for( k in 1:length(d_s)){
  ind = d_s[[k]]
  c.mat = x[,ind]
  hr_sum[,k] = rowSums(c.mat)
}

hr = cbind(ID = act10[,1],hr_sum)

library(dplyr)
hr = as.data.frame(hr %>% group_by(ID) %>% summarise_each(funs = funs(mean(.,na.rm = T))))

rm(list = setdiff(ls(),"hr"))

#2. get the activity data with full covariance data
#rm(list = ls())
load("D:/Dropbox/NHANES Data/Processed/cov.rda")

cov = mc
cov = subset(cov, select = c(SEQN,Age, Gender, BMI,Diabetes, CHF, Stroke, Cancer, CHD,
                             wave,permth_exm,mortstat,SDMVPSU,SDMVSTRA))
cov = na.omit(cov)
cov$mort_yr  = cov$permth_exm/12
cov$yr5_mort = as.integer(ifelse(cov$mort_yr <= 5 & cov$mortstat == 1, 1,0))
library(nhanesdata)

keep = intersect(cov$SEQN,hr$ID)
cov = subset(cov, SEQN %in% keep)
hr = subset(hr, ID %in% keep)
id50 = cov$SEQN[which(cov$Age >= 50)]
hr50 = subset(hr, ID %in% id50)
cov50 = subset(cov, SEQN %in% id50)
cov50 = reweight_accel(cov50)

names(cov)[1] = "ID"
names(cov50)[1] = "ID"

cov50$permth_exm = NULL
cov50$mort_yr = NULL
cov50$wtmec2yr_adj = cov50$wtmec4yr_adj = NULL


setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/data/")
# save(cov, file = "cov.rda")
save(cov50, file = "cov50.rda")
# save(hr, file = "hr.rda")
save(hr50, file = "hr50.rda")