###################################################################
#       Data processing for Activity and Covariance Data         ##
##                             04/05/2018                        ##
###################################################################

#1. created half-hourly level data
rm(list = ls())
load("~/Dropbox/Junrui Di/fragmentation/Review AJE/data/processed/act/act_c.rda")
load("~/Dropbox/Junrui Di/fragmentation/Review AJE/data/processed/act/act_d.rda")

act_c = act_c[,-c(2:4)]
act_d = act_d[,-c(2:4)]

act = rbind(act_c, act_d)

x = act[,-1]
x = act[,c(421:1380)]


d = c(1:960)
d_s = split(d, ceiling(seq_along(d)/30))

hr_sum = data.frame(matrix(NA,nrow = nrow(act),ncol = 32))
names(hr_sum) = paste0("HR",c(1:32))
for( k in 1:length(d_s)){
  ind = d_s[[k]]
  c.mat = x[,ind]
  hr_sum[,k] = rowSums(c.mat)
}

hr = cbind(ID = act[,1],hr_sum)

library(dplyr)
hr = as.data.frame(hr %>% group_by(ID) %>% summarise_each(funs = funs(mean(.,na.rm = T))))

#2. get the activity data with full covariance data
#rm(list = ls())
load("~/Dropbox/Junrui Di/fragmentation/Review AJE/data/processed/covariates/mc_c.rda")
load("~/Dropbox/Junrui Di/fragmentation/Review AJE/data/processed/covariates/mc_d.rda")
mc_c$yr = 34
mc_d$yr = 56

cov = rbind(mc_c,mc_d)
rm(list = c("mc_c","mc_d"))

cov$age = cov$RIDAGEMN/12
cov = subset(cov, select = c(SEQN,age, Male, BMXBMI,diabetes, CHF, stroke, cancer, CHD,
                                MobilityProblem, Mexican, HispanicOther, Black, OtherRace, 
                                LessThanHS, HighSchool, MissingEduc, CurrentSmoker, FormerSmoker, 
                                FormerDrinker, ModerateDrinker, HeavyDrinker, MissingAlcohol,yr,permth_exm,mortstat,SDMVPSU,SDMVSTRA))
cov = na.omit(cov)
cov$include = 1
cov$mort_yr  = cov$permth_exm/12
cov$yr5_mort = as.integer(ifelse(cov$mort_yr <= 5 & cov$mortstat == 1, 1,0))

keep = intersect(cov$SEQN,hr$ID)
cov = subset(cov, SEQN %in% keep)
hr = subset(hr, ID %in% keep)
id50 = cov$SEQN[which(cov$age > 50)]
hr50 = subset(hr, ID %in% id50)
cov50 = subset(cov, SEQN %in% id50)
names(cov)[1] = "ID"
names(cov50)[1] = "ID"

cov50$permth_exm = NULL

# reweight for the subjects
library(nhanesaccel)
include.column = which(names(cov50)=="include")
c1 = subset(cov50,yr==34)
d1 = subset(cov50,yr==56)

c2 = nhanes.accel.reweight(acceldata = c1, wave = 1, seqn.column = 1, 
                           include.column = include.column)
d2 = nhanes.accel.reweight(acceldata = d1, wave = 2, seqn.column = 1, 
                           include.column = include.column)
cov50$wt4yr = c(c2$wtmec2yr_adj/2,d2$wtmec2yr_adj/2)
normalize=function(v){return(v/sum(v))}
cov50$wt4yr_norm= normalize(cov50$wt4yr)*nrow(cov50)
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/data/")
# save(cov, file = "cov.rda")
save(cov50, file = "cov50.rda")
# save(hr, file = "hr.rda")
save(hr50, file = "hr50.rda")

