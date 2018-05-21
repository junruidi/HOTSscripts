###################################################################
##                  Main Analytical Pipelines                    ##
##                        04/05/2018                             ##
###################################################################
rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("Data/hr50.rda")
source("HOTSscripts/scripts.R")
library(MASS)

# 1. Standardization ---------------------------------------------------
Y = scale(hr50[,-1],center = T, scale = F)
u = svd(Y)$u
v = svd(Y)$v
s = diag(svd(Y)$d)
w = v %*% solve(s)

Yp = Y %*% w
normalization = function(x){return(x/sqrt(sum(x^2)))}


# 2. PCA (both covariance and correlation) --------------------------------
pca = prcomp(Y, center = F, scale. = F)
phi2 = pca$rotation
sc2 = pca$x
summary(pca)

pca.sd = prcomp(Y, center = F, scale. = T)
phi2.sd = pca.sd$rotation
sc2sd = pca.sd$x
summary(pca.sd)


# 3. third order decomposition --------------------------------------------
# 3.1 HOSVD
# decomposition
moment3 = MGT3(Yp)
hosvd3_v = hoevd(moment3,rank = 32)
phi3 = hosvd3_v$u

# rank
core3_v = hosvd3_v$z
unfold_3_1_v = k_unfold(as.tensor(core3_v),2)@data
svd_3_v = svd(unfold_3_1_v)$d
pct_ev3_v = cumsum(svd_3_v)/sum(svd_3_v)  

rank_ev3_v = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev3_v < i) + 1, i)
  rank_ev3_v = rbind(rank_ev3_v,row.i)
}
names(rank_ev3_v) = c("rank","lambda")
# rank lambda
# 1    2    0.1
# 2    4    0.2
# 3    6    0.3
# 4    8    0.4
# 5   11    0.5
# 6   14    0.6
# 7   18    0.7
# 8   22    0.8
# 9   26    0.9


# convert back to gain interpretability
convert_phi3 = ginv(t(w)) %*% phi3
convert_phi3_norm = apply(convert_phi3, 2, normalization)

# scores
sc3 = Yp %*% phi3 
convert_sc3 = Y %*% convert_phi3_norm


# 3.2 RTPM
# decomposition
order3_yp = Rob_TPM(Yp, L = 30, N = 30, order = 3, p = 32)
phi3_rob = order3_yp$eigenv
ev3_rob = order3_yp$eigenl

# rank
pct_ev3_rob = cumsum(ev3_rob)/sum(ev3_rob)  
rank_ev3_rob = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev3_rob < i) + 1, i)
  rank_ev3_rob = rbind(rank_ev3_rob,row.i)
}
names(rank_ev3_rob) = c("rank","lambda")

# rank lambda
# 1    2    0.1
# 2    4    0.2
# 3    6    0.3
# 4    8    0.4
# 5   11    0.5
# 6   14    0.6
# 7   18   0.7
# 8   22    0.8
# 9   27    0.9


# convert back to gain interpretability
convert_phi3_rob = ginv(t(w)) %*% phi3_rob
convert_phi3_rob_norm = apply(convert_phi3_rob, 2, normalization)

# scores
sc3_rob = Yp %*% phi3_rob 
convert_sc3_rob = Y %*% convert_phi3_rob_norm



# 4. fourth order ---------------------------------------------------------
# 4.1 HOSVD
# decomposition
moment4 = MGT4(Yp)
hosvd4_v = hoevd(moment4,rank = 32)
phi4 = hosvd4_v$u

# rank
core4_v = hosvd4_v$z
unfold_4_1_v = k_unfold(as.tensor(core4_v),1)@data
svd_4_v = svd(unfold_4_1_v)$d
pct_ev4_v = cumsum(svd_4_v)/sum(svd_4_v) 
rank_ev4_v = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev4_v < i) + 1, i)
  rank_ev4_v = rbind(rank_ev4_v,row.i)
}
names(rank_ev4_v) = c("rank","lambda")

# rank lambda
# 1    1    0.1
# 2    2    0.2
# 3    4    0.3
# 4    6    0.4
# 5    8    0.5
# 6   11    0.6
# 7   15    0.7
# 8   20    0.8
# 9   25    0.9



# convert back to gain interpretability
convert_phi4 = ginv(t(w)) %*% phi4
convert_phi4_norm = apply(convert_phi4, 2, normalization)

# scores
sc4 = Yp %*% phi4  
convert_sc4 = Y %*% convert_phi4_norm

# 4.2 RTPM
# decomposition
order4_yp = Rob_TPM(Yp, L = 30, N = 30, order = 4, p = 32)
phi4_rob = order4_yp$eigenv
ev4_rob = order4_yp$eigenl

# rank
pct_ev4_rob= cumsum(ev4_rob)/sum(ev4_rob) 
rank_ev4_rob = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev4_rob < i) + 1, i)
  rank_ev4_rob = rbind(rank_ev4_rob,row.i)
}
names(rank_ev4_rob) = c("rank","lambda")
# rank lambda
# 1    1    0.1
# 2    2    0.2
# 3    4    0.3
# 4    6    0.4
# 5    8    0.5
# 6   11    0.6
# 7   15    0.7
# 8   19    0.8
# 9   25    0.9

# convert back to gain interpretability
convert_phi4_rob = ginv(t(w)) %*% phi4_rob
convert_phi4_rob_norm = apply(convert_phi4_rob, 2, normalization)

# scores 
sc4_rob = Yp %*% phi4_rob 
convert_sc4_rob = Y %*% convert_phi4_rob_norm


# 5. Assemble data --------------------------------------------------------

# 5.1 hosvd phi
hosvd_phi = as.data.frame(cbind(phi2, phi3, phi4))
names(hosvd_phi) = paste0("hosvd_phi",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

# 5.2 hosvd phi  convert and norm
hosvd_phi_convert_norm = as.data.frame(cbind(phi2, convert_phi3_norm, convert_phi4_norm))
names(hosvd_phi_convert_norm) = paste0("hosvd_phi_convert_norm",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

# 5.3 robust power tensor iteration phi
rob_phi = as.data.frame(cbind(phi2, phi3_rob,phi4_rob))
names(rob_phi) = paste0("rob_phi",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

# 5.4 robust power tensor iteration phi convert norm
rob_phi_convert_norm = as.data.frame(cbind(phi2, convert_phi3_rob_norm,convert_phi4_rob_norm))
names(rob_phi_convert_norm) = paste0("rob_phi_convert_norm",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

# 5.5 pc score hosvd
hosvd_sc = as.data.frame(cbind(sc2, sc3, sc4))
names(hosvd_sc) = paste0("hosvd_sc",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))
names(hosvd_sc)[1:32] = paste0("sc2","_",c(1:32))

# 5.6 pc score hosvd convert y
hosvd_sc_y_convert = as.data.frame(cbind(sc2,convert_sc3, convert_sc4))
names(hosvd_sc_y_convert) = paste0("hosvd_sc_y_convert",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))
names(hosvd_sc_y_convert)[1:32] = paste0("sc2","_",c(1:32))

# 5.7 pc score robust power tensor iteration
rob_sc = as.data.frame(cbind(sc2,sc3_rob,sc4_rob))
names(rob_sc) = paste0("rob_sc",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))
names(rob_sc)[1:32] = paste0("sc2","_",c(1:32))

# 5.8 pc score rob convert y
rob_sc_y_convert = as.data.frame(cbind(sc2,convert_sc3_rob, convert_sc4_rob))
names(rob_sc_y_convert) = paste0("rob_sc_y_convert",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))
names(rob_sc_y_convert)[1:32] = paste0("sc2","_",c(1:32))

# # 6. Exploration on scores correlations -----------------------------------
# 
# pdf(file = "results/exploration/hosvd_sc.pdf", width = 28, height = 28)
# par(mar = c(4,5,9,6))
# par(oma = c(1,0,1,0))
# corrplot::corrplot(cor(hosvd_sc[,c(1:96)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
# dev.off()
# pdf(file = "results/exploration/hosvd_sc_convert.pdf", width = 28, height = 28)
# par(mar = c(4,5,9,6))
# par(oma = c(1,0,1,0))
# corrplot::corrplot(cor(hosvd_sc_y_convert[,c(1:96)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
# dev.off()
# pdf(file = "results/exploration/rob_sc.pdf", width = 28, height = 28)
# par(mar = c(4,5,9,6))
# par(oma = c(1,0,1,0))
# corrplot::corrplot(cor(rob_sc[,c(1:96)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
# dev.off()
# pdf(file = "results/exploration/rob_sc_convert.pdf", width = 28, height = 28)
# par(mar = c(4,5,9,6))
# par(oma = c(1,0,1,0))
# corrplot::corrplot(cor(rob_sc_y_convert[,c(1:96)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
# dev.off()
# 
# 
# 
# # 7. Explore eigen vectors ------------------------------------------------
# library(qdap)
# library(timeDate)
# library(lubridate)
# TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
# TIME = beg2char(TIME,":",2)
# 
# pdf("results/exploration/phi.pdf",width = 10,height = 10)
# par(mfrow = c(3,1))
# for(i in 1:32){
#   plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
#   axis(1, at = c(seq(1,32,2)))
#   abline(h = 0,lty = 3)
#   plot(phi3[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
#   lines(phi3_rob[,i], col = "red")
#   axis(1, at = c(seq(1,32,2)))
#   abline(h = 0,lty = 3)
#   plot(phi4[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
#   lines(phi4_rob[,i], col = "red")
#   axis(1, at = c(seq(1,32,2)))
#   abline(h = 0,lty = 3)
# }
# dev.off()
# 
# pdf("results/exploration/phi_convert_norm.pdf",width = 10,height = 10)
# par(mfrow = c(3,1))
# for(i in 1:32){
#   plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
#   axis(1, at = c(seq(1,32,2)),labels = TIME)
#   abline(h = 0,lty = 3)
#   plot(convert_phi3_norm[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
#   lines(convert_phi3_rob_norm[,i], col = "red")
#   axis(1, at = c(seq(1,32,2)),labels = TIME)
#   abline(h = 0,lty = 3)
#   plot(convert_phi4_norm[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
#   lines(convert_phi4_rob_norm[,i], col = "red")
#   axis(1, at = c(seq(1,32,2)),labels = TIME)
#   abline(h = 0,lty = 3)
# }
# dev.off()
# 
save(rob_phi_convert_norm, file = "data/rtpm_eig.rda")
save(hosvd_phi_convert_norm, file = "data/hosvd_eig.rda")


# 8. Create datasets for prediction model ---------------------------------
all = cbind(hosvd_sc,hosvd_sc_y_convert[,-c(1:32)],rob_sc[,-c(1:32)],rob_sc_y_convert[,-c(1:32)])

# load("data/survivaforall.rda")
# all = survival_all[,c(8:295)]
load("data/cov50.rda")
cov50 = subset(cov50, select = c(yr5_mort, Diabetes, Cancer, CHF, CHD, wave, SDMVPSU, SDMVSTRA,NormWts))
survival_all = cbind(cov50,all)


save(survival_all, file = "data/survivaforall.rda")
