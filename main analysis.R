rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("Data/hr50.rda")
source("HOTSscripts/scripts.R")
library(MASS)
Y = scale(hr50[,-1],center = T, scale = F)
u = svd(Y)$u
v = svd(Y)$v
s = diag(svd(Y)$d)
w = v %*% solve(s)

Yp = Y %*% w
normalization = function(x){return(x/sqrt(sum(x^2)))}

# 1. 2nd order PC
pca = prcomp(Y, center = F, scale. = F)
phi2 = pca$rotation
sc2 = pca$x
summary(pca)

pca.sd = prcomp(Y, center = F, scale. = T)
phi2.sd = pca.sd$rotation
sc2sd = pca.sd$x
summary(pca.sd)


##########################################################################
# 2. 3rd order

# 2.1 HOSVD

# 2.1.1 PC score
moment3 = MGT3(Yp)
hosvd3_v = hoevd(moment3,rank = 32)
phi3 = hosvd3_v$u

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

#11 50%+, 14 60%+, 17 70%+, 22 80%+, 26 90%+

convert_phi3 = ginv(t(w)) %*% phi3
convert_phi3_norm = apply(convert_phi3, 2, normalization)

sc3 = Yp %*% phi3 
convert_sc3 = Yp %*% convert_phi3
convert_sc3Y = Y %*% convert_phi3

# 2.1.2 MDS score
hosvd3_u = Gram3_hosvd(Yp)


# 2.2  RPTI

# 2.2.1 PC score
order3_yp = Rob_TSM(Yp, L = 30, N = 30, order = 3, p = 32)
phi3_rob = order3_yp$eigenv
ev3_rob = order3_yp$eigenl
pct_ev3_rob = cumsum(ev3_rob)/sum(ev3_rob)  # 1 50%, 14 60%, 18 70%, 

rank_ev3_rob = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev3_rob < i) + 1, i)
  rank_ev3_rob = rbind(rank_ev3_rob,row.i)
}
names(rank_ev3_rob) = c("rank","lambda")


convert_phi3_rob = ginv(t(w)) %*% phi3_rob
convert_phi3_rob_norm = apply(convert_phi3_rob, 2, normalization)

sc3_rob = Yp %*% phi3_rob 
convert_sc3_rob = Yp %*% convert_phi3_rob
convert_sc3Y_rob = Y %*% convert_phi3_rob


# 2.2.2 MDS score
order3_yp_mds = Rob_TSM(t(Yp), L = 30, N = 30, order = 3, p = 32)
rob_mds3 = order3_yp_mds$eigenv


##########################################################################
# 3. 4th order

# 3.1 HOSVD

# 3.1.1 PC score
moment4 = MGT4(Yp)
hosvd4_v = hoevd(moment4,rank = 32)
phi4 = hosvd4_v$u

core4_v = hosvd4_v$z
unfold_4_1_v = k_unfold(as.tensor(core4_v),1)@data
svd_4_v = svd(unfold_4_1_v)$d
pct_ev4_v = cumsum(svd_4_v)/sum(svd_4_v) # 8 50%+, 11 60%+, 15 70%+, 19 80%+, 25 90%+
rank_ev4_v = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev4_v < i) + 1, i)
  rank_ev4_v = rbind(rank_ev4_v,row.i)
}
names(rank_ev4_v) = c("rank","lambda")


convert_phi4 = ginv(t(w)) %*% phi4
convert_phi4_norm = apply(convert_phi4, 2, normalization)

sc4 = Yp %*% phi4  
convert_sc4 = Yp %*% convert_phi4
convert_sc4Y = Y %*% convert_phi4

# 3.1.2 MDS score
hosvd4_u = Gram4_hosvd(Yp)


# 3.2  RPTI

# 3.2.1 PC score
order4_yp = Rob_TSM(Yp, L = 30, N = 30, order = 4, p = 32)
phi4_rob = order4_yp$eigenv
ev4_rob = order4_yp$eigenl

pct_ev4_rob= cumsum(ev4_rob)/sum(ev4_rob) # 8 50%+, 11 60%+, 15 70%+, 19 80%+, 25 90%+
rank_ev4_rob = data.frame()
for(i in seq(0.1,0.9,0.1)){
  row.i = c(sum(pct_ev4_rob < i) + 1, i)
  rank_ev4_rob = rbind(rank_ev4_rob,row.i)
}
names(rank_ev4_rob) = c("rank","lambda")



convert_phi4_rob = ginv(t(w)) %*% phi4_rob
convert_phi4_rob_norm = apply(convert_phi4_rob, 2, normalization)

sc4_rob = Yp %*% phi4_rob 
convert_sc4_rob = Yp %*% convert_phi4_rob
convert_sc4Y_rob = Y %*% convert_phi4_rob


# 3.2.2 MDS score
order4_yp_mds = Rob_TSM(t(Yp), L = 30, N = 30, order = 4, p = 32)
rob_mds4 = order4_yp_mds$eigenv


############################################################################

# 4. Assemble data

# 4.1 hosvd phi

hosvd_phi = as.data.frame(cbind(phi2, phi3, phi4))
names(hosvd_phi) = paste0("hosvd_phi",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# 4.2 hosvd phi convert
hosvd_phi_convert = as.data.frame(cbind(phi2, convert_phi3, convert_phi4))
names(hosvd_phi_convert) = paste0("hosvd_phi_convert",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# 4.3 hosvd phi  convert and norm
hosvd_phi_convert_norm = as.data.frame(cbind(phi2, convert_phi3_norm, convert_phi4_norm))
names(hosvd_phi_convert_norm) = paste0("hosvd_phi_convert_norm",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# 4.4 robust power tensor iteration phi
rob_phi = as.data.frame(cbind(phi2, phi3_rob,phi4_rob))
names(rob_phi) = paste0("rob_phi",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# 4.5 robust power tensor iteration phi convert
rob_phi_convert = as.data.frame(cbind(phi2, convert_phi3_rob,convert_phi4_rob))
names(rob_phi_convert) = paste0("rob_phi_convert",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# 4.6 robust power tensor iteration phi convert norm
rob_phi_convert_norm = as.data.frame(cbind(phi2, convert_phi3_rob_norm,convert_phi4_rob_norm))
names(rob_phi_convert_norm) = paste0("rob_phi_convert_norm",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# 4.7 pc score hosvd
hosvd_sc = as.data.frame(cbind(sc2, sc3, sc4))
names(hosvd_sc) = paste0("hosvd_sc",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

# 4.8 pc score hosvd convert yp
hosvd_sc_yp_convert = as.data.frame(cbind(convert_sc3, convert_sc3))
names(hosvd_sc_yp_convert) = paste0("hosvd_sc_yp_convert",rep(c(3,4),each = 32),"_",rep(c(1:32),2))

# 4.9 pc score hosvd convert y
hosvd_sc_y_convert = as.data.frame(cbind(convert_sc3Y, convert_sc4Y))
names(hosvd_sc_y_convert) = paste0("hosvd_sc_y_convert",rep(c(3,4),each = 32),"_",rep(c(1:32),2))

# 4.10 pc score robust power tensor iteration
rob_sc = as.data.frame(cbind(sc3_rob,sc4_rob))
names(rob_sc) = paste0("rob_sc",rep(c(3,4),each = 32),"_",rep(c(1:32),2))

# 4.11 pc score rob convert yp
rob_sc_yp_convert = as.data.frame(cbind(convert_sc3_rob, convert_sc4_rob))
names(rob_sc_yp_convert) = paste0("rob_sc_yp_convert",rep(c(3,4),each = 32),"_",rep(c(1:32),2))

# 4.12 pc score rob convert y
rob_sc_y_convert = as.data.frame(cbind(convert_sc3Y_rob, convert_sc4Y_rob))
names(rob_sc_y_convert) = paste0("rob_sc_y_convert",rep(c(3,4),each = 32),"_",rep(c(1:32),2))


# 4.13 u score hosvd
hosvd_u = as.data.frame(cbind(u,hosvd3_u$u,hosvd4_u$u))
names(hosvd_u) = paste0("hosvd_u",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

# 4.14 u score robust power tensor iteration
rob_u = as.data.frame(cbind(u,rob_mds3,rob_mds4))
names(rob_u) = paste0("rob_u",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

hosvd_pcscores = cbind(hosvd_sc,hosvd_sc_yp_convert,hosvd_sc_y_convert)
rob_pcscores = cbind(rob_sc,rob_sc_yp_convert,rob_sc_y_convert)
############################################################################
# 5. correlation plot

pdf(file = "Write Up/mar19/correlation/hosvd_sc.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(hosvd_pcscores[,c(1:96)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/mar19/correlation/hosvd_sc_yp_convert.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(hosvd_pcscores[,c(1:32,97:160)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/mar19/correlation/hosvd_sc_y_convert.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(hosvd_pcscores[,c(1:32,161:224)]),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()


x = cbind(hosvd_pcscores[,1:32],rob_pcscores[,c(1:64)])
pdf(file = "Write Up/mar19/correlation/rob_sc.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(x),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

x = cbind(hosvd_pcscores[,1:32],rob_pcscores[,c(65:128)])
pdf(file = "Write Up/mar19/correlation/rob_sc_yp_convert.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(x),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

x = cbind(hosvd_pcscores[,1:32],rob_pcscores[,c(129:192)])
pdf(file = "Write Up/mar19/correlation/rob_sc_y_convert.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(x),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/mar19/correlation/hosvd_u.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(hosvd_u),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/mar19/correlation/rob_u.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(rob_u),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off

################################################################################################

# 6. plot the eigen vectors
library(qdap)
library(timeDate)
library(lubridate)
TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

pdf("Write Up/mar19/vectors/hosvd_phi.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(phi2.sd[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(phi4[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
}
dev.off()

pdf("Write Up/mar19/vectors/hosvd_phi_convert_norm.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(convert_phi3_norm[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(convert_phi4_norm[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
}
dev.off()

pdf("Write Up/mar19/vectors/rob_phi.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(phi3_rob[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(phi4_rob[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
}
dev.off()

pdf("Write Up/mar19/vectors/rob_phi_convert_norm.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(convert_phi3_rob_norm[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(convert_phi4_rob_norm[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
}
dev.off()

save(phi2,convert_phi3_rob_norm,convert_phi4_rob_norm, file = "Data/rtpm_eig.rda")
save(phi2,convert_phi3_norm,convert_phi4_norm, file = "Data/hosvd_eig.rda")
################################################################################################\

# 8. create data for analysis
all = cbind(hosvd_sc,hosvd_sc_yp_convert,rob_sc,rob_sc_yp_convert,hosvd_u,rob_u[,33:96])

load("Data/cov50.rda")
survival_all = cbind(cov50,all)

save(survival_all, file = "Data/survivaforall.rda")
