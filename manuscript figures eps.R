###################################################################
##              Generate figures for the manuscripts             ##
##                            eps                                ##
##                        04/05/2018                             ##
###################################################################


# 1. Univariate characteristics -------------------------------------------
rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/hr50.rda")
library(qdap)
library(lubridate)
library(timeDate)
library(e1071)

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("23:00"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

Y = scale(hr50[,-1],center = T, scale = F)
moment_p = function(x,p){
  return(mean(x^p))
}
Y2nd = apply(Y,2,moment_p,p = 2)
Y3rd = apply(Y,2,moment_p,p = 3)
Y4th = apply(Y,2,moment_p,p = 4)

Y = hr50[,-1]

std2y = apply(scale(hr50[,-1],center = T, scale = T), 2, var)
skwY = apply(Y,2,skewness)
kurY = apply(Y,2,kurtosis) + 3


postscript(file = "results/manuscript/2nd version/univariate.eps",
           width = 11 ,height = 14,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow = c(3,1))
par(mar=c(5,5,4,5))

plot(Y2nd,main = "2nd Order Properties",type = "l",xaxt = "n",ylab = "",yaxt = "n",lwd = 2, xlab = "Time",cex.main = 1.5,cex.lab = 1.7)
axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
a = min(Y2nd)
b = mean(c(min(Y2nd),max(Y2nd)))
c =  max(Y2nd)
ab = mean(c(a,b))
bc = mean(c(b,c))
lab = formatC(c(a,ab,b,bc,c), format = "e", digits = 2)
axis(2, at = c(a,ab,b,bc,c), labels = lab,cex.axis = 1.3)
mtext(side = 2, line = 3, 'Central Moment',cex =  1.1)
par(new = T)
plot(std2y, lty = 4,type = "l",xaxt = "n", ylab = "",yaxt = "n", lwd = 2, xlab = "", ylim = c(0,2))
legend(6,2,legend = c("Standardized Moment","Central Moment (Variance)"), lty = c(4,1), lwd = 2, bty =  "n", horiz = T,cex = 1.7)
axis(side = 4,at = 1, cex.axis = 1.3)
mtext(side = 4, line = 2.9, 'Standardized Moment',cex =  1.2)

plot(Y3rd,main = "3rd Order Properties",type = "l",xaxt = "n", ylab = "", yaxt = "n",lwd = 2, xlab = "Time", cex.lab = 1.7, cex.main = 1.5)
axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
a = min(Y3rd)
b = mean(c(min(Y3rd),max(Y3rd)))
c =  max(Y3rd)
ab = mean(c(a,b))
bc = mean(c(b,c))
lab = formatC(c(a,ab,b,bc,c), format = "e", digits = 2)
axis(2, at = c(a,ab,b,bc,c), labels = lab,cex.axis = 1.3)
mtext(side = 2, line = 3, 'Central Moment',cex =  1.1)
par(new = T)
plot(skwY, lty = 4,type = "l",xaxt = "n", ylab = "",yaxt = "n", lwd = 2, xlab = "", ylim = c(1.9,6))
legend(6,6.1,legend = c("Standardized Moment (Skewness)","Central Moment"), lty = c(4,1), lwd = 2, bty =  "n", horiz = T,cex = 1.7)
axis(side = 4,at = c(2,3,4,5,6), cex.axis = 1.3)
mtext(side = 4, line = 2.9, 'Standardized Moment',cex =  1.2)

plot(Y4th,main = "4th Order Properties",type = "l",xaxt = "n", ylab = "", yaxt = "n",lwd = 2, xlab = "Time", cex.lab = 1.7, cex.main = 1.5)
axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
a = min(Y4th)
b = mean(c(min(Y4th),max(Y4th)))
c =  max(Y4th)
ab = mean(c(a,b))
bc = mean(c(b,c))
lab = formatC(c(a,ab,b,bc,c), format = "e", digits = 2)
axis(2, at = c(a,ab,b,bc,c), labels = lab,cex.axis = 1.3)
mtext(side = 2, line = 3, 'Central Moment',cex =  1.1)
par(new = T)
plot(kurY, lty = 4,type = "l",xaxt = "n", ylab = "",yaxt = "n", lwd = 2, xlab = "", ylim = c(10,75))
legend(6,75,legend = c("Standardized Moment (Kurtosis)","Central Moment"), lty = c(4,1), lwd = 2, bty =  "n", horiz = T,cex = 1.7)
axis(side = 4,at = c(10,30,50,70), cex.axis = 1.3)
mtext(side = 4, line = 2.9, 'Standardized Moment',cex =  1.2)

dev.off()


# 2. Plot eigen vectors ---------------------------------------------------
rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/")
library(timeDate)
library(qdap)
library(lubridate)
TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("23:00"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)
load("data/hosvd_eig.rda")
load("data/rtpm_eig.rda")

# 2.1 HOSVD top4
col = gray.colors(n=3, start = 0.2, end = 0.8, gamma = 2.2, alpha = NULL)
postscript(file = "results/manuscript/2nd version/blackwhite/hosvd_eigen_4.eps",
           width = 15 ,height = 13,horizontal = FALSE, onefile = FALSE, paper = "special")

par(oma = c(4, 1, 1, 1))
par(mfrow = c(2,2))
for(i in 1:4){
  pc.i = cbind(hosvd_phi_convert_norm[,i],hosvd_phi_convert_norm[,i+32],hosvd_phi_convert_norm[,i+64])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1]
    pc.i[,2] = -pc.i[,2]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  if(i == 3){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  if(i == 4){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.5,col = col[3],
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2],lty = 2, lwd = 2,col = col[2])
  lines(pc.i[,3], lty = 4, lwd = 2.5,col = col[1])
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(tilde(Phi)^(3)),expression(tilde(Phi)^(4))), 
       xpd = TRUE, col = rev(col), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = c(1,2,4))
dev.off()

postscript(file = "results/manuscript/2nd version/color/hosvd_eigen_4.eps",
           width = 15 ,height = 13,horizontal = FALSE, onefile = FALSE, paper = "special")

par(oma = c(4, 1, 1, 1))
par(mfrow = c(2,2))
for(i in 1:4){
  pc.i = cbind(hosvd_phi_convert_norm[,i],hosvd_phi_convert_norm[,i+32],hosvd_phi_convert_norm[,i+64])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1]
    pc.i[,2] = -pc.i[,2]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  if(i == 3){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  if(i == 4){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.5,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2],lty = 1, lwd = 2,col = "red")
  lines(pc.i[,3], lty = 1, lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(tilde(Phi)^(3)),expression(tilde(Phi)^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()


# 2.2 RTPM top4
col = gray.colors(n=3, start = 0.2, end = 0.8, gamma = 2.2, alpha = NULL)
postscript(file = "results/manuscript/2nd version/blackwhite/rtpm_eigen_4.eps",
           width = 15 ,height = 13,horizontal = FALSE, onefile = FALSE, paper = "special")

par(oma = c(4, 1, 1, 1))
par(mfrow = c(2,2))
for(i in 1:4){
  pc.i = cbind(rob_phi_convert_norm[,i],rob_phi_convert_norm[,i+32],rob_phi_convert_norm[,i+64])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  if(i == 3){
    pc.i[,1] = -pc.i[,1]
  }
  
  if(i == 4){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.5,col = col[3],
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2],lty = 2, lwd = 2,col = col[2])
  lines(pc.i[,3], lty = 4, lwd = 2.5,col = col[1])
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(tilde(Phi)^(3)),expression(tilde(Phi)^(4))), 
       xpd = TRUE, col = rev(col), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = c(1,2,4))
dev.off()


postscript(file = "results/manuscript/2nd version/color/rtpm_eigen_4.eps",
           width = 15 ,height = 13,horizontal = FALSE, onefile = FALSE, paper = "special")

par(oma = c(4, 1, 1, 1))
par(mfrow = c(2,2))
for(i in 1:4){
  pc.i = cbind(hosvd_phi_convert_norm[,i],hosvd_phi_convert_norm[,i+32],hosvd_phi_convert_norm[,i+64])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1]
    pc.i[,2] = -pc.i[,2]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  if(i == 3){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  if(i == 4){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.5,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2],lty = 1, lwd = 2,col = "red")
  lines(pc.i[,3], lty = 1, lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(tilde(Phi)^(3)),expression(tilde(Phi)^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()

# 3. Prediction Models Comparison -----------------------------------------
rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/finalmodel0518.rda")
col = gray.colors(n=3, start = 0.2, end = 0.5, gamma = 2.2, alpha = NULL)

postscript(file = "results/manuscript/2nd version/blackwhite/model_selection.eps",
           width = 7 ,height = 10,horizontal = FALSE, onefile = FALSE, paper = "special")

par(mfrow = c(3,1))
for(i in c(1,4,7)){
  if(i == 1){
    main.i = "5-Year Mortality"
  }
  if(i == 4){
    main.i = "Diabetes"
  }
  if(i == 7){
    main.i = "Cancer"
  }
  
  x1 = fwd.selection[[i]]$AUC
  x2 = fwd.selection[[i+1]]$AUC
  x3 = fwd.selection[[i+2]]$AUC
  
  lowlim = floor(min(c(x1,x2,x3)) / 0.005) * 0.005
  uplim = ceiling(max(c(x1,x2,x3)) / 0.005) * 0.005
  
  plot(x1,xlim = c(1,max(length(x1),length(x2),length(x3))),ylim = c(lowlim,uplim),type="b",pch=NA,xaxt = "n",col = col[3],
       ylab = "AUC", xlab = "Number of Predictors", main =main.i, cex.main = 1.8, cex.lab = 1.5)
  axis(1, at = c(1:max(length(x1),length(x2),length(x3))), cex = 1.5)
  text(c(1:length(x1)), x1,labels=c("Model1"),cex=0.9,col=col[3])
  lines(x2, type = "b",col = col[2], pch = NA)
  text(c(1:length(x2)), x2,labels=c("Model2"),cex=0.9,col=col[2])
  lines(x3, type = "b", col = col[1], pch = NA)
  text(c(1:length(x3)), x3,labels=c("Model3"),cex=0.9,col=col[1])
}


dev.off()


rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/HOTS/")
load("data/finalmodel0518.rda")



postscript(file = "results/manuscript/2nd version/color/model_selection.eps",
           width = 10 ,height = 10,horizontal = FALSE, onefile = FALSE, paper = "special")

par(mfrow = c(2,1))
for(i in c(4,7)){

  if(i == 4){
    main.i = "Diabetes"
  }
  if(i == 7){
    main.i = "Cancer"
  }
  
  x1 = fwd.selection[[i]]$AUC
  x2 = fwd.selection[[i+1]]$AUC
  x3 = fwd.selection[[i+2]]$AUC
  
  lowlim = floor(min(c(x1,x2,x3)) / 0.005) * 0.005
  uplim = ceiling(max(c(x1,x2,x3)) / 0.005) * 0.005
  
  plot(x1,xlim = c(1,max(length(x1),length(x2),length(x3))),ylim = c(lowlim,uplim),type="b",pch=NA,xaxt = "n",
       ylab = "AUC", xlab = "Number of Predictors", main =main.i, cex.main = 1.8, cex.lab = 1.5)
  axis(1, at = c(1:max(length(x1),length(x2),length(x3))), cex = 1.5)
  text(c(1:length(x1)), x1,labels=c("Model1"),cex=0.9,col="black")
  lines(x2, type = "b",col = "red", pch = NA)
  text(c(1:length(x2)), x2,labels=c("Model2"),cex=0.9,col="red")
  lines(x3, type = "b", col = "blue", pch = NA)
  text(c(1:length(x3)), x3,labels=c("Model3"),cex=0.9,col="blue")
}


dev.off()


