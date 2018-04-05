rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
library(timeDate)
library(qdap)
library(lubridate)


#1. hosvd
rm(list = ls())

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("23:00"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)
load("Data/hosvd_eig.rda")

# 1.1 12
png(file = "Write Up/manuscript/hosvd_eigen_12.png", 
    units = "in",width = 17 ,height = 16,res = 300)
par(oma = c(4, 1, 1, 1))
par(mfrow = c(4,3))
for(i in 1:12){
pc.i = cbind(phi2[,i],convert_phi3_norm[,i],convert_phi4_norm[,i])

if(i == 1){
  pc.i[,1] = -pc.i[,1] 
  pc.i[,3] = -pc.i[,3] 
}

min.i = min(c(pc.i))
max.i = max(c(pc.i))

plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
     type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.6,
     xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
abline(h = 0,lty = 3)
lines(pc.i[,2], type = "l",lwd = 2,col = "red")
lines(pc.i[,3], type = "l",lwd = 2,col = "blue")

}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(Phi^(3)),expression(Phi^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()



# 1.2 all
pdf(file = "Write Up/manuscript/hosvd_eigen_all.pdf", 
    width = 24 ,height = 16)
par(oma = c(4, 1, 1, 1))
par(mfrow = c(4,4))
for(i in 1:32){
  pc.i = cbind(phi2[,i],convert_phi3_norm[,i],convert_phi4_norm[,i])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1] 
    pc.i[,3] = -pc.i[,3] 
  }
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.6,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2], type = "l",lwd = 2,col = "red")
  lines(pc.i[,3], type = "l",lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(Phi^(3)),expression(Phi^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()


# 1.3 first 2
png(file = "Write Up/manuscript/hosvd_eigen_2.png", 
    units = "in",width = 18 ,height = 7,res = 300)
par(oma = c(4, 1, 1, 1))
par(mfrow = c(1,2))
for(i in 1:2){
  pc.i = cbind(phi2[,i],convert_phi3_norm[,i],convert_phi4_norm[,i])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1]
    pc.i[,3] = -pc.i[,3]
  }
  
  if(i == 2){
    pc.i[,3] = -pc.i[,3]
  }

  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.6,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2], type = "l",lwd = 2,col = "red")
  lines(pc.i[,3], type = "l",lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(tilde(Phi)^(3)),expression(tilde(Phi)^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()



#2. rtpm
rm(list = ls())

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("23:00"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)
load("Data/rtpm_eig.rda")

# 1.1 12
png(file = "Write Up/manuscript/rptm_eigen_12.png", 
    units = "in",width = 17 ,height = 16,res = 300)
par(oma = c(4, 1, 1, 1))
par(mfrow = c(4,3))
for(i in 1:12){
  pc.i = cbind(phi2[,i],convert_phi3_rob_norm[,i],convert_phi4_rob_norm[,i])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1] 
  }
  
  if(i == 3){
    pc.i[,3] = -pc.i[,3] 
  }
  
  
  if(i == 4){
    pc.i[,3] = -pc.i[,3] 
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.6,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2], type = "l",lwd = 2,col = "red")
  lines(pc.i[,3], type = "l",lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(Phi^(3)),expression(Phi^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()

# 2.2 all
pdf(file = "Write Up/manuscript/rtpm_eigen_all.pdf", 
    width = 24 ,height = 16)
par(oma = c(4, 1, 1, 1))
par(mfrow = c(4,4))
for(i in 1:32){
  pc.i = cbind(phi2[,i],convert_phi3_rob_norm[,i],convert_phi4_rob_norm[,i])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1] 
  }
  
  if(i == 3){
    pc.i[,3] = -pc.i[,3] 
  }
  
  
  if(i == 4){
    pc.i[,3] = -pc.i[,3] 
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.6,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2], type = "l",lwd = 2,col = "red")
  lines(pc.i[,3], type = "l",lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(Phi^(3)),expression(Phi^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()


# 1.3 first 2
png(file = "Write Up/manuscript/rtpm_eigen_2.png", 
    units = "in",width = 18 ,height = 7,res = 300)
par(oma = c(4, 1, 1, 1))
par(mfrow = c(1,2))
for(i in 1:2){
  pc.i = cbind(phi2[,i],convert_phi3_rob_norm[,i],convert_phi4_rob_norm[,i])
  
  if(i == 1){
    pc.i[,1] = -pc.i[,1] 
  }
  
  if(i == 3){
    pc.i[,3] = -pc.i[,3] 
  }
  
  
  if(i == 4){
    pc.i[,3] = -pc.i[,3] 
  }
  
  
  min.i = min(c(pc.i))
  max.i = max(c(pc.i))
  
  plot(pc.i[,1], main = paste("Component - ",i,sep = ""),cex.lab = 1.5,
       type = "l",xaxt = "n",ylab = "",lwd = 2,cex.main = 1.6,
       xlab = "Time",ylim = c(min.i,max.i),cex.axis = 1.5)
  axis(1, at = c(seq(1,32,4),33),labels = TIME[seq(1,17,2)],cex.axis = 1.5)
  abline(h = 0,lty = 3)
  lines(pc.i[,2], type = "l",lwd = 2,col = "red")
  lines(pc.i[,3], type = "l",lwd = 2,col = "blue")
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c(expression(Phi^(2)),expression(Phi^(3)),expression(Phi^(4))), 
       xpd = TRUE, col = c("black","red","blue"), cex = 2,
       horiz = TRUE, inset = c(0,0), bty = "n", lwd = 2, lty = 1)
dev.off()
