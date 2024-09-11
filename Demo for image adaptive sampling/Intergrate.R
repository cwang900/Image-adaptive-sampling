library(MASS)
library(RSpectra)

source('model1_simulation.R')
source('model2_simulation.R')


#generate fig1.a from table A1 
MFPCA001<-output1[seq(1,35,by = 2),1]
FPCA001<-output1[seq(1,35,by = 2),4]
FPCARS<-output1[seq(1,35,by = 2),7]
MFULL1<-full1[seq(3,13,by = 2)]


text = c("MFPCA fully observed", expression("MFPCA "*Delta*"= 0.01"), expression("FPCA "*Delta*"= 0.01"),"FPCA RS")
# r = 10
par(mfcol = c(3,1))
par(mar=c(2,2,2,1))
par(oma = c(2.5,2.5,1,1))
tick<-c(0.4,0.8,1.2,1.6,2.0,3.0)
plot(x = tick,y = MFPCA001[1:6],pch = 5,xlab = "",ylab = "",font.main = 1, main = paste("Model ", as.roman(1), "with r = 10"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.5,cex.lab=1.5, cex.axis=1.5,col = rainbow(12)[6])
axis(side = 1, at = tick,cex.axis = 1.5)
lines(x = tick,y = MFPCA001[1:6],lty = 1,col = rainbow(12)[6])
points(x = tick,y = FPCA001[1:6],cex = 1.2,pch = 15,col = rainbow(12)[8])
lines(x = tick,y = FPCA001[1:6],lty = 3,col = rainbow(12)[8])
points(x = tick,y = FPCARS[1:6],cex = 1.2, pch = 19)
lines(x = tick,y = FPCARS[1:6],lty = 5)
lines(x = tick, y = MFULL1,lty = 1,col = 'orange',pch = 5)
points(x = tick, y = MFULL1,pch = 5, cex = 1.2,col = 'orange')
legend(2.0, 195, legend= text,text.width = strwidth(text)[1]*1.4,
       col=c('orange',rainbow(12)[6],rainbow(12)[8],'black'), lty = c(1,1,3,5), pch = c(5,5,15,19), cex=1.2)


# r = 6
plot(x = tick,y = MFPCA001[7:12],pch = 5,xlab = "",ylab = "",font.main = 1,main = paste("Model ", as.roman(1), "with r = 7"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.5,cex.lab=1.5, cex.axis=1.5,col = rainbow(12)[6])
axis(side = 1, at = tick,cex.axis = 1.5)
lines(x = tick,y = MFPCA001[7:12],lty = 1,col = rainbow(12)[6])
points(x = tick,y = FPCA001[7:12],cex = 1.2,pch = 15,col = rainbow(12)[8])
lines(x = tick,y = FPCA001[7:12],lty = 3,col = rainbow(12)[8])
points(x = tick,y = FPCARS[7:12],cex = 1.2, pch = 19)
lines(x = tick,y = FPCARS[7:12],lty = 5)
lines(x = tick, y = MFULL1,lty = 1,col = 'orange',pch = 5)
points(x = tick, y = MFULL1,pch = 5, cex = 1.2,col = 'orange')


# r = 4
plot(x = tick,y = MFPCA001[13:18],pch = 5,xlab = "",ylab = "",font.main = 1, main = paste("Model ", as.roman(1), "with r = 4"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.5,cex.lab=1.5, cex.axis=1.5,col = rainbow(12)[6])
axis(side = 1, at = tick,cex.axis = 1.5)
lines(x = tick,y = MFPCA001[13:18],lty = 1,col = rainbow(12)[6])
points(x = tick,y = FPCA001[13:18],cex = 1.2,pch = 15,col = rainbow(12)[8])
lines(x = tick,y = FPCA001[13:18],lty = 3,col = rainbow(12)[8])
points(x = tick,y = FPCARS[13:18],cex = 1.2, pch = 19)
lines(x = tick,y = FPCARS[13:18],lty = 5)
lines(x = tick, y = MFULL1,lty = 1,col = 'orange',pch = 5)
points(x = tick, y = MFULL1,pch = 5, cex = 1.2,col = 'orange')
mtext(expression("Shift magnitude "*delta), side = 1, line = 0.8, outer = TRUE, cex = 1.1 )
mtext(expression('ARL'[1]),side = 2,line = 0.8,outer = TRUE,cex = 1.1,las = 0)

output2 #display table A2 

#generate plot from table A2

MFCA001<-output2[seq(1,35,by = 2),1]
FCA001<-output2[seq(1,35,by = 2),4]
FCARS<-output2[seq(1,35,by = 2),7]
MFULL2<-full2[seq(3,13,by = 2)]


text = c("MFPCA fully observed",expression("MFPCA "*Delta*"= 0.01"), expression("FPCA "*Delta*"= 0.01"),"FPCA RS")
par(mfcol = c(3,1))
par(mar=c(2,2,2,1))
par(oma = c(2.5,2.5,1,1))
tick<-c(0.8,1.2,1.6,2.0,3.0,4.0)

#r = 8
plot(x = tick,y = MFCA001[1:6],pch = 5,xlab = "",ylab = "", font.main = 1,main = paste("Model ", as.roman(2), "with r = 10"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.5,cex.lab=1.5, cex.axis=1.5,col = rainbow(12)[6])
axis(side = 1, at = tick,cex.axis = 1.5)
lines(x = tick,y = MFCA001[1:6],lty = 1,col = rainbow(12)[6])
points(x = tick,y = FCA001[1:6],cex = 1.2,pch = 15,col = rainbow(12)[8])
lines(x = tick,y = FCA001[1:6],lty = 3,col = rainbow(12)[8])
points(x = tick,y = FCARS[1:6],cex = 1.2, pch = 19)
lines(x = tick,y = FCARS[1:6],lty = 5)
lines(x = tick, y = MFULL2,lty = 1,col = 'orange',pch = 5)
points(x = tick, y = MFULL2,pch = 5, cex = 1.2,col = 'orange')
legend(2.8, 195, legend= text,text.width = strwidth(text)[1]*1.4,
       col=c('orange',rainbow(12)[6],rainbow(12)[8],'black'), lty = c(1,1,3,5), pch = c(5,5,15,19), cex=1.2)


# r = 6
plot(x = tick,y = MFCA001[7:12],pch = 5,xlab = "",ylab = "", font.main = 1, main = paste("Model ", as.roman(2), "with r = 7"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.5,cex.lab=1.5, cex.axis=1.5,col = rainbow(12)[6])
axis(side = 1, at = tick,cex.axis = 1.5)
lines(x = tick,y = MFCA001[7:12],lty = 1,col = rainbow(12)[6])
points(x = tick,y = FCA001[7:12],cex = 1.2,pch = 15,col = rainbow(12)[8])
lines(x = tick,y = FCA001[7:12],lty = 3,col = rainbow(12)[8])
points(x = tick,y = FCARS[7:12],cex = 1.2, pch = 19)
lines(x = tick,y = FCARS[7:12],lty = 5)
lines(x = tick, y = MFULL2,lty = 1,col = 'orange',pch = 5)
points(x = tick, y = MFULL2,pch = 5, cex = 1.2,col = 'orange')


# r = 4
plot(x = tick,y = MFCA001[13:18],pch = 5,xlab = "",ylab = "", font.main = 1, main = paste("Model ", as.roman(2), "with r = 4"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.5,cex.lab=1.5, cex.axis=1.5,col = rainbow(12)[6])
axis(side = 1, at = tick,cex.axis = 1.5)
lines(x = tick,y = MFCA001[13:18],lty = 1,col = rainbow(12)[6])
points(x = tick,y = FCA001[13:18],cex = 1.2,pch = 15,col = rainbow(12)[8])
lines(x = tick,y = FCA001[13:18],lty = 3,col = rainbow(12)[8])
points(x = tick,y = FCARS[13:18],cex = 1.2, pch = 19)
lines(x = tick,y = FCARS[13:18],lty = 5)
lines(x = tick, y = MFULL2,lty = 1,col = 'orange')
points(x = tick, y = MFULL2,pch = 5, cex = 1.2,col = 'orange')
mtext(expression("Shift magnitude "*delta), side = 1, line = 0.8, outer = TRUE, cex = 1.1)
mtext(expression('ARL'[1]),side = 2,line = 0.8,outer = TRUE,cex = 1.1,las = 0)



