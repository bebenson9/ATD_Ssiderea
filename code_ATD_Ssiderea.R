library(arm)
library(MuMIn)
library(lme4)
library(tidyverse)
library(ggplot2)
library(plotrix)

setwd("/Users/brookebenson/ATD/")
data <- as.data.frame(read.csv('Data_ATD_Ssiderea'))

# MODEL AVERAGING
## ATDot == temp 5 year values
## ATDop = par 5 year values
model <- lmer(ATDot~ tissue*meanLE5yr*E1+corelength+(1|reefzone/site), data, na.action = na.fail)
std.model <- standardize(model, standardize.y=F)
modelset <- dredge(std.model)
top <-get.models(modelset, cumsum(weight)<=0.95)
average <- model.avg(top)
summary(average)

# GET MODEL PREDICTIONS
data$MATDpar <-predict(average,se.fit = F, type = 'response', backtransform = F)
data$MATDsst <-predict(average,se.fit = F, type = 'response', backtransform = F) 
sstmod <- lm(MATDsst~oatd, data)
parmod <- lm(MATDpar~oatd, data)

# PLOT MODEL PREDICTIONS
plot.new()
plot(MATDsst~oatd, data, xlim = c(0,13), ylim = c(2,11), pch = 19, xlab = 'Observed SST-derived ATD (months)', ylab = 'Modeled ATD (months)', col = 'red')
points(MATDpar~oatd, data, pch = 17, col ='darkgoldenrod2')
abline(sstmod, col = 'red')
abline(parmod, lty = 'dashed', col = 'darkgoldenrod2')

##EFFECT SIZE PLOT
#coefficients from model output
full.coef <- data.frame(Coef = c(), SE = c())
confInt <- cbind(1:8,full.coef, apply(full.coef,1,function(x) x[1]+qnorm(0.025)*x[2]), apply(full.coef,1,function(x) x[1]+qnorm(0.975)*x[2]), apply(full.coef,1,function(x) x[1]+qnorm(0.1)*x[2]), apply(full.coef,1,function(x) x[1]+qnorm(0.9)*x[2]))
colnames(confInt) <- c('N','Coef','SE','-95%','+95%','-80%','+80%')
E1 <- expression('E'[1])
E1ext <- (expression('E'[1]* phantom(x)* 'x Linear Extension'))
E1tis <- (expression('E'[1]* phantom(x)* 'x Tissue Thickness'))
E1lintis <- (expression('E'[1]* 'Linear Extension x \n x Tissue Thickness'))
param.names<-c(E1,'Linear Extension', 'Tissue Thickness',E1ext, E1tis,'Linear Extension \nx Tissue Thickness', E1lintis, 'Core Length')
#expand left margin
par(mar=c(5.0,10,4.1,2.1))
#change bar ends to square
par(lend=2)  
#set up plot
plot.new()
plot(range(confInt[,4:7]), c(1,2), type='n', axes=F, ylab='', xlim=c(-15.0,15.0), ylim=c(.75,8.25), xlab='Effect Size', cex.lab=1.3)
axis(1,cex.axis=1.3)
axis(2,at=1:8, labels=param.names, cex.axis=1, las=2)
abline(v=0, col='grey60', lty=2)
apply(confInt, 1, function(x) segments(x[4], x[1], x[5], x[1], lwd=7, col='grey70'))
apply(confInt, 1, function(x) segments(x[6], x[1], x[7], x[1], lwd=9, col='grey40'))
apply(confInt, 1, function(x) points(x[2], x[1], pch=16, col='white', cex=2))
apply(confInt, 1, function(x) points(x[2], x[1], pch=1, cex=2))
box()
effectSize <- recordPlot()
