#!/usr/bin/Rscript
# 
#   Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
# 
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or (at
#   your option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
#
#
#
#
#


library(coda)
library(gam)

# if TRUE, then calibration curve will be fit using a GAM instead of a linear regression
use.gam = FALSE

# if TRUE, then the calibration curve will be weighted based on the number of cells from
# the validation dataset that fall into each bin
use.weights = TRUE

validation = readRDS("dat/rawData.rds")$valid
load("results/posteriors.rdata")

makePred = function(bn, data) {
	with(data, plogis(bn[1] + bn[2]*ddeg + bn[3]*ddeg^2 + bn[4]*ddeg^3 + 
			bn[5]*an_prcp + bn[6]*an_prcp^2 + bn[7]*pToPET + bn[8]*pToPET^2))
}


presDat = with(validation, PresObs + ifelse(is.na(pres), 0, pres))

# this function controls the 'binning' into nclass categories ranging from 0 to 1
# returns a vector of the same length as x, with values as the midpoint of each bin
classify = function(x, nclass) ifelse(x==1, nclass-1, floor(x*nclass))/nclass

ncats=20
xvals = seq(0.05, 0.95, length.out=ncats)

# x values for drawing the curve
predict.x = data.frame(xvals=seq(0,1, length.out=200))

# data structures to hold the y-values for each posterior sample
curve.y.naive = matrix(NA, nrow=nrow(naivePosterior), ncol=nrow(predict.x))
curve.y.int = matrix(NA, nrow=nrow(intPresPosterior), ncol=nrow(predict.x))

for(i in 1:nrow(naivePosterior))
{
	preds.naive = data.frame(pres=presDat, predict=makePred(naivePosterior[i,], validation))
	preds.int = data.frame(pres=presDat, predict=makePred(intPresPosterior[i,], validation))
	# classify predictions
	preds.naive$class = classify(preds.naive$predict, ncats)
	preds.int$class = classify(preds.int$predict, ncats)
	nclass.naive = length(unique(preds.naive$class))
	nclass.int = length(unique(preds.int$class))
	freqs.naive = freqs.int = nn = ni = rep(NA, ncats)
	freqs.naive[1:nclass.naive] = as.numeric(with(preds.naive, table(pres, class)[2,] / table(class)))
	freqs.int[1:nclass.int] = as.numeric(with(preds.int, table(pres, class)[2,] / table(class)))
	
	if(use.weights)
	{
		nn[1:nclass.naive] = as.numeric(with(preds.naive, table(class)))
		ni[1:nclass.int] = as.numeric(with(preds.int, table(class)))
		if(use.gam)
		{
			naive.mod = gam(freqs.naive~s(xvals, 3), weights=nn)
			int.mod = gam(freqs.int~s(xvals, 3), weights=ni)
		} else {
			naive.mod = glm(freqs.naive~xvals, weights=nn)
			int.mod = glm(freqs.int~xvals, weights=ni)
		}
	} else {
		if(use.gam) {
			naive.mod = gam(freqs.naive~s(xvals, 3))
			int.mod = gam(freqs.int~s(xvals, 3))
		} else {
			naive.mod = glm(freqs.naive~xvals)
			int.mod = glm(freqs.int~xvals)
		}
	}
	curve.y.naive[i,] = predict(naive.mod, newdata=predict.x)
	curve.y.int[i,] = predict(int.mod, newdata=predict.x)
}

col.curve='#00ff33'

plot.int = function(x, y, col, tr='33')
{
	upper = apply(y, 2, quantile, 0.975, na.rm=T)
	lower = apply(y, 2, quantile, 0.025, na.rm=T)
	polygon(c(x, rev(x)), c(upper, rev(lower)), col=paste(col, tr, sep=""), border=NA)	
}

paperwidth = 8
dpi = 600
hToWRatio = 0.45
width = as.integer(dpi*paperwidth)
height = as.integer(width * hToWRatio)
fontsize = 12
png(w=width, h=height, file="ex2_calib.png", pointsize=fontsize, res = dpi)

cal.lims = c(0,0.7)
par(mfrow=c(1,2), mar=c(3,3,0.3, 3), mgp=c(2,0.5,0), tck=-0.01, bty='n', cex=0.8)
plot(predict.x$xvals, colMeans(curve.y.naive, na.rm=T), col=col.curve, type='l', 
		 xlim=cal.lims, ylim=cal.lims, xlab="Predicted Probability of Presence", ylab="Observed Proportion Present")
plot.int(predict.x$xvals, curve.y.naive, col.curve)
lines(cal.lims, cal.lims, lty=2)
text(0.2, max(cal.lims), "Naive", cex=0.8)


plot(predict.x$xvals, colMeans(curve.y.int, na.rm=T), col=col.curve, type='l', 
		 xlim=cal.lims, ylim=cal.lims, xlab="Predicted Probability of Presence", ylab="Observed Proportion Present")
plot.int(predict.x$xvals, curve.y.int, col.curve)
lines(cal.lims, cal.lims, lty=2)
text(0.2, max(cal.lims), "Integrated-Present", cex=0.8)

dev.off()
