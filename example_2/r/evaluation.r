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
library(pROC)
validation = readRDS("dat/rawData.rds")$valid
calib = readRDS("dat/rawData.rds")$calib
load("results/posteriors.rdata")

makeAUC = function(bn, bi, data) {
	rocDat = with(data,
	{
		data.frame(
			naivePred = plogis(bn[1] + bn[2]*ddeg + bn[3]*ddeg^2 + bn[4]*an_prcp + bn[5]*an_prcp^2 + bn[6]*pToPET + bn[7]*pToPET^2),
			intPred = plogis(bi[1] + bi[2]*ddeg + bi[3]*ddeg^2 + bi[4]*an_prcp + bi[5]*an_prcp^2 + bi[6]*pToPET + bi[7]*pToPET^2),
			PresObs = PresObs
		)
	})

	c(naive1=roc(response = rocDat$PresObs, predictor = rocDat$naivePred)$auc,
		int1=roc(response = rocDat$PresObs, predictor = rocDat$intPred)$auc)
}


auc = matrix(NA, nrow=nrow(naivePosterior), ncol=2)
for(i in 1:nrow(naivePosterior))
{
	auc[i,] = makeAUC(naivePosterior[i,], intPresPosterior[i,], validation)
}

adj = 2
lwd=1.5
col.n = '#377eb8'
col.i = '#4daf4a'
nd = density(auc[,1], adjust=adj)
id = density(auc[,2], adjust=adj)
nd$y = nd$y / max(nd$y)
id$y = id$y / max(id$y)

paperwidth = 5.5
dpi = 600
hToWRatio = 0.85
width = as.integer(dpi*paperwidth)
height = as.integer(width * hToWRatio)
fontsize = 12
png(w=width, h=height, file="ex2_auc.png", pointsize=fontsize, res = dpi)
layout(matrix(c(2,1), nrow=2), heights=c(.2,1))
par(tcl=-0.2, mgp=c(2,0.5,0), cex=0.8, mar=c(3, 3, 0, 0.5), bty='n')
plot(nd, main="", sub="", xlab="AUC", ylab="Relative Density", col=col.n, lwd=lwd)
lines(id, col=col.i, lwd=lwd)
legend(0.82, 1, legend=c('Naive', 'Integrated'), col=c(col.n, col.i), lwd=lwd, bty='n')
par(mar=c(0, 3, 0.5, 0.5))
boxplot(auc, pch='|', cex=0.4, xaxt='n', yaxt='n', col=c(col.n, col.i), horizontal=T, lty=1)
dev.off()

saveRDS(auc, 'results/auc.rds')

