# ex2_makeFigures.r
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
# produce the figures for the second example

load("results/predictions.rdata")
load("dat/maple.rdata")
source("ex2_Functions.r")

pdfFileName <- "ex2.pdf"
pdfWidth <- 9
pdfHeight <- 5.5
titleCEX <- 0.7

library(fields)

# plot the predictions
pdf(file=pdfFileName, width=pdfWidth,height=pdfHeight)
plotLayout <- matrix(c(
	1,2,11,11,
	1,2,12,12,
	3,4,7,8,
	5,6,9,10),
  byrow=T, nrow=4)
  
layout(plotLayout, heights=c(0.5, 0.5, 1, 1))
par(mar=c(1,1,1,1))


meanCols <- c(3,7,11,15)
SECols <- meanCols + 1
ciCols <- SECols + 1
useSE <- TRUE	# should we plot the SE or the credible interval


## pick the data to display
latitude <- predictions[,2]
longitude <- predictions[,1]
meanPredictions <- cbind(maple$Phenofit_CRU, maple$Phenofit_HadA2, predictions[,meanCols])
if(useSE) {
	errorPredictions <- predictions[,SECols]
} else {
	errorPredictions <- predictions[,ciCols + 1] - predictions[,ciCols]
}

## plotting settings
meanZlims <- c(0,1)
meanColors <- bpy.colors()
meanTitles <- c("A. Phenofit: Present", "B. Phenofit: HadA2", "C. Naive SDM: Present",  "D. Naive SDM: HadA2", "G. Integrated: Present", "H. Integrated: HadA2")
errorTitles <- c("E. SDM: Present",  "F. SDM: HadA2", "I. Integrated: Present", "J. Integrated: HadA2")
meanScaleTitle <- "Suitability"

if(useSE) {
	errorZlims <- range(errorPredictions)
	errorColors <- colorRampPalette(c('#000000', '#FF0000', '#FFFF00'), interpolate='spline', bias=3)(200)
	errorTitles <- paste(errorTitles, "SE", sep=" ")
	errorScaleTitle <- "Posterior Standard Error"
} else {
	errorZlims <- c(0,1)
	errorColors <- colorRampPalette(c('#000000', '#FF0000', '#FFFF00'), interpolate='spline', bias=2)(100)
	errorTitles <- paste(errorTitles, "CI Width", sep=" ")
	errorScaleTitle <- "Posterior Credible Interval Width"
}


# produce the plot

for(i in 1:ncol(meanPredictions)) {
	quilt.plot(longitude, latitude, meanPredictions[,i], col=meanColors, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
	mtext(meanTitles[i], adj=0, cex = titleCEX)
}

for(i in 1:ncol(errorPredictions)) {
	quilt.plot(longitude, latitude, errorPredictions[,i], col=errorColors, zlim=errorZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
	mtext(errorTitles[i], adj=0, cex = titleCEX)
}

par(mar=c(3,1,3,1))
image(x=seq(meanZlims[1],meanZlims[2],length.out=101), z=matrix(seq(meanZlims[1], meanZlims[2],length.out=100), nrow=100, ncol=1), zlim=meanZlims, col=meanColors, yaxt='n', xlab='', ylab='', useRaster=T)
mtext(meanScaleTitle, line=0.5, cex = titleCEX)
image(x=seq(errorZlims[1],errorZlims[2],length.out=101), z=matrix(seq(errorZlims[1], errorZlims[2],length.out=100), nrow=100, ncol=1), zlim=errorZlims, col=errorColors, yaxt='n', xlab='', ylab='', useRaster=T)
mtext(errorScaleTitle, line=0.5, cex = titleCEX)

dev.off()