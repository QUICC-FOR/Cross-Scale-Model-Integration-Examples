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
# produce the figure for the second example



load("results/predictions.rdata")
load("dat/maple.rdata")
source("ex2_Functions.r")
library(fields)
library(rgdal)

in.range = function(x, lims) x >= lims[1] & x <= lims[2]

if(!exists('ocean')) ocean = readOGR(dsn="dat/ne_110m_ocean", layer="ne_110m_ocean")

pdfFileName = "ex2.pdf"
pdfWidth = 7.5
pdfHeight = 6
titleCEX = 0.7

meanCols = c(3,5,7,11)
SECols = meanCols + 1

## pick the data to display
latLimits = c(27,65)
longLimits = c(-105, -55)
predictionsSub = subset(predictions, in.range(lat, latLimits) & in.range(long, longLimits))
mapleAllSub = subset(mapleAll, in.range(lat, latLimits) & in.range(long, longLimits))

latitude = predictionsSub[,2]
longitude = predictionsSub[,1]
meanPredictions = cbind(mapleAllSub$Phenofit_CRU, mapleAllSub$Phenofit_HadA2, predictionsSub[,meanCols])
errorPredictions = predictionsSub[,SECols]

## plotting settings
meanZlims = c(0,1)
# meanColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d"), interpolate='spline', bias=1, space="rgb")(200)
meanColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), interpolate='spline', bias=1, space="rgb")(200)
meanTitles = c("A. Phenofit: Present", "B. Phenofit: HadA2", "C. Naive SDM: Present",  "D. Naive SDM: HadA2", "G. Integrated: Present", "H. Integrated: HadA2")
errorTitles = c("E. SDM: Present SE",  "F. SDM: HadA2 SE", "I. Integrated: Present SE", "J. Integrated: HadA2 SE")
meanScaleTitle = "Suitability"

errorZlims = range(errorPredictions)
errorColors = colorRampPalette(c("#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"), interpolate='spline', space="rgb", bias=1.3)(200)
errorTitles = paste(errorTitles, "SE", sep=" ")
errorScaleTitle = "Posterior Standard Error"


# produce the plot
pdf(file=pdfFileName, width=pdfWidth,height=pdfHeight)
plotLayout = matrix(c(
	1,2,11,11,
	1,2,12,12,
	3,4,7,8,
	5,6,9,10),
  byrow=T, nrow=4)
  
layout(plotLayout, heights=c(0.5, 0.5, 1, 1))
par(mar=c(1,1,1,1))

for(i in 1:ncol(meanPredictions)) {
	quilt.plot(longitude, latitude, meanPredictions[,i], col=meanColors, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
	plot(ocean, col="white", add=T)
	mtext(meanTitles[i], adj=0, cex = titleCEX)
}

for(i in 1:ncol(errorPredictions)) {
	quilt.plot(longitude, latitude, errorPredictions[,i], col=errorColors, zlim=errorZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
	plot(ocean, col="white", add=T)
	mtext(errorTitles[i], adj=0, cex = titleCEX)
}

## scale bars
par(mar=c(3,1,3,1))
image(x=seq(meanZlims[1],meanZlims[2],length.out=101), z=matrix(seq(meanZlims[1], meanZlims[2],length.out=100), 
	nrow=100, ncol=1), zlim=meanZlims, col=meanColors, yaxt='n', xlab='', ylab='', useRaster=T)
mtext(meanScaleTitle, line=0.5, cex = titleCEX)
image(x=seq(errorZlims[1],errorZlims[2],length.out=101), z=matrix(seq(errorZlims[1], errorZlims[2],length.out=100), 
	nrow=100, ncol=1), zlim=errorZlims, col=errorColors, yaxt='n', xlab='', ylab='', useRaster=T)
mtext(errorScaleTitle, line=0.5, cex = titleCEX)
