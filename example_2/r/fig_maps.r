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


library(fields)
library(rgdal)
library(coda)

in.range = function(x, lims) x >= lims[1] & x <= lims[2]

rawDat = readRDS("dat/rawData.rds")
load("results/posteriors.rdata")
useR = TRUE

ocean = readOGR(dsn="dat/figure/ne_50m_ocean", layer="ne_50m_ocean")
lakes = readOGR(dsn="dat/figure/ne_50m_lakes", layer="ne_50m_lakes")
mapleRange = readOGR(dsn="dat/figure/acersacr", layer="acersacr")
# grab specific lakes
lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]


meanColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), interpolate='spline', bias=1, space="rgb")(200)
errorColors = colorRampPalette(c("#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"), interpolate='spline', space="rgb", bias=1.4)(200)
errDiffCols = colorRampPalette(c("blue", "white", "orange"), interpolate='spline', space="rgb", bias=1)(200)
meanZlims = c(0,1)
rangeBorder = '#FF3333bb'
rangeLwd = 1.2
rangeCol = "#66666600"
xlim=c(-105, -55)
ylim=c(30,65)
titleCEX = 0.7

plotbg = function(txt="")
{
	plot(ocean, col="white", add=T)
	plot(mapleRange[c(1,3,47,43),], border=rangeBorder, col=rangeCol, add=T, lwd=rangeLwd)
	plot(grLakes, col="white", add=T)
	mtext(txt, adj=0, cex = titleCEX)
}

pdfWidth = 6.5
pdfHeight = 2.75
pdf(file = "ex2_pres_map.pdf", height=pdfHeight, width=pdfWidth)
plotLayout = matrix(c(
	1,2,3,4,
	5,5,6,6),
  byrow=T, nrow=2)
layoutHeights = c(1,0.3)
layout(plotLayout, heights=layoutHeights)
mapMar = c(1,1,1,1)
par(mar=mapMar)
scaleMar = c(2.25,0.75,1.75,1.75)

statsSub = which(in.range(naivePredictions$long, xlim) & in.range(naivePredictions$lat, ylim))
errorZlims = 
{
	vals = c(naivePredictions[statsSub, 'presSE'], intPresPredictions[statsSub, 'presSE'])
	c(min(vals), max(vals))
}	
meanTitles = c("A. Phenofit", "B. Naive", "C. Integrated-Present")
errorTitles = c("D. Uncertainty")
meanScaleTitle = 'Probability of presence'
errScaleTitle = expression('Integrated SE - Naive SE')
quilt.plot(rawDat$all$long, rawDat$all$lat, rawDat$all$Phenofit_CRU, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(meanTitles[1])

quilt.plot(naivePredictions$long, naivePredictions$lat, naivePredictions$pres, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(meanTitles[2])

quilt.plot(intPresPredictions$long, intPresPredictions$lat, intPresPredictions$pres, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(meanTitles[3])

errDiffZlims = c(-0.0075, 0.0075)
quilt.plot(naivePredictions$long, naivePredictions$lat, intPresPredictions$presSE - 
		naivePredictions$presSE, col=errDiffCols, xlim=xlim, ylim=ylim, zlim=errDiffZlims, 
		add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(errorTitles[1])

## scale bars
par(mar=scaleMar)
image(x=seq(meanZlims[1], meanZlims[2],length.out=101), z=matrix(seq(meanZlims[1], meanZlims[2],length.out=100), 
	nrow=100, ncol=1), zlim=meanZlims, col=meanColors, yaxt='n', xlab='', ylab='', useRaster=useR)
mtext(meanScaleTitle, line=0.5, cex = titleCEX)
image(x=seq(errDiffZlims[1],errDiffZlims[2],length.out=101), 
	z=matrix(seq(errDiffZlims[1], errDiffZlims[2],length.out=100), 
	nrow=100, ncol=1), zlim=errDiffZlims, col=errDiffCols, yaxt='n', xaxt='n', xlab='', ylab='', useRaster=useR)
axis(side=1, at=seq(-0.006,0.006,0.003))
mtext(errScaleTitle, line=0.5, cex = titleCEX)
dev.off()

pdf(file = "ex2_fut_map.pdf", height=pdfHeight, width=pdfWidth)
layout(plotLayout, heights=layoutHeights)
par(mar=mapMar)

statsSub = which(in.range(naivePredictions$long, xlim) & in.range(naivePredictions$lat, ylim))
errorZlims = 
{
	vals = c(naivePredictions[statsSub, 'futSE'], intFutPredictions[statsSub, 'futSE'])
	c(min(vals), max(vals))
}	
meanTitles = c("A. Phenofit", "B. Naive", "C. Integrated-Future")
errorTitles = c("D. Uncertainty")

quilt.plot(rawDat$all$long, rawDat$all$lat, rawDat$all$Phenofit_HadA2, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(meanTitles[1])

quilt.plot(naivePredictions$long, naivePredictions$lat, naivePredictions$fut, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(meanTitles[2])

quilt.plot(intFutPredictions$long, intFutPredictions$lat, intFutPredictions$fut, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(meanTitles[3])


errDiffZlims = c(-0.015, 0.015)
quilt.plot(naivePredictions$long, naivePredictions$lat, intFutPredictions$futSE - 
		naivePredictions$futSE, col=errDiffCols, xlim=xlim, ylim=ylim, zlim=errDiffZlims, 
		add.legend=F, xaxt='n', yaxt='n', useRaster=useR)
plotbg(errorTitles[1])


## scale bars
par(mar=scaleMar)
image(x=seq(meanZlims[1], meanZlims[2],length.out=101), z=matrix(seq(meanZlims[1], meanZlims[2],length.out=100), 
	nrow=100, ncol=1), zlim=meanZlims, col=meanColors, yaxt='n', xlab='', ylab='', useRaster=useR)
mtext(meanScaleTitle, line=0.5, cex = titleCEX)
image(x=seq(errDiffZlims[1],errDiffZlims[2],length.out=101), 
	z=matrix(seq(errDiffZlims[1], errDiffZlims[2],length.out=100), 
	nrow=100, ncol=1), zlim=errDiffZlims, col=errDiffCols, yaxt='n',xaxt='n', xlab='', ylab='', useRaster=useR)
axis(side=1, at=seq(-0.01,0.01,0.005))
mtext(errScaleTitle, line=0.5, cex = titleCEX)
dev.off()

