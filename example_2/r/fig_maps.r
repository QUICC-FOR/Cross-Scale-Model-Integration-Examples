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
predictionDat = read.csv("dat/mcmc/predictionData.csv", header=FALSE)
naiveStats = readRDS("results/naiveStats.rds")
intFutStats = readRDS("results/intFutStats.rds")
intPresStats = readRDS("results/intPresStats.rds")


ocean = readOGR(dsn="dat/figure/ne_50m_ocean", layer="ne_50m_ocean")
lakes = readOGR(dsn="dat/figure/ne_50m_lakes", layer="ne_50m_lakes")
mapleRange = readOGR(dsn="dat/figure/acersacr", layer="acersacr")
# grab specific lakes
lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]


meanColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), interpolate='spline', bias=1, space="rgb")(200)
errorColors = colorRampPalette(c("#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"), interpolate='spline', space="rgb", bias=1.4)(200)
meanZlims = c(0,1)
rangeBorder = '#FF3333bb'
rangeLwd = 1.2
rangeCol = "#66666600"
xlim=c(-105, -55)
ylim=c(27,65)
titleCEX = 0.7

plotbg = function(txt="")
{
	plot(ocean, col="white", add=T)
	plot(mapleRange[c(1,3,47,43),], border=rangeBorder, col=rangeCol, add=T, lwd=rangeLwd)
	plot(grLakes, col="white", add=T)
	mtext(txt, adj=0, cex = titleCEX)
}

pdfWidth = 7.5
pdfHeight = 5.25
pdf(file = "ex2_pres_map.pdf", height=pdfHeight, width=pdfWidth)
plotLayout = matrix(c(
	1,1,2,3,
	6,7,4,5),
  byrow=T, nrow=2)
layoutWidths = c(0.5,0.5,1,1)
layout(plotLayout, widths=layoutWidths)
mapMar = c(1,1,1,1)
par(mar=mapMar)
scaleMar = c(1.75,3.75,1.75,3.75)

statsSub = which(in.range(naiveStats$long, xlim) & in.range(naiveStats$lat, ylim))
errorZlims = 
{
	vals = c(naiveStats[statsSub, 'pres_SE'], intPresStats[statsSub, 'pres_SE'])
	c(min(vals), max(vals))
}	
meanTitles = c("A. Phenofit", "B. Naive", "C. Integrated-Present")
errorTitles = c("D. Standard error (Naive)", "E. Standard error (Integrated-Present)")

quilt.plot(rawDat$all$long, rawDat$all$lat, rawDat$all$Phenofit_CRU, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(meanTitles[1])

quilt.plot(naiveStats$long, naiveStats$lat, naiveStats$pres_mean, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(meanTitles[2])

quilt.plot(intPresStats$long, intPresStats$lat, intPresStats$pres_mean, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(meanTitles[3])


quilt.plot(naiveStats$long, naiveStats$lat, naiveStats$pres_SE, col=errorColors, 
		xlim=xlim, ylim=ylim, zlim=errorZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(errorTitles[1])

quilt.plot(intPresStats$long, intPresStats$lat, intPresStats$pres_SE, col=errorColors, 
		xlim=xlim, ylim=ylim, zlim=errorZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(errorTitles[2])

## scale bars
par(mar=scaleMar)
image(y=seq(meanZlims[1], meanZlims[2],length.out=101), z=matrix(seq(meanZlims[1], meanZlims[2],length.out=100), 
	ncol=100, nrow=1), zlim=meanZlims, col=meanColors, xaxt='n', xlab='', ylab='', useRaster=T)
# mtext(meanScaleTitle, line=0.5, cex = titleCEX)
image(y=seq(errorZlims[1],errorZlims[2],length.out=101), z=matrix(seq(errorZlims[1], errorZlims[2],length.out=100), 
	ncol=100, nrow=1), zlim=errorZlims, col=errorColors, xaxt='n', xlab='', ylab='', useRaster=T)
dev.off()









pdf(file = "ex2_fut_map.pdf", height=pdfHeight, width=pdfWidth)
layout(plotLayout, widths=layoutWidths)
par(mar=mapMar)

statsSub = which(in.range(naiveStats$long, xlim) & in.range(naiveStats$lat, ylim))
errorZlims = 
{
	vals = c(naiveStats[statsSub, 'fut_SE'], intPresStats[statsSub, 'fut_SE'])
	c(min(vals), max(vals))
}	
meanTitles = c("A. Phenofit", "B. Naive", "C. Integrated-Future")
errorTitles = c("D. Standard error (Naive)", "E. Standard error (Integrated-Future)")

quilt.plot(rawDat$all$long, rawDat$all$lat, rawDat$all$Phenofit_HadA2, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(meanTitles[1])

quilt.plot(naiveStats$long, naiveStats$lat, naiveStats$fut_mean, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(meanTitles[2])

quilt.plot(intPresStats$long, intPresStats$lat, intPresStats$fut_mean, col=meanColors, 
		xlim=xlim, ylim=ylim, zlim=meanZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(meanTitles[3])


quilt.plot(naiveStats$long, naiveStats$lat, naiveStats$fut_SE, col=errorColors, 
		xlim=xlim, ylim=ylim, zlim=errorZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(errorTitles[1])

quilt.plot(intPresStats$long, intPresStats$lat, intPresStats$fut_SE, col=errorColors, 
		xlim=xlim, ylim=ylim, zlim=errorZlims, add.legend=F, xaxt='n', yaxt='n', useRaster=T)
plotbg(errorTitles[2])

## scale bars
par(mar=scaleMar)
image(y=seq(meanZlims[1], meanZlims[2],length.out=101), z=matrix(seq(meanZlims[1], meanZlims[2],length.out=100), 
	ncol=100, nrow=1), zlim=meanZlims, col=meanColors, xaxt='n', xlab='', ylab='', useRaster=T)
# mtext(meanScaleTitle, line=0.5, cex = titleCEX)
image(y=seq(errorZlims[1],errorZlims[2],length.out=101), z=matrix(seq(errorZlims[1], errorZlims[2],length.out=100), 
	ncol=100, nrow=1), zlim=errorZlims, col=errorColors, xaxt='n', xlab='', ylab='', useRaster=T)
dev.off()

