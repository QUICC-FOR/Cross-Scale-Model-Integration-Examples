# ex1_makeMapFig.r
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
# produces a figure showing a map of model integration results

## Rscript usage
## rscript ex1_results.r filename.pdf
## will output plots to pdf

## setup
library(rjags)

height.map <- 6		# this plot gets a different height since it is 2x2
xlab.line <- 2
ylab.line <- 1

# define color palettes for the two types of image plots
presence.palette=colorRampPalette(c("#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d"),space="rgb")
uncertainty.palette=colorRampPalette(c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"), space="rgb")

main <- function() {

	filename <- commandArgs(TRUE)[1]
	source("ex1_globals.r")
	pdf(file=filename, width=width, height=height.map)

	layout(matrix(c(1,2,3,4,5,6,7,8), nrow=2, byrow=TRUE), widths=c(4,1,4,1), heights=c(4,4))
	par(xpd=FALSE, bty=bty, mgp=mgp, tcl=tcl)

	# get data from disk
	load("dat/ex1_m1.rdata")
	load("dat/ex1_m2.rdata")
	load("dat/ex1_mm.rdata")
	
	# set up matrices for plotting
	m1Matrix <- posterior_map(m1Predictions$map$mean, m1Predictions$mapDomain)
	m1Uncertainty <- posterior_map(m1Predictions$map$SD, m1Predictions$mapDomain)
	mmMatrix <- posterior_map(mmPredictions$map$mean, mmPredictions$mapDomain)
	mmUncertainty <- posterior_map(mmPredictions$map$SD, mmPredictions$mapDomain)


	zlSDM <- range(c(m1Matrix, mmMatrix))
	zlUn <- range(c(m1Uncertainty, mmUncertainty))
	
	with(m1Predictions$mapDomain, 
		sdm_plot(m1Matrix, unique(precip), unique(temp), zlim=zlSDM, title="(a) Naive model, probability of presence"))
	with(m1Predictions$mapDomain, 
		sdm_plot(m1Uncertainty, unique(precip), unique(temp), zlim=zlUn, col.palette=uncertainty.palette, title="(b) Naive model, standard deviation"))
	with(mmPredictions$mapDomain, 
		sdm_plot(mmMatrix, unique(precip), unique(temp), zlim=zlSDM, title="(c) Integrated model, probability of presence"))
	with(mmPredictions$mapDomain, 
		sdm_plot(mmUncertainty, unique(precip), unique(temp), zlim=zlUn, col.palette=uncertainty.palette, title="(d) Integrated model, standard deviation"))

	dev.off()
}


sdm_plot <- function(mapdata, precip, temp, title="", zlim=NA,  col.palette=presence.palette, n.color=25) {
	if(!is.vector(zlim)) zlim <- range(mapdata)
	par(mar=c(4,2.5,3,0.5), xpd=TRUE)
	image(x=precip, y=temp, z=mapdata, xaxt='n', yaxt='n', col=col.palette(25), xlab="", ylab="", zlim=zlim)
	mtext("Temperature", side=2, line=ylab.line, cex = cex.lab)
	mtext(title, side=3, line=title.line, cex=cex.title, adj=0)
	sdm_xaxis(precipRegimes)
	
	par(mar=c(5,1,3,3))
	image.scale(mapdata, col=col.palette(n.color), axis.pos=4, zlim=zlim)
}

posterior_map <- function(predictions, domain) {
	posteriorMatrix <- matrix(predictions, nrow=length(unique(domain$precip)), ncol=length(unique(domain$temp)), byrow=T)
	return(posteriorMatrix)
}

sdm_xaxis <- function(precip.range) {
	hist.yloc <- -1.21
	hist.line.yloc <- -1.14
	hist.xloc <- precip.range$current[2] + (precip.range$current[2] - precip.range$current[1])/11 #11 is the number of cells in the raster
	fut.yloc <- -1.25
	fut.line.yloc <- -1.18
	fut.xloc <- precip.range$future[1] - (precip.range$current[2] - precip.range$current[1])/11 #11 is the number of cells in the raster
	hist.pos <- 2
	fut.pos <- 4
	lines(c(precip.range$current[1],precip.range$current[2]), rep(hist.line.yloc,2), lwd=2)
	lines(c(precip.range$future[1], precip.range$future[2]), rep(fut.line.yloc,2), lwd=2, col='blue')
	text(hist.xloc, hist.yloc, "Historical Range", pos=hist.pos, cex=0.8)
	text(fut.xloc, fut.yloc, "Projected Range", pos=fut.pos, col='blue', cex=0.8)
	mtext("Precipitation", side=1, line=xlab.line, cex = cex.lab)
}






## SHAMELESSLY BORROWED FROM:
# http://menugget.blogspot.ca/2011/08/adding-scale-to-image-plot.html#more
# all credit and blame to the authors

# This function creates a color scale for use with the image()
# function. Input parameters should be consistent with those
# used in the corresponding image plot. The "axis.pos" argument
# defines the side of the axis. The "add.axis" argument defines
# whether the axis is added (default: TRUE) or not (FALSE).
image.scale <- function(z, zlim, col = heat.colors(12), breaks, axis.pos=1, add.axis=TRUE, ...){
	if(!missing(breaks)){
		if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
	}
	if(missing(breaks) & !missing(zlim)){
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
	}
	if(missing(breaks) & missing(zlim)){
		zlim <- range(z, na.rm=TRUE)
		zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
		zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
	}
	poly <- vector(mode="list", length(col))
	for(i in seq(poly)){
		poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
	}
	if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
	if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
	plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
	for(i in seq(poly)){
		if(axis.pos %in% c(1,3)){
			polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
		}
		if(axis.pos %in% c(2,4)){
			polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
		}
	}
	box()
	if(add.axis) {axis(axis.pos)}
}



main()