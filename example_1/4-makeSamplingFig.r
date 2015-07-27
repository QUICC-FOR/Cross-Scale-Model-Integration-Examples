# ex1_makeSamplingFig.r
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
# produces a figure showing how generated points were sampled for the first example



## Rscript usage
## rscript 4-makeSamplingFig.r filename.pdf
## will output plots to pdf


model2_predict_start <- 4 # the column index in the posterior of model 2 that marks the start of the predictions

sdmPoints.col <- function(x) {
	# x: vector of presence/absence coded as 1 or 0
	return(rep("black", length(x)))
}
sdmPoints.sym <- function(x) {
	return(ifelse(x, 16, 4))
}

height <- 3

filename <- commandArgs(TRUE)[1]
source("ex1_globals.r")
pdf(file=filename, width=width, height=height)

layout(matrix(c(1,2), nrow=1, byrow=TRUE), widths=c(1,1), heights=c(1))
par(mar=margins, xpd=FALSE, bty=bty, mgp=mgp, tcl=tcl)

# get data from disk
load("dat/ex1_m1.rdata")

## first plot: simulation data
plot(m1SimData[,1], m1SimData[,2], col=sdmPoints.col(m1SimData[,3]), 
	pch=sdmPoints.sym(m1SimData[,3]), xlab="", ylab="", xaxt='n', yaxt='n', ylim=c(-1,1), xlim=c(0,1))
make_axes(xlab="Precipitation", ylab="Temperature", main="(a) Simulated presence/absence samples")


# get data for second model
load("dat/ex1_m2.rdata")

## second plot
## experimental data and results
with(m2SimData, {
	yl=c(min(r)-max(SE), max(r)+max(SE))
	plot(precip, r, ylim=yl, pch=sdmPoints.sym(TRUE),  xlab="", ylab="", xaxt='n', yaxt='n')
	abline(h=0, lty=2)
	invisible(sapply(1:length(precip), function(i) lines(rep(precip[i],2), r[i] + SE[i]*c(-1,1))))
	make_axes(xlab="Precipitation", ylab="Growth rate", main="(b) Hypothetical experimental results")
})

dev.off()


