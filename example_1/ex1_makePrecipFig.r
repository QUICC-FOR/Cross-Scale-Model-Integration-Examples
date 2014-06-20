# ex1_makePrecipFig.r
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
# produces a figure showing results of integration as a function of precipitation


## Rscript usage
## rscript ex1_results.r filename.pdf
## will output plots to pdf

## setup
library(rjags)
height <- 2.5

precip_plot <- function(x, y, ylim=c(0,1), xlab="Precipitation", ylab=expression(psi), main="") {
	# x is a vector giving the precip values for plotting
	# y is a data frame from the predict_psi function, giving mean, SD and quantiles of the predictions for each x value
	
	plot(x, y$mean, type='l', ylim=ylim, xaxt='n', yaxt='n', xlab='', ylab='')
	make_axes(xlab=xlab, ylab=ylab, main=main)
	lines(x, y$lower, lty=2)
	lines(x, y$upper, lty=2) 
}


filename <- commandArgs(TRUE)[1]
source("ex1_globals.r")
pdf(file=filename, width=width, height=height)


par(mfrow=c(1,3))
par(mar=margins, xpd=FALSE, bty=bty, mgp=mgp, tcl=tcl)

# get data from disk
load("dat/ex1_m1.rdata")
load("dat/ex1_m2.rdata")
load("dat/ex1_mm.rdata")

yl <- c(
	min(c(m1Predictions$precipPredict$lower, m2Predictions$lower, mmPredictions$precipPredict$lower)), 
	max(c(m1Predictions$precipPredict$upper, m2Predictions$upper, mmPredictions$precipPredict$upper)))
yl <- round(yl,2)
with(m1Predictions, precip_plot(precipDomain$precip, precipPredict, ylim=yl, main='(a) Naive model'))
lines(c(0,1), c(0.99,0.99), col='blue')
lines(c(0,0), c(0.98, 1), col='blue')
lines(c(1,1), c(0.98, 1), col='blue')

precip_plot(m2Predictions$precip, m2Predictions, ylim=yl, main='(b) Mechanistic model')
with(mmPredictions, precip_plot(precipDomain$precip, precipPredict, ylim=yl, main='(c) Integrated Model'))

dev.off()