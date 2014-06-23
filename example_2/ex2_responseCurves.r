# ex2_responseCurves.r
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
# produces a set of response curves for the climate variables in example 2
# presently these are not included as part of the manuscript


load("results/naiveModelResults.rdata")
load("results/integratedModelResults.rdata")
source("ex2_Functions.r")
library(rjags)


make_response_curve_plot <- function(intPosterior, naivePosterior, climDat, varName, title=varName, legend=T) {
	xvals <- seq(min(climDat[[varName]]), max(climDat[[varName]]), 0.02)
	ycols <- c(1, which(grepl(varName, colnames(intPosterior))))
	
	yInt <- response_curve(intPosterior[,ycols], xvals)
	yNaive <- response_curve(naivePosterior[,ycols], xvals)
	
	plot(xvals, yInt$estimate, col='red', lwd=2, ylim=c(0,1), xlab=varName, ylab="Suitability", type='l')
	polygon(c(xvals, rev(xvals)), c(yInt$lower, rev(yInt$upper)), col="#FF444444", border=NA)
	lines(xvals, yNaive$estimate, col='blue', lwd=2)
	polygon(c(xvals, rev(xvals)), c(yNaive$lower, rev(yNaive$upper)), col="#4444FF44", border=NA)

	if(legend) legend(4.5, 0.25, c("Integrated", "Naive"), col=c('red', 'blue'), lwd=2, cex=0.7, bty='n')
}

response_curve <- function(posterior, X, responseScale=TRUE, quantiles=c(0.05, 0.95)) {
	# constructs a response curve with credible intervals, with certain assumptions about the variable posterior:
	# posterior[,1] is the intercept
	# all subsequent columns represent increasing powers of the X variable, starting at 1 and skipping none
	# if nrow(posterior) is 1, returns only the estimate, with no credible intervals
	# this allows the function to be used recursively to calculate credible intervals
	
	postMeans <- sapply(1:ncol(posterior), function(i) mean(posterior[,i]))
	
	resultValue <- postMeans[1]
	for(i in 2:length(postMeans))
		resultValue <- resultValue + postMeans[i] * X^(i-1)
	
	if(responseScale) resultValue <- inv_logit(resultValue)
	
	
	# calculate credible interval
	if(nrow(posterior) > 1) {
		resultPosterior <- sapply(1:nrow(posterior), function(i) response_curve(posterior[i,,drop=FALSE], X, responseScale=F))
		resultQuant <- t(sapply(1:nrow(resultPosterior), function(i) quantile(resultPosterior[i,], quantiles)))
		resultValue <- list( estimate = resultValue, lower=resultQuant[,1], upper=resultQuant[,2])
		if(responseScale) {
			resultValue$lower <- inv_logit(resultValue$lower)
			resultValue$upper <- inv_logit(resultValue$upper)
		}
	}
	
	
	return(resultValue)
}


par(mfrow=c(2,2))
for(varName in unique(integratedModel$variables$varNames)[-1])
	make_response_curve_plot( integratedModel$posteriorSamples[[1]], naiveModel$posteriorSamples[[1]], integratedModel$mcmcData, varName, legend=(varName == "an_prcp"))
