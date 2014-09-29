# 5-make_figure.r
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
# produce the figure for the second example
#
#




main = function()
{
	require(fields)
	require(rgdal)
	require(coda)
	
	load("dat/naive_model.rdata")
	intPredictions = read.csv("results/integratedStats.csv")
	load("results/integratedModel.rdata")
	
	predictions = process_predictions(maple, naiveModel, intPredictions, integratedModel)
	ocean = readOGR(dsn="dat/ne_110m_ocean", layer="ne_110m_ocean")

	pdfFileName = "ex2.pdf"
	pdfWidth = 7.5
	pdfHeight = 6
	titleCEX = 0.7
	
	# which columns to loop thru for mean and errors
	meanCols = c(3,5,7,11)
	SECols = meanCols + 1

	## pick the data to display
	latLimits = c(27,65)
	longLimits = c(-105, -55)
	predictionsSub = subset(predictions$all, in.range(lat, latLimits) & in.range(long, longLimits))
	mapleAllSub = subset(maple$all, in.range(lat, latLimits) & in.range(long, longLimits))
	latitude = predictionsSub[,2]
	longitude = predictionsSub[,1]
	meanPredictions = cbind(mapleAllSub$Phenofit_CRU, mapleAllSub$Phenofit_HadA2, predictionsSub[,meanCols])
	errorPredictions = predictionsSub[,SECols]

	## plotting settings
	meanZlims = c(0,1)
	meanColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), interpolate='spline', bias=1, space="rgb")(200)
	meanTitles = c("A. Phenofit: Present", "B. Phenofit: HadA2", "C. Naive SDM: Present",  "D. Naive SDM: HadA2", "G. Integrated: Present", "H. Integrated: HadA2")
	errorTitles = c("E. SDM: Present SE",  "F. SDM: HadA2 SE", "I. Integrated: Present SE", "J. Integrated: HadA2 SE")
	meanScaleTitle = "Suitability"

	errorZlims = range(errorPredictions)
	errorColors = colorRampPalette(c("#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"), interpolate='spline', space="rgb", bias=1.3)(200)
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
	dev.off()
}


process_predictions = function(maple, naiveModel, intPredictions, intPosterior)
{
#
#	builds an object containing the predictions of all models, both for all of the
#	original data points (in result$all) and for the validation points (result$valid)
#
#	maple: the same maple object saved in step 1
#	naiveModel: the glm object representing the naive model
#	intPredictions: a data frame of predictions (produced in step 4)
#	intPosterior: the posterior samples from the integrated model, produced by step 3
#
	require(pROC)

	# extract only the relevant climate variables
	presClimate = maple$all[,which(colnames(maple$calib) %in% 
		unique(naiveModel$variables$varNames))]
	futClimate = maple$all[,which(substr(colnames(maple$calib),5, nchar(colnames(maple$calib)))
		%in% unique(naiveModel$variables$varNames))]
	validationClimate = maple$valid[,which(colnames(maple$calib) %in% 
		unique(naiveModel$variables$varNames))]
	colnames(futClimate) = colnames(presClimate)

	# apply the transformations used in the calibration data to the projection datasets
	presClimate = as.data.frame(sapply(names(presClimate), function(name) {
		maple$transformations[[name]]$forward(presClimate[,name])}))
	futClimate = as.data.frame(sapply(names(futClimate), function(name) {
		maple$transformations[[name]]$forward(futClimate[,name])}))
	validationClimate = as.data.frame(sapply(names(validationClimate), function(name) {
		maple$transformations[[name]]$forward(validationClimate[,name])}))

	# produce predictions for the naive model
	naivePresPred = predict(naiveModel$model, newdata=presClimate, type='response', se.fit=TRUE)
	naiveFutPred = predict(naiveModel$model, newdata=futClimate, type='response', se.fit=TRUE)
	naiveValidPred = predict(naiveModel$model, newdata=validationClimate, type='response', se.fit=FALSE)

	# wrap everything into a data frame
	predictions = cbind(maple$all[,1:2], naivePresPred$fit, naivePresPred$se.fit, 
		naiveFutPred$fit, naiveFutPred$se.fit, intPredictions)
	colnames(predictions) = c("long", "lat", 
		'naivePresent', 'naivePresentSE', 'naiveFuture', 'naiveFutureSE',
		'intPresent', 'intPresentSE', 'intPresentLower', 'intPresentUpper',
		'intFuture', 'intFutureSE', 'intFutureLower', 'intFutureUpper')

	# validation
	intValidPred = predict_mcmc(intPosterior, validationClimate)
	validation = list(data = cbind(maple$valid[,1:5], naiveValidPred, intValidPred))
	colnames(validation$data) = c('long', 'lat', 'presence', 'phenofitPres', 'phenofitFut', 'naive', 'integrated')
	validation$naiveROC = roc(validation$data$presence, validation$data$naive)
	validation$phenofitROC = roc(validation$data$presence, validation$data$phenofitPres)
	validation$integratedROC = roc(validation$data$presence, validation$data$integrated)
	
	return(list(all=predictions, valid=validation))
}


predict_mcmc = function(posteriorSample, newData) 
{
#
#	Works like predict(), but for mcmc samples
#	Returns the average prediction for each point in the new dataset
#	
#	posteriorSample: mcmc object; the posterior sample from which to produce predictions
#	newData: data frame of predictors where predictions should be produced
#
	parPrefix = 2
	parEnd = 1	
	
	# predict the prob of presence for each point in newData
	# assumes the FIRST parameter in posteriorSample is the intercept
	predictedPrPresenceMean = mean(posteriorSample[,1])

	currentCol = 2
	while(currentCol <= ncol(posteriorSample)) {
		ndColumn = which(colnames(newData) == substr(colnames(posteriorSample)[currentCol], parPrefix + 1, nchar(colnames(posteriorSample)[currentCol]) - parEnd))
		exponent = as.integer(substr(colnames(posteriorSample)[currentCol], nchar(colnames(posteriorSample)[currentCol]), nchar(colnames(posteriorSample)[currentCol])))
		predictedPrPresenceMean = predictedPrPresenceMean + mean(posteriorSample[,currentCol]) * newData[,ndColumn]^exponent		
		currentCol = currentCol + 1
	}
	predictedPrPresenceMean = inv_logit(predictedPrPresenceMean)	
	return(predictedPrPresenceMean)
}


inv_logit = function(x) {
	inf.ind = which(x == Inf)		# negative infinity works normally, but Inf returns NaN
	result = exp(x) / (1+exp(x))
	result[inf.ind] = 1
	return(result)
}

in.range = function(x, lims) x >= lims[1] & x <= lims[2]


main()