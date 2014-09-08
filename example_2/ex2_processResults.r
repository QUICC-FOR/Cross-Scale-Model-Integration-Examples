# ex2_processResults.r
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
# produce posterior predictions and standard errors suitable for mapping from the 
# posterior parameter estimates of the meta-model
# note that this step is very memory intensive; expect multiple GBs of memory usage, with lots of swapping
# runtime can be quite long, especially if there isn't enough system RAM and lots of swapping is needed


load("dat/maple.rdata")
load("results/naiveModel.rdata")
load("results/integratedModel.rdata")
source("ex2_Functions.r")


# extract only the relevant climate variables
presClimate = mapleAll[,which(colnames(maple) %in% unique(variables$varNames))]
futClimate = mapleAll[,which(substr(colnames(maple),5, nchar(colnames(maple))) %in% unique(variables$varNames))]
validationClimate = mapleValidation[,which(colnames(maple) %in% unique(variables$varNames))]
validationFutClimate = mapleValidation[,which(substr(colnames(maple),5, nchar(colnames(maple))) %in% unique(variables$varNames))]
colnames(futClimate) = colnames(validationFutClimate) = colnames(presClimate)

# apply the transformations used in the calibration data to the projection datasets
presClimate = as.data.frame(sapply(names(presClimate), function(name) transformations[[name]]$forward(presClimate[,name])))
futClimate = as.data.frame(sapply(names(futClimate), function(name) transformations[[name]]$forward(futClimate[,name])))
validationClimate = as.data.frame(sapply(names(validationClimate), function(name) transformations[[name]]$forward(validationClimate[,name])))
validationFutClimate = as.data.frame(sapply(names(validationFutClimate), function(name) transformations[[name]]$forward(validationFutClimate[,name])))

# produce predictions for the naive model
naivePresPred = predict(naiveModel, newdata=presClimate, type='response', se.fit=TRUE)
naiveFutPred = predict(naiveModel, newdata=futClimate, type='response', se.fit=TRUE)
naiveValidPred = predict(naiveModel, newdata=validationClimate, type='response', se.fit=FALSE)

# integrated predictions
intPosterior = integratedModel
intPresPred = process_output(intPosterior, newData=presClimate)
intFutPred = process_output(intPosterior, newData=futClimate)

# validation
intValidPred = process_output(intPosterior, newData = validationClimate, SE=FALSE, credInterval=FALSE)
intValidPredFut = process_output(intPosterior, newData = validationFutClimate, SE=FALSE, credInterval=FALSE)

predictions = cbind(mapleAll[,1:2], naivePresPred$fit, naivePresPred$se.fit, naiveFutPred$fit, naiveFutPred$se.fit, intPresPred, intFutPred)

colnames(predictions) = c("long", "lat", 
	'naivePresent', 'naivePresentSE', 'naiveFuture', 'naiveFutureSE',
	'intPresent', 'intPresentSE', 'intPresentLower', 'intPresentUpper',
	'intFuture', 'intFutureSE', 'intFutureLower', 'intFutureUpper')

validation = list(data = cbind(mapleValidation[,1:5], naiveValidPred, intValidPred))
colnames(validation$data) = c('long', 'lat', 'presence', 'phenofitPres', 'phenofitFut', 'naive', 'integrated')

validation$integratedPresR2 = cor(validation$data$phenofitPres, validation$data$integrated)^2
validation$integratedPresR2 = cor(validation$data$phenofitPres, validation$data$integrated)^2
# this fails on older versions of R, so wrapped in try to finish the script
try({
	validation$naiveROC = ROC(validation$data$presence, validation$data$naive)
	validation$integratedROC = ROC(validation$data$presence, validation$data$integrated)
})

save(predictions, validation, file="results/predictions.rdata")