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
load("results/naiveModelResults.rdata")
load("results/integratedModelResults.rdata")
source("ex2_Functions.r")


presClimate = mapleAll[,which(colnames(maple) %in% unique(naiveModel$variables$varNames))]
futClimate = mapleAll[,which(substr(colnames(maple),5, nchar(colnames(maple))) %in% unique(naiveModel$variables$varNames))]
validationClimate = mapleValidation[,which(colnames(maple) %in% unique(naiveModel$variables$varNames))]
colnames(futClimate) = colnames(presClimate)



# naive predictions
naivePresPred = process_output(naiveModel$posteriorSamples[[1]], transformations, newData=presClimate, do.transform=TRUE)
naiveFutPred = process_output(naiveModel$posteriorSamples[[1]], transformations, newData=futClimate, do.transform=TRUE)
naiveValidPred = process_output(naiveModel$posteriorSamples[[1]], transformations, newData = validationClimate, do.transform=TRUE, SE=FALSE, credInterval=FALSE)

# integrated predictions
# drop integration terms from the posterior
intPosterior = integratedModel$posteriorSamples[[1]][,colnames(integratedModel$posteriorSamples[[1]]) %in% integratedModel$variables$parameter]
intPresPred = process_output(intPosterior, transformations, newData=presClimate, do.transform=TRUE)
intFutPred = process_output(intPosterior, transformations, newData=futClimate, do.transform=TRUE)
intValidPred = process_output(intPosterior, transformations, newData = validationClimate, do.transform=TRUE, SE=FALSE, credInterval=FALSE)

predictions = cbind(maple[,1:2], naivePresPred, naiveFutPred, intPresPred, intFutPred)

colnames(predictions) = c("long", "lat", 
	'naivePresent', 'naivePresentSE', 'naivePresentLower', 'naivePresentUpper',
	'naiveFuture', 'naiveFutureSE', 'naiveFutureLower', 'naiveFutureUpper',
	'intPresent', 'intPresentSE', 'intPresentLower', 'intPresentUpper',
	'intFuture', 'intFutureSE', 'intFutureLower', 'intFutureUpper')

validation = list(data = cbind(mapleValidation[,1:3], naiveValidPred, intValidPred))
colnames(validation$data) = c('long', 'lat', 'presence', 'naive', 'integrated')
validation$naiveROC = ROC(validation$data$presence, validation$data$naive)
validation$integratedROC = ROC(validation$data$presence, validation$data$integrated)

save(predictions, validation, file="results/predictions.rdata")