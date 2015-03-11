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

# set the seed for the random number generator, for reproducibility
# comment this line out for random pseudo absence draws
# this seed was chosen by a single call to:
# runif(1,0,.Machine$integer.max)
set.seed(626786234)

# determine which predictors to use
predictors=c("ddeg", "sum_prcp", "pToPET")


library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-v", "--validation", default=1/3, type="double", help="proportion of presences to use for validation")
parser$add_argument("-d", "--datafile", default = "dat/raw/AceSac.csv", help="input file name")
parser$add_argument("-o", "--outfile", default = "dat/mapleDat_processed.rds", help="RDS output file name")
parser$add_argument("-p", "--prior", default='dat/mcmc/naivePriors.csv', help="output file for naive priors")
parser$add_argument("-n", "--inits", default='dat/mcmc/naiveInits.csv', help="output file for naive inits")
argList = parser$parse_args()

rescale = function(x, return.functions=FALSE, stdev = 1)
{
	# centers data on 0 with a standard deviation of stdev
	# if return.functions is TRUE, returns a list of two functions giving forward and backward transformations
	# if not, just returns the transformed data
	xbar = mean(x)
	xsd = sd(x)
	tr = list(
		forward = (function(x) stdev * (x - xbar) / xsd),
		backward = (function(x) (x * xsd / stdev) + xbar))
	if(return.functions) {
		return(tr)
	} else {
		return(tr$forward(x))
	}
}

smithson_transform = function(x, N = length(x), s = 0.5) 
{
	# transformation from supplemental material of Smithson and Verkuilen 2006
	# takes data on interval [0,1] and transforms to (0,1)
	# s is the scaling factor given in the same reference
	val = (x * (N-1) + s)/N
	return(val)
}

rawData = read.csv(argList$datafile)
load("dat/raw/currentClim.rdata")
load("dat/raw/futureClim.rdata")

# make the input data spatial and use a spatial join to relate them
library(sp)
rawData = SpatialPointsDataFrame(rawData[,1:2], rawData)
current_clim = SpatialPointsDataFrame(current_clim[,1:2],current_clim)
future_clim = SpatialPointsDataFrame(future_clim[,1:2],future_clim)
rawData = cbind(rawData, over(rawData, current_clim), over(rawData, future_clim))

## drop NAs & extra latlong columns
rawData = rawData[complete.cases(rawData),c(1:5, 10:15, 18:23)]

# add the ratio of annual precip to pet
rawData = as.data.frame(append(rawData, list(pToPET = with(rawData, an_prcp/pet)), after=11))
rawData = within(rawData, {fut_pToPET = fut_an_prcp/fut_pet})


# select only the predictors we're going to be using, along with squared and cubed terms
rawData = rawData[,c(1:5, sapply(predictors, grep, colnames(rawData)))]
# predictors2 = predictors3 = rawData[,-(1:5)]
# predictors2 = predictors2^2
# # predictors3 = predictors3^3
# colnames(predictors2) = paste(colnames(predictors2), "2", sep='')
# # colnames(predictors3) = paste(colnames(predictors3), "3", sep='')
# # rawData = cbind(rawData, predictors2, predictors3)
# rawData = cbind(rawData, predictors2)

## split data into calibration and validation sets
## we do this by drawing 2/3 of the presences for calibration, plus enough absences to keep
## the same ratio as the original dataset. For validation, we take the remaining presences
## plus an equal number of absences. "Extra" absences are discarded
presences = rawData[rawData$PresObs==1,]
absences = rawData[rawData$PresObs==0,]
calibPresN = as.integer((1-argList$validation) * nrow(presences))
calibAbsN = calibPresN * as.integer(nrow(absences) / nrow(presences))
calibPr = sample(nrow(presences), calibPresN)
calibAb = sample(nrow(absences), calibAbsN)

allData = list(
		calib = rbind(presences[calibPr,], absences[calibAb,]),
		valid = presences[-calibPr,],
		all = rawData)
validAb = sample(nrow(absences[-calibAb,]), nrow(allData$valid))
allData$valid = rbind(allData$valid, absences[-calibAb,][validAb,])

# center and scale the data
# we scale all parameters to a mean of 0 and sd of 0.5 following Gelman et al 2008
# we apply the scaling to the calibration dataset, then save it to apply the same scaling
# to all projection and validation data later

# the transformations object that gets appended to the data list includes functions
# for performing the forward and backward transformations
# finally, we also use a transformation to move the phenofit observations from the closed 
# [0,1] interval to the open (0,1) interval to allow some variance
prNames = colnames(allData$calib)[6:ncol(allData$calib)]
prNames = prNames[substr(prNames,1,4) != "fut_"]
allData$transformations = lapply(allData$calib[,prNames], rescale, return.functions=TRUE, stdev=1)
for(nm in names(allData$transformations))
{
	nmfut = paste("fut_", nm, sep="")
	allData$calib[,nm] = allData$transformations[[nm]]$forward(allData$calib[,nm])
	allData$calib[,nmfut] = allData$transformations[[nm]]$forward(allData$calib[,nmfut])
	allData$calib[,paste(nm,"2",sep="")] = allData$calib[,nm]^2
	allData$calib[,paste(nmfut,"2",sep="")] = allData$calib[,nmfut]^2
	
}

allData$calib$Phenofit_CRU = smithson_transform(allData$calib$Phenofit_CRU)
allData$calib$Phenofit_HadA2 = smithson_transform(allData$calib$Phenofit_HadA2)

# weight data so that the total weight of presences and absences is the same
wght = as.integer(sum(allData$calib$PresObs == 0) / sum(allData$calib$PresObs == 1))
allData$calib$weightedPresence = wght*allData$calib$PresObs
allData$calib$weightedN = rep(1, nrow(allData$calib))
allData$calib$weightedN[allData$calib$PresObs == 1] = wght

# save variable names in order for convenience
allData$variables = as.vector(sapply(predictors, function(x) colnames(allData$calib)[grep(paste("^", x, "[1-9]?", sep=""), colnames(allData$calib))]))

saveRDS(allData, argList$outfile)





## prepare data for MCMC
naivePriors = with(allData, data.frame(
		mean = rep(0, length(variables) + 1),
		sd = c(10, rep(2.5, length(variables))),
		dist = rep(1, length(variables) + 1),
		row.names = c("intercept", variables)))
write.csv(naivePriors, file=argList$prior, row.names = FALSE)

# use the prior distribution to draw random starts for the parameters
# naiveInits = sapply(1:nrow(naivePriors), function(i)
# {
# 	with(naivePriors[i,],
# 	{
# 		func = ifelse(dist == 0, rnorm, rcauchy)
# 		func(1, mean, sd)
# 	})
# })

# new idea
# draw starts from a fairly conservative gaussian
naiveInits = rnorm(nrow(naivePriors), 0, 1)
write.csv(naiveInits, file=argList$inits, row.names = FALSE)


presPredictors = allData$variables
futPredictors = paste("fut_", presPredictors, sep="")


naiveData_weighted = naiveData_unweighted = intData_Pres = allData$calib[,presPredictors]
intData_Fut = allData$calib[,futPredictors]
naiveData_weighted$response = allData$calib$weightedPresence
naiveData_unweighted$response = allData$calib$PresObs
intData_Pres$response = allData$calib$Phenofit_CRU
intData_Fut$response = allData$calib$Phenofit_HadA2
naiveData_weighted$weights = allData$calib$weightedN
naiveData_unweighted$weights = intData_Pres$weights = intData_Fut$weights = rep(1, nrow(allData$calib))

write.csv(intData_Fut, file='dat/mcmc/integratedData_Fut.csv', row.names = FALSE)
write.csv(intData_Pres, file='dat/mcmc/integratedData_Pres.csv', row.names = FALSE)
write.csv(naiveData_weighted, file='dat/mcmc/naiveData_weighted.csv', row.names = FALSE)
write.csv(naiveData_unweighted, file='dat/mcmc/naiveData_unweighted.csv', row.names = FALSE)

predictionDat = allData$all[,6:ncol(allData$all)]
predictionDat$ddeg = allData$transformations$ddeg$forward(predictionDat$ddeg)
predictionDat$sum_prcp = allData$transformations$sum_prcp$forward(predictionDat$sum_prcp)
predictionDat$pToPET = allData$transformations$pToPET$forward(predictionDat$pToPET)
predictionDat$fut_ddeg = allData$transformations$ddeg$forward(predictionDat$fut_ddeg)
predictionDat$fut_sum_prcp = allData$transformations$sum_prcp$forward(predictionDat$fut_sum_prcp)
predictionDat$fut_pToPET = allData$transformations$pToPET$forward(predictionDat$fut_pToPET)
write.table(predictionDat, "dat/mcmc/predictionData.csv", sep=",", col.names=FALSE, row.names=FALSE)

