# ex2_mcmcIntegrated.r
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
# this script sets up a model of the same form as the one chosen by the stepwise procedure
# with the addition of conditioning the model on the predictions for the present from phenofit
# then runs it in JAGS to estimate posterior distributions for all parameters
# this is the integrated model

source("ex2_Functions.r")
load("results/stepResults.rdata")
load("dat/maple.rdata")
library(rjags)
settings <- set_mcmc_settings()

# set up an object to hold everything related to this model, including posterior samples
# parameters, data, convergence checks, etc
integratedModel <- list(
	variables = variables,
	startingValues = starting_values,
	allData = maple,
	modelFilename="integrated_model.jags"
)

# set up the model file
integratedModel$modelText = make_model_text(integratedModel$variables, linearPredictor = integrated_linear_predictor, prior = integrated_prior)
cat(integratedModel$modelText, file=integratedModel$modelFilename)

# set up the data
# jags wants only data that are used in the analysis, so use the variables$varNames (minus the intercept) to determine which
integratedModel$mcmcData <- lapply(as.character(unique(variables$varNames)[-1]), function(x) integratedModel$allData[[x]])
names(integratedModel$mcmcData) <- unique(variables$varNames)[-1]
integratedModel$mcmcData$weightedPresence <- integratedModel$allData$weightedPresence
integratedModel$mcmcData$weightedN <- integratedModel$allData$weightedN
integratedModel$mcmcData$phenofit <- integratedModel$allData$Phenofit_CRU
integratedModel$mcmcData$N <- nrow(integratedModel$allData)


# set up the jags model
# start model and run initial burnin
# for some reason, inits were breaking this, so we are letting jags pick
load.module("glm")
#integratedModel$jagsModel <- jags.model(integratedModel$modelFilename, data = integratedModel$mcmcData, n.chains=settings$n.chains, n.adapt=settings$n.adapt, inits = integratedModel$startingValues)
integratedModel$jagsModel <- jags.model(integratedModel$modelFilename, data = integratedModel$mcmcData, n.chains=settings$n.chains, n.adapt=settings$n.adapt)
update(integratedModel$jagsModel, settings$burninMin)
integratedModel$burninLength <- settings$burninMin
integratedModel$converged <- FALSE

# monitor additional terms introduced by the integration
integratedModel$paramsToMonitor <- c(integratedModel$variables$parameter, 'phi', 'b_Ph0', 'b_Ph1')

#### try to converge the model
done <- FALSE
tries <- 0
sampleSize <- settings$startingSampleSize
while(!done) {
	integratedModel$gdSamples <- coda.samples(integratedModel$jagsModel, integratedModel$paramsToMonitor, sampleSize, thin=settings$thin)
	tries <- tries + 1
	integratedModel$burninLength <- integratedModel$burninLength + sampleSize
	integratedModel$gd <- gelman.diag(integratedModel$gdSamples)
	print(paste("Tried an additional ", sampleSize, " samples; total burnin now = ", integratedModel$burninLength, sep=""))
	print(integratedModel$gd)
	flush.console()

	if(integratedModel$gd$mpsrf <= settings$mpsrfThreshold) {  
		done <- TRUE
		integratedModel$converged <- TRUE
		msg <- "Converged; taking final posterior samples"
	} else if(integratedModel$burninLength > settings$burninMax) {
		done <- TRUE
		msg <- "convergence failure; some estimates may be unstable"
	} else if(tries == 5) {
			sampleSize <- as.integer(1.2 * sampleSize)
			tries <- 0
	}
}

print(msg)
flush.console()


# set up new single-chain model
# get posterior samples
integratedModel$jagsModel <- jags.model(integratedModel$modelFilename, data = integratedModel$mcmcData, n.chains=1, n.adapt=settings$n.adapt)

# burnin is whatever we got for burnin last time, times 1.5, rounded to the nearest 10,000
update(integratedModel$jagsModel, round(integratedModel$burninLength * 1.5, -4))

integratedModel$posteriorSamples <- coda.samples(integratedModel$jagsModel, integratedModel$paramsToMonitor, n.iter=settings$finalSampleSize, thin=settings$thin)

save(integratedModel, file="results/integratedModelResults.rdata")