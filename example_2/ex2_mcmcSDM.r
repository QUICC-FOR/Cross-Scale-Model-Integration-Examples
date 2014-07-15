# ex2_mcmcSDM.r
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
# run the naive model (i.e., the meta-model without any integration with phenofit)
# this script sets up a model of the same form as the one chosen by the stepwise procedure
# then runs it in JAGS to estimate posterior distributions for all parameters

source("ex2_Functions.r")
load("results/stepResults.rdata")
load("dat/maple.rdata")
library(rjags)
settings <- set_mcmc_settings()

# set up an object to hold everything related to this model, including posterior samples
# parameters, data, convergence checks, etc
naiveModel <- list(
	variables = variables,
	startingValues = starting_values,
	allData = maple,
	modelFilename="naive_model.jags"
)

# set up the model file
naiveModel$modelText = make_model_text(naiveModel$variables)
cat(naiveModel$modelText, file=naiveModel$modelFilename)

# set up the data
# jags wants only data that are used in the analysis, so use the variables$varNames (minus the intercept) to determine which
naiveModel$mcmcData <- lapply(as.character(unique(variables$varNames)[-1]), function(x) naiveModel$allData[[x]])
names(naiveModel$mcmcData) <- unique(variables$varNames)[-1]
naiveModel$mcmcData$PresObs <- naiveModel$allData$PresObs
naiveModel$mcmcData$N <- nrow(naiveModel$allData)


# set up the jags model
# start model and run initial burnin
load.module("glm")
# for some reason, inits were breaking this, so we are letting jags pick
#naiveModel$jagsModel <- jags.model(naiveModel$modelFilename, data = naiveModel$mcmcData, n.chains=settings$n.chains, n.adapt=settings$n.adapt, inits = naiveModel$startingValues)
naiveModel$jagsModel <- jags.model(naiveModel$modelFilename, data = naiveModel$mcmcData, n.chains=settings$n.chains, n.adapt=settings$n.adapt)
update(naiveModel$jagsModel, settings$burninMin)
naiveModel$burninLength <- settings$burninMin
naiveModel$converged <- FALSE


#### try to converge the model
done <- FALSE
tries <- 0
sampleSize <- settings$startingSampleSize
while(!done) {
	naiveModel$gdSamples <- coda.samples(naiveModel$jagsModel, naiveModel$variables$parameter, sampleSize, thin=settings$thin)
	tries <- tries + 1
	naiveModel$burninLength <- naiveModel$burninLength + sampleSize
	naiveModel$gd <- gelman.diag(naiveModel$gdSamples)
	print(paste("Tried an additional ", sampleSize, " samples; total burnin now = ", naiveModel$burninLength, sep=""))
	print(naiveModel$gd)
	flush.console()

	if(naiveModel$gd$mpsrf <= settings$mpsrfThreshold) {  
		done <- TRUE
		naiveModel$converged <- TRUE
		msg <- "Converged; taking final posterior samples"
	} else if(naiveModel$burninLength > settings$burninMax) {
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
naiveModel$jagsModel <- jags.model(naiveModel$modelFilename, data = naiveModel$mcmcData, n.chains=1, n.adapt=settings$n.adapt)

# burnin is whatever we got for burnin last time, times 1.5, rounded to the nearest 10,000
update(naiveModel$jagsModel, round(naiveModel$burninLength * 1.5, -4))

naiveModel$posteriorSamples <- coda.samples(naiveModel$jagsModel, naiveModel$variables$parameter, n.iter=settings$finalSampleSize, thin=settings$thin)

save(naiveModel, file="results/naiveModelResults.rdata")