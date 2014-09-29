# 1-setup_naive.r
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
# Sets up the naive model, including some data processing as well as model fitting
# Uses results of naive model to set up the input data for the integrated model
#  usage: Rscript 1-setup_naive.r
# results are automatically written to dat/naive_model.rdata
#
# see main() for an overview of the steps
# prep_data describes the preparation of the data
# fit_naive_model has the procedure for the stepwise regression
#
#


main = function()
{
## this program gets executed in global scope when running/sourcing this script
	# source("ex2_Functions.r")

	# set the seed for the random number generator, for reproducibility
	# comment this line out for random pseudo absence draws
	# this seed was chosen by a single call to:
	# runif(1,0,.Machine$integer.max)
	set.seed(626786234)
	
	maple = prep_data()
	
	# fit the naive model
	# we use the weighted presences as the response
	# we have selected 3 of 6 predictors due to collinearities in the predictors
	naiveModel = fit_naive_model(maple$calib, response=c(20,21), predictors=c(6,8,12))
	
	# finally, prep the information for the integrated model
	priors = data.frame(mean=c(naiveModel$model$coefficients), sd = summary(
		naiveModel$model)$coefficients[,2])
	# put the priors into the same order as the other pieces
	row.names(priors) = sapply(row.names(priors), function(x) {
		naiveModel$variables$parameter[naiveModel$variables$coefName == x]})
	priors = priors[c(1, order(row.names(priors)))[-(nrow(priors)+1)],]	

	inits = data.frame(inits = unlist(naiveModel$starts()))

	# unfortunately, this part must be done manually, and must match the model description
	# the order of the columns of intData must matchs the order in priors and inits
	# finally, the response (phenofit predictions) must come last
	intData = data.frame(
		ddeg1 = maple$calib$fut_ddeg,
		ddeg2 = maple$calib$fut_ddeg^2,
		ddeg3 = maple$calib$fut_ddeg^3,
		pToPET1 = maple$calib$fut_pToPET,
		pToPET2 = maple$calib$fut_pToPET^2,
		sum_prcp1 = maple$calib$fut_sum_prcp,
		sum_prcp2 = maple$calib$fut_sum_prcp^2,
		sum_prcp3 = maple$calib$fut_sum_prcp^3,
		phenofit = maple$calib$Phenofit_HadA2
	)

	save(maple, naiveModel, file="dat/naive_model.rdata")
	write.csv(priors, file='dat/integratedPriors.csv', row.names = FALSE)
	write.csv(inits, file='dat/integratedInits.csv', row.names = FALSE)
	write.csv(intData, file='dat/integratedData.csv', row.names = FALSE)
}



###################
#
#	Main helper functions
#
###################


fit_naive_model = function(dat, response=c(20,21), predictors=c(6,8,12))
{
#	fits the naive model using stepwise regression
#	dat: the data frame to be used for the analysis
#	response: the columns containing the weighted response data
#	predictors: columns containing the predictors
	require(glm2)
	
	stepScope = list(
		lower = cbind(weightedPresence, weightedN) ~ 1,
		upper = generate_formula(dat, response, predictors)
	)
	mod = step( glm2(stepScope$lower, family=binomial, data=dat), 
		scope=stepScope, direction="both", k = log(nrow(dat)))
	vars = variable_names(mod)
	startingVals = setup_starting_values(vars, mod)
	list(model=mod, variables=vars, starts = startingVals)


}


prep_data = function()
{
	# settings
	validationSize = 1/3

	mapleAll = read.csv("dat/AceSac.csv")
	load("dat/currentClim.rdata")
	load("dat/futureClim.rdata")

	require(sp)
	# make spatial objects and perform a spatial join
	mapleAll = SpatialPointsDataFrame(mapleAll[,1:2], mapleAll)
	current_clim = SpatialPointsDataFrame(current_clim[,1:2],current_clim)
	future_clim = SpatialPointsDataFrame(future_clim[,1:2],future_clim)
	mapleAll = cbind(mapleAll, over(mapleAll, current_clim), over(mapleAll, future_clim))

	## drop NAs & extra latlong columns
	dropRows = which( (apply(mapleAll, 1, function(x) sum(is.na(x)))) > 0)	# find rows with at least 1 NA
	mapleAll = mapleAll[-dropRows, c(1:5, 10:15, 18:23)]

	# add the ratio of annual precip to pet
	mapleAll = as.data.frame(append(mapleAll, list(pToPET = with(mapleAll, an_prcp/pet)), after=11))
	mapleAll$fut_pToPET = with(mapleAll, fut_an_prcp/fut_pet)

	## split data into calibration (maple$calib) and validation (maple$valid) sets
	## we do this by drawing 2/3 of the presences for calibration, plus enough absences to keep
	## the same ratio as the original dataset. For validation, we take the remaining presences
	## plus an equal number of absences. "Extra" absences are discarded
	mapleAllPres = mapleAll[mapleAll$PresObs==1,]
	mapleAllAbs = mapleAll[mapleAll$PresObs==0,]
	calibPresN = as.integer((1-validationSize) * nrow(mapleAllPres))
	calibAbsN = calibPresN * as.integer(nrow(mapleAllAbs) / nrow(mapleAllPres))
	calibPr = sample(nrow(mapleAllPres), calibPresN)
	calibAb = sample(nrow(mapleAllAbs), calibAbsN)
	maple = list(
		calib = rbind(mapleAllPres[calibPr,], mapleAllAbs[calibAb,]),
		valid = mapleAllPres[-calibPr,],
		all = mapleAll)
	validAb = sample(nrow(mapleAllAbs[-calibAb,]), nrow(maple$valid))
	maple$valid = rbind(maple$valid, mapleAllAbs[-calibAb,][validAb,])

	## transform data
	## all climate data gets zero-centered (based on the PRESENT climatic conditions, including for future climate)
	## phenofit predictions get "squeezed" to fit on the (0,1) interval instead of [0,1]
	maple$transformations = lapply(maple$calib[,6:12], rescale, return.functions=TRUE)
	for(nm in names(maple$transformations)) {
		nmfut = paste("fut_", nm, sep="")
		maple$calib[,nm] = maple$transformations[[nm]]$forward(maple$calib[,nm])
		maple$calib[,nmfut] = maple$transformations[[nm]]$forward(maple$calib[,nmfut])
	}

	maple$calib$Phenofit_CRU = smithson_transform(maple$calib$Phenofit_CRU)
	maple$calib$Phenofit_HadA2 = smithson_transform(maple$calib$Phenofit_HadA2)


	# weight data so that the total weight of presences and absences is the same
	wght = as.integer(sum(maple$calib$PresObs == 0) / sum(maple$calib$PresObs == 1))
	maple$calib$weightedPresence = wght*maple$calib$PresObs
	maple$calib$weightedN = rep(1, nrow(maple$calib))
	maple$calib$weightedN[maple$calib$PresObs == 1] = wght

	return(maple)
}




###################
#
#	Additional utility functions
#
###################

rescale = function(x, return.functions=FALSE)
{
	# centers data on 0 with a standard deviation of 1
	# if return.functions is TRUE, returns a list of two functions giving forward and backward transformations
	# if not, just returns the transformed data
	
	xbar = mean(x)
	xsd = sd(x)
	
	tr = list(
		forward = (function(x) (x - xbar) / xsd),
		backward = (function(x) (x * xsd) + xbar))
	
	if(return.functions) {
		return(tr)
	} else {
		return(tr$forward(x))
	}
}


generate_formula = function(data, response.column=1, predictor.columns=NULL, degree=3)
{
	## danger! There is no guarantee that this will work with non-numeric input data
	#
	# This function builds a formula from the input data for use in glm and step type models
	# data: data frame to be used for building the matrix
	# response.column: single value or a pair of values (for binomial models) giving the
	#	index of the column(s) to be used for the response variable
	# predictor.columns: vector of column indices to be used as predictors. if omitted
	#	all columns except the response will be used
	# degree: the maximum degree of predictors to be used in the model
	#
	# the formula does not produce interactions among predictors  

	expandedFormula = c()
	if(is.null(predictor.columns)) {
		allCols = 1:(dim(data)[2])
		predictor.columns = allCols[-which(allCols == response.column)]
	}
	
	varNames = names(data)[predictor.columns]
	
	for(variable in varNames) {
		currentFormula = variable
		if(degree > 1) {
			for(i in 2:degree)
				currentFormula = paste(currentFormula, " + I(", variable, "^", i, ")", sep="")
		}
		expandedFormula = c(expandedFormula, currentFormula)
	}
	
	if(length(response.column) == 1) {
		response = names(data)[response.column]
	} else {
		response = paste("cbind(", names(data)[response.column[1]], ",", names(data)[response.column[2]], ")", sep="")
	}
	
	collapsedFormula = paste(response, " ~ ", paste(expandedFormula, collapse = " + "))
	
	# turn the text string of the formula into an actual formula object
	collapsedFormula = eval(parse(text=collapsedFormula))
	return(collapsedFormula)
	
}


smithson_transform = function(x, N = length(x), s = 0.5) 
{
	# transformation from supplemental material of Smithson and Verkuilen 2006
	# takes data on interval [0,1] and transforms to (0,1)
	# s is the scaling factor given in the same reference
	val = (x * (N-1) + s)/N
	return(val)
}

variable_names = function(model)
{
#
#	Simple utility function; given a glm-type model, will produce an object listing
#	variable names, corresponding parameter names, and other bits of information that
#	will be useful for the model-fitting process
#
	coefNames = names(model$coefficients)
	modelTerms = unname(sapply(coefNames, function(x) {
		if(substr(x,1,2) == "I(") {
			return(substr(x, 3, nchar(x) - 1)) 
		} else {
			return(x)
		}
	}))
	
	variables = data.frame(
		powers = unname(sapply(modelTerms, function(x) {
			if(substr(x, nchar(x)-1, nchar(x)-1) == "^")
				return(as.integer(substr(x, nchar(x), nchar(x))))
			return(1)
		})),
		varNames = unname(sapply(modelTerms, function(x) {
			if(substr(x, nchar(x)-1, nchar(x)-1) == "^")
				return(substr(x, 1, nchar(x)-2))
			return(x)
		})),
		coefName = coefNames,
		stringsAsFactors = FALSE
	)
	
	variables$parameter = paste("b_", variables$varNames, variables$powers, sep="")

	## fix a bug where the intercept isn't sorted first on some OSs
	interceptRow = variables[1,]
	interceptRow$parameter[1] = "b0"
	variables = variables[-1,]
	
	# sort the variables without the intercept
	variables = variables[order(variables$varNames, variables$powers),]
	
	# return the intercept
	variables = rbind(interceptRow, variables)

	return(variables)
}


setup_starting_values = function(vars, model)
{
#
#	returns a function that will draw random starting values for the mcmc using the given
#	parameter names
#	vars: an object returned from the variable_names function
#	model: a glm-type model
#
	coefs = model$coefficients
	ses = summary(model)$coefficients[,2]
	valText = c()
	for(i in 1:nrow(vars)) {
		modInd = which(names(coefs) == vars$coefName[i])
		valText = c(valText, paste(vars$parameter[i], " = rnorm(1,", coefs[modInd], ",", ses[modInd], ")", sep=""))
	}
	valText = paste(valText, collapse=', ')
	valText = paste("function() list(", valText, ")", sep="")
	return(eval(parse(text=valText)))
}



main()