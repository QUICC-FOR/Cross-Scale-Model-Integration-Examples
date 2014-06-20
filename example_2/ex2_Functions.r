# ex2_Functions.r
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
# Supporting functions for the second example


set_mcmc_settings <- function() {
	# returns a list containing constants for the mcmc analysis
	list(
		n.chains = 4,
		n.adapt = 1000,
		burninMin = 50000,
		startingSampleSize = 10000,
		thin = 5,
		burninMax = 300000,
		finalSampleSize = 100000,
		mpsrfThreshold = 1.1
		)
}


process_output <- function(posteriorSample, transformations, newData, SE = TRUE, credInterval = TRUE, allreps = FALSE) {
	# posterior sample: a sample from a posterior distribution, e.g., from coda.samples, class is mcmc
	# transformations: a named list of lists; names correspond to variables, each list is a list of 2 functions, $forward and $backward that will apply data transformations
	# newData: a new dataset with names matching those in transformations
	
	# the number of characters in the prefix and end of the parameter names in posteriorSample
	# such that 
	# substr(colnames(posteriorSample), parPrefix + 1, nchar(colnames(posteriorSample)) - parEnd) == names(transformation) == colnames(newData)
	parPrefix <- 2
	parEnd <- 1


	# transform predictors
	newData <- sapply(names(newData), function(varName) transformations[[varName]]$forward(newData[,varName]))
	

	# predict the prob of presence for each point in newData
	# assumes the FIRST parameter in posteriorSample is the intercept
	predictedPrPresenceMean <- mean(posteriorSample[,1])
	predictedPrPresenceReps <- matrix( posteriorSample[,1], nrow = nrow(newData), ncol=nrow(posteriorSample) )
	currentCol <- 2
	while(currentCol <= ncol(posteriorSample)) {
		ndColumn <- which(colnames(newData) == substr(colnames(posteriorSample)[currentCol], parPrefix + 1, nchar(colnames(posteriorSample)[currentCol]) - parEnd))
		exponent <- as.integer(substr(colnames(posteriorSample)[currentCol], nchar(colnames(posteriorSample)[currentCol]), nchar(colnames(posteriorSample)[currentCol])))
		predictedPrPresenceMean <- predictedPrPresenceMean + mean(posteriorSample[,currentCol]) * newData[,ndColumn]^exponent


		# to calculate the standard error, we need to determine the predicted prob of presence at EACH point for each posterior replicate
		## memory hungry line here
		predictedPrPresenceReps <- predictedPrPresenceReps + sapply(posteriorSample[,currentCol], function(x) x * newData[,ndColumn]^exponent)
		
		currentCol <- currentCol + 1
	}
	predictedPrPresenceMean <- inv_logit(predictedPrPresenceMean)
	predictedPrPresenceReps <- inv_logit(predictedPrPresenceReps)
	
	result <- data.frame(prediction = predictedPrPresenceMean)
	if(SE) {
		predictedPrPresenceSE <- sapply(1:nrow(predictedPrPresenceReps), function(i) sd(predictedPrPresenceReps[i,]))
		result$SE <- predictedPrPresenceSE
	}
	if(credInterval) {
		quants <- t(sapply(1:nrow(predictedPrPresenceReps), function(i) quantile(predictedPrPresenceReps[i,], c(0.05, 0.95))))
		result$lower <- quants[,1]
		result$upper <- quants[,2]
	}
	if(allreps) {
		result <- list(statistics = result, posteriorPredictions = predictedPrPresenceReps)
	}
	
	return(result)
}


make_model_text <- function(vars, linearPredictor = naive_linear_predictor, prior = naive_prior, responseName = "PresObs") {
	mt <- c(
		"model {",
		"  for(i in 1:N) {",
		paste("    ", responseName, "[i] ~ dbern(pr[i])", sep=""),
		paste("    ", linearPredictor(vars), sep=""),
		"  }",
		prior(vars),
		"}")
	mt <- paste(mt, collapse="\n")
	return(mt)
}


naive_prior <- function(vars) {
	pr <- c()
	for(p in vars$parameter) {
		pstr <- paste("  ", p, "~dnorm(0,0.001)", sep="")
		pr <- c(pr, pstr)
	}
	return(paste(pr, collapse="\n"))
}

naive_linear_predictor <- function(vars) {
	mod <- ("logit(pr[i]) <- b0")
	if(nrow(vars) > 1) {
		for(i in 2:nrow(vars)) {
			mod <- paste(mod, " + ", vars$parameter[i], " * ", vars$varNames[i], "[i]", sep="")
			if(vars$powers[i] > 1) mod <- paste(mod, "^", vars$powers[i], sep="")
		}
	}
	return(mod)
}

integrated_linear_predictor <- function(vars, phenofitPredictionName = "phenofit") {
	mod <- naive_linear_predictor(vars)
	mod <- paste(mod, "\n", 
		phenofitPredictionName, "[i] ~ dbeta(p[i],q[i])\n",
		"p[i] <- pr[i] * phi\n",
		"q[i] <- (1-pr[i]) * phi\n",
		sep="")
		return(mod)
}

integrated_prior <- function(vars) {
	pr <- naive_prior(vars)
	pr <- paste(pr,
		"## uninformative prior for phi from Gelman 2006\n",
		"phi <- U^2\n",
		"U ~ dunif(0,50)", sep="")
	return(pr)
}

generate_formula <- function(data, response.column=1, predictor.columns=NULL, degree=3) {
	## danger! There is no guarantee that this will work with non-numeric input data

	expandedFormula <- c()
	if(is.null(predictor.columns)) {
		allCols <- 1:(dim(data)[2])
		predictor.columns <- allCols[-which(allCols == response.column)]
	}
	
	varNames <- names(data)[predictor.columns]
	
	for(variable in varNames) {
		currentFormula <- variable
		if(degree > 1) {
			for(i in 2:degree)
				currentFormula <- paste(currentFormula, " + I(", variable, "^", i, ")", sep="")
		}
		expandedFormula <- c(expandedFormula, currentFormula)
	}
	
	collapsedFormula <- paste(names(data)[response.column], " ~ ", paste(expandedFormula, collapse = " + "))
	
	# turn the text string of the formula into an actual formula object
	collapsedFormula <- eval(parse(text=collapsedFormula))
	return(collapsedFormula)
	
}


setup_starting_values <- function(vars, model) {
	# returns a function that will draw random starting values for the mcmc using the given parameter names
	coefs <- model$coefficients
	ses <- summary(model)$coefficients[,2]
	valText <- c()
	for(i in 1:nrow(vars)) {
		modInd <- which(names(coefs) == vars$coefName[i])
		valText <- c(valText, paste(vars$parameter[i], " = rnorm(1,", coefs[modInd], ",", ses[modInd], ")", sep=""))
	}
	valText <- paste(valText, collapse=', ')
	valText <- paste("function() list(", valText, ")", sep="")
	return(eval(parse(text=valText)))
}

variable_names <- function(model) {
	coefNames <- names(model$coefficients)
	modelTerms <- unname(sapply(coefNames, function(x) {
		if(substr(x,1,2) == "I(") {
			return(substr(x, 3, nchar(x) - 1)) 
		} else {
			return(x)
		}
	}))
	
	variables <- data.frame(
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
	
	variables$parameter <- paste("b_", variables$varNames, variables$power, sep="")

	## fix a bug where the intercept isn't sorted first on some OSs
	interceptRow <- variables[1,]
	interceptRow$parameter[1] <- "b0"
	variables <- variables[-1,]
	
	# sort the variables without the intercept
	variables <- variables[order(variables$varNames, variables$power),]
	
	# return the intercept
	variables <- rbind(interceptRow, variables)

	return(variables)
}

bpy.colors <- function (n = 100, cutoff.tails = 0.1, alpha = 1) 
{	# "borrowed" from package sp
    n <- as.integer(n[1])
    if (n <= 0) 
        return(character(0))
    if (cutoff.tails >= 1 || cutoff.tails < 0) 
        stop("cutoff.tails should be in [0, 1]")
    i = seq(0.5 * cutoff.tails, 1 - 0.5 * cutoff.tails, length = n)
    r = ifelse(i < 0.25, 0, ifelse(i < 0.57, i/0.32 - 0.78125, 
        1))
    g = ifelse(i < 0.42, 0, ifelse(i < 0.92, 2 * i - 0.84, 1))
    b = ifelse(i < 0.25, 4 * i, ifelse(i < 0.42, 1, ifelse(i < 
        0.92, -2 * i + 1.84, i/0.08 - 11.5)))
    rgb(r, g, b, alpha)
}

inv_logit <- function(x) {
	inf.ind <- which(x == Inf)		# negative infinity works normally, but Inf returns NaN
	result <- exp(x) / (1+exp(x))
	result[inf.ind] <- 1
	return(result)
}