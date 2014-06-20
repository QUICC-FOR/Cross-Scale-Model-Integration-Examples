# ex1_m1.r
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
# run the meta-model with no integration


## set the RNG seed for reproducibility
## disable this for a true random run
seed <- 2158798
set.seed(seed)

library(sp)
source("ex1_globals.r")

# determine the "true" probability of presence given a set of temperature and precipitation
# values. parameters of the model are arbitrary, but they give reasonable results when T and P are on [-1,1]
# assumes that temp and precip are the only things affecting presence
presence <- function(T, P) {
	b <- c(2,1.5,-2,0.5,-3.5)
	psi <- logit_inv(b[1] + b[2] * T + b[3] * T^2 + b[4] * P + b[5] * P^2)
}

sampling <- function(N, rangeP, clustering = 0.25, ...) {
	# clustering: integer between 0 and 1 controlling the amount of clustering in the sample
	#   clustering is the ratio of clusters to sample size, so the amount of clustering increases as the value decreases
	#   when clustering == 1, the sample is completely spatially random
	#	when clustering == 0, all points will be drawn from one randomly located cluster
  
    sgrid = SpatialPoints(expand.grid(precip = seq(rangeP[1], rangeP[2], 0.1), temp = seq(-1,1,0.1)))

	if(clustering < 1) {
		nc <- N * clustering
		if(nc <= 0) nc <- 1
		type <- 'clustered'
	} else {
		type <- 'random'
		nc <- 1
	}
    sampledPts = spsample(sgrid,n=N, type=type, nclusters=as.integer(nc), ...)
    temp = as.data.frame(sampledPts)
    colnames(temp) = c("precip", "temp")
    return(temp)
    
}


# constants
N <- 100 # number of data points for the SDM

# generate some simulated data, representing a presence-absence sampling technique
m1SimData <- sampling(N, precipRegimes$current, clustering = 0.2, iter=10)
m1SimData$presence <- rbinom(nrow(m1SimData), size=1, p=presence(m1SimData$temp, m1SimData$precip))



# run a simple bayesian logistic regression on the simulated data

bPrior <- matrix(c(
	0.0, 1.0E-4,	# mean and tau for b0 (intercept)
	0.0, 1.0E-4,	# b1 (first temp parameter)
	0.0, 1.0E-4,	# b2 (temp^2 parameter)
	0.0, 1.0E-4,	# b3 (precip)
	0.0, 1.0E-4),	# b4 (precip^2)
	byrow=TRUE, ncol=2)

m1Data <- list(precip=m1SimData$precip, temp = m1SimData$temp, presence = m1SimData$presence, bPrior=bPrior)
m1Pars <- c('b0', 'b1', 'b2', 'b3', 'b4')


# create the model
model1 <- jags.model('ex1_metamodel.jags', data = m1Data, n.chains=settings$chains, n.adapt=settings$tuning )
update(model1, settings$burnin)
if(settings$diagnostics) {
	test <- coda.samples(model1, m1Pars, settings$samples, thin=settings$thin)
	print(gelman.diag(test))
}

m1Results <- coda.samples(model1, m1Pars, settings$samples, thin=settings$thin)
m1Predictions <- list(
	map = predict_psi(m1Results, predictionMap$temp, predictionMap$precip, quantiles=settings$quantiles),
	mapDomain = predictionMap,
	precipPredict = predict_psi(m1Results, rep(0, length(precipPrediction)), precipPrediction, quantiles=settings$quantiles),
	precipDomain = data.frame(temp =  rep(0, length(precipPrediction)), precip=precipPrediction))

save(m1SimData, m1Results, m1Predictions, file="dat/ex1_m1.rdata")