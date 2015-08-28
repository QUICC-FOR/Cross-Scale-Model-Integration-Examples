# ex1_m2.r
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
# run the integrated model


## set the RNG seed for reproducibility
## disable this for a true random run
seed <- 2158798
set.seed(seed)

source("ex1_globals.r")
load('dat/ex1_m1.rdata')
load('dat/ex1_m2.rdata')

# use model2 to compute priors for the metamodel
m2Stats <- as.data.frame(summary(m2Results)$statistics)
m2Stats$tau <- 1/(m2Stats$SD^2)

# factor by which to reduce the precision of the prior; is allows the user to decrease the confidence
# in the prior estimates
precisionReduction <- 20

bPrior <- matrix(c(
	m2Stats$Mean[1], m2Stats$tau[1]/precisionReduction,	# mean and tau for b0 (intercept)
	0.0, 1.0E-4,										# b1 (first temp parameter)
	0.0, 1.0E-4,										# b2 (temp^2 parameter) - this and the one above are uninformative because we only have info on precipitation
	m2Stats$Mean[2],m2Stats$tau[2]/precisionReduction,	# b3 (precip)
	m2Stats$Mean[3],m2Stats$tau[3]/precisionReduction),	# b4 (precip^2)
	byrow=TRUE, ncol=2)

## compute the integrated meta-model
mmData <- list(
	presence = m1SimData$presence,
	temp = m1SimData$temp,
	precip = m1SimData$precip,
	bPrior = bPrior)
mmPars <- c('b0', 'b1', 'b2', 'b3', 'b4')

metamodel <- jags.model("ex1_metamodel.jags", data = mmData,  n.chains=settings$chains, n.adapt=settings$tuning)

# burn in the model
update(metamodel, settings$burnin)
if(settings$diagnostics) {
	test <- coda.samples(metamodel, mmPars, settings$samples, thin=settings$thin)
	print(gelman.diag(test))
}
mmResults <- coda.samples(metamodel, mmPars, settings$samples, thin=settings$thin)

mmPredictions <- list(
	map = predict_psi(mmResults, predictionMap$temp, predictionMap$precip, quantiles=settings$quantiles),
	mapDomain = predictionMap,
	precipPredict = predict_psi(mmResults, rep(0, length(precipPrediction)), precipPrediction, quantiles=settings$quantiles),
	precipDomain = data.frame(temp =  rep(0, length(precipPrediction)), precip=precipPrediction))

save(mmResults, mmPredictions, file="dat/ex1_mm.rdata")
