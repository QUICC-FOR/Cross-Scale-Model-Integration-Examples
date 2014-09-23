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
# run the model relating the second model to the metamodel

source("ex1_globals.r")
NPredict <- 50	# how many data points for prediction/plotting


## the "data" from an experiment relating manipulated precipitation to growth rate
m2SimData <- data.frame(
	precip = c(-1, -0.5, 0, 0.5, 1),
	r      = c(-0.10, 0.08, 0.23, 0.26, 0.13),
	SE     = c(0.27, 0.29, 0.31, 0.25, 0.27))

m2SimData$psi <- 1-pnorm(0,m2SimData$r, m2SimData$SE)

# calculate the upscaled parameters from model 2
m2Data <- list(psi=m2SimData$psi, precip=m2SimData$precip)
m2Pars <- c('a0', 'a1', 'a2', 'phi')
model2 <- jags.model("ex1_model2.jags", data = m2Data, n.chains=settings$chains, n.adapt=settings$tuning)
# burn in the model
update(model2, settings$burnin)
if(settings$diagnostics) {
	test <- coda.samples(model2, m2Pars, settings$samples, thin=settings$thin)
	print(gelman.diag(test))
}
m2Results <- coda.samples(model2, m2Pars, settings$samples, thin=settings$thin)


precipPredictions <- seq(min(m2Data$precip), max(m2Data$precip), length.out=NPredict)
m2Predictions <- cbind(precipPredictions, predict_psi(m2Results, precipPredictions))
colnames(m2Predictions)[1] <- 'precip'

save(m2SimData, m2Results, m2Predictions, file="dat/ex1_m2.rdata")




