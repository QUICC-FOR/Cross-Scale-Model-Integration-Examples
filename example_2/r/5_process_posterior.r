#!/usr/bin/Rscript
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
#
#
#

library(coda)
library(abind)

predict.mcmc = function(b, newdat)
{
	pres = with(newdat, plogis(b[1] + b[2]*ddeg + b[3]*ddeg^2 + b[4]*ddeg^3 + 
			b[5]*an_prcp + b[6]*an_prcp^2 + b[7]*pToPET + b[8]*pToPET^2))
	fut = with(newdat, plogis(b[1] + b[2]*fut_ddeg + b[3]*fut_ddeg^2 + b[4]*fut_ddeg^3 + 
			b[5]*fut_an_prcp + b[6]*fut_an_prcp^2 + b[7]*fut_pToPET + b[8]*fut_pToPET^2))
	data.frame(pres=pres, fut=fut)
}

coef.names = c('intercept', 'ddeg', 'ddeg2', 'ddeg3', 'an_prcp', 'an_prcp2', 'pToPET', 
		'pToPET2')

predictionGrid = readRDS("dat/predictionData.rds")

naive = read.csv("results/mcmc/naivePosterior.csv", h=F)
colnames(naive) = coef.names
naive = mcmc(naive, thin=50, start=50*20001)

intPres = read.csv("results/mcmc/integratedPresent.csv", h=F)
colnames(intPres) = coef.names
intPres = mcmc(intPres, thin=50, start=50*20001)

intFut = read.csv("results/mcmc/integratedFuture.csv", h=F)
colnames(intFut) = coef.names
intFut = mcmc(intFut, thin=50, start=50*20001)

naivePr = predict.mcmc(summary(naive)$statistics[,1], predictionGrid)
intPrePr = predict.mcmc(summary(intPres)$statistics[,1], predictionGrid)
intFutPr = predict.mcmc(summary(intFut)$statistics[,1], predictionGrid)

naivePosPr = apply(naive, 1, predict.mcmc, newdat=predictionGrid)
intPresPosPr = apply(intPres, 1, predict.mcmc, newdat=predictionGrid)
intFutPosPr = apply(intFut, 1, predict.mcmc, newdat=predictionGrid)

naivePosPr = abind(naivePosPr, along=3)
naiveSE = apply(naivePosPr, 1:2, sd)
naiveCI = apply(naivePosPr, 1:2, quantile, c(0.025, 0.975))
intPresPosPr = abind(intPresPosPr, along=3)
intPresSE = apply(intPresPosPr, 1:2, sd)
intPresCI = apply(intPresPosPr, 1:2, quantile, c(0.025, 0.975))
intFutPosPr = abind(intFutPosPr, along=3)
intFutSE = apply(intFutPosPr, 1:2, sd)
intFutCI = apply(intFutPosPr, 1:2, quantile, c(0.025, 0.975))

naivePosterior = naive
naivePredictions = cbind(naivePr, naiveSE, naiveCI['2.5%',,'pres'], 
		naiveCI['97.5%',,'pres'], naiveCI['2.5%',,'fut'], naiveCI['97.5%',,'fut'], 
		predictionGrid[,c('long', 'lat')])
colnames(naivePredictions) = c('pres', 'fut', 'presSE', 'futSE', 'pres2.5', 'pres97.5', 
	'fut2.5', 'fut97.5', 'long', 'lat')
intPresPosterior = intPres
intPresPredictions = cbind(intPrePr, intPresSE, intPresCI['2.5%',,'pres'], 
		intPresCI['97.5%',,'pres'], intPresCI['2.5%',,'fut'], intPresCI['97.5%',,'fut'], 
		predictionGrid[,c('long', 'lat')])
colnames(intPresPredictions) = colnames(naivePredictions)
intFutPosterior = intFut
intFutPredictions = cbind(intFutPr, intFutSE, intFutCI['2.5%',,'pres'], 
		intFutCI['97.5%',,'pres'], intFutCI['2.5%',,'fut'], intFutCI['97.5%',,'fut'], 
		predictionGrid[,c('long', 'lat')])
colnames(intFutPredictions) = colnames(naivePredictions)

save(naivePosterior, naivePredictions, intPresPosterior, intPresPredictions, 
		intFutPosterior, intFutPredictions, file="results/posteriors.rdata")
