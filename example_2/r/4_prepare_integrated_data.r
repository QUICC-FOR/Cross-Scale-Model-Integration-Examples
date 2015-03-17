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

library(coda)
naivePosterior = readRDS("results/naivePosterior.rds")
rawData = readRDS("dat/rawData.rds")
naiveSummary = summary(naivePosterior)

intPriors = data.frame(
	mean = naiveSummary$statistics[,1],
	sd = naiveSummary$statistics[,2],
	dist = rep(0, length(naiveSummary$statistics[,1])),
	row.names=rownames(naiveSummary$statistics))
write.csv(intPriors, file='dat/mcmc/integratedPriors.csv', row.names = FALSE)

intInits = sapply(1:nrow(intPriors), function(x) rnorm(1, intPriors$mean[x], intPriors$sd[x]))
write.csv(intInits, file='dat/mcmc/integratedInits.csv', row.names = FALSE)

integratedPresData = with(rawData$calib, data.frame(
	ddeg1 = ddeg,
	ddeg2 = ddeg^2,
	sum_prcp1 = sum_prcp,
	sum_prcp2 = sum_prcp^2,
	pToPET1 = pToPET,
	pToPET2 = pToPET^2,
	pres = Phenofit_CRU,
	count = rep(1, length(Phenofit_CRU))))
write.csv(integratedPresData, file="dat/mcmc/integratedPresData.csv", row.names = FALSE)

integratedFutureData = with(rawData$calib, data.frame(
	ddeg1 = fut_ddeg,
	ddeg2 = fut_ddeg^2,
	sum_prcp1 = fut_sum_prcp,
	sum_prcp2 = fut_sum_prcp^2,
	pToPET1 = fut_pToPET,
	pToPET2 = fut_pToPET^2,
	pres = Phenofit_HadA2,
	count = rep(1, length(Phenofit_HadA2))))
write.csv(integratedFutureData, file="dat/mcmc/integratedFutureData.csv", row.names = FALSE)

