# ex2_prepIntegrated.r
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
# prepares the data for the integrated model

load('results/naiveModel.rdata')
load('dat/maple.rdata')

if(interactive()) summary(naiveModel)

priors = data.frame(mean=c(naiveModel$coefficients), sd = summary(naiveModel)$coefficients[,2])
# do some wizardry to get things from the naiveModel (which is in whatever order step felt
# like putting it into) into the same order as the initial values etc:
row.names(priors) = sapply(row.names(priors), function(x) variables$parameter[variables$coefName == x])
priors = priors[c(1, order(row.names(priors)))[-(nrow(priors)+1)],]
write.csv(priors, file='dat/integratedPriors.csv', row.names = FALSE)

inits = data.frame(inits = unlist(starting_values()))
write.csv(inits, file='dat/integradtedModelInits.csv', row.names = FALSE)

## unfortunately, this part must be done manually, and must match the model description
## it is very important that the order of the columns matches the order in priors and inits
## finally, the response (phenofit predictions) must come last
intData = data.frame(
	ddeg1 = maple$ddeg,
	ddeg2 = maple$ddeg^2,
	ddeg3 = maple$ddeg^3,
	pToPET1 = maple$pToPET,
	pToPET2 = maple$pToPET^2,
	sum_prcp1 = maple$sum_prcp,
	sum_prcp2 = maple$sum_prcp^2,
	sum_prcp3 = maple$sum_prcp^3,
	phenofit = maple$Phenofit_HadA2
)
write.csv(intData, file='dat/integratedModelData.csv', row.names = FALSE)
