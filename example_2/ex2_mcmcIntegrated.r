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
# takes the raw mcmc output from the python program and turns it into something
# we can easily use in R
# requires coda package

library(coda)
load('results/naiveModel.rdata')

thinLength = 50
burnin = 500000


thin = function(x, n) {
	ind = seq(1, nrow(x),  n)
	return(x[ind,])
}

integratedModel = read.csv("results/integratedModel.csv", header=FALSE, stringsAsFactors=FALSE)
#for(i in 1:ncol(integratedModel)) integratedModel[,i] = as.numeric(integratedModel[,i])
colnames(integratedModel) = variables$parameter
startVal = burnin + 1
endVal = nrow(integratedModel)
integratedModel = integratedModel[(burnin+1):nrow(integratedModel),]
integratedModel = thin(integratedModel, thinLength)
write.table(integratedModel, "results/integratedModelThinned.csv", sep=',', col.names=FALSE, row.names=FALSE)
integratedModel = mcmc(integratedModel, start = startVal, end = endVal, thin = thinLength)
save(integratedModel, file="results/integratedModel.rdata")


## prep the data to calculate the predictions
load('dat/maple.rdata')
# maplePredict = mapleAll[,c(1,2,4,5,6,6,6,12,12,8,8,8,13,13,13,19,19,15,15,15)]
maplePredict = mapleAll[,c(6,12,8,13,19,15)]
maplePredict$ddeg = transformations$ddeg$forward(maplePredict$ddeg)
# maplePredict$ddeg.1 = maplePredict$ddeg^2
# maplePredict$ddeg.2 = maplePredict$ddeg^3
maplePredict$sum_prcp = transformations$sum_prcp$forward(maplePredict$sum_prcp)
# maplePredict$sum_prcp.1 = maplePredict$sum_prcp^2
# maplePredict$sum_prcp.2 = maplePredict$sum_prcp^3
maplePredict$pToPET = transformations$pToPET$forward(maplePredict$pToPET)
# maplePredict$pToPET.1 = maplePredict$pToPET^2
maplePredict$fut_ddeg = transformations$ddeg$forward(maplePredict$fut_ddeg)
# maplePredict$fut_ddeg.1 = maplePredict$fut_ddeg^2
# maplePredict$fut_ddeg.2 = maplePredict$fut_ddeg^3
maplePredict$fut_sum_prcp = transformations$sum_prcp$forward(maplePredict$fut_sum_prcp)
# maplePredict$fut_sum_prcp.1 = maplePredict$fut_sum_prcp^2
# maplePredict$fut_sum_prcp.2 = maplePredict$fut_sum_prcp^3
maplePredict$fut_pToPET = transformations$pToPET$forward(maplePredict$fut_pToPET)
# maplePredict$fut_pToPET.1 = maplePredict$fut_pToPET^2
write.table(maplePredict, "dat/predictionData.csv", sep=",", col.names=FALSE, row.names=FALSE)


