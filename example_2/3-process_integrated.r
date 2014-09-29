# 3-process_integrated.r
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
# takes the raw mcmc output from the python program and turns it into something
# we can easily use in R
# requires coda package
#
#

library(coda)
load('dat/naive_model.rdata')

thinLength = 50
burnin = 500000


thin = function(x, n) {
	ind = seq(1, nrow(x),  n)
	return(x[ind,])
}

integratedModel = read.csv("results/integratedModel.csv", header=FALSE, stringsAsFactors=FALSE)
colnames(integratedModel) = naiveModel$variables$parameter
startVal = burnin + 1
endVal = nrow(integratedModel)
integratedModel = integratedModel[(burnin+1):nrow(integratedModel),]
integratedModel = thin(integratedModel, thinLength)
integratedModel = mcmc(integratedModel, start = startVal, end = endVal, thin = thinLength)


## prep the data to calculate the predictions
maplePredict = maple$all[,c(6,12,8,13,19,15)]
maplePredict$ddeg = maple$transformations$ddeg$forward(maplePredict$ddeg)
maplePredict$sum_prcp = maple$transformations$sum_prcp$forward(maplePredict$sum_prcp)
maplePredict$pToPET = maple$transformations$pToPET$forward(maplePredict$pToPET)
maplePredict$fut_ddeg = maple$transformations$ddeg$forward(maplePredict$fut_ddeg)
maplePredict$fut_sum_prcp = maple$transformations$sum_prcp$forward(maplePredict$fut_sum_prcp)
maplePredict$fut_pToPET = maple$transformations$pToPET$forward(maplePredict$fut_pToPET)


save(integratedModel, file="results/integratedModel.rdata")
write.table(integratedModel, "results/integratedModelThinned.csv", sep=',', col.names=FALSE, row.names=FALSE)
write.table(maplePredict, "dat/predictionData.csv", sep=",", col.names=FALSE, row.names=FALSE)


