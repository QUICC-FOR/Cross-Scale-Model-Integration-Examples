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
# takes the raw mcmc output and turns it into something
# we can easily use in R
# requires coda package
#
#


THIN_LEN = 50
BURN_IN = 500000


thin = function(x, n) {
	ind = seq(1, nrow(x),  n)
	return(x[ind,])
}


prep_integrated = function(im, naiveModel, burnin, thinInterval) {
	require(coda, quietly = TRUE)
	colnames(im) = naiveModel$variables$parameter
	startVal = burnin + 1
	endVal = nrow(im)
	im = im[(burnin+1):nrow(im),]
	im = thin(im, thinInterval)
	im = mcmc(im, start = startVal, end = endVal, thin = thinInterval)
}




## prep the predictors
load('dat/naive_model.rdata')
maplePredict = maple$all[,c(6,12,8,13,19,15)]
maplePredict$ddeg = maple$transformations$ddeg$forward(maplePredict$ddeg)
maplePredict$sum_prcp = maple$transformations$sum_prcp$forward(maplePredict$sum_prcp)
maplePredict$pToPET = maple$transformations$pToPET$forward(maplePredict$pToPET)
maplePredict$fut_ddeg = maple$transformations$ddeg$forward(maplePredict$fut_ddeg)
maplePredict$fut_sum_prcp = maple$transformations$sum_prcp$forward(maplePredict$fut_sum_prcp)
maplePredict$fut_pToPET = maple$transformations$pToPET$forward(maplePredict$fut_pToPET)
write.table(maplePredict, "dat/predictionData.csv", sep=",", col.names=FALSE, row.names=FALSE)



# try to read integrated model from disk
integratedModel_Pres = tryCatch(
	read.csv("results/integratedModel_Pres.csv", header=FALSE, stringsAsFactors=FALSE),
	error = function(e) {
		warning(e)
		return(NULL)
	}
)
integratedModel_Fut = tryCatch(
	read.csv("results/integratedModel_Fut.csv", header=FALSE, stringsAsFactors=FALSE),
	error = function(e) {
		warning(e)
		return(NULL)
	}
)


# prep the results from the present
if(!is.null(integratedModel_Pres)) {
	integratedModel_Pres = prep_integrated(integratedModel_Pres, naiveModel, BURN_IN, THIN_LEN)
	save(integratedModel_Pres, file="results/integratedModel_Pres.rdata")
	write.table(integratedModel_Pres, "results/integratedModel_Pres_Thinned.csv", sep=',', col.names=FALSE, row.names=FALSE)
}


# prep the results from the future
if(!is.null(integratedModel_Fut)) {
	integratedModel_Fut = prep_integrated(integratedModel_Fut, naiveModel, BURN_IN, THIN_LEN)
	save(integratedModel_Fut, file="results/integratedModel_Fut.rdata")
	write.table(integratedModel_Fut, "results/integratedModel_Fut_Thinned.csv", sep=',', col.names=FALSE, row.names=FALSE)
}





