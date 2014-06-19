# ex2_setUpSDM.r
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
# This script will perform a stepwise regression (similar to the biomod2 package)
# this is used to help quickly determine the variables to include in the SDM, as well as 
# to determine starting values to give to the MCMC algorithm

library(glm2)
library(fields)

load("dat/maplePA.rdata")
source("ex2_Functions.r")

responseColumn <- 3
predictorColumns <- 6:11
stepScope <- list( lower = PresObs ~ 1, upper = generate_formula(maplePA, responseColumn, predictorColumns))

# run a stepwise regression, using BIC because it is a bit more parsimonious than AIC
stepModel <- step( glm2(stepScope$lower, family=binomial, data=maplePA), scope=stepScope, direction="both", k = log(nrow(maplePA)))

# update the model to include all lower-order terms when higher-order terms are included
variables <- variable_names(stepModel)
for(nm in unique(variables$varNames)) {
	pows <- 1:(max(variables$power[variables$varNames==nm]))
	for(p in pows) {
		if(!(p  %in% variables$power[variables$varNames==nm])) {
			if(p > 1) {
				term <- paste("~. + I(", nm, "^", p, ")", sep="")
			} else {
				term <- paste("~. + ", nm, sep="")
			}
			stepModel <- update(stepModel, eval(parse(text=term)))
		}
	}
}

## view the predictions of the best model next to phenofit predictions; load maple data and apply transformations
load("dat/maple.rdata")
maple_fut <- as.data.frame(sapply(names(transformations), function(nm) transformations[[nm]]$forward(maple[,paste("fut_",nm,sep="")])))
maple_pres <- as.data.frame(sapply(names(transformations), function(nm) transformations[[nm]]$forward(maple[,nm])))

predictions <- data.frame(present=predict(stepModel, newdata=maple_pres, type="response"), 
	future=predict(stepModel, newdata=maple_fut, type="response"),
	pheno_present <- maple$Phenofit_CRU,
	pheno_future <- maple$Phenofit_HadA2)


## produce a plot of the results of the stepwise regression; uncomment to examine the intermediate results
# pdf("img/SDM_stepwiseTest.pdf", width=9,height=8.5)
# layout(matrix(c(1,2,5,3,4,5), byrow=T, nrow=2), widths=c(1,1,0.2))
# par(mar=c(1,1,1,1))
# titles <- c('sdm present', 'sdm HadA2 2100', 'Phenofit present', 'Phenofit HadA2 2100')
# for(i in c(1:4))
# 	quilt.plot(maple[,1], maple[,2], predictions[,i], col=bpy.colors(), zlim=c(0,1), main=titles[i], add.legend=F, xaxt='n', yaxt='n')
# plot.new()
# quilt.plot(maple[,1], maple[,2], predictions[,i], legend.only=T, col=bpy.colors(), zlim=c(0,1), add=F, smallplot=c(0.2,0.4, 0.05,0.95), bigplot=c(0,0,0,0))	
# dev.off()



## construct a list of variable names, parameter names, and starting values based on the results of the stepwise regression
variables <- variable_names(stepModel)
starting_values <- setup_starting_values(variables, stepModel)

save(stepModel, variables, starting_values, file="results/stepResults.rdata")