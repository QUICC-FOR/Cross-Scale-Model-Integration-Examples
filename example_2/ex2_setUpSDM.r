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

load("dat/maple.rdata")
source("ex2_Functions.r")

responseColumns = 18:19
predictorColumns <- 6:10 #### NOTE: Dropping pet because it is highly correlated with ddeg and min_temp
futureColumns <- predictorColumns+6 
stepScope <- list( lower = cbind(weightedPresence, weightedN) ~ 1, upper = generate_formula(maple, responseColumns, predictorColumns))

# run a stepwise regression, using BIC because it is a bit more parsimonious than AIC
stepModel <- step( glm2(stepScope$lower, family=binomial, data=maple), scope=stepScope, direction="both", k = log(nrow(maple)))

# update the model to include all lower-order terms when higher-order terms are included
variables <- variable_names(stepModel)
for(nm in unique(variables$varNames)) {
	pows <- 1:(max(variables$powers[variables$varNames==nm]))
	for(p in pows) {
		if(!(p  %in% variables$powers[variables$varNames==nm])) {
			if(p > 1) {
				term <- paste("~. + I(", nm, "^", p, ")", sep="")
			} else {
				term <- paste("~. + ", nm, sep="")
			}
			stepModel <- update(stepModel, eval(parse(text=term)))
		}
	}
}

## update the list of variable and parameter names, and get starting values based on the results of the stepwise regression
variables <- variable_names(stepModel)
starting_values <- setup_starting_values(variables, stepModel)


## view the predictions of the best model
## produce a plot of the results of the stepwise regression; uncomment to examine the intermediate results
# old code, will help when I split data into calibration/validation sets
# maple_fut = as.data.frame(sapply(colnames(mapleAll)[predictorColumns], function(nm) transformations[[nm]]$forward(mapleAll[,paste("fut_",nm,sep="")])))
# maple_pres = as.data.frame(sapply(colnames(mapleAll)[predictorColumns], function(nm) transformations[[nm]]$forward(mapleAll[,nm])))
# 
# predictions <- data.frame(
# 	present=predict(stepModel, newdata=maple_pres, type="response"), 
# 	future=predict(stepModel, newdata=maple_fut, type="response"))
# 
# # pdf("img/SDM_stepwiseTest.pdf", width=9,height=8.5)
# quartz(width=6,height=9)
# layout(matrix(c(1,2,3,3), byrow=F, nrow=2), widths=c(1, 0.15))
# par(mar=c(1,1,1,1))
# titles <- c('sdm present', paste('sdm HadA2 2080, ', paste(colnames(mapleAll)[predictorColumns], collapse=','), sep=""))
# quilt.plot(mapleAll[,1], mapleAll[,2], predictions$present, col=bpy.colors(), zlim=c(0,1), main=titles[1], add.legend=F, xaxt='n', yaxt='n')
# quilt.plot(mapleAll[,1], mapleAll[,2], predictions$future, col=bpy.colors(), zlim=c(0,1), main=titles[2], add.legend=F, xaxt='n', yaxt='n')
# plot.new()
# quilt.plot(mapleAll[,1], mapleAll[,2], predictions[,1], legend.only=T, col=bpy.colors(), zlim=c(0,1), add=F, smallplot=c(0.2,0.4, 0.05,0.95), bigplot=c(0,0,0,0))	
# dev.off()




save(stepModel, variables, starting_values, file="results/stepResults.rdata")