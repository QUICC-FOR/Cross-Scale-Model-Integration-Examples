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

library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-i", "--infile", default="dat/mapleDat_processed.rds", help="RDS input file name")
parser$add_argument("-o", "--outfile", default="results/tmp/naive_glm.rds", help="RDS output file name")
argList = parser$parse_args()


# fit the naive model - this is an initial step for selecting predictors. we will
# later re-fit the model using MCMC
# we use the weighted presences as the response
# we have selected 3 of 6 predictors due to collinearities in the predictors
dat = readRDS(argList$infile)$calib

# column numbers for response and predictor variables
response="cbind(weightedPresence, weightedN)"
predictors=c("ddeg", "sum_prcp", "pToPET")
#predictors=c(6,8,12)

generate_formula = function(x, response, predictors)
{
	allPredictors = colnames(x)[sapply(paste("^", predictors, sep=""), grep, colnames(x))]
	eval(paste(response, "~", paste(allPredictors, collapse=" + "), sep=" "))
}

stepScope = list(
	lower = cbind(weightedPresence, weightedN) ~ 1,
	upper = generate_formula(dat, response, predictors)
)

mod = step( glm(stepScope$lower, family=binomial, data=dat), scope=stepScope, 
		direction="both", k = log(nrow(dat)))
saveRDS(mod, file=argList$outfile)


