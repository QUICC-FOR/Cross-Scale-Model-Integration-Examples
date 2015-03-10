#!/usr/bin/Rscript

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


library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-t", "--thin", default=50, type="integer", help="thinning interval")
parser$add_argument("-b", "--burnin", default=500000, type="integer", help="burn-in period")
parser$add_argument("-i", "--infile", help="input file name - csv file from MCMC")
parser$add_argument("-d", "--datafile", default = "dat/mapleDat_processed.rds", help="processed data in rds format")
parser$add_argument("-o", "--outfile", help="output file name (without extension)")

argList = parser$parse_args()

if(is.null(argList$infile)) stop("Input file must be specified with -i")
if(is.null(argList$outfile)) argList$outfile = argList$infile


do_thin = function(x, n) {
	ind = seq(1, nrow(x),  n)
	return(x[ind,])
}


prep_posterior = function(samples, variables, burnin, thinInterval) {
	require(coda, quietly = TRUE)
	colnames(samples) = variables
	startVal = burnin + 1
	endVal = nrow(samples)
	samples = samples[(burnin+1):nrow(samples),]
	samples = do_thin(samples, thinInterval)
	posterior = mcmc(samples, start = startVal, end = endVal, thin = thinInterval)
	return(posterior)
}

dat = readRDS(argList$datafile)
model = read.csv(argList$infile, header=FALSE, stringsAsFactors=FALSE)
modelPosterior = prep_posterior(model, c("intercept", dat$variables), argList$burnin, argList$thin)
saveRDS(modelPosterior, paste(argList$outfile, ".rds", sep=""))
write.table(modelPosterior, paste(argList$outfile, ".csv", sep=""), sep=',', row.names=FALSE, col.names=FALSE)
