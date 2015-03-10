library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-m", "--model", default="results/tmp/naive_glm.rds", help="fitted glm model to use as basis for the MCMC")
parser$add_argument("-d", "--datafile", default="dat/mapleDat_processed.rds", help="data file")
parser$add_argument("-p", "--prior", default='dat/mcmc/naivePriors.csv', help="output file for naive priors")
parser$add_argument("-n", "--inits", default='dat/mcmc/naiveInits.csv', help="output file for naive inits")
# parser$add_argument("-o", "--outfile", default="results/tmp/naive_glm.rds", help="RDS output file name")
argList = parser$parse_args()

mod = readRDS(argList$model)
dat = readRDS(argList$datafile)

naivePriors = with(mod, data.frame(
		mean = rep(0, length(coefficients)),
		sd = c(10, rep(2.5, length(coefficients)-1)),
		dist = rep(1, length(coefficients)),
		row.names = sort(names(coefficients))))
write.csv(naivePriors, file=argList$prior, row.names = FALSE)


# use the prior distribution to draw random starts for the parameters
naiveInits = sapply(1:nrow(naivePriors), function(i)
{
	with(naivePriors[i,],
	{
		func = ifelse(dist == 0, rnorm, rcauchy)
		func(1, mean, sd)
	})
})
write.csv(naiveInits, file=argList$inits, row.names = FALSE)
		

presPredictors = sort(names(mod$coefficients))[-1]
futPredictors = paste("fut_", presPredictors, sep="")

naiveData_weighted = naiveData_unweighted = intData_Pres = dat$calib[,presPredictors]
intData_Fut = dat$calib[,futPredictors]
naiveData_weighted$response = dat$calib$weightedPresence
naiveData_unweighted$response = dat$calib$PresObs
intData_Pres$response = dat$calib$Phenofit_CRU
intData_Fut$response = dat$calib$Phenofit_HadA2
naiveData_weighted$weights = dat$calib$weightedN
naiveData_unweighted$weights = intData_Pres$weights = intData_Fut$weights = rep(1, nrow(dat$calib))

write.csv(intData_Fut, file='dat/mcmc/integratedData_Fut.csv', row.names = FALSE)
write.csv(intData_Pres, file='dat/mcmc/integratedData_Pres.csv', row.names = FALSE)
write.csv(naiveData_weighted, file='dat/mcmc/naiveData_weighted.csv', row.names = FALSE)
write.csv(naiveData_unweighted, file='dat/mcmc/naiveData_unweighted.csv', row.names = FALSE)

predictionDat = dat$all[,c(6,12,8,13,19,15)]
predictionDat$ddeg = dat$transformations$ddeg$forward(predictionDat$ddeg)
predictionDat$sum_prcp = dat$transformations$sum_prcp$forward(predictionDat$sum_prcp)
predictionDat$pToPET = dat$transformations$pToPET$forward(predictionDat$pToPET)
predictionDat$fut_ddeg = dat$transformations$ddeg$forward(predictionDat$fut_ddeg)
predictionDat$fut_sum_prcp = dat$transformations$sum_prcp$forward(predictionDat$fut_sum_prcp)
predictionDat$fut_pToPET = dat$transformations$pToPET$forward(predictionDat$fut_pToPET)
write.table(predictionDat, "dat/mcmc/predictionData.csv", sep=",", col.names=FALSE, row.names=FALSE)

