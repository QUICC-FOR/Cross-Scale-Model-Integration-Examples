library(coda)
library(pROC)
validation = readRDS("dat/rawData.rds")$valid
validation = validation[complete.cases(validation),]
naive = readRDS("results/naivePosterior.rds")
integrated = readRDS("results/integratedPresentPosterior.rds")

bn = colMeans(naive)
bi = colMeans(integrated)

# transform validation dataset into something usable with pROC
validation = within(validation,
{
	naivePred = plogis(bn[1] + bn[2]*ddeg + bn[3]*ddeg^2 + bn[4]*sum_prcp + bn[5]*sum_prcp^2 + bn[6]*pToPET + bn[7]*pToPET^2)
	intPred = plogis(bi[1] + bi[2]*ddeg + bi[3]*ddeg^2 + bi[4]*sum_prcp + bi[5]*sum_prcp^2 + bi[6]*pToPET + bi[7]*pToPET^2)
	nZeros = count - pres
})

zeroInd = rep(1:nrow(validation), times=validation$nZeros)
oneInd =  rep(1:nrow(validation), times=validation$pres)
pres = c(rep(0, length(zeroInd)), rep(1, length(oneInd)))
rocDat = cbind(pres, validation[c(zeroInd, oneInd), c('intPred', 'naivePred')])

naiveROC = roc(response = rocDat$pres, predictor = rocDat$naivePred)
intROC = roc(response = rocDat$pres, predictor = rocDat$intPred)

cat(paste("Naive AUC: ", naiveROC$auc, "\nIntegrated AUC: ", intROC$auc, "\n", sep=""), file="auc.txt")
