setwd("/Users/mtalluto/Documents/git_projects/Cross-Scale-Model-Integration-Examples/example_2/src")
# TEST_LEN = 200
# 
# inv_logit = function(x) exp(x) / (1 + exp(x))
# 
# testDat = data.frame(x1 = rnorm(TEST_LEN), x2 = rpois(TEST_LEN, 4))
# testDat$x3 = testDat$x1^2
# # starting with only 0 or 1 probabilities for now so we get 100% consistent datasets
# testDat$probs = rbinom(TEST_LEN, 1, prob = inv_logit(0.5 + 2.3 * testDat$x1 - 0.050*testDat$x3 + 0.795*testDat$x2))
# 
# # make a model to show us the "right" answers
# summary(glm(probs ~ x1+x2+x3, family=binomial, data=testDat))
# 
# priors = data.frame(mean=c(0,0,0.8,-0.05), sd = c(100000, 100000, 0.2, 0.03))
# 
# write.csv(testDat, file='testData.csv', row.names = FALSE)
# write.csv(priors, file='testPriors.csv', row.names = FALSE)
# 
# 
# thin = function(x, n) {
# 	ind = seq(1, nrow(x),  n)
# 	return(x[ind,])
# }
# 
# 
# ## check out the mcmc results
# library(coda)
# mcmcDat = read.csv('mcmc_output.csv', header=FALSE, stringsAsFactors=FALSE)
# for(i in 1:ncol(mcmcDat)) mcmcDat[,i] = as.numeric(mcmcDat[,i])
# colnames(mcmcDat) = c('b0', 'b1', 'b2', 'b3')
# mcmcDat = thin(mcmcDat, 10)
# mcmcDat = mcmc(mcmcDat, start=33000-5000, end=38000, thin=10)
# plot(mcmcDat)
# summary(mcmcDat)

load('../dat/maple.rdata')
load('../results/naiveModel.rdata')
# subset
mapleSubset = maple[sample(1:nrow(maple), 1000),]

# make a dummy response dataset for starting values
dumY = rbinom(nrow(mapleSubset), 1, mapleSubset$Phenofit_HadA2)

# model for starting values
coefs = glm(dumY ~ ddeg + I(ddeg^2) + I(ddeg^3) + sum_prcp + I(sum_prcp^2) + I(sum_prcp^3) + pToPET + I(pToPET^2), data = mapleSubset, family=binomial)$coefficients
write.csv(coefs, file='inits.csv', row.names = FALSE)

priors = data.frame(mean=rep(0,length(coefs)), sd=rep(100000, length(coefs)))
write.csv(priors, file='testPriors.csv', row.names = FALSE)

testDat = data.frame(ddeg = mapleSubset$ddeg, ddeg2 = mapleSubset$ddeg^2, ddeg3 = mapleSubset$ddeg^3,
	sum_prcp = mapleSubset$sum_prcp, sum_prcp2 = mapleSubset$sum_prcp^2, sum_prcp3 = mapleSubset$sum_prcp^3,
	pToPET = mapleSubset$pToPET, pToPET2 = mapleSubset$pToPET^2, phenofit = mapleSubset$Phenofit_HadA2)
write.csv(testDat, file='testData.csv', row.names = FALSE)
