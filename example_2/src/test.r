setwd('/Users/mtalluto/Documents/git_projects/Cross-Scale-Model-Integration-Examples/example_2/src')

datLen = 5000

ilogit <- function(x) exp(x)/(1+exp(x))
dat <- data.frame(x1=rnorm(datLen,0,1),x2=rnorm(datLen,0,1))
dat$y = ilogit(3.1 + 2.5 * dat$x1 - 2*dat$x2 + rnorm(datLen))
dat$y = rbinom(nrow(dat), 1, dat$y)
priors = data.frame(mu=c(2.3,3.5,-1), sd=c(1,1,1))
inits = data.frame(x=c(2.3,3.5,-1))

write.csv(dat, "integratedData.csv", row.names=FALSE)
write.csv(priors, "integratedPriors.csv", row.names=FALSE)
write.csv(inits, "integratedInits.csv", row.names=FALSE)




#####

library(coda)

thinLength = 50
burnin = 9000


thin = function(x, n) {
	ind = seq(1, nrow(x),  n)
	return(x[ind,])
}

integratedModel = read.csv("test.csv", header=FALSE, stringsAsFactors=FALSE)
for(i in 1:ncol(integratedModel)) integratedModel[,i] = as.numeric(integratedModel[,i])
colnames(integratedModel) = c('b1','b2','b3')
startVal = burnin + 1
endVal = nrow(integratedModel)
integratedModel = integratedModel[(burnin+1):nrow(integratedModel),]
integratedModel = thin(integratedModel, thinLength)
integratedModel = mcmc(integratedModel, start = startVal, end = endVal, thin = thinLength)

plot(integratedModel)
summary(integratedModel)
summary(glm(y ~ x1+x2, data=dat, family=binomial))