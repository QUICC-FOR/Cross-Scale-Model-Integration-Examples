setwd('/Users/mtalluto/Documents/git_projects/Cross-Scale-Model-Integration-Examples/example_2/src')

datLen = 5000

ilogit <- function(x) exp(x)/(1+exp(x))
dat <- data.frame(x1=rnorm(datLen,0,1),x2=rnorm(datLen,0,1))
dat$y = ilogit(3.1 + 2.5 * dat$x1 - 2*dat$x2 + rnorm(datLen))
priors = data.frame(mu=c(2.3,3.5,-1), sd=c(1,1,1))
inits = data.frame(x=c(2.3,3.5,-1))

write.csv(dat, "integratedData.csv", row.names=FALSE)
write.csv(priors, "integratedPriors.csv", row.names=FALSE)
write.csv(inits, "integratedInits.csv", row.names=FALSE)

