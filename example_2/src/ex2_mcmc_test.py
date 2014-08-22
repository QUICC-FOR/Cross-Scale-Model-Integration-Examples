import numpy as np
import mcmc

# read in data from csv files
priors = np.genfromtxt('testPriors.csv', delimiter=',', skip_header=1)
rawData = np.genfromtxt('testData.csv', delimiter=',', skip_header=1)
predictors = rawData[:,0]
predictors = np.reshape(predictors, (len(predictors), 1))
phenofit = rawData[:,1]
inits = np.array([1.,4.])

mcmc.do_mcmc(priors, phenofit, predictors, 1000, inits)