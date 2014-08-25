import numpy as np
import mcmc

# read in data from csv files
priors = np.genfromtxt('testPriors.csv', delimiter=',', skip_header=1)
rawData = np.genfromtxt('testData.csv', delimiter=',', skip_header=1)
predictors = rawData[:,0]
predictors = np.reshape(predictors, (len(predictors), 1))
phenofit = rawData[:,1]
inits = np.array([1.,4.])
sampler = mcmc.Sampler(priors, mcmc.make_response_function(phenofit), predictors, inits, verbose=False, autoAdapt = True)
sampler._outputIncrement = 1000
sampler.run_sampler(20000)