import numpy as np
import mcmc

# read in data from csv files
priors = np.genfromtxt('testPriors.csv', delimiter=',', skip_header=1)
rawData = np.genfromtxt('testData.csv', delimiter=',', skip_header=1)
inits = np.genfromtxt('inits.csv', delimiter=',', skip_header=1)
predictors = rawData[:,0:8]
phenofit = rawData[:,8]

sampler = mcmc.Sampler(priors, mcmc.make_response_function(phenofit), predictors, inits, verbose=True, autoAdapt = True)
sampler.verbose = False
sampler._outputIncrement = 5000
sampler.run_sampler(50000)

