""" integrated_model.py: script for invoking the integrated model and taking posterior
	samples. For details on implementation, invoke:
		import mcmc
		help(mcmc)
		help(mcmc.Sampler)

    License:
    Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
  
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.
    
    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

"""

import numpy as np
import mcmc
import pickle

# read in data from csv files
priors = np.genfromtxt('dat/integratedPriors.csv', delimiter=',', skip_header=1)
rawData = np.genfromtxt('dat/integratedModelData.csv', delimiter=',', skip_header=1)
inits = np.genfromtxt('dat/integradtedModelInits.csv', delimiter=',', skip_header=1)
predictors = rawData[:,0:8]
phenofit = rawData[:,8]

sampler = mcmc.Sampler(priors, mcmc.make_response_function(phenofit), predictors, inits, 
  outputFile = "results/integratedModel.csv", verbose=True)
sampler.verbose = False
sampler.run_sampler(1000000)

# save the sampler object for future use; because the current state is preserved,
# after reloading it can be resumed with sampler.run_sampler() to take additional samples
pFile = 'sampler.pkl'
with open(pFile, 'wb') as f:
	pickle.dump(sampler, f, pickle.HIGHEST_PROTOCOL)
	
# to restore:
# with open(pFile, 'rb') as f:
# 	sampler = pickle.load(f)
	