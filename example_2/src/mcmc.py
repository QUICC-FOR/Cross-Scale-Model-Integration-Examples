#  mcmc.py
#  
#    Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
#  
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#    
#    This program is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#  
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#  
# 	Primary implementation of Metropolis within gibbs sampling for the second example


## things to do:
# read in data and model spec
# some adaptation; figure for an acceptance prob of about 0.2 to 0.4
# some printing for debugging
# testing

import numpy as np
import scipy.stats


def do_mcmc(priors, phenofitPredictions, predictors, numIterations = 100, inits = None):
    nParams = priors.shape[0]
    # create a 2-D list to hold the posterior samples
    # first dimension is SAMPLES, second is PARAMETERS
    posteriorSamples = np.empty(shape=(numIterations, nParams))
    posteriorSamples.fill(np.NAN)
    currentVals = initial_values(inits, priors)
    posteriorSamples[0,:] = currentVals   # COPY the values to the posterior samples

    # all tuning values initially set to one; this is appropriate for a gaussian sampler
    tuning = np.ones(nParams)

    for i in range(1,numIterations):    # starting at one because of initial values
        X = generate_data(phenofitPredictions)
        # update the parameters one at a time
        # randomize the order in which parameters are updated
        indices = range(nParams)
        np.random.shuffle(indices)
        for k in indices:
            proposedVal = propose_parameter(currentVals[k], tuning[k])
            currentVals[k] = evaluate_parameter(proposedVal, k, X, predictors, currentVals, priors)
        posteriorSamples[i,:] = currentVals
        print currentVals
    return(posteriorSamples)


def model_linear_predictor(x, theta):
    '''this function must change with model specification, though it is quite general
    currently assumes that every predictor in x is simply multiplied by its parameter
    to specify powers, just add a column in x corresponding to whatever power is desired'''
    assert len(theta) == len(x) + 1
    psi = theta[0]
    for xx,th in zip(x, theta[1:]):
        psi += xx*th
    return psi
    

def log_posterior_prob(X, parameters, predictors, prior, k):
    """ calculates the full conditional log likelihood for the kth parameter
    WILL NEED TO CHANGE WITH MODEL SPECIFICATION (but not with number of params or prior vals)
    X: the data, a vector
    parameters: the vector of parameters for the model
    predictors: predictors for the linear model
    k: the index in parameters of the current parameter being evaluated
    """    
    # calculate conditional likelihood
    sumlogl = 0.0
    for i in range(len(X)):
        y = predictors[i,:]
        x = X[i]
        p = inv_logit(model_linear_predictor(y, parameters))
        sumlogl += x * np.log(np.power(p,x)) + (1-x) * np.log(1-p)

    # incorporate the prior for the kth parameter (all others are constant)
    # assumes a normal prior for all parameters
    sumlogl += np.log(scipy.stats.norm(prior[k,0], prior[k,1]).pdf(parameters[k]))
    return sumlogl


def initial_values(inits, prior):
    """Define initial values for all parameters
    WILL NEED TO CHANGE WITH MODEL SPECIFICATION (assumes all priors are normal)"""
    if inits is None or inits.shape[0] != prior.shape[0]:
        inits = np.zeros(prior.shape[0])
        for i in range(prior.shape[0]):
            inits[i] = np.random.normal(prior[i,0], prior[i,1],1)
    return inits


def inv_logit(x):
    # some logic here to avoid overflows
    if x > 0:
        psi = 1.0 / (1.0 + np.exp(-x))
    elif x <= 0:
        psi = np.exp(x) / (1.0 + np.exp(x))
    else:
        raise ValueError
    return(psi)


    

def generate_data(inputProbabilities):
    outputData = np.zeros(shape=inputProbabilities.shape)
    for i in range(inputProbabilities.shape[0]):
        outputData[i] = np.random.binomial(1,inputProbabilities[i])
    return(outputData)
        


def propose_parameter(prevVal, A):
    """prevVal is the previous value
    A is the tuning parameter; here it is the standard deviation of a normal distribution"""
    return np.random.normal(prevVal, A)
    


def evaluate_parameter(proposedVal, k, X, predictors, currentVals, prior):
    proposedParams = currentVals
    proposedParams[k] = proposedVal
    acceptanceProb = np.exp(log_posterior_prob(X, proposedParams, predictors, prior, k) - log_posterior_prob(X, currentVals, predictors, prior, k))
    testVal = np.random.uniform(0,1,1)
    if testVal < acceptanceProb:
        result = proposedVal
    else:
        result = currentVals[k]
    return result
    

