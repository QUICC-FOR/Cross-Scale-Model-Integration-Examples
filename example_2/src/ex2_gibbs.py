#  ex2_gibbs.py
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
# 	Primary implementation of gibbs sampling for the second example


## things to do:
# change to numpy.array instead of a list
# set up priors
# write data generation
# read in data and model spec
# set up initial vals
# some adaptation figure for an acceptance prob of about 0.2
# some printing for debugging
# testing

import numpy
import scipy.stats

NUM_PARAMS = 5
MARKOV_ITERATIONS = 100

# create a 2-D list to hold the posterior samples
# first dimension is SAMPLES, second is PARAMETERS
posteriorSamples = list()
currentVals = initial_values()
posteriorSamples.append(currentVals)

# all tuning values initially set to one; this is appropriate for a gaussian sampler
tuning = numpy.ones(NUM_PARAMS)

while(len(posteriorSamples) < MARKOV_ITERATIONS):
    X = generate_data(phenofitPredictions)
    newVals = list()
    # update the parameters one at a time
    for k in range(NUM_PARAMS):
        proposedVal = propose_parameter(currentVals[k], tuning[k])
        currentVals[k] = evaluate_parameter(proposedVal, k, X, predictors, currentVals)
    posteriorSamples.append(currentVals)


def model_linear_predictor(x, theta)
    '''this function must change with model specification, though it is quite general
    currently assumes that every predictor in x is simply multiplied by its parameter
    to specify powers, just add a column in x corresponding to whatever power is desired'''
    assert len(theta) == len(x) + 1
    psi = theta[0]
    for xx,th in zip(x, theta[1:])
        psi += xx*th
    return psi
    

def log_posterior_prob(X, parameters, predictors, k):
    """ calculates the full conditional log likelihood for the kth parameter
    WILL NEED TO CHANGE WITH MODEL SPECIFICATION (but not with number of params or prior vals)
    X: the data, a vector
    parameters: the vector of parameters for the model
    predictors: predictors for the linear model
    k: the index in parameters of the current parameter being evaluated
    """    
    # calculate conditional likelihood
    sumlogl = 0.0
    for x in X:
        p = inv_logit(model_linear_predictor(predictors, parameters))
        sumlogl += x * numpy.log(p^x) + (1-x) * numpy.log(1-p)

    # incorporate the prior for the kth parameter (all others are constant)
    sumlogl += numpy.log(scipy.stats.norm(prior[k][0], prior[k][1]).pdf(parameters[k]))
    return sumlogl


def inv_logit(x):
    # some logic here to avoid overflows
    if x > 0:
        psi = 1.0 / (1.0 + numpy.exp(-x))
    elif x <= 0:
        psi = numpy.exp(x) / (1.0 + numpy.exp(x))
    else:
        raise ValueError
    return(psi)


def initial_values():
    pass
    

def generate_data():
    pass


def propose_parameter(prevVal, A):
    """prevVal is the previous value
    A is the tuning parameter; here it is the standard deviation of a normal distribution"""
    return numpy.random.normal(prevVal, A)
    


def evaluate_parameters(proposedVal, k, X, predictors, currentVals):
    proposedParams = currentVals
    proposedParams[k] = proposedVal
    acceptanceProb = numpy.exp(log_posterior_prob(X, proposedParams, predictors) - log_posterior_prob(X, currentVals, predictors))
    testVal = numpy.random.uniform(0,1,1)
    if testVal < acceptanceProb:
        result = proposedVal
    else:
        result = currentVals[k]
    return result
    
