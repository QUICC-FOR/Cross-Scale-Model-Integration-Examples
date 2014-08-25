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
#     Primary implementation of Metropolis within gibbs sampling for the second example


## things to do:
# test mcmc
# documentation

import numpy as np
import scipy.stats


def make_response_function(respData):
    """returns a function that, each time called, returns a randomly generated
    response dataset drawn from a binomial distribution using the input probabilities"""
    def resp_generator(respData = respData):
        outputData = np.zeros(shape=respData.shape)
        for i in range(respData.shape[0]):
            outputData[i] = np.random.binomial(1,respData[i])
        return outputData
    return resp_generator


class Sampler(object):
    def __init__(self, priors, response, predictors, initialValues = None, autoAdapt = True, 
      verbose = False, outputFile = 'mcmc_output.csv'):
        self.verbose = verbose
        self.priors = priors
        self.response = response
        self.predictors = predictors
        self.samplesTaken = 0
        self.outputFileName = outputFile
        self.posteriorSamples = self._currentState = None
        self._nParams = self.priors.shape[0]
        self._tuning = np.ones(self._nParams)
        self._completedIterations = 0
        self._retainPreAdaptationSamples = False
        self._autoAdaptIncrement = 1000
        self._targetAcceptanceRate = (0.27, 0.33)
        self._adapted = False
        self._adaptationRate = 1.1
        self._maxAdaptation = 50000
        self._flushOnWrite = False
        self._firstWriteComplete = False
        self._outputIncrement = 100000
        self.set_initial_values(initialValues)
        if autoAdapt:
            self.auto_adapt()
            
    def run_sampler(self, numSamples):
        numCompleted = 0
        while numCompleted < numSamples:
            samplesToTake = (numSamples - numCompleted) if (numSamples - numCompleted < 
              self._outputIncrement) else self._outputIncrement
            newSamples = self._do_sample(samplesToTake, verbose = self.verbose)
            self._add_samples(newSamples)
            numCompleted += samplesToTake
            print("MCMC Iteration " + str(self.samplesTaken) + "; current job completed " + 
              str(numCompleted) + " of " + str(numSamples))
            self._write_samples()
    
    def set_initial_values(self, inits = None):
        """Define initial values for all parameters
        WILL NEED TO CHANGE WITH MODEL SPECIFICATION (assumes all priors are normal)"""
        if inits is None or inits.shape[0] != self.priors.shape[0]:
            print("Automatically generating initial values")
            inits = np.zeros(self.priors.shape[0])
            for i in range(self.priors.shape[0]):
                inits[i] = np.random.normal(self.priors[i,0], self.priors[i,1],1)
        self._currentState = inits.copy()
        self.posteriorSamples = np.empty((0, self._nParams))

    def auto_adapt(self):
        """High level interface for automatically adapting the sampler.
        Runs until all samplers are running within an internally-specified efficiency
        interval"""
        adaptAttempts = 0
        print("Starting automatic adaptation...")
        while not self._adapted and not adaptAttempts >= self._maxAdaptation:
            samples, acceptanceRates = self._do_sample(self._autoAdaptIncrement, adapt=True)
            adaptAttempts += self._autoAdaptIncrement
            self._adapted = True
            for k in range(self._nParams):
                if acceptanceRates[k] < self._targetAcceptanceRate[0]:
                    self._tuning[k] /= self._adaptationRate
                    self._adapted = False
                elif acceptanceRates[k] > self._targetAcceptanceRate[1]:
                    self._tuning[k] *= self._adaptationRate
                    self._adapted = False
            if self.verbose: print("Adapting: acceptance rates: " + str(acceptanceRates) +
              " with tuning parameters: " + str(self._tuning))
            if not self._retainPreAdaptationSamples:
                samples = np.reshape(self._currentState.copy(), (1,self._nParams))
            self._add_samples(samples)
        if self._adapted:
            print("Adaptation completed successfully")
        else:
            raise RuntimeError("Automatic adaptation failed after " + str(adaptAttempts) + " iterations")
        
    def _do_sample(self, numIterations, adapt=False, verbose=False):
        samples = np.empty(shape=(numIterations, self._nParams))
        samples.fill(np.NAN)    # technically unneeded; it is here for debugging only
        nAccepted = np.zeros((self._nParams), int)
        indices = range(self._nParams)
        for i in range(numIterations):
            Y = self.response() # designed to be used with a dynamically generated response
            # update the parameters one at a time
            # randomize the order in which parameters are updated
            np.random.shuffle(indices)
            for k in indices:
                proposedVal = self._propose_parameter(k)
                nAccepted[k] += self._choose_value(proposedVal, k, Y, samples)
            samples[i,:] = self._currentState
            if verbose: print self._currentState
        self.samplesTaken += numIterations
        if adapt:
            acceptanceRates = nAccepted / float(numIterations)
            return (samples, acceptanceRates)
        else:
            return samples

    def _propose_parameter(self, k):
        return np.random.normal(self._currentState[k], self._tuning[k])
    
    def _choose_value(self, proposedVal, k, Y, samples):
        proposedParams = self._currentState.copy()
        proposedParams[k] = proposedVal
        acceptanceProb = np.exp(self.log_posterior_prob(Y, proposedParams, k) - 
           self.log_posterior_prob(Y, self._currentState, k))
        testVal = np.random.uniform(0,1,1)
        if testVal < acceptanceProb:
            self._currentState[k] = proposedVal
            return True
        else:
            return False
    
    def _add_samples(self, samples):
        if self.posteriorSamples is None:
            self.posteriorSamples = samples
        else:
            self.posteriorSamples = np.concatenate((self.posteriorSamples, samples))
            
    def log_posterior_prob(self, Y, parameters, k):
        """ calculates the full conditional log likelihood for the kth parameter
        WILL NEED TO CHANGE WITH MODEL SPECIFICATION (but not with number of params or prior vals)
        X: the data, a vector
        parameters: the vector of parameters for the model
        predictors: predictors for the linear model
        k: the index in parameters of the current parameter being evaluated
        """    
        # calculate conditional likelihood
        sumlogl = 0.0
        for i in range(len(Y)):
            x = self.predictors[i,:]
            y = Y[i]
            p = inv_logit(self.model_linear_predictor(x, parameters))
            sumlogl += y * np.log(np.power(p,y)) + (1-y) * np.log(1-p) # binomial density
        # incorporate the prior for the kth parameter (all others are constant)
        # assumes a normal prior for all parameters
        sumlogl += np.log(scipy.stats.norm(self.priors[k,0], self.priors[k,1]).pdf(parameters[k]))
        return sumlogl

    def model_linear_predictor(self, x, theta):
        '''this function must change with model specification, though it is quite general
        currently assumes that every predictor in x is simply multiplied by its parameter
        to specify powers, just add a column in x corresponding to whatever power is desired'''
        assert len(theta) == len(x) + 1
        psi = theta[0]
        for xx,th in zip(x, theta[1:]):
            psi += xx*th
        return psi
        
    def _write_samples(self):
        if self._flushOnWrite and self._firstWriteComplete:
            mode = 'a'
        else:
            mode = 'w'
        with open(self.outputFileName, mode) as file:
            np.savetxt(file, self.posteriorSamples, delimiter=',')
        if self._flushOnWrite:
            self.posteriorSamples = None
        if not self._firstWriteComplete:
            self._firstWriteComplete = True
        


def inv_logit(x):
    # some logic here to avoid overflows
    x = np.longdouble(x)
    if x > 0:
        psi = 1.0 / (1.0 + np.exp(-x))
    elif x <= 0:
        psi = np.exp(x) / (1.0 + np.exp(x))
    else:
        raise ValueError
    return(float(psi))