""" mcmc.py: Implementation of Metropolis-Hastings sampler for the phenofit integration

Provides two convenience functions and one class

make_response_function(respData):
    This returns a function that can be passed to the sampler class as the response
    object. The function returned uses the input vector of probabilities to generate
    a random simulated dataset

Sampler(priors, response, predictors, initialValues = None, autoAdapt = True, 
      verbose = False, outputFile = 'mcmc_output.csv')
    Sampler class for performing all the mcmc work. See help(Sampler) for details.

inv_logit(x):
    Inverse of the logistic transformation


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
import scipy.stats


def make_response_function(respData):
    """Returns a closure for generating random response datasets.
    
    respData: a vector of probabilities
    Value: a function that, when invoked, generates a random simulated dataset of 1's and
      0's based on the probabilities in respData
    """
    def resp_generator(respData = respData):
        outputData = np.zeros(shape=respData.shape)
        for i in range(respData.shape[0]):
            outputData[i] = np.random.binomial(1,respData[i])
        return outputData
    return resp_generator


class Sampler(object):
    """Class for implementing the MCMC algorithm.
    
    This class makes extensive use of numpy features. In particular, numpy.ndarray
    is used throughout, and users should be familiar with with numpy arrays. Unless
    otherwise noted, when an array is expected, it is a numpy array, and passing other
    array-like structures may produce errors or undefined behavior.
    The sampler class provides for setting a number of options at initialization:
    priors: 2-D array of priors; each row is a parameter; the first column is the mean and
      the second is the standard deviation (NOT the precision)
    response: a callable object that returns a vector of presence-absence observations.
      The utility function mcmc.make_response_function is provided to generate such an
      object using a vector of probabilities
    predictors: a 2-D array of predictor variables, one per column; same number of rows
      as the object returned by response()
    initialValues: an optional 1-D array of values to be used to initialize the sampler.
      If omitted, the sampler will attempt to choose initial values from the prior
      distribution. With uninformative priors, this will likely result in very bad guesses
      and subsequent long convergence times.
    outputFileName: Default is 'mcmc_output.csv.' By default, the sampler will save all
      samples to this file every 100,000 iterations. Samples are not flushed from memory
      on write by default (see methods below to change this behavior); for very long
      runs, this can cause the memory footprint of the sampler to grow quite large and
      can result in long writes.
    autoAdapt: if True (the default), the sampler will automatically adapt itself to
      provide acceptance rates in the interval [0.27, 0.34]. Will raise a warning if this 
      takes more than 50,000 iterations. In this case, it is fine to try another round of 
      adaptation with Sampler.auto_adapt(). Alternatively, try different initial values.
    verbose: Set to True if you want to see each step in the Markov chain at the cost of
      greatly slowing down the algorithm. False is the default.
      
    Additional attributes available by default:
    posteriorSamples: 2-D array; each row is a single MCMC replicate, columns store
      each parameter.
      
    The following high-level methods are available in this class. For implementing
    custom samplers based on this class, it will be necessary to override some of them
    (or to subclass Sampler); pay particular attention to functions focused on calculating
    the model likelihood to functions called by those functions.

    run_sampler(numSamples):
      By default, this is the only method that the user needs to invoke after 
      initialization, assuming that adaptation completes successfully. This will run
      the sampler for numSamples iterations. Output is written to disk periodically,
      and at completion, all samples are written to disk.

    set_initial_values(inits = None):
      This function is run by default at initialization, and ordinarily does not need to
      be invoked by the user. However, in the case of bad or forgotten inits, it can be
      invoked manually. When inits = None (the default) or if the value passed in inits
      does not match the number of parameters, will set initial values by drawing from the
      prior distribution. Otherwise, Inits should be a 1-D array with one value per
      parameter in the model. Will need to be overridden to change from a gaussian prior.
    
    auto_adapt():
      By default, runs automatically at initialization. Will proceed until the sampler
      is adapted for all parameters or until it has tried 50,000 iterations. If adaptation
      fails, it can be invoked again to try an additional 50k iterations.
      
    set_output_behavior(file = None, increment = None, flush = None):
      Change how the sampler writes its output to disk. The default values (None) result
      in no changes being made. Set individual values to change behavior:
      file: a string; the name of the file to write to
      increment: an integer; after how many iterations should output be written
        (output is always written when run_sampler() terminates)
      flush: boolean; if True, when output is written to disk, it will be removed from
        memory, reducing the memory footprint for long runs. Subsequent writes will append
        to the file rather than overwriting it.
    
    log_posterior_prob(Y, parameters, k):
      The full conditional probability for parameter k given the other parameters and the
      response data. It is ordinarily never necessary for the user to invoke this method
      directly. Override this method to change the model specification.
    
    model_linear_predictor(x, theta):
      The returns the linear predictor of the model with parameters in theta and
      predictor variables in x. It is ordinarily never necessary for the user to invoke 
      this method directly. Override this method to change the model specification.
      
     
    """
    def __init__(self, priors, response, predictors, initialValues = None, outputFileName = 
    'mcmc_output.csv', autoAdapt = True, verbose = False):
        self.priors = priors
        self.response = response
        self.predictors = predictors
        self.verbose = verbose
        self.samplesTaken = 0
        self.outputFileName = outputFileName
        self.posteriorSamples = self._currentState = None
        self._nParams = self.priors.shape[0]
        self._tuning = np.ones(self._nParams)
        self._retainPreAdaptationSamples = False
        self._autoAdaptIncrement = 1000
        self._targetAcceptanceRate = (0.27, 0.34)
        self._adapted = False
        self._adaptationRate = 1.1
        self._maxAdaptation = 50000
        self._flushOnWrite = False
        self._firstWriteComplete = False
        self._outputIncrement = 100000
        self._guardFittedZerosOnes = True
        self.set_initial_values(initialValues)
        if autoAdapt:
            self.auto_adapt()
            
    def run_sampler(self, numSamples):
        """ Run the MCMC sampler for numSamples.
        
        By default, this is the only method that the user needs to invoke after 
        initialization, assuming that adaptation completes successfully. This will run
        the sampler for numSamples iterations. Output is written to disk periodically,
        and at completion, all samples are written to disk."""
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
        """Define initial values for all parameters.
        
        This function is run by default at initialization, and ordinarily does not need to
        be invoked by the user. However, in the case of bad or forgotten inits, it can be
        invoked manually. When inits = None (the default) or if the value passed in inits
        does not match the number of parameters, will set initial values by drawing from the
        prior distribution. Otherwise, Inits should be a 1-D array with one value per
        parameter in the model. Will need to be overridden to change from a gaussian prior."""
        if inits is None or inits.shape[0] != self.priors.shape[0]:
            print("Automatically generating initial values")
            inits = np.zeros(self.priors.shape[0])
            for i in range(self.priors.shape[0]):
                inits[i] = np.random.normal(self.priors[i,0], self.priors[i,1],1)
        self._currentState = inits.copy()
        self.posteriorSamples = np.empty((0, self._nParams))

    def auto_adapt(self):
        """High level interface for automatically adapting the sampler.
        
        By default, runs automatically at initialization. Will proceed until the sampler
        is adapted for all parameters or until it has tried 50,000 iterations. If adaptation
        fails, it can be invoked again to try an additional 50k iterations."""
        adaptAttempts = 0
        print("Starting automatic adaptation...")
        while not self._adapted and not adaptAttempts >= self._maxAdaptation:
            samples, acceptanceRates = self._do_sample(self._autoAdaptIncrement, adapt=True)
            adaptAttempts += self._autoAdaptIncrement
            self._adapted = True
            if self.verbose: print("Adapting: acceptance rates: " + str(acceptanceRates) +
              " with tuning parameters: " + str(self._tuning))
            for k in range(self._nParams):
                if acceptanceRates[k] < self._targetAcceptanceRate[0]/2.0:
                    self._tuning[k] /= 2.0*self._adaptationRate
                    self._adapted = False                    
                elif acceptanceRates[k] < self._targetAcceptanceRate[0]:
                    self._tuning[k] /= self._adaptationRate
                    self._adapted = False
                elif acceptanceRates[k] > self._targetAcceptanceRate[1]*2.0:
                    self._tuning[k] *= 2.0*self._adaptationRate
                    self._adapted = False
                elif acceptanceRates[k] > self._targetAcceptanceRate[1]:
                    self._tuning[k] *= self._adaptationRate
                    self._adapted = False
            if not self._retainPreAdaptationSamples:
                samples = np.reshape(self._currentState.copy(), (1,self._nParams))
            self._add_samples(samples)
        if self._adapted:
            print("Adaptation completed successfully")
        else:
            raise RuntimeWarning("Automatic adaptation failed after " + str(adaptAttempts) + 
              " iterations with tuning parameters " + str(self._tuning))
    
    def set_output_behavior(file = None, increment = None, flush = None):
        """Change how the sampler writes its output to disk.
        
        The default values (None) result in no changes being made. Set individual values 
        to change behavior:
        file: a string; the name of the file to write to
        increment: an integer; after how many iterations should output be written
          (output is always written when run_sampler() terminates)
        flush: boolean; if True, when output is written to disk, it will be removed from
          memory, reducing the memory footprint for long runs. Subsequent writes will append
          to the file rather than overwriting it."""
        if file is not None:
            self.outputFileName = file
        if increment is not None:
            self._outputIncrement = increment
        if flush is not None:
            self._flushOnWrite = flush
        
    def log_posterior_prob(self, Y, parameters, k):
        """ calculates the full conditional log likelihood for the kth parameter

        Override this method to change the model specification.
        X: the data, a vector
        parameters: the vector of parameters for the model
        predictors: predictors for the linear model
        k: the index in parameters of the current parameter being evaluated
        
        if _guardFittedZerosOnes is True (the default), propagation of inf or nan in the
        binomial density will be avoided when prob == 1 or prob == 0 by setting the 
        probability to the largest representable floating point value that is less than 1
        (in the case of p == 1) or the smallest representable value greater than 0
        (for p == 0)
        if _guardFittedZerosOnes is False, inf and nan resulting from fitted zeros and ones
        will be allowed to propagate, resulting in an inf or nan likelihood and a rejection
        of the parameter combination.
        """    
        # calculate conditional likelihood
        sumlogl = np.longdouble(0.0)
        for i in range(len(Y)):
            x = self.predictors[i,:]
            y = Y[i]
            p = inv_logit(self.model_linear_predictor(x, parameters))
            if self._guardFittedZerosOnes and (p == 0 or p == 1):
                p = np.nextafter(p, np.abs(1 - p))  # prevents nan or inf with 0 or 1 probs
            sumlogl += y * np.log(np.power(p,y)) + (1-y) * np.log(1-p) # binomial density
        # incorporate the prior for the kth parameter (all others are constant)
        # assumes a normal prior for all parameters
        sumlogl += np.log(scipy.stats.norm(self.priors[k,0], self.priors[k,1]).pdf(parameters[k]))
        return float(sumlogl)

    def model_linear_predictor(self, x, theta):
        """Returns the linear predictor for the model.
        
        x: a vector of predictor variables
        theta: a vector of parameters; should be one longer than x
        It is ordinarily never necessary for the user to invoke this method directly. 
        Override this method to change the model specification.
        Currently assumes that every predictor in x is simply multiplied by its parameter
        to specify powers, just add a column in x corresponding to whatever power is desired"""
        assert len(theta) == len(x) + 1
        psi = theta[0]
        for xx,th in zip(x, theta[1:]):
            psi += xx*th
        return psi

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
    """Return the inverse of the logit function"""
    # some logic here to avoid overflows
    x = np.longdouble(x)
    if x > 0:
        psi = 1.0 / (1.0 + np.exp(-x))
    elif x <= 0:
        psi = np.exp(x) / (1.0 + np.exp(x))
    else:
        raise ValueError
    return psi