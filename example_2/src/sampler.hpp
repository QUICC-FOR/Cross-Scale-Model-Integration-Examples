#ifndef SAMPLER_H
#define SAMPLER_H

/*
Model integration example 2: sampler.cpp
	
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
	
	


	Primary implementation of the Metropoolis-within-Gibbs sampler

	The sampler is implemented as a single class (named Sampler) with two public
	member functions.
	
	Constructor:
	Sampler::Sampler(std::vector<std::vector<double> > priors, 
			std::vector<double> response,
			std::vector<std::vector<double> > predictors, 
			std::vector<double> initialValues=std::vector<double>(), 
			std::vector<double> tuningParameters=std::vector<double>(), 
			int verbose=0, 
			bool autoAdapt=true)
		priors: a 2-D vector
		  the first dimension is of length k and contains one entry per parameter
		  priors[i][0] is the mean, priors[i][1] is the standard deviation for parameter i
		response: the vector of the response variable; this is a double distributed from
		  0 to 1; it is the probability of occurrence in simulated datasets used as the
		  pseudo-response in the likelihood function
		predictors: a 2-D vector of predictor variables
		  the first dimension contains one vector per observation/site
		  the second dimension contains individual variables
		  predictors[i][j] gives variable j at site i
		initialValues: an optional vector of initial values for the parameters
		  if omitted, the model will attempt to choose values using the prior
		  however, it is not very smart about it, and often the selected values will
		  be quite bad. Specifying inits is highly recommended
		tuningParameters: an optional vector of tuning parameters for the sampler
		  in normal sampling, given a value for parameter i at time t = x[i][t], the 
		  sampler will propose a new value for x[i][t+1] using:
		  		x[i][t] + rnorm(mean = 0, sd = tuningParameters[i])
		  If automatic adaptation is turned on (highly recommended), the sampler will 
		  adjust these values to provide an acceptance rate of ca. 0.3 per parameter.
		  If approximate values are known (perhaps as a result of preliminary runs),
		  they can be specified to speed up the adaptation process. Otherwise, the model
		  will default to a value of 1 for every parameter.
		verbose: determine how much output to send to the console via stderr
		  0: only basic status messages will be logged
		  1: also report tuning values during adaptation
		  2: also report the value of parameters at each step of sampling
		autoAdapt: boolean; if true, perform adaptation, otherwise go directly to sampling
		
	Sampler::run(const size_t n)
		Runs the sampler for n iterations. All samples are automatically sent to stdout
		on completion. Note that the Sampler automagically maintains its own internal
		state, so that:
			mySampler.run(n1);
			mySampler.run(n2);
		is functionally equivalent to
			mySampler.run(n1+n2);
		This can be used to provide simple breakpoints in the code, where the writing
		over long runs is broken into smaller pieces. Assuming stdout is redirected to a
		file, this provides a simple guard against crashes or other disruptions.

*/


// This turns on some additional debugging and logging to the console. Normally it is not
// needed
// #define SAMPLER_DEBUG

#include <vector>
#include <cstdlib>
#include <string>

class Sampler {
  public:
  	void run(const size_t n);
	Sampler(std::vector<std::vector<double> > priors, 
			std::vector<std::string> priorDistros, 
			std::vector<double> response,
			std::vector<int> weights,
			std::vector<std::vector<double> > predictors, 
			std::vector<double> initialValues=std::vector<double>(), 
			std::vector<double> tuningParameters=std::vector<double>(), 
			bool simResponse = false, int verbose=0, bool autoAdapt=true);

  private:
  	// methods
  	void auto_adapt();
  	void set_initial_values(std::vector<double> vals);
  	std::vector<double> do_sample(std::vector<std::vector<double> > &dest, size_t n); // returns the acceptance rates
  	void add_samples(const std::vector<std::vector<double> > &samples);
	std::vector<int> make_simulated_response() const;
	double propose_parameter(const size_t i) const;
	int choose_parameter(const long double proposal, const size_t i, 
			const std::vector<int> &Y, const std::vector<int> &N);
	long double log_posterior_prob(const std::vector<int> &Y, const std::vector<int> &N, 
			const std::vector<double> &params, const size_t i) const;
	long double log_prior(const long double & par, int index) const;
	long double model_linear_predictor(const std::vector<double> &x, 
			const std::vector<double> &params) const;
	void output();
  
  	// data
	std::vector<std::vector<double> > priors;
	std::vector<std::string> priorDist;
	std::vector<double> response;
	std::vector<int> weight;	// size parameter of the binomial distribution
	std::vector<std::vector<double> > predictors;
	std::vector<std::vector<double> > posteriorSamples;
	
	// state variables
	std::vector<double> currentState;
	std::vector<double> tuningParameters;
	size_t samplesTaken;
	size_t nParams;
	bool adapted;	
	
	// settings
	int verbose;
	bool retainPreAdaptationSamples;
	size_t autoAdaptIncrement;
	double targetAcceptanceRateInterval [2];
	float adaptationRate;
	size_t maxAdaptation;
	bool flushOnWrite;
	size_t outputIncrement;
	bool allowFittedExtremes;
	bool simulateResponse;
};

#endif