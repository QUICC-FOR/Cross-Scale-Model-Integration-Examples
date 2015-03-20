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
	See the file sampler.hpp for the public interface to this implementation
	
	This program takes advantage of the openMP library to speed execution on multi-core
	systems. See the makefile and associated documentation for instructions on building
	for multi-threaded usage. The global constant S_NUM_THREADS defines how many threads
	to use; for ideal performance, this should be 1 or 2 less than the number of physical
	CPU cores on your machine; the value used here (6) was tested for optimum performance
	on an 8-core Mac Pro.

*/


#include <iostream>
#include <ctime>
#include <cassert>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "sampler.hpp"

using std::vector;
using std::cout;
using std::cerr;

// total number of threads to be used if the program is built for multi-threaded use
#define S_NUM_THREADS 8


namespace {
	gsl_rng * set_up_rng() {
		gsl_rng * newRNG = gsl_rng_alloc(gsl_rng_mt19937);
		assert(newRNG);
		gsl_rng_set(newRNG, (int) time(NULL));
		return newRNG;
	}

	gsl_rng * rng = set_up_rng();
	
	std::string vec_to_str(const vector<double> &data)
	{
		std::stringstream out;
		out << std::setprecision(3);
		for(vector<double>::const_iterator i = data.begin(); i != data.end() - 1; ++i)
			out << *i << ", ";
		out << data.back();
		return out.str();
	}
	
	long double inv_logit(long double val)
	{
		if(val > 0)
			return 1.0 / (1.0 + exp(-val));
		else
			return exp(val) / (1.0 + exp(val));		
	}
}


void Sampler::run(const size_t n)
{
	size_t burninCompleted = 0; 
	size_t samplesCompleted = 0;
	while(samplesCompleted < n) {
		size_t samplesToTake;
		if(burninCompleted < burnin)
			samplesToTake = ( (burnin - burninCompleted < outputIncrement) ? (burnin - burninCompleted) : outputIncrement);
		else
			samplesToTake = ( (n - samplesCompleted < outputIncrement) ? (n - samplesCompleted) : outputIncrement);
		vector<vector<double> > newSamples;
		newSamples.reserve(samplesToTake);
		do_sample(newSamples, samplesToTake);

		time_t rawtime;
		time(&rawtime);
		struct tm * timeinfo = localtime(&rawtime);
		char fmtTime [20];
		strftime(fmtTime, 20, "%F %T", timeinfo);

		if(burninCompleted < burnin)
		{
			burninCompleted += samplesToTake;
			cerr << fmtTime << "   MCMC Iteration " << samplesTaken << "; burnin sample " << burninCompleted << " of " << burnin << '\n';
		}
		else
		{
			add_samples(newSamples);
			samplesCompleted += samplesToTake;
			cerr << fmtTime << "   MCMC Iteration " << samplesTaken << "; current job completed " << samplesCompleted << " of " << n << '\n';
			output();
		}		
	}
}

void Sampler::output()
{
	for(vector<vector<double> >::const_iterator row = posteriorSamples.begin(); row != posteriorSamples.end(); row++) {
		for(vector<double>::const_iterator val = (*row).begin(); val != (*row).end() - 1; val++)
			cout << *val << ',';
		// once more to get the end
		cout << (*row).back() << '\n';
	}
	cout.flush();
	if(flushOnWrite)
		posteriorSamples.clear();
}


Sampler::Sampler(vector<vector<double> > priors, std::vector<std::string> priorDistros, 
vector<double> response,
 vector<int> weights, vector<vector<double> > predictors, vector<double> initialValues, 
vector<double> tuningParameters, size_t thin, size_t burn, bool simResponse, int verbose, 
bool autoAdapt) :

// initializers for data/settings via the parameter list
priors(priors), priorDist(priorDistros), response(response), weight(weights), 
predictors(predictors), tuning(tuningParameters), verbose(verbose), 
simulateResponse(simResponse), thinning(thin), burnin(burn),

// magic numbers here are default values that have no support for initialization via parameters
retainPreAdaptationSamples(false), autoAdaptIncrement(5000), targetAcceptanceRateInterval {0.27, 0.34},
adaptationRate(1.1), maxAdaptation(100000), flushOnWrite(true), outputIncrement(50000),
allowFittedExtremes(false)
{

	// initialize the model state variables; these are only defined AFTER settings are initialized
	adapted = false;
	nParams = priors.size();
	size_t samplesTaken = 0;
	set_initial_values(initialValues);	// this will update the model state as well
	
	// fill in tuning parameters to make sure we have as many tuning values as parameters
	while(tuning.size() < nParams)
		tuning.push_back(1);
		
	if(autoAdapt)
		auto_adapt();	
}


void Sampler::auto_adapt()
{
	size_t adaptationSamplesTaken = 0;
	std::cerr << "Starting automatic adaptation...\n";
	while( !adapted && !(adaptationSamplesTaken >= maxAdaptation) ) {
		vector<vector<double> > newSamples;
		vector<double> acceptanceRates = do_sample(newSamples, autoAdaptIncrement);
		adaptationSamplesTaken += autoAdaptIncrement;
		adapted = true;
		if(verbose)
			cerr << "Adapting: acceptance rates: " << vec_to_str(acceptanceRates) << "; with tuning params: " << vec_to_str(tuning) << '\n';
		for(size_t k = 0; k < nParams; k++) {
			if(acceptanceRates[k] < targetAcceptanceRateInterval[0]/2.0) {
				tuning[k] /= 2.0*adaptationRate;
				adapted = false;
			}
			else if(acceptanceRates[k] < targetAcceptanceRateInterval[0]){
				tuning[k] /= adaptationRate;
				adapted = false;
			}
			else if(acceptanceRates[k] > targetAcceptanceRateInterval[1]*2.0){
				tuning[k] *= 2.0*adaptationRate;
				adapted = false;
			}
			else if(acceptanceRates[k] > targetAcceptanceRateInterval[1]){
				tuning[k] *= adaptationRate;
				adapted = false;
			}
		}
		if(retainPreAdaptationSamples)
			add_samples(newSamples);
	}
// 	if(!retainPreAdaptationSamples)
// 		posteriorSamples.push_back(currentState);	// save only the new "starting value" if we are not saving the adaptation samples

	if(adapted)
		cerr << "Adaptation completed successfully\n";
	else
		throw std::runtime_error("Automatic adaptation failed; try setting tuning parameters manually or using different inits");
}


vector<int> Sampler::make_simulated_response() const
{
	vector<int> result;
	result.reserve(response.size());
	if(simulateResponse)
	{
		for(size_t i = 0; i < response.size(); i++) 
			result.push_back(gsl_ran_binomial(rng, response[i], weight[i]));
	}
	else
	{
		for(size_t i = 0; i < response.size(); i++) 
			result.push_back(int(response[i]));
	}
	
	return(result);
}


double Sampler::propose_parameter(const size_t i) const
{
	return currentState[i] + gsl_ran_gaussian(rng, tuning[i]);
}


int Sampler::choose_parameter(const long double proposal, const size_t i, 
		const vector<int> &Y, const vector<int> &N)
{
	// returns 1 if proposal is accepted, 0 otherwise
	vector<double> proposedParameters = currentState;
	proposedParameters[i] = proposal;
	long double proposalLL = log_posterior_prob(Y, N, proposedParameters, i);
	long double currentLL = log_posterior_prob(Y, N, currentState, i);
	long double acceptanceProb = exp(proposalLL - currentLL);

	/*
		reject any parameter combinations that cause overflow or underflow
		this is perhaps excessively paranoid; it is possible that overflows could be
		handled more elgantly. but this cruder solution will help protect against runaway
		parameters
		
		temporarily disabled because we found some (hopefully) more appropriate solutions
	*/
// 	if(not std::isfinite(acceptanceProb))
// 	{
// 		std::cerr << "    proposalLL: " << proposalLL << "\n";
// 		std::cerr << "    currentLL: " << currentLL << "\n";
// 		std::cerr << "    acceptance probability: " << acceptanceProb << "\n";
// 		acceptanceProb = 0;
// 	}

	// 	check for nan -- right now this is not being handled, but it should be
	if(std::isnan(acceptanceProb))
		acceptanceProb = 0;
//		throw std::runtime_error("NaN detected in likelihood");
				
	double testVal = gsl_rng_uniform(rng);
	if(testVal < acceptanceProb) {
		currentState[i] = proposal;
		return 1;
	}
	else
	{
		return 0;
	}
}


long double Sampler::model_linear_predictor(const vector<double> &x, const vector<double> &params) const
{
	#ifdef SAMPLER_DEBUG
		assert(x.size() == params.size() - 1);
	#endif
	
	long double psi = params[0];
	for(size_t i = 1; i < params.size(); i++)
		psi += x[i-1] * params[i];
	return psi;
}



long double Sampler::log_posterior_prob(const vector<int> &Y, const vector<int> &N, 
		const vector<double> &params, const size_t index) const
{
	long double sumlogl = 0;

	const int nDataPoints = Y.size();
	#pragma omp parallel num_threads(S_NUM_THREADS)
	{
	#pragma omp for reduction(+:sumlogl)
	for(int i = 0; i < nDataPoints; i++) {
		long double p = inv_logit(model_linear_predictor(predictors[i], params));

		/* 
			if allowFittedExtremes is true, we detect overflow/underflow (due to inv logit
			of very large numbers or divide by very large number) and use nextafter to
			nudge those values back to the nearest representative value
		*/
		
		if(allowFittedExtremes && (p == 0.0 || p == 1.0))
		{
			std::cerr << "Warning: over/underflow detected in logit function is being ignored due to allowFittedExtremes = true\n";
			p = nextafter(p, abs(1.0 - p));
		}
			
//		long double logl = Y[i] * log(p) + (1-Y[i])*log(1-p); // binomial density
		long double logl;
		if(p == 0.0 || p == 1.0)	// penalize all fitted zeroes or ones
			logl = 0;
		else
			logl = std::log(gsl_ran_binomial_pdf(Y[i], p, N[i]));
		
		sumlogl += logl;
	}
	}
	
	// we only need to compute the prior prob of the parameter being evaluated
	// this is because the acceptance prob is a ratio of the probabilities, so all
	// other parameters, which are constant, will cancel
	sumlogl += log_prior(params[index], index);	
	return sumlogl;
}


long double Sampler::log_prior(const long double & par, int index) const
{
	// currently only gaussian and cauchy priors are supported
	// for cauchy, the SD is interpreted as the scale parameter
	long double value;
	if(priorDist[index] == "cauchy")
	{
		value = gsl_ran_cauchy_pdf(par - priors[index][0], priors[index][1]);
	}
	else if(priorDist[index] == "gaussian" or priorDist[index] == "normal")
	{
		value = gsl_ran_gaussian_pdf(par - priors[index][0], priors[index][1]);
	}
	else
	{
		std::string err = "Invalid distribution: " + priorDist[index] + " for parameter " + 
				std::to_string(index);
		throw std::runtime_error(err);
	}
	return std::log(value);
}


vector<double> Sampler::do_sample(vector<vector<double> > &dest, size_t n)
{
	dest.reserve(dest.size() + n);
	vector<size_t> nAccepted = vector<size_t>(nParams, 0);
	vector<double> acceptanceRates = vector<double>(nParams, 0);

	// we will shuffle the order each time we sample; this sets up an array of indices to do so
	size_t indices [nParams];
	for(size_t i = 0; i < nParams; i++) indices[i] = i;
	
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < thinning; j++) {
			vector<int> Y = make_simulated_response();	// this does nothing if simulateResponse is False
		
			gsl_ran_shuffle(rng, indices, nParams, sizeof(size_t));
			for(size_t j = 0; j < nParams; j++) {
				size_t k = indices[j];
				long double proposedVal = propose_parameter(k);
				nAccepted[k] += choose_parameter(proposedVal, k, Y, weight);
			}
		}
		dest.push_back(currentState);
		
		#ifdef SAMPLER_DEBUG
			if(verbose > 1)
				cerr << vec_to_str(currentState) << '\n';
		#endif
	}
	
	samplesTaken += n;
	for(size_t i = 0; i < nParams; i++)
		acceptanceRates[i] = double(nAccepted[i]) / n;
	
	return acceptanceRates;
} 


void Sampler::add_samples(const vector<vector<double> > &samples)
{
	posteriorSamples.insert(posteriorSamples.end(), samples.begin(), samples.end());
}


void Sampler::set_initial_values(vector<double> vals)
{
	if(vals.size() == 0 || vals.size() != priors.size()) {
		cerr << "Warning: bad or missing initial values; setting default values\n";
		vals = vector<double> (0, nParams);
		for(int i = 0; i < nParams; i++) {
			vals[i] = priors[i][0] + gsl_ran_gaussian(rng, priors[i][1]);
		}
	}
	currentState = vals;
}

