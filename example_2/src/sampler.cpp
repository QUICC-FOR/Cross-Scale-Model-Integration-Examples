#include <iostream>
#include <ctime>
#include <cassert>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include "sampler.hpp"

using std::vector;
using std::cout;
using std::cerr;


// unnamed namespace for file-global objects
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
	size_t numCompleted = 0;
	while(numCompleted < n) {
		size_t samplesToTake = ( (n - numCompleted < outputIncrement) ? (n - numCompleted) : outputIncrement);
		vector<vector<double> > newSamples;
		newSamples.reserve(samplesToTake);
		do_sample(newSamples, samplesToTake);
		add_samples(newSamples);
		numCompleted += samplesToTake;
		cerr << "MCMC Iteration " << samplesTaken << "; current job completed " << numCompleted << " of " << n << '\n';
		output();
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


Sampler::Sampler(vector<vector<double> > priors, vector<double> response,
vector<vector<double> > predictors, vector<double> initialValues, vector<double> tuningParameters,
int verbose, bool autoAdapt) :

// initializers for data/settings via the parameter list
priors(priors), response(response), predictors(predictors), tuningParameters(tuningParameters), verbose(verbose),

// magic numbers here are default values that have no support for initialization via parameters
retainPreAdaptationSamples(true), autoAdaptIncrement(1000), targetAcceptanceRateInterval {0.27, 0.34},
adaptationRate(1.1), maxAdaptation(50000), flushOnWrite(true), outputIncrement(10000),
preventFittedZeroesOnes(true)
{

	// initialize the model state variables; these are only defined AFTER settings are initialized
	adapted = false;
	nParams = priors.size();
	size_t samplesTaken = 0;
	set_initial_values(initialValues);	// this will update the model state as well
	if(tuningParameters.size() == 0)
		this->tuningParameters = vector<double> (nParams, 1);	// initialize all tuning parameters to 1 by default

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
			cerr << "Adapting: acceptance rates: " << vec_to_str(acceptanceRates) << "; with tuning params: " << vec_to_str(tuningParameters) << '\n';
		for(size_t k = 0; k < nParams; k++) {
			if(acceptanceRates[k] < targetAcceptanceRateInterval[0]/2.0) {
				tuningParameters[k] /= 2.0*adaptationRate;
				adapted = false;
			}
			else if(acceptanceRates[k] < targetAcceptanceRateInterval[0]){
				tuningParameters[k] /= adaptationRate;
				adapted = false;
			}
			else if(acceptanceRates[k] > targetAcceptanceRateInterval[1]*2.0){
				tuningParameters[k] *= adaptationRate;
				adapted = false;
			}
			else if(acceptanceRates[k] > targetAcceptanceRateInterval[1]){
				tuningParameters[k] *= adaptationRate;
				adapted = false;
			}
		}
		if(retainPreAdaptationSamples)
			add_samples(newSamples);
	}
	if(!retainPreAdaptationSamples)
		posteriorSamples.push_back(currentState);	// safe only the new "starting value" if we are not saving the adaptation samples

	if(adapted)
		cerr << "Adaptation completed successfully\n";
	else
		throw std::runtime_error("Automatic adaptation failed; try setting tuning parameters manually or using different inits");
}


vector<int> Sampler::make_simulated_response() const
{
	vector<int> result;
	result.reserve(response.size());
	for(size_t i = 0; i < response.size(); i++) 
		result.push_back(gsl_ran_bernoulli(rng, response[i]));

	return(result);
}


double Sampler::propose_parameter(const size_t i) const
{
	return currentState[i] + gsl_ran_gaussian(rng, tuningParameters[i]);
}


int Sampler::choose_parameter(const long double proposal, const size_t i, const vector<int> &Y)
{
	// returns 1 if proposal is accepted, 0 otherwise
	vector<double> proposedParameters = currentState;
	proposedParameters[i] = proposal;
// std::cerr << "getting acceptance prob...\n";
// std::cerr << "Current likelihood: " << log_posterior_prob(Y, currentState, i) << "\n";
	long double acceptanceProb = exp( log_posterior_prob(Y, proposedParameters, i) - log_posterior_prob(Y, currentState, i));
	double testVal = gsl_rng_uniform(rng);
// std::cerr << " test value = " << testVal << "\n";
	if(testVal < acceptanceProb) {
		currentState[i] = proposal;
		return 1;
	}
	else
		return 0;
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



long double Sampler::log_posterior_prob(const vector<int> &Y, const vector<double> &params, const size_t i) const
{
	long double sumlogl = 0;

	const int nDataPoints = Y.size();
	#pragma omp parallel num_threads(3)
	{
	#pragma omp for reduction(+:sumlogl)
	for(int i = 0; i < nDataPoints; i++) {
		long double p = inv_logit(model_linear_predictor(predictors[i], params));

		/* in this case, it is ok if the model fits a 0 or 1; we still want to count those
		    parameters. A zero or one is an overflow, not an error (the value isn't actually
		    zero or one, it just exceeds the machine's floating point precision. So to avoid
		    these values resuling in inf or nan (thus poisoning the whole likelihood), we
		    set p to the closest representable number greater than 0 (for a fitted zero) or
		    less than one (for a fitted one) */
		if(preventFittedZeroesOnes && (p == 0.0 || p == 1.0))
			p = nextafter(p, abs(1.0 - p));
			
		long double logl = Y[i] * log(p) + (1-Y[i])*log(1-p); // binomial density
		sumlogl += logl;
	}
	}
	
	sumlogl += log(gsl_ran_gaussian_pdf(params[i] - priors[i][0], priors[i][1]));	
	return sumlogl;
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
		vector<int> Y = make_simulated_response();
		gsl_ran_shuffle(rng, indices, nParams, sizeof(size_t));
		for(size_t j = 0; j < nParams; j++) {
			size_t k = indices[j];
			long double proposedVal = propose_parameter(k);
			nAccepted[k] += choose_parameter(proposedVal, k, Y);
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

