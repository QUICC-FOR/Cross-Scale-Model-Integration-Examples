#ifndef SAMPLER_H
#define SAMPLER_H

#define SAMPLER_DEBUG

#include <vector>
#include <cstdlib>

class Sampler {
  public:
  	void run(const size_t n);
	Sampler(std::vector<std::vector<double> > priors, std::vector<double> response,
		std::vector<std::vector<double> > predictors, std::vector<double> initialValues=std::vector<double>(), 
		std::vector<double> tuningParameters=std::vector<double>(), int verbose=0, bool autoAdapt=true);

  private:
  	// methods
  	void auto_adapt();
  	void set_initial_values(std::vector<double> vals);
  	std::vector<double> do_sample(std::vector<std::vector<double> > &dest, size_t n); // returns the acceptance rates
  	void add_samples(const std::vector<std::vector<double> > &samples);
	std::vector<int> make_simulated_response() const;
	double propose_parameter(const size_t i) const;
	int choose_parameter(const long double proposal, const size_t i, const std::vector<int> &Y);
	long double log_posterior_prob(const std::vector<int> &Y, const std::vector<double> &params, const size_t i) const;
	long double model_linear_predictor(const std::vector<double> &x, const std::vector<double> &params) const;
	void output();
  
  	// data
	std::vector<std::vector<double> > priors;
	std::vector<double> response;
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
	bool preventFittedZeroesOnes;
};

#endif