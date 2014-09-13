#include "sampler.hpp"
#include "csv.hpp"
#include <vector>
#include <exception>
#include <iostream>

using std::vector;

// defines the column of the dataset that has the response variable
#define RESP_COL 8

int main(void) 
{
	// data structures for input data
	CSV rawData;
	vector<vector<double> > priors;
	vector<double> inits;
	
	try {
		priors = CSV("dat/integratedPriors.csv", 1).data();
		rawData = CSV("dat/integratedData.csv", 1);
		inits = CSV("dat/integratedInits.csv", 1).columns(0);
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
		abort();
	}
	
	vector<vector<double> > predictors = rawData.columns(0,RESP_COL);
	vector<double> phenofit = rawData.columns(RESP_COL);
// 	vector<double> tuning = vector<double>((predictors.at(0)).size() + 1, 1.000);
//	vector<double> tuning {0.188, 0.342, 0.376, 0.233, 0.467, 1.000, 0.310, 0.282, 0.0388};
	vector<double> tuning {0.150, 0.210, 0.151, 0.049, 0.552, 0.610, 0.141, 0.213, 0.0310};

	size_t verbose = 1;
	Sampler sampler = Sampler(priors, phenofit, predictors, inits, tuning, verbose);
	sampler.run(5000000);

	return 0;
}

