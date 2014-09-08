#include "sampler.hpp"
#include "csv.hpp"
#include <vector>

// #include <sstream>
// #include <iomanip>

using std::vector;

#define RESP_COL 2

// std::string vec_to_str(const vector<double> &data)
// {
// 	std::stringstream out;
// 	out << std::setprecision(3);
// 	for(vector<double>::const_iterator i = data.begin(); i != data.end() - 1; ++i)
// 		out << *i << ", ";
// 	out << data.back();
// 	return out.str();
// }


int main(void) 
{
	// data structures for input data
	vector<vector<double> > priors = CSV("integratedPriors.csv", 1).data();
	CSV rawData = CSV("integratedData.csv", 1);
	vector<double> inits = CSV("integratedInits.csv", 1).columns(0);
	vector<vector<double> > predictors = rawData.columns(0,RESP_COL);
	vector<double> phenofit = rawData.columns(RESP_COL);
	vector<double> tuning = vector<double>(); // fill this in with vals from python
	
// std::cerr << "Priors: \n";
// for(int i = 0; i < priors.size(); i++)
// 	std::cerr << " [" << vec_to_str(priors[i]) << "]\n";
// std::cerr << "\nInits:\n [" << vec_to_str(inits) << "]\n";
// std::cerr << "\nResponse:\n [" << vec_to_str(phenofit) << "]\n";
// std::cerr << "Predictors: \n";
// for(int i = 0; i < predictors.size(); i++)
// 	std::cerr << " [" << vec_to_str(predictors[i]) << "]\n";
// 


int bug = 0;
	size_t verbose = 0;
	Sampler sampler = Sampler(priors, phenofit, predictors, inits, tuning, verbose);
	sampler.run(10000);


	return 0;
}

