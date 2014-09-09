#include "csv.hpp"
#include <vector>
#include <exception>
#include <iostream>
#include <cmath>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>

using std::vector;
using std::cout;
using std::cerr;

#define RESP_COL 8
#define CI_LO 0.05
#define CI_HI 0.95

long double inv_logit(long double val)
{
	if(val > 0)
		return 1.0 / (1.0 + exp(-val));
	else
		return exp(val) / (1.0 + exp(val));		
}
	
int main(void) 
{
	// input data
	vector<vector<double> > mcmcDat;
	vector<vector<double> > predictors;
	try {
		mcmcDat = CSV("../results/integratedModelThinned.csv", 0).data();
		predictors = CSV("../dat/integratedData.csv", 1).columns(0,RESP_COL);
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
		abort();
	}
	
	const int nmcmc = mcmcDat.size();
	const int ndata = predictors.size();
	const int nparams = mcmcDat[0].size();
	
	// do some old school arrays here, to make sure we don't overflow later
//	double * means = new double[ndata];
//	double * sds = new double[ndata];
	double * dpPredictions = new double [nmcmc];
	 
	cout << "mean,SE,lower,upper\n";
	
	for(int datapoint = 0; datapoint < ndata; datapoint++) {
		const vector<double> & currentData = predictors.at(datapoint);
		for(int mcmcrep = 0; mcmcrep < nmcmc; mcmcrep++) {
			const vector<double> & currentRep = mcmcDat.at(mcmcrep);
			dpPredictions[mcmcrep] = currentRep.at(0);
			for(int par = 1; par < nparams; par++) {
				dpPredictions[mcmcrep] += currentRep.at(par) * currentData.at(par-1);
			}
			dpPredictions[mcmcrep] = (double) inv_logit(dpPredictions[mcmcrep]);
		}
//		means[datapoint] = gsl_stats_mean(dpPredictions, sizeof(double), nmcmc);
//		sds[datapoint] = gsl_stats_sd_m(dpPredictions, sizeof(double), nmcmc, means[datapoint];
		// storage is great and all, but not strictly required, since I can just print the values straightaway;
		double mean = gsl_stats_mean(dpPredictions, sizeof(double), nmcmc);
		double sd = gsl_stats_sd_m(dpPredictions, sizeof(double), nmcmc, mean);
		gsl_sort(dpPredictions, sizeof(double), nmcmc);
		double lower = gsl_stats_quantile_from_sorted_data(dpPredictions,sizeof(double), nmcmc, CI_LO);
		double upper = gsl_stats_quantile_from_sorted_data(dpPredictions,sizeof(double), nmcmc, CI_HI);
		cout << mean << ',' << sd << ',' << lower << ',' << upper << '\n';
		if(mean > 1) {
			cerr << "Problem on line " << datapoint << "\n";
			for(int i = 0; i < nmcmc; i++)
				cerr << dpPredictions[i] << ',';
 			cerr << "  " << mean << ',' << sd << ',' << lower << ',' << upper << '\n';

			break;
		}
	}


//	delete [] means;
//	delete [] sds;
//	delete [] dpPredictions;
	
	return 0;
}

