/*
Model integration example 2: statistics.cpp
	
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
	
	


	This program (which is compiled as a separate binary from the main program) computes
	statistics on the posterior samples.
	In particular, it is designed to be given a set of predictors (not necessarily the
	ones used for estimation), on which it will compute the predicted probability of
	presence for each data point for each posterior sample.
	As written, it expects a set of present climate variables and a set of future
	climate variables as input
	Output is the mean probability of presence, as well as the standard error and 5% and 
	95% quantiles (i.e., the 90% credible interval) for each of the present and future 
	datasets.
	
	This functionality is hard-coded; see the code below to change the way predictions
	are calculated.
*/


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
	vector<double> ddeg, pToPET, sum_prcp, fut_ddeg, fut_pToPET, fut_sum_prcp;
	
	CSV rawData;
	try {
		mcmcDat = CSV("results/integratedModelThinned.csv", 0).data();
		rawData = CSV("dat/predictionData.csv", 0);
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
		abort();
	}
	
	ddeg = rawData.columns(0);
	pToPET = rawData.columns(1);
	sum_prcp = rawData.columns(2);
	fut_ddeg = rawData.columns(3);
	fut_pToPET = rawData.columns(4);
	fut_sum_prcp = rawData.columns(5);


	const int nmcmc = mcmcDat.size();
	const int ndata = ddeg.size();
	cerr << "Computing statistics for " << nmcmc << " MCMC reps and " << ndata << " data points\n";
	
	double * presPredictions = new double [nmcmc];
	double * futPredictions = new double [nmcmc];
	 
	cout << "pres_mean,pres_SE,pres_lower,pres_upper,fut_mean,fut_SE,fut_lower,fut_upper\n";
	
	for(int datapoint = 0; datapoint < ndata; datapoint++) {
		for(int mcmcrep = 0; mcmcrep < nmcmc; mcmcrep++) {
			const vector<double> & b = mcmcDat.at(mcmcrep);
			presPredictions[mcmcrep] = b[0] + b[1]*ddeg[datapoint] + b[2]*pow(ddeg[datapoint],2) + 
				b[3]*pow(ddeg[datapoint],3) + b[4]*pToPET[datapoint] + b[5]*pow(pToPET[datapoint],2) + 
				b[6]*sum_prcp[datapoint] + b[7]*pow(sum_prcp[datapoint],2) + b[8]*pow(sum_prcp[datapoint],3);
			futPredictions[mcmcrep] = b[0] + b[1]*fut_ddeg[datapoint] + b[2]*pow(fut_ddeg[datapoint],2) + 
				b[3]*pow(fut_ddeg[datapoint],3) + b[4]*fut_pToPET[datapoint] + b[5]*pow(fut_pToPET[datapoint],2) + 
				b[6]*fut_sum_prcp[datapoint] + b[7]*pow(fut_sum_prcp[datapoint],2) + b[8]*pow(fut_sum_prcp[datapoint],3);
			presPredictions[mcmcrep] = (double) inv_logit(presPredictions[mcmcrep]);
			futPredictions[mcmcrep] = (double) inv_logit(futPredictions[mcmcrep]);
		}
		double mean [2];
		mean[0] = gsl_stats_mean(presPredictions, 1, nmcmc);
		mean[1] = gsl_stats_mean(futPredictions, 1, nmcmc);
		double sd [2];
		sd[0] = gsl_stats_sd_m(presPredictions, 1, nmcmc, mean[0]);
		sd[1] = gsl_stats_sd_m(futPredictions, 1, nmcmc, mean[1]);
		gsl_sort(presPredictions, 1, nmcmc);
		gsl_sort(futPredictions, 1, nmcmc);
		double lower [2];
		lower[0] = gsl_stats_quantile_from_sorted_data(presPredictions,1, nmcmc, CI_LO);
		lower[1] = gsl_stats_quantile_from_sorted_data(futPredictions,1, nmcmc, CI_LO);
		double upper [2];
		upper[0] = gsl_stats_quantile_from_sorted_data(presPredictions,1, nmcmc, CI_HI);
		upper[0] = gsl_stats_quantile_from_sorted_data(futPredictions,1, nmcmc, CI_HI);
		cout << mean[0] << ',' << sd[0] << ',' << lower[0] << ',' << upper[0] << ',';
		cout << mean[1] << ',' << sd[1] << ',' << lower[1] << ',' << upper[1] << '\n';
	}
	
	return 0;
}

