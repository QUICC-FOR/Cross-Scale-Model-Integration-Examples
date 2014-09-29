/*
Model integration example 2: main.cpp
	
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
	
	


	Main program block for the second example. See the SI Appendix 2 text and the included
	makefile for building instructions


*/


// some global model settings follow; they can be changed to alter overall model behavior
// verbosity; a value of 0 suppresses all but basic status messages
// 1 prints the status of adaptation/tuning
// 2 (not recommended except for testing) prints parameter values at every step
#define VERBOSE_LEVEL 1

// defines the column of the dataset that has the response variable
#define RESP_COL 8

// number of MCMC replicates
#define MCMC_REPS 50000000

// input file names
static const char * PRIOR_FILE = "dat/integratedPriors.csv";
static const char * DATA_FILE = "dat/integratedData.csv";
static const char * INIT_FILE = "dat/integratedInits.csv";

#include "sampler.hpp"
#include "csv.hpp"
#include <vector>
#include <exception>
#include <iostream>

using std::vector;

int main(void) 
{
	CSV rawData;
	vector<vector<double> > priors;
	vector<double> inits;
	
	try {
		priors = CSV(PRIOR_FILE, 1).data();
		rawData = CSV(DATA_FILE, 1);
		inits = CSV(INIT_FILE, 1).columns(0);
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
		abort();
	}
	
	vector<vector<double> > predictors = rawData.columns(0,RESP_COL);
	vector<double> phenofit = rawData.columns(RESP_COL);

	// uncomment the first to use no a priori tuning values
	// the second are values that are pretty close (developed from previous testing;
	// using them will allow for faster tuning
//	vector<double> tuning (9,0);
	vector<double> tuning {0.150, 0.210, 0.151, 0.049, 0.552, 0.610, 0.141, 0.213, 0.0310};

	Sampler sampler = Sampler(priors, phenofit, predictors, inits, tuning, VERBOSE_LEVEL);
	sampler.run(MCMC_REPS);

	return 0;
}
