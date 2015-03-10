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
#define WEIGHT_COL 9

// number of MCMC replicates
static int MCMC_REPS = 1000000;

// input file names
static const char * PRIOR_FILE = "dat/integratedPriors.csv";
static const char * DATA_FILE = "dat/integratedData.csv";
static const char * INIT_FILE = "dat/integratedInits.csv";
static bool SIM_RESPONSE = true;

#include "sampler.hpp"
#include "csv.hpp"
#include <vector>
#include <exception>
#include <iostream>
#include <unistd.h> // for getopt
#include <cstdlib>
#include <string>

using std::vector;
using std::string;

void parse_args(int argc, char **argv);
void print_help();

int main(int argc, char **argv) 
{
	parse_args(argc, argv);
	
	CSV rawData;
	vector<vector<double> > rawPriors;
	vector<vector<double> > priors;
	vector<string> priorDists;
	vector<double> inits;
	
	try {
		rawPriors = CSV(PRIOR_FILE, 1).data();
		rawData = CSV(DATA_FILE, 1);
		inits = CSV(INIT_FILE, 1).columns(0);
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
		abort();
	}
	
	// for the priors, we interpret a 1 in the distribution as cauchy, a 0 as gaussian
	for(int i = 0; i < rawPriors.size(); i++)
	{
		priors.push_back(vector<double> {rawPriors[i][0], rawPriors[i][1]});
		if(rawPriors[i][2] == 0)
			priorDists.push_back("gaussian");
		else if(rawPriors[i][2] == 1)
			priorDists.push_back("cauchy");
		else
			priorDists.push_back("unknown");
	}
	
	vector<vector<double> > predictors = rawData.columns(0,RESP_COL);
	vector<double> response = rawData.columns(RESP_COL);
	vector<int> weights;
	for(const auto & w : rawData.columns(WEIGHT_COL))
		weights.push_back(int(w));


	// tuning parameters; using predefined values will allow for faster tuning
	// no a priori tuning values
	// vector<double> tuning (9,0);

	// values that are pretty close (developed from previous testing) for the FUTURE
	// vector<double> tuning {0.150, 0.210, 0.151, 0.049, 0.552, 0.610, 0.141, 0.213, 0.0310};

	// values for the present
	// vector<double> tuning {0.150, 0.174, 0.0938, 0.0368, 0.377, 0.738, 0.155, 0.0994, 0.0192};

	// unweighted naive model
	vector<double> tuning {0.266, 2.967, 2.924, 1.482, 2.1, 7, 3.782, 3.98, 10.643};

	Sampler sampler = Sampler(priors, priorDists, response, weights, predictors, inits, 
			tuning, SIM_RESPONSE, VERBOSE_LEVEL);
	sampler.run(MCMC_REPS);

	return 0;
}

void parse_args(int argc, char **argv)
{
	int thearg;
	while((thearg = getopt(argc, argv, "hnd:p:i:r:")) != -1)
	{
		switch(thearg)
		{
			case 'h':
				print_help();
				break;
			case 'n':
				SIM_RESPONSE = false;
			case 'd':			
				DATA_FILE = optarg;
				break;
			case 'p':
				PRIOR_FILE = optarg;
				break;
			case 'i':
				INIT_FILE = optarg;
				break;
			case 'r':
				MCMC_REPS = atoi(optarg);
				break;
			case '?':
				print_help();
		}
	}
}

void print_help()
{
	std::cerr << "Command line options:\n";
	std::cerr << "    -h:             display this help\n";
	std::cerr << "    -n:             do NOT simulate the response (simulation is the default -- see guide)\n";
	std::cerr << "    -d <filname>:   use the file <filename> for the main data file\n";
	std::cerr << "    -i <filname>:   use the file <filename> for the initial parameter values\n";
	std::cerr << "    -p <filname>:   use the file <filename> for the prior distribution data\n";
	std::cerr << "    -r <integer>:   specify the number of mcmc replicates\n";
	exit(1);
}
