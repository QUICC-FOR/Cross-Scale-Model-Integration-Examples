#### individual scripts for example 2:
## order of operations to run the full analysis:
# 1: Rscript 1-setup_naive_model.r
# 2: ./bin/integrated_model >results/integratedModel.csv
# 3: Rscript 3-process_integrated.r
# 4: ./bin/stats >results/integratedStats.csv
# 5: Rscript 5-make_figure.r


SI: ELEtallutoSA2.pdf
model: bin/integrated_model bin/stats

### build the SI appendix pdf
ELEtallutoSA2.pdf: tex/ELEtallutoSA2.tex
	cd tex; pdflatex ELEtallutoSA2
	mv tex/ELEtallutoSA2.pdf ELEtallutoSA2.pdf




## step 1: prepare the data
dat/mapleDat_processed.rds: r/1_prep_data.r dat/raw/AceSac.csv dat/raw/currentClim.rdata \
dat/raw/futureClim.rdata
	Rscript r/1_prep_data.r -d dat/raw/AceSac.csv -o dat/mapleDat_processed.rds

dat/mcmc/naivePriors.csv: dat/mapleDat_processed.rds
dat/mcmc/naiveInits.csv: dat/mcmc/naivePriors.csv
dat/mcmc/integratedData_Fut.csv: dat/mcmc/naivePriors.csv
dat/mcmc/integratedData_Pres.csv: dat/mcmc/naivePriors.csv
dat/mcmc/naiveData_weighted.csv: dat/mcmc/naivePriors.csv
dat/mcmc/naiveData_unweighted.csv: dat/mcmc/naivePriors.csv
dat/mcmc/predictionData.csv: dat/mcmc/naivePriors.csv


## step 2: run naive model
results/mcmc/naivePosteriorRaw_unweighted.csv: bin/integrated_model \
dat/mcmc/naiveData_unweighted.csv dat/mcmc/naiveInits.csv dat/mcmc/naivePriors.csv
	bin/integrated_model -n -r 5000000 -d dat/mcmc/naiveData_unweighted.csv \
	-i dat/mcmc/naiveInits.csv -p dat/mcmc/naivePriors.csv \
	>results/mcmc/naivePosteriorRaw_unweighted.csv

results/mcmc/naivePosteriorRaw_weighted.csv: bin/integrated_model \
dat/mcmc/naiveData_weighted.csv dat/mcmc/naiveInits.csv dat/mcmc/naivePriors.csv
	bin/integrated_model -n -r 5000000 -d dat/mcmc/naiveData_weighted.csv \
	-i dat/mcmc/naiveInits.csv -p dat/mcmc/naivePriors.csv \
	>results/mcmc/naivePosteriorRaw_weighted.csv


## step 3: process the posterior distributions
results/naiveModel_uw.rds: r/3_process_posterior.r dat/mapleDat_processed.rds \
results/mcmc/naivePosteriorRaw_unweighted.csv 
	Rscript r/3_process_posterior.r -i results/mcmc/naivePosteriorRaw_unweighted.csv \
	-o results/naiveModel_uw

results/naiveModel_w.rds: r/3_process_posterior.r dat/mapleDat_processed.rds \
results/mcmc/naivePosteriorRaw_weighted.csv 
	Rscript r/3_process_posterior.r -i results/mcmc/naivePosteriorRaw_weighted.csv \
	-o results/naiveModel_w








##### pre-cleanup scripts below






#### step 6: run integrated model future



#### step 7: compute statistics (mean and sd) on PREDICTIONS of the 3 models
results/naiveStats.csv: bin/stats results/naiveModel.csv dat/predictionData.csv
	./bin/stats -d results/naiveModel.csv >results/naiveStats.csv


#### step 8: save the resulting stats as .rds files to facilitate quick loading later



#### step 9: make figures





## really old stuff that will definitely need updating 




#### step 2a: run the future integrated model
# requires the integrated_model binary, which must be built before running these steps
# see the end of this makefile for this build procedure
# note that this is very CPU intensive
results/integratedModel_Fut.csv: bin/integrated_model dat/integratedData_Fut.csv \
dat/integratedPriors.csv dat/integratedInits.csv
	mkdir -p results
	./bin/integrated_model -r 5000000 -d dat/integratedData_Fut.csv \
	>results/integratedModel_Fut.csv
	
#### step 2b: run the present integrated model
results/integratedModel_Pres.csv: bin/integrated_model dat/integratedData_Pres.csv \
dat/integratedPriors.csv dat/integratedInits.csv
	mkdir -p results
	./bin/integrated_model -r 5000000 -d dat/integratedData_Pres.csv \
	>results/integratedModel_Pres.csv

#### step 3: process the results (e.g., throw away burnin samples, thin results)
results/integratedModel_Pres.rdata: 3-process_integrated.r results/integratedModel_Fut.csv \
results/integratedModel_Pres.csv dat/naive_model.rdata
	Rscript 3-process_integrated.r
	
#### step 4: compute statistics (mean & sd) on integrated model predictions
stats: results/integratedStats_Pres.csv results/integratedStats_Fut.csv

results/integratedStats_Pres.csv: bin/stats results/integratedModel_Pres_Thinned.csv \
dat/predictionData.csv
	./bin/stats -d results/integratedModel_Pres_Thinned.csv >results/integratedStats_Pres.csv

results/integratedStats_Fut.csv: bin/stats results/integratedModel_Fut_Thinned.csv \
dat/predictionData.csv
	./bin/stats -d results/integratedModel_Fut_Thinned.csv >results/integratedStats_Fut.csv

#### step 5: make the figure
ex2.pdf: 5-make_figure.r dat/naive_model.rdata results/integratedStats_Pres.csv \
results/integratedStats_Fut.csv results/integratedModel_Pres.rdata \
results/integratedModel_Fut.rdata dat/ne_110m_ocean/ne_110m_ocean.shp
	Rscript 5-make_figure.r


## extra stuff just to make sure all dependencies are resolved	    
## these are created by the scripts above, so all are no-operation
#dat/integratedPriors.csv: dat/naive_model.rdata

#dat/integratedInits.csv: dat/naive_model.rdata


#results/integratedModel_Fut.rdata: results/integratedModel_Pres.rdata

#results/integratedModel_Pres_Thinned.csv: results/integratedModel_Pres.rdata

#results/integratedModel_Fut_Thinned.csv: results/integratedModel_Fut.rdata
	

################
#
#  C++ build procedure for the integrated model
#
################

CC=/opt/local/bin/c++-mp-4.9
CF=-std=c++11

# for compiling with openMP, use the first
# otherwise, use the second
CO=$(CF) -fopenmp
#CO=$(CF)

bin/integrated_model: bin/main.o bin/sampler.o bin/csv.o
	$(CC) $(CO) -lgsl -o bin/integrated_model bin/main.o bin/sampler.o bin/csv.o

bin/main.o: src/main.cpp src/sampler.hpp src/csv.hpp
	mkdir -p bin
	$(CC) $(CF) -c -o bin/main.o src/main.cpp

bin/sampler.o: src/sampler.cpp src/sampler.hpp
	mkdir -p bin
	$(CC) $(CO) -c -o bin/sampler.o src/sampler.cpp

bin/csv.o: src/csv.cpp src/csv.hpp
	mkdir -p bin
	$(CC) $(CF) -c -o bin/csv.o src/csv.cpp

bin/stats: bin/stats.o bin/csv.o
	$(CC) $(CF) -lgsl -o bin/stats bin/stats.o bin/csv.o

bin/stats.o: src/statistics.cpp src/csv.hpp
	mkdir -p bin
	$(CC) $(CF) -c -o bin/stats.o src/statistics.cpp




