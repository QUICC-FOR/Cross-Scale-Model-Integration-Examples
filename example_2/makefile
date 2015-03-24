

## step 1: prepare the data
dat/rawData.rds: r/1_prepare_naive_data.r dat/raw/currentClim.rdata \
dat/raw/futureClim.rdata dat/raw/AceSac_pres.rds dat/raw/AceSac.csv
	Rscript r/1_prepare_naive_data.r

dat/mcmc/naiveData.csv: dat/rawData.rds
dat/mcmc/naiveData.rds: dat/rawData.rds
dat/mcmc/naivePriors.csv: dat/rawData.rds
dat/mcmc/naiveInits.csv: dat/rawData.rds
dat/predictionData.rds: dat/rawData.rds
dat/parameterScaling.rds: dat/rawData.rds

## step 2: run naive model
results/mcmc/naivePosterior.csv: bin/integrated_model dat/mcmc/naiveData.csv \
dat/mcmc/naiveInits.csv dat/mcmc/naivePriors.csv
	./bin/integrated_model -n -d dat/mcmc/naiveData.csv -i dat/mcmc/naiveInits.csv -p \
	dat/mcmc/naivePriors.csv -r 25000 -t 50 -b 20000 >results/mcmc/naivePosterior.csv

## step 3: prep integrated model
dat/mcmc/integratedPriors.csv: r/4_prepare_integrated_data.r results/naivePosterior.rds dat/rawData.rds
	Rscript r/4_prepare_integrated_data.r

dat/mcmc/integratedInits.csv: dat/mcmc/integratedPriors.csv
dat/mcmc/integratedPresData.csv: dat/mcmc/integratedPriors.csv
dat/mcmc/integratedFutureData.csv: dat/mcmc/integratedPriors.csv


## step 4: run integrated models on present phenofit data
results/mcmc/integratedPresent.csv: bin/integrated_model dat/mcmc/integratedPresData.csv \
dat/mcmc/integratedInits.csv dat/mcmc/integratedPriors.csv
	./bin/integrated_model -d dat/mcmc/integratedPresData.csv -i dat/mcmc/integratedInits.csv -p \
	dat/mcmc/integratedPriors.csv -r 25000 -t 50 -b 20000 >results/mcmc/integratedPresent.csv

results/mcmc/integratedFuture.csv: bin/integrated_model dat/mcmc/integratedFutureData.csv \
dat/mcmc/integratedInits.csv dat/mcmc/integratedPriors.csv
	./bin/integrated_model -d dat/mcmc/integratedFutureData.csv -i dat/mcmc/integratedInits.csv -p \
	dat/mcmc/integratedPriors.csv -r 25000 -t 50 -b 20000 >results/mcmc/integratedFuture.csv

## step 5: process results
results/posteriors.rdata: dat/predictionData.rds results/mcmc/naivePosterior.csv \
results/mcmc/integratedPresent.csv results/mcmc/integratedFuture.csv
	Rscript r/5_process_posterior.r


## make figures

figs: ex2_params.pdf ex2_response.pdf ex2_pres_map.pdf ex2_fut_map.pdf

ex2_params.pdf: r/fig_params.r results/posteriors.rdat
	Rscript r/fig_params.r

ex2_response.pdf: r/fig_response_curves.r results/posteriors.rdata \
dat/parameterScaling.rds dat/rawData.rds
	Rscript r/fig_response_curves.r

ex2_pres_map.pdf: r/fig_maps.r dat/figure/ne_50m_ocean dat/figure/ne_50m_lakes \
dat/figure/acersacr dat/rawData.rds results/naiveStats.rds results/intPresStats.rds \
results/intFutStats.rds
	Rscript r/fig_maps.r

ex2_fut_map.pdf: ex2_pres_map.pdf



## Validation
auc.txt: r/validation.r dat/rawData.rds results/naivePosterior.rds results/integratedPresentPosterior.rds
	Rscript r/validation.r


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





