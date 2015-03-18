

## step 1: prepare the data
dat/rawData.rds: r/1_prepare_naive_data.r dat/raw/currentClim.rdata \
dat/raw/futureClim.rdata dat/raw/AceSac_pres.rds dat/raw/AceSac.csv
	Rscript r/1_prepare_naive_data.r

dat/mcmc/naiveData.csv: dat/rawData.rds
dat/mcmc/naiveData.rds: dat/rawData.rds
dat/mcmc/naivePriors.csv: dat/rawData.rds
dat/mcmc/naiveInits.csv: dat/rawData.rds
dat/mcmc/predictionData.csv: dat/rawData.rds
dat/parameterScaling.rds: dat/rawData.rds

## step 2: run naive model
# disabled so it doesn't get re-run now that I have decent results
# results/mcmc/naive.csv: bin/integrated_model dat/mcmc/naiveData.csv \
#dat/mcmc/naiveInits.csv dat/mcmc/naivePriors.csv
#	./bin/integrated_model -n -d dat/mcmc/naiveData.csv -i dat/mcmc/naiveInits.csv -p \
#	dat/mcmc/naivePriors.csv -r 5000000 >results/mcmc/naive.csv

## step 3: process the posterior distribution
# results/naivePosterior.rds: r/3_process_posterior.r results/mcmc/naive.csv dat/mcmc/naiveData.rds
#	Rscript r/3_process_posterior.r -i results/mcmc/naive.csv -d dat/mcmc/naiveData.rds -o results/naivePosterior

# results/naivePosterior.csv: results/naivePosterior.rds

## step 4: prep integrated model
dat/mcmc/integratedPriors.csv: r/4_prepare_integrated_data.r results/naivePosterior.rds dat/rawData.rds
	Rscript r/4_prepare_integrated_data.r

dat/mcmc/integratedInits.csv: dat/mcmc/integratedPriors.csv
dat/mcmc/integratedPresData.csv: dat/mcmc/integratedPriors.csv
dat/mcmc/integratedFutureData.csv: dat/mcmc/integratedPriors.csv

## step 5: run integrated model on present phenofit data
results/mcmc/integratedPresent.csv: bin/integrated_model dat/mcmc/integratedPresData.csv \
dat/mcmc/integratedInits.csv dat/mcmc/integratedPriors.csv
	./bin/integrated_model -d dat/mcmc/integratedPresData.csv -i dat/mcmc/integratedInits.csv -p \
	dat/mcmc/integratedPriors.csv -r 5000000 >results/mcmc/integratedPresent.csv

## step 6: run integrated model on future phenofit data
results/mcmc/integratedFuture.csv: bin/integrated_model dat/mcmc/integratedFutureData.csv \
dat/mcmc/integratedInits.csv dat/mcmc/integratedPriors.csv
	./bin/integrated_model -d dat/mcmc/integratedFutureData.csv -i dat/mcmc/integratedInits.csv -p \
	dat/mcmc/integratedPriors.csv -r 5000000 >results/mcmc/integratedFuture.csv

## step 7: prepare integrated posterior distributions
results/integratedPresentPosterior.rds: r/3_process_posterior.r results/mcmc/integratedPresent.csv dat/mcmc/naiveData.rds
	Rscript r/3_process_posterior.r -i results/mcmc/integratedPresent.csv -d dat/mcmc/naiveData.rds -o results/integratedPresentPosterior

results/integratedPresentPosterior.csv: results/integratedPresentPosterior.rds

results/integratedFuturePosterior.rds: r/3_process_posterior.r results/mcmc/integratedFuture.csv dat/mcmc/naiveData.rds
	Rscript r/3_process_posterior.r -i results/mcmc/integratedFuture.csv -d dat/mcmc/naiveData.rds -o results/integratedFuturePosterior

results/integratedFuturePosterior.csv: results/integratedFuturePosterior.rds


## step 8: compute posterior statistics (mean, sd, cred. intervals) on the predictions of the 3 models and save as .rds files
stats: results/naiveStats.rds results/intPresStats.rds results/intFutStats.rds

results/naiveStats.csv: bin/stats results/naivePosterior.csv dat/mcmc/predictionData.csv
	./bin/stats -d results/naivePosterior.csv >results/naiveStats.csv

results/intPresStats.csv: bin/stats results/integratedPresentPosterior.csv dat/mcmc/predictionData.csv
	./bin/stats -d results/integratedPresentPosterior.csv >results/intPresStats.csv

results/intFutStats.csv: bin/stats results/integratedFuturePosterior.csv dat/mcmc/predictionData.csv
	./bin/stats -d results/integratedFuturePosterior.csv >results/intFutStats.csv

## make RDS files
results/naiveStats.rds: r/8_rds_stats.r results/naiveStats.csv
	Rscript r/8_rds_stats.r -d results/naiveStats.csv

results/intPresStats.rds: r/8_rds_stats.r results/intPresStats.csv
	Rscript r/8_rds_stats.r -d results/intPresStats.csv

results/intFutStats.rds: r/8_rds_stats.r results/intFutStats.csv
	Rscript r/8_rds_stats.r -d results/intFutStats.csv

## step 9+: make figures

figs: ex2_params.pdf ex2_response.pdf ex2_pres_map.pdf ex2_fut_map.pdf

ex2_params.pdf: r/fig_params.r results/integratedPresentPosterior.rds results/integratedFuturePosterior.rds results/naivePosterior.rds
	Rscript r/fig_params.r

ex2_response.pdf: r/fig_response_curves.r results/integratedPresentPosterior.rds results/integratedFuturePosterior.rds results/naivePosterior.rds dat/parameterScaling.rds dat/rawData.rds
	Rscript r/fig_response_curves.r

ex2_pres_map.pdf: r/fig_maps.r dat/figure/ne_50m_ocean dat/figure/ne_50m_lakes \
dat/figure/acersacr dat/rawData.rds results/naiveStats.rds
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

bin/stats: bin/stats.o bin/csv.o
	$(CC) $(CO) -lgsl -o bin/stats bin/stats.o bin/csv.o

bin/stats.o: src/statistics.cpp src/csv.hpp
	mkdir -p bin
	$(CC) $(CO) -c -o bin/stats.o src/statistics.cpp




