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

#### step 1: set up the models
dat/naive_model.rdata: 1-setup_models.r dat/AceSac.csv dat/currentClim.rdata \
dat/futureClim.rdata
	Rscript 1-setup_models.r

#### step 2: run the integrated model
# requires the integrated_model binary, which must be built before running these steps
# see the end of this makefile for this build procedure
# note that this is very CPU intensive
results/integratedModel.csv: bin/integrated_model dat/integratedData.csv \
dat/integratedPriors.csv dat/integratedInits.csv
	mkdir -p results
	./bin/integrated_model >results/integratedModel.csv
	
#### step 3: process the results (e.g., throw away burnin samples, thin results)
results/integratedModel.rdata: 3-process_integrated.r results/integratedModel.csv \
dat/naive_model.rdata
	Rscript 3-process_integrated.r
	
#### step 4: compute statistics (mean & sd) on integrated model predictions
results/integratedStats.csv: bin/stats results/integratedModelThinned.csv \
dat/predictionData.csv
	./bin/stats >results/integratedStats.csv

#### step 5: make the figure
ex2.pdf: 5-make_figure.r dat/naive_model.rdata results/integratedStats.csv \
results/integratedModel.rdata dat/ne_110m_ocean/ne_110m_ocean.shp
	Rscript 5-make_figure.r


## extra stuff just to make sure all dependencies are resolved	    
## these are created by the scripts above, so all are no-operation
dat/integratedPriors.csv: dat/naive_model.rdata
	;

dat/integratedInits.csv: dat/naive_model.rdata
	;

dat/integratedData.csv: dat/naive_model.rdata
	;
	
results/integratedModelThinned.csv: results/integratedModel.rdata
	;
	
dat/predictionData.csv: results/integratedModel.rdata
	;

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




