all: example1 example2

example1: example_1/ex1_Sampling.pdf example_1/ex1_precip.pdf example_1/ex1_map.pdf

example2: example_2/ex2.pdf


example_1/ex1_Sampling.pdf: example_1/ex1_makeSamplingFig.r example_1/ex1_globals.r \
example_1/dat/ex1_m1.rdata example_1/dat/ex1_m2.rdata
	cd example_1; Rscript ex1_makeSamplingFig.r ex1_Sampling.pdf

example_1/ex1_precip.pdf: example_1/ex1_makePrecipFig.r example_1/ex1_globals.r \
example_1/dat/ex1_m1.rdata example_1/dat/ex1_m2.rdata example_1/dat/ex1_mm.rdata 
	cd example_1; Rscript ex1_makePrecipFig.r ex1_precip.pdf

example_1/ex1_map.pdf: example_1/ex1_makeMapFig.r example_1/ex1_globals.r \
example_1/dat/ex1_m1.rdata example_1/dat/ex1_m2.rdata example_1/dat/ex1_mm.rdata
	cd example_1; Rscript ex1_makeMapFig.r ex1_map.pdf

example_1/dat/ex1_m1.rdata: example_1/ex1_m1.r example_1/ex1_globals.r 
	cd example_1; Rscript ex1_m1.r

example_1/dat/ex1_m2.rdata: example_1/ex1_m2.r example_1/ex1_globals.r 
	cd example_1; Rscript ex1_m2.r

example_1/dat/ex1_mm.rdata: example_1/ex1_mm.r example_1/ex1_globals.r \
example_1/dat/ex1_m2.rdata example_1/dat/ex1_m1.rdata 
	cd example_1; Rscript ex1_mm.r



#### individual scripts for example 2:
## order of operations to run the full analysis:
# 1. ex2_prepMapleData.r
# 2. ex2_naiveModel.r
# 3. ex2_mcmcIntegrated.r
# 4. ex2_processResults.r		
# 5. ex2_makeFigures.r

#### step 1: prepare the data from the raw data
example_2/dat/maple.rdata: example_2/ex2_prepMapleData.r example_2/dat/ex2_AceSac.csv \
example_2/dat/ex2_currentClim.rdata example_2/dat/ex2_futureClim.rdata example_2/ex2_Functions.r
	cd example_2; Rscript ex2_prepMapleData.r

#### step 2: use a stepwise regression to build the naive model
example_2/results/naiveModel.rdata: example_2/ex2_naiveModel.r example_2/ex2_Functions.r \
example_2/dat/maple.rdata
	cd example_2; Rscript ex2_naiveModel.r

### step 3: run the integrated model
### requires python, numpy, and scipy
### note that this is very CPU intensive; expect hours to days of runtime on a typical desktop computer
example2/dat/integratedModelData.csv: example_2/ex2_prepIntegrated.r results/naiveModel.rdata \
dat/maple.rdata
	cd example_2; Rscript ex2_prepIntegrated.r

example2/results/integratedModel.csv: example2/integratedModel.py example2/mcmc.py \
example2/dat/integratedModelData.csv
	cd example_2; python integrated_model.py

example_2/results/integratedModel.rdata: example_2/ex2_mcmcIntegrated.r \
example2/results/integratedModel.csv example_2/results/naiveModel.rdata
	cd example_2; Rscript ex2_mcmcIntegrated.r




# step 4: process results
# note that this step is very memory intensive; expect multiple GBs of memory usage, with lots of swapping
# runtime can be quite long, especially if there isn't enough system RAM and lots of swapping is needed
example_2/results/predictions.rdata: example_2/ex2_processResults.r example_2/dat/maple.rdata \
example_2/results/naiveModel.rdata example_2/results/integratedModel.rdata example_2/ex2_Functions.r
	cd example_2; Rscript ex2_processResults.r

# step 5: make the figures
example_2/ex2.pdf: example_2/ex2_makeFigures.r example_2/results/predictions.rdata \
example_2/dat/maple.rdata
	cd example_2; Rscript ex2_makeFigures.r