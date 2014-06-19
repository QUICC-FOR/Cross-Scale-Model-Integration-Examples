all: example1 example2

example1:

example2: example_2/ex2.pdf



#### individual scripts for example 2:
## order of operations to run the full analysis:
# 1. ex2_prepMapleData.r
# 2. ex2_drawPseudoAbsences.r
# 3. ex2_setUpSDM.r
# 4. ex2_mcmcSDM.r
# 5. ex2_mcmcIntegrated.r
# 6. ex2_processResults.r		
# 7. ex2_makeFigures.r

#### step 1: prepare the data from the raw data
example_2/dat/maple.rdata: example_2/ex2_prepMapleData.r example_2/dat/ex2_AceSac.csv \
example_2/dat/ex2_currentClim.rdata example_2/dat/ex2_futureClim.rdata
	cd example_2; Rscript ex2_prepMapleData.r

#### step 2: draw pseudo-absences & transform variables for analysis
example_2/dat/maplePA.rdata: example_2/ex2_drawPseudoAbsences.r example_2/dat/maple.rdata \
example_2/ex2_Transformations.r
	cd example_2; Rscript ex2_drawPseudoAbsences.r

#### step 3: use a stepwise regression to choose form & starting values for the SDM
example_2/results/stepResults.rdata: example_2/ex2_setUpSDM.r example_2/ex2_Functions.r \
example_2/dat/maplePA.rdata
	cd example_2; Rscript ex2_setUpSDM.r

### steps 4 and 5 can be run concurrently; these will run the models then process the posterior samples
### step 4: run the naive sdm via jags 
### note that they are very CPU intensive; expect hours to days of runtime on a typical desktop computer
example_2/results/naiveModelResults.rdata: example_2/ex2_mcmcSDM.r \
example_2/ex2_Functions.r example_2/results/stepResults.rdata example_2/dat/maplePA.rdata
	cd example_2; Rscript ex2_mcmcSDM.r

example_2/results/integratedModelResults.rdata: example_2/ex2_mcmcIntegrated.r \
example_2/ex2_Functions.r example_2/results/stepResults.rdata example_2/dat/maplePA.rdata
	cd example_2; Rscript ex2_mcmcIntegrated.r

# step 6: process results
# note that this step is very memory intensive; expect multiple GBs of memory usage, with lots of swapping
# runtime can be quite long, especially if there isn't enough system RAM and lots of swapping is needed
example_2/results/predictions.rdata: example_2/ex2_processResults.r \
example_2/dat/maple.rdata example_2/dat/maplePA.rdata example_2/results/naiveModelResults.rdata \
example_2/results/integratedModelResults.rdata example_2/ex2_Functions.r
	cd example_2; Rscript ex2_processResults.r

# step 7: make the figures
example_2/ex2.pdf: example_2/ex2_makeFigures.r example_2/results/predictions.rdata \
example_2/dat/maple.rdata
	cd example_2; Rscript ex2_makeFigures.r