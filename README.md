Cross-Scale-Model-Integration-Examples
======================================

These scripts are provided to reproduce the examples given in the paper "Cross-scale 
integration of data and knowledge for predicting species ranges." All source code is
licensed under the GPL3, details are provided in the LICENSE file.

For convenience, a makefile is provided to demonstrate how to build the examples and
provide an automated way to generate all necessary files. The command `make all` will 
build both examples, while `make example1` and `make example2` will build each example 
individually. Note that these will likely take a long time and consume a lot of memory, 
so plan your usage accordingly.

Dependencies:
-------------
* R version 3.x
* JAGS
* Rscript (optional, but recommended. Included with most R installations)

R packages:
-----------
* sp
* glm2
* fields
* rjags

Running the models:
-------------------
Both examples are provided as a collection of scripts. Order matters (more so for example 
2, but to some extent for example 1 as well). A list of all scripts and the proper order
follows. The examples are easiest to run by invoking them with `Rscript` at the command line, as follows:

	cd path/to/desired/example
    Rscript file.r

Alternatively, you can run the scripts from the console by using `source()`, or by pasting
one line at a time. Just make sure you set your working directory to the appropriate
place, either `example_1` or `example_2`.

Example 1:
----------
1. `Rscript ex1_m1.r` -- Run the naive model
1. `Rscript ex1_m2.r` -- Run the mechanistic submodel
3. `Rscript ex1_mm.r` -- Run the metamodel
4. `Rscript ex1_makeSamplingFig.r ex1_Sampling.pdf` -- Make figure 2 from the manuscript
4. `Rscript ex1_makePrecipFig.r ex1_precip.pdf` -- Make figure 3 from the manuscript
4. `Rscript ex1_makeMapFig.r ex1_map.pdf` -- Make figure 4 from the manuscript

Example 2:
----------
1. `Rscript ex2_prepMapleData.r` -- Loads data from original sources and formats it for the analysis
2. `Rscript ex2_naiveModel.r` -- Use conventional methods to select the form for the naive metamodel
3. `Rscript ex2_mcmcIntegrated.r` -- Run the mcmc for the integrated model
4. `Rscript ex2_processResults.r` -- Produce predictions from the posterior parameter distributions
5. `Rscript ex2_makeFigure.r` -- Produce figure 5 from the manuscript
