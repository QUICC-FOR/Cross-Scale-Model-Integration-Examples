
example1: example_1/ex1_Sampling.pdf example_1/ex1_precip.pdf example_1/ex1_map.pdf

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


