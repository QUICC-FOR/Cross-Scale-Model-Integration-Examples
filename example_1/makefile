#### individual scripts for example 1:
## order of operations to run the full analysis:
# 1: Rscript 1-m1.r
# 2: Rscript 2-m2.r
# 3: Rscript 3-mm.r
# 4: Rscript 4-makeSamplingFig.r ex1_Sampling.pdf
# 5: Rscript 5-makePrecipFig.r ex1_precip.pdf
# 6: Rscript 6-makeMapFig.r ex1_map.pdf

SI: ELEtallutoSA1.pdf

figures: ex1_Sampling.pdf ex1_precip.pdf ex1_map.pdf

analysis: dat/ex1_mm.rdata

### build the SI appendix pdf
ELEtallutoSA1.pdf: tex/ELEtallutoSA1.tex
	cd tex; pdflatex ELEtallutoSA1
	mv tex/ELEtallutoSA1.pdf ELEtallutoSA1.pdf


# scripts for performing the analysis
dat/ex1_m1.rdata: 1-m1.r ex1_globals.r 
	Rscript 1-m1.r

dat/ex1_m2.rdata: 2-m2.r ex1_globals.r 
	Rscript 2-m2.r

dat/ex1_mm.rdata: 3-mm.r ex1_globals.r dat/ex1_m2.rdata dat/ex1_m1.rdata 
	Rscript 3-mm.r


# Scripts for building the figures
ex1_Sampling.pdf: 4-makeSamplingFig.r ex1_globals.r dat/ex1_m1.rdata dat/ex1_m2.rdata
	Rscript 4-makeSamplingFig.r ex1_Sampling.pdf

ex1_precip.pdf: 5-makePrecipFig.r ex1_globals.r dat/ex1_m1.rdata dat/ex1_m2.rdata \
dat/ex1_mm.rdata 
	Rscript 5-makePrecipFig.r ex1_precip.pdf

ex1_map.pdf: 6-makeMapFig.r ex1_globals.r dat/ex1_m1.rdata dat/ex1_m2.rdata \
dat/ex1_mm.rdata
	Rscript 6-makeMapFig.r ex1_map.pdf



