# ex2_drawPseudoAbsences.r
# 
#   Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
# 
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or (at
#   your option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
# 
# 
# this script produces a pseudoabsence draw required for the naive SDM portion of the 
# meta-model

load("dat/maple.rdata")
source("ex2_Transformations.r")


# set the seed for the random number generator, for reproducibility
# comment this line out for random pseudo absence draws
# this seed was chosen by a single call to:
# runif(1,0,.Machine$integer.max)
set.seed(626786234)


prevalence <- 0.33					## prevalence of ABSENCES in the post-draw dataset
outFile <- "dat/maplePA.rdata"
presence.column <- 'PresObs'		## name or index of the column containing presence/absence data


dat.presence <- maple[maple[,presence.column]==1,]
dat.absence <- maple[maple[,presence.column]==0,]

# decide how many pseudo-absences to draw
n.draw <- as.integer(nrow(dat.presence)/(1-prevalence)) - nrow(dat.presence)

# make sure we have enough data
if(n.draw > nrow(dat.absence)) stop(paste(nrow(dat.absence),  " rows in pseudo absence data, 
	input prevalence of ", prevalence, " requires ", n.draw, " pseudo absences.", sep=""))

maplePA <- rbind(dat.presence, dat.absence[sample(1:nrow(dat.absence), n.draw, replace=F),])

## transform data
## all climate data gets zero-centered (based on the PRESENT climatic conditions, including for future climate)
## phenofit predictions get "squeezed" to fit on the (0,1) interval instead of [0,1]
transformations <- lapply(maplePA[,6:11], rescale, return.functions=TRUE)
for(nm in names(transformations)) {
	nmfut <- paste("fut_", nm, sep="")
	maplePA[,nm] <- transformations[[nm]]$forward(maplePA[,nm])
	maplePA[,nmfut] <- transformations[[nm]]$forward(maplePA[,nmfut])
}

maplePA$Phenofit_CRU <- smithson_transform(maplePA$Phenofit_CRU)
maplePA$Phenofit_HadA2 <- smithson_transform(maplePA$Phenofit_HadA2)


save(maplePA, transformations, file=outFile)