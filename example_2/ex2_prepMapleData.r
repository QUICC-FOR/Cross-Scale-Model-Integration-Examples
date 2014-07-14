# ex2_prepMapleData.r
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
# this script extracts the data, combines it with the climate, prepares it for analysis
# and exports it to an .rdata file


library(sp)

maple <- read.csv("dat/ex2_AceSac.csv")
load("dat/ex2_currentClim.rdata")
load("dat/ex2_futureClim.rdata")
source("ex2_Transformations.r")

# set the seed for the random number generator, for reproducibility
# comment this line out for random pseudo absence draws
# this seed was chosen by a single call to:
# runif(1,0,.Machine$integer.max)
set.seed(626786234)


# make them spatial objects
maple <- SpatialPointsDataFrame(maple[,1:2], maple)
current_clim <- SpatialPointsDataFrame(current_clim[,1:2],current_clim)
future_clim <- SpatialPointsDataFrame(future_clim[,1:2],future_clim)

## perform a spatial join of the features
maple <- cbind(maple, over(maple, current_clim), over(maple, future_clim))

## drop NAs & extra latlong columns
dropRows <- which( (apply(maple, 1, function(x) sum(is.na(x)))) > 0)	# find rows with at least 1 NA
maple <- maple[-dropRows, c(1:5, 10:15, 18:23)]


## code below will be useful for splitting off a validation dataset
## dat.presence <- maple[maple[,presence.column]==1,]
## dat.absence <- maple[maple[,presence.column]==0,]


## transform data
## all climate data gets zero-centered (based on the PRESENT climatic conditions, including for future climate)
## phenofit predictions get "squeezed" to fit on the (0,1) interval instead of [0,1]
transformations <- lapply(maple[,6:11], rescale, return.functions=TRUE)
for(nm in names(transformations)) {
	nmfut <- paste("fut_", nm, sep="")
	maple[,nm] <- transformations[[nm]]$forward(maple[,nm])
	maple[,nmfut] <- transformations[[nm]]$forward(maple[,nmfut])
}

maple$Phenofit_CRU <- smithson_transform(maple$Phenofit_CRU)
maple$Phenofit_HadA2 <- smithson_transform(maple$Phenofit_HadA2)


# weight data so that the total weight of presences and absences is the same
numAbs = sum(maple$PresObs == 0)
numPres = sum(maple$PresObs == 1)
wght = as.integer(numAbs / numPres)
maple$weightedPresence = wght*maple$PresObs
maple$weightedN = rep(1, nrow(maple))
maple$weightedN[maple$PresObs == 1] = wght

save(maple, transformations, file="dat/maple.rdata")