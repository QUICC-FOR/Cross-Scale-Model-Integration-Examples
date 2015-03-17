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

library(sp)
library(raster)

load('dat/raw/currentClim.rdata')
load('dat/raw/futureClim.rdata')
presDat = readRDS("dat/raw/AceSac_pres.rds")
phenDat = read.csv("dat/raw/AceSac.csv")

climBaseCols = c("long", "lat", "Phenofit_CRU", "Phenofit_HadA2")
climPredictorCols = c("ddeg", "sum_prcp", "pToPET")
climCols= c(climBaseCols, climPredictorCols)

coordinates(phenDat) = c('long', 'lat')
coordinates(current_clim) = c("X", "Y")
coordinates(future_clim) = c("X", "Y")

# spatial merge all data frames
# note that this returns a normal data frame, not a spatial data frame
climDat = cbind(phenDat, over(phenDat, current_clim), over(phenDat, future_clim))

# add the ratio of annual precip to pet
climDat = within(climDat,
{
	pToPET <- an_prcp/pet
	fut_pToPET <- fut_an_prcp/fut_pet
})

## drop NAs and any extra columns
climDat = climDat[complete.cases(climDat), unlist(sapply(climCols, grep, colnames(climDat)))]

# center and scale all climate data, and transform phenofit observations
scales = list()
for(v in climPredictorCols)
{
	v.sc = scale(climDat[,v])
	climDat[,v] = v.sc
	fv = paste("fut_", v, sep="")
	climDat[,fv] = scale(climDat[,fv], center = attr(v.sc, "scaled:center"), scale = attr(v.sc, "scaled:scale"))
	scales[[v]] = v.sc
}
saveRDS(scales, "dat/parameterScaling.rds")

# transform phenofit observations to be on the open interval, instead of [0,1]
# transformation from supplemental material of Smithson and Verkuilen 2006
# s is the scaling factor given in the same reference
st = function(x, N = length(x), s = 0.5) (x * (N-1) + s)/N
climDat = within(climDat,
{
	Phenofit_CRU <- st(Phenofit_CRU)
	Phenofit_HadA2 <- st(Phenofit_HadA2)
})



# add a grid ID column to facilitate easy aggregation with the point-based presence data
climDat$grID = seq(1, nrow(climDat))
climPixels = climDat
coordinates(climPixels) = c("long", "lat")
gridded(climPixels) = TRUE
climRas = raster(climPixels, layer=which(names(climPixels) == 'grID'))

# now aggregate the presence-absence data within each grid cell at the climate data scale
# need to do it once for each of calibration and validation sets
coordinates(presDat) = c('lon', 'lat')
presDat$grID = extract(climRas, presDat)

set.seed(626786234)
validRows = sample(nrow(presDat), as.integer(1/3 * nrow(presDat)))
presList = list(
	valid = presDat[validRows,],
	calib = presDat[-validRows,],
	all = presDat)
	
prSum = lapply(presList, function(x) {
	res = aggregate(x$presence, by=list(x$grID), FUN=sum)
	colnames(res) = c("grID", "pres")
	res})
prCount = lapply(presList, function(x) {
	res = as.data.frame(table(x$grID))
	colnames(res) = c("grID", "count")
	res})
presence = lapply(1:length(prSum), function(i) merge(prSum[[i]], prCount[[i]], by="grID"))
baseData = lapply(presence, function(x) merge(climDat, x, by="grID", all.x=TRUE))
names(baseData) = names(presList)
saveRDS(baseData, file="dat/rawData.rds")



# prepare data for naive model
naiveDat = with(baseData$calib[complete.cases(baseData$calib),], data.frame(
	ddeg1 = ddeg,
	ddeg2 = ddeg^2,
	sum_prcp1 = sum_prcp,
	sum_prcp2 = sum_prcp^2,
	pToPET1 = pToPET,
	pToPET2 = pToPET^2,
	pres = pres,
	count = count)
)
write.csv(naiveDat, file="dat/mcmc/naiveData.csv", row.names = FALSE)
saveRDS(naiveDat, file="dat/mcmc/naiveData.rds") # to save us from having to load the csv later

naivePriors = data.frame(
		mean = rep(0, ncol(naiveDat) - 1),
		sd = c(10, rep(2.5, ncol(naiveDat) - 2)),
		dist = rep(1, ncol(naiveDat) - 1),
		row.names = c("intercept", colnames(naiveDat)[1:(ncol(naiveDat)-2)]))
write.csv(naivePriors, file='dat/mcmc/naivePriors.csv', row.names = FALSE)

# draw starts from a fairly conservative gaussian
naiveInits = rnorm(nrow(naivePriors), 0, 1)
write.csv(naiveInits, file='dat/mcmc/naiveInits.csv', row.names = FALSE)




# prepare prediction data
predictionDat = baseData$all[,unlist(sapply(c(climPredictorCols, "long", "lat"), function(v) grep(v, colnames(baseData$all))))]
write.table(predictionDat, "dat/mcmc/predictionData.csv", sep=",", col.names=FALSE, row.names=FALSE)

