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

climBaseCols = c("long", "lat", "PresObs", "Phenofit_CRU", "Phenofit_HadA2")
climPredictorCols = c("ddeg", "an_prcp", "pToPET")
climCols= c(climBaseCols, climPredictorCols)

# subset the data using lat/long
# we don't need to fit the model for Alaska and Greenland, for example
xlims = c(-130,-50)
ylims = c(20,65)
climInd = with(current_clim, which(X >= xlims[1] & X <= xlims[2] & Y >= ylims[1] & Y <= ylims[2]))
current_clim = current_clim[climInd,]
future_clim = future_clim[climInd,]

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
coordinates(presDat) = c('lon', 'lat')
presDat$grID = extract(climRas, presDat)
prSum = aggregate(presDat$presence, by=list(presDat$grID), FUN=sum)
colnames(prSum) = c("grID", "pres")
prCount = as.data.frame(table(presDat$grID))
colnames(prCount) = c("grID", "count")
presence = merge(prSum, prCount, by="grID")
baseData = merge(climDat, presence, by="grID", all.x=TRUE)
baseData = within(baseData, 
{
	countAug = 1 + ifelse(is.na(count), 0, count)
	presAug = PresObs + ifelse(is.na(pres), 0, pres)
})

set.seed(626786234)
validRows = sample(nrow(baseData), as.integer(1/3 * nrow(baseData)))
baseData = list(
	valid = baseData[validRows,],
	calib = baseData[-validRows,],
	all = baseData)

saveRDS(baseData, file="dat/rawData.rds")



# prepare data for naive model
naiveDat = with(baseData$calib[complete.cases(baseData$calib),], data.frame(
	ddeg1 = ddeg,
	ddeg2 = ddeg^2,
	ddeg3 = ddeg^3,
	an_prcp1 = an_prcp,
	an_prcp2 = an_prcp^2,
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
saveRDS(predictionDat, file="dat/predictionData.rds")

