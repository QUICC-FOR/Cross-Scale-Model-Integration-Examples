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


# make them spatial objects
maple <- SpatialPointsDataFrame(maple[,1:2], maple)
current_clim <- SpatialPointsDataFrame(current_clim[,1:2],current_clim)
future_clim <- SpatialPointsDataFrame(future_clim[,1:2],future_clim)

## perform a spatial join of the features
maple <- cbind(maple, over(maple, current_clim), over(maple, future_clim))

## drop NAs & extra latlong columns
dropRows <- which( (apply(maple, 1, function(x) sum(is.na(x)))) > 0)	# find rows with at least 1 NA
maple <- maple[-dropRows, c(1:5, 10:15, 18:23)]

save(maple, file="dat/maple.rdata")