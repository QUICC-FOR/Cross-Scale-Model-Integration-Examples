# ex2_Transformations.r
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
# set of functions for transforming data for the second example


smithson_transform <- function(x, N = length(x), s = 0.5) {
	# transformation from supplemental material of Smithson and Verkuilen 2006
	# takes data on interval [0,1] and transforms to (0,1)
	# s is the scaling factor given in the same reference
	val <- (x * (N-1) + s)/N
	return(val)
}

smithson_reverse <- function(val, N = length(x), s = 0.5) {
	x <- (val * N - s)/(N-1)
	return(x)
}

rescale <- function(x, return.functions=FALSE) {
	# centers data on 0 with a standard deviation of 1
	# if return.functions is TRUE, returns a list of two functions giving forward and backward transformations
	# if not, just returns the transformed data
	
	xbar <- mean(x)
	xsd <- sd(x)
	
	transformations <- list(
		forward = (function(x) (x - xbar) / xsd),
		backward = (function(x) (x * xsd) + xbar))
	
	if(return.functions) {
		return(transformations)
	} else {
		return(transformations$forward(x))
	}
	
	return(result)
}

