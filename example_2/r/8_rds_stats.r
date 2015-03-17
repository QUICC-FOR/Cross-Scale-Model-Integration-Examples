#!/usr/bin/Rscript

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


library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-d", "--datafile", default = "dat/mcmc/naiveData.rds", help="posterior stats to convert to rds format")

argList = parser$parse_args()

if(is.null(argList$datafile)) stop("Input file must be specified with -d")
outfile = paste(substr(argList$datafile, 1, nchar(argList$datafile) - 3), "rds", sep="")
dat = read.csv(argList$datafile, stringsAsFactors=FALSE)
saveRDS(dat, outfile)
