# ex1_makeMapFig.r
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
# produce density plots of the parameter posterior distributions


library(rjags)
load("dat/ex1_m1.rdata")
load("dat/ex1_m2.rdata")
load("dat/ex1_mm.rdata")




adj=1.5
m1p = do.call(rbind, m1Results)
m1d = apply(m1p, 2, density, adjust=adj)
m1d = lapply(m1d, function(x) {x$y = x$y / max(x$y); x})

m2p = do.call(rbind, m2Results)
m2d = apply(m2p, 2, density, adjust=adj)
m2d = lapply(m2d, function(x) {x$y = x$y / max(x$y); x})

mmp = do.call(rbind, mmResults)
mmd = apply(mmp, 2, density, adjust=adj)
mmd = lapply(mmd, function(x) {x$y = x$y / max(x$y); x})


lwd=1.5
col.m1 = '#1b9e77'
col.m2 = '#d95f02'
col.mm = '#7570b3'
j = c(1,0,0,2,3)
labs = c(expression(b[0]), expression(b[1]), expression(b[2]), expression(b[3]), expression(b[4]), expression(b[5]))

paperwidth = 5.5
dpi = 600
hToWRatio = 0.65
width = as.integer(dpi*paperwidth)
height = as.integer(width * hToWRatio)
fontsize = 10
png(w=width, h=height, file="ex1_params.png", pointsize=fontsize, res = dpi)

layout(matrix(c(2,4,6,1,3,5,8,10,11,7,9,11), nrow=4, byrow=T), heights=c(.2,1,.2,1))
for(i in 1:5)
{
	par(tcl=-0.2, mgp=c(2,0.5,0), cex=0.8, mar=c(3, 3, 0, 0.5), bty='n')
	plot(m1d[[i]], main="", sub="", xlab=labs[i], ylab="Relative Density", lwd=lwd, col=col.m1)
	if(i %in% c(1,4,5))
		lines(m2d[[j[i]]], lwd=lwd, col=col.m2)
	lines(mmd[[i]], lwd=lwd, col=col.mm)

	par(mar=c(0, 3, 0.5, 0.5))
	if(i %in% c(1,4,5))
	{
		boxplot(cbind(m1p[,i], m2p[,j[i]], mmp[,i]), pch='|', cex=0.4, xaxt='n', yaxt='n', col=c(col.m1, col.m2, col.mm), horizontal=T, lty=1, range=0)
	} else 
	{
		boxplot(cbind(m1p[,i], mmp[,i]), pch='|', cex=0.4, xaxt='n', yaxt='n', col=c(col.m1, col.mm), horizontal=T, lty=1, range=0)
	}
}
plot(0,0,xlab='', ylab='', type='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1))
legend(0, 1, legend=c("Naive Model", "Sub-model", "Integrated Model"), col=c(col.m1, col.m2, col.mm), lwd=lwd, bty='n', cex=0.8)


dev.off()

