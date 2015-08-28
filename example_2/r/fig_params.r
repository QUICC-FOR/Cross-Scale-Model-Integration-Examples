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

library(coda)
load("results/posteriors.rdata")

intPreStats = summary(intPresPosterior)
naiveStats = summary(naivePosterior)
intFutStats = summary(intFutPosterior)

parameters = list(naive = naiveStats$statistics[,1], intPre = intPreStats$statistics[,1],
		intFut = intFutStats$statistics[,1])
lower = list(naive = naiveStats$quantiles[,1], intPre = intPreStats$quantiles[,1], 
		intFut = intFutStats$quantiles[,1])
upper = list(naive = naiveStats$quantiles[,5], intPre = intPreStats$quantiles[,5], 
		intFut = intFutStats$quantiles[,5])


# plot settings
pdfFileName = "ex2_params.pdf"
pdfWidth = 3.5
pdfHeight = 3.5
yle = 2.5
allVals = rbind(naiveStats, intPreStats, intFutStats)
ylims = c(-10,10)
pch=21
cols=c('#66c2a5', '#fc8d62', '#8da0cb', 'black', 'black', 'black')
outline="#888888"
oset = 0.5
xx = seq(1,nrow(naiveStats$statistics)*3,3)
textYPos = apply(cbind(upper$naive, upper$intPre, upper$intFut), 1, max) + 0.2
bty='n'
texcex = 0.5
margin = c(0,1.75,0,0) + 0.5
cex.axis = 0.6
cex.lab = 0.7
mgp = c(1.25,0.5,0)
tcl=-0.3
cex.legend=0.55
cex=0.9
lwd.axis=0.75

labs = rownames(naiveStats$statistics)
ones = grep("1$", labs)
labs[ones] = substr(labs[ones], 1, nchar(labs[ones]) - 1)
twos = grep("[2-9]$", labs)
labs[twos] = paste(substr(labs[twos], 1, nchar(labs[twos])-1), "^", substr(labs[twos], nchar(labs[twos]), nchar(labs[twos])), sep="")

# make plot
pdf(file=pdfFileName, width=pdfWidth,height=pdfHeight)

par(mar=margin, mgp=mgp, tcl=tcl)
plot(0,0, ylim=ylims, type='n', pch=pch, xlim=c(min(xx)-oset, max(xx)+oset), 
		bty='n', ylab="Parameter value", xaxt='n', yaxt='n', xlab='', cex.lab = cex.lab)
axis(side=1, pos=0, labels=F, lwd.ticks=0, lty=2, at=c(min(xx), max(xx)+1), col='grey', outer=T)
axis(side=2, labels=T, lwd.ticks=lwd.axis, lwd=lwd.axis, cex.axis=cex.axis)
points(xx-oset, parameters$naive, pch=pch, bg=cols[1], col=outline, cex=cex)
points(xx, parameters$intPre, pch=pch, bg=cols[2], col=outline, cex=cex)
points(xx+oset, parameters$intFut, pch=pch, bg=cols[3], col=outline, cex=cex)
segments(xx-oset, lower$naive, xx-oset, upper$naive, col=cols[4])
segments(xx, lower$intPre, xx, upper$intPre, col=cols[5])
segments(xx+oset, lower$intFut, xx+oset, upper$intFut, col=cols[6])
text(xx, textYPos, parse(text=labs), pos=3, cex=texcex)
legend(max(xx)-9, -7, legend=c("Naive", "Integrated (Present)", "Integrated (Future)"),
		bty='o', pch=pch, pt.bg=cols[1:3], col=outline, cex=cex.legend, box.lwd=lwd.axis)
dev.off()