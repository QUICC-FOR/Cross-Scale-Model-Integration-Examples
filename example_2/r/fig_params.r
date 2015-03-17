library(coda)
intPres = readRDS("results/integratedPresentPosterior.rds")
intFut = readRDS("results/integratedFuturePosterior.rds")
naive = readRDS("results/naivePosterior.rds")

intPreStats = summary(intPres)
naiveStats = summary(naive)
intFutStats = summary(intFut)

parameters = list(naive = naiveStats$statistics[,1], intPre = intPreStats$statistics[,1],
		intFut = intFutStats$statistics[,1])
lower = list(naive = naiveStats$quantiles[,1], intPre = intPreStats$quantiles[,1], 
		intFut = intFutStats$quantiles[,1])
upper = list(naive = naiveStats$quantiles[,5], intPre = intPreStats$quantiles[,5], 
		intFut = intFutStats$quantiles[,5])


# plot settings
pdfFileName = "ex2_params.pdf"
pdfWidth = 5.5
pdfHeight = 5.5
yle = 2.5
allVals = rbind(naiveStats, intPreStats, intFutStats)
ylims = c(-10,10)
pch=21
cols=c('#66c2a5', '#fc8d62', '#8da0cb', 'black', 'black', 'black')
outline="#888888"
oset = 0.2
xx = 1:nrow(naiveStats$statistics)
textYPos = apply(cbind(upper$naive, upper$intPre, upper$intFut), 1, max) + 0.2
bty='n'
texcex = 0.6
margin = c(0,3,0,0) + 0.5
cex.axis = 0.8
cex.lab = 0.8
mgp = c(2,0.5,0)
tcl=-0.3
cex.legend=0.8
cex=0.9

labs = rownames(naiveStats$statistics)
ones = grep("1$", labs)
labs[ones] = substr(labs[ones], 1, nchar(labs[ones]) - 1)
twos = grep("[2-9]$", labs)
labs[twos] = paste(substr(labs[twos], 1, nchar(labs[twos])-1), "^", substr(labs[twos], nchar(labs[twos]), nchar(labs[twos])), sep="")

# make plot
pdf(file=pdfFileName, width=pdfWidth,height=pdfHeight)

par(mar=margin, mgp=mgp, tcl=tcl)
plot(0,0, ylim=ylims, pch=pch, xlim=c(min(xx)-oset, max(xx)+oset), 
		bty='n', ylab="Parameter value", xaxt='n', xlab='', cex.axis=cex.axis, 
		cex.lab = cex.lab)
axis(side=1, pos=0, labels=F, lwd.ticks=0, lty=2, at=c(min(xx), max(xx)+1), col='grey', outer=T)
points(xx-oset, parameters$naive, pch=pch, bg=cols[1], col=outline, cex=cex)
points(xx, parameters$intPre, pch=pch, bg=cols[2], col=outline, cex=cex)
points(xx+oset, parameters$intFut, pch=pch, bg=cols[3], col=outline, cex=cex)
segments(xx-oset, lower$naive, xx-oset, upper$naive, col=cols[4])
segments(xx, lower$intPre, xx, upper$intPre, col=cols[5])
segments(xx+oset, lower$intFut, xx+oset, upper$intFut, col=cols[6])
text(xx, textYPos, parse(text=labs), pos=3, cex=texcex)
legend(max(xx)-2, 0.9*ylims[2], legend=c("Naive", "Integrated (Present)", 
		"Integrated (Future)"), pch=pch, pt.bg=cols[1:3], col=outline, cex=cex.legend)
dev.off()