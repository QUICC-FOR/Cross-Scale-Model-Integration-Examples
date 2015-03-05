library(coda)
setwd("/Users/mtalluto/Dropbox/work/projects/Cross-Scale-Model-Integration-Examples_git/example_2")

load("results/integratedModel_Pres.rdata")
load("results/integratedModel_Fut.rdata")
load("dat/naive_model.rdata")
	

integratedStats_Pres = summary(integratedModel_Pres)$statistics
naiveStats = summary(naiveModel$model)$coefficients
integratedStats_Fut = summary(integratedModel_Fut)$statistics

# fix names and row orders
rownames(naiveStats) = sub("I?\\((.+)\\)", "\\1", rownames(naiveStats))
rownames(integratedStats_Pres)[1] = "Intercept"
rownames(integratedStats_Pres) = 
		sub("b_([a-zA-Z_]+)([2-3])", "\\1\\^\\2", rownames(integratedStats_Pres))
rownames(integratedStats_Pres) = 
		sub("b_([a-zA-Z_]+)(1)", "\\1", rownames(integratedStats_Pres))
naiveStats = naiveStats[rownames(integratedStats_Pres),,drop=FALSE]
rownames(integratedStats_Fut) = rownames(integratedStats_Pres)


# plot settings
pdfFileName = "ex2_params.pdf"
pdfWidth = 5.5
pdfHeight = 5.5
yle = 2.5
allVals = rbind(naiveStats, integratedStats_Pres, integratedStats_Fut)
ylims = range(c(allVals[,1] + yle*allVals[,2], allVals[,1] - yle*allVals[,2]))
pch=16
cols=c('cyan', 'blue', 'orange', 'black', 'black', 'black')
oset = 0.2
xx = 1:nrow(naiveStats)
textYPos = apply(cbind(naiveStats[,1], integratedStats_Pres[,1], integratedStats_Fut[,1]), 1, max) + naiveStats[,2]
bty='n'
texcex = 0.6
margin = c(0,3,0,0) + 0.5
cex.axis = 0.8
cex.lab = 0.8
mgp = c(2,0.5,0)
tcl=-0.3
cex.legend=0.7

# make plot
pdf(file=pdfFileName, width=pdfWidth,height=pdfHeight)

par(mar=margin, mgp=mgp, tcl=tcl)
plot(xx-oset, naiveStats[,1], ylim=ylims, pch=pch, col=cols[1], xlim=c(min(xx)-oset, max(xx)+oset), 
		bty='n', ylab="Parameter value", xaxt='n', xlab='', cex.axis=cex.axis, 
		cex.lab = cex.lab)
axis(side=1, pos=0, labels=F, lwd.ticks=0, lty=2, at=c(min(xx), max(xx)+1), col='grey', outer=T)
points(xx, integratedStats_Pres[,1], pch=pch, col=cols[2])
points(xx+oset, integratedStats_Fut[,1], pch=pch, col=cols[3])
segments(xx-oset, naiveStats[,1]+naiveStats[,2], xx-oset, naiveStats[,1]-naiveStats[,2], col=cols[4])
segments(xx, integratedStats_Pres[,1]+integratedStats_Pres[,2], xx, 
		integratedStats_Pres[,1]-integratedStats_Pres[,2], col=cols[5])
segments(xx+oset, integratedStats_Fut[,1]+integratedStats_Fut[,2], xx+oset, 
		integratedStats_Fut[,1]-integratedStats_Fut[,2], col=cols[6])
text(xx, textYPos, parse(text=rownames(naiveStats)), pos=3, cex=texcex)
legend(max(xx)-2, 0.6*ylims[1], legend=c("Naive", "Integrated (Present)", 
		"Integrated (Future)"), pch=pch, col=cols[1:3], cex=cex.legend)
dev.off()