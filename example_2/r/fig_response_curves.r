library(coda)
intPres = readRDS("results/integratedPresentPosterior.rds")
intFut = readRDS("results/integratedFuturePosterior.rds")
naive = readRDS("results/naivePosterior.rds")
scaling = readRDS("dat/parameterScaling.rds")
rawDat = readRDS("dat/rawData.rds")

unscale = function(x, scale.fit)
{
	cen = attr(scale.fit, "scaled:center")
	sc = attr(scale.fit, "scaled:scale")
	(x * sc) + cen
}
predict.mc = function(b0, b1, b2, newdat)
{
	ylog = sapply(newdat, function(x) b0 + b1*x + b2 * x^2)
	plogis(ylog)
}

response_curve = function(v, v.name, lims, calib.range, draw.legend=FALSE)
{
#	color = c('#a6cee3', '#1f78b4', '#b2df8a')
	color = c('#66c2a5', '#fc8d62', '#8da0cb')
	calCols = c('#aaaaaa44', "#444444")
	bcol = paste(color, "88", sep="")
	pcol = paste(color, "66", sep="")
	v.unscaled = seq(lims[1], lims[2], length.out=lims[3])
	v.scaled = scale(v.unscaled, center = attr(scaling[[v]], "scaled:center"), scale = attr(scaling[[v]], "scaled:scale"))
	calib.range.unscaled = unscale(calib.range, scaling[[v]])

	postNames = paste(v, c('1', '2'), sep="")
	postPreds = list(
		naive = predict.mc(naive[,1], naive[,postNames[1]], naive[,postNames[2]], v.scaled),
		intPres = predict.mc(intPres[,1], intPres[,postNames[1]], intPres[,postNames[2]], v.scaled),
		intFut = predict.mc(intFut[,1], intFut[,postNames[1]], intFut[,postNames[2]], v.scaled))

	yy = lapply(postPreds, colMeans)
	quant = lapply(postPreds, function(x) t(apply(x, 2, quantile, c(0.025, 0.975))))
	plot(0,0, type='n', ylim=c(0,1), col=color[1], ylab = "Probability of presence", xlab=v.name, bty='n', xlim=range(v.unscaled))
	polygon(c(calib.range.unscaled, rev(calib.range.unscaled)), c(0,0,1,1), col=calCols[1], border=calCols[2], lwd=0.5)
	lines(v.unscaled, yy[[1]], col=color[1])
	lines(v.unscaled, yy[[2]], col=color[2])
	lines(v.unscaled, yy[[3]], col=color[3])
	polygon(c(v.unscaled, rev(v.unscaled)), c(quant[[1]][,1], rev(quant[[1]][,2])), col=pcol[1], border=bcol[1])
	polygon(c(v.unscaled, rev(v.unscaled)), c(quant[[2]][,1], rev(quant[[2]][,2])), col=pcol[2], border=bcol[2])
	polygon(c(v.unscaled, rev(v.unscaled)), c(quant[[3]][,1], rev(quant[[3]][,2])), col=pcol[3], border=bcol[3])
	if(draw.legend) 
		legend(600, 0.97, bty='o', legend=c("Naive", "Integrated-Present", "Integrated-Future", "Naive Calibration Range"), fill=c(color, calCols[1]), bg="#FFFFFF")
}

# figure out calibration ranges for the variables
ddeg.calib = with(rawDat$calib[complete.cases(rawDat$calib),], range(ddeg))
sum_prcp.calib = with(rawDat$calib[complete.cases(rawDat$calib),], range(sum_prcp))
pToPET.calib = with(rawDat$calib[complete.cases(rawDat$calib),], range(pToPET))

pdf(file="ex2_response.pdf", w=7, h=2.75)
par(mfrow=c(1,3), mar=c(4,4,0,0))
sm = 30 # controls how smooth the curves are; higher is smoother but slower
response_curve('ddeg', 'Degree days', c(0, 6000, sm), ddeg.calib, TRUE)
response_curve('sum_prcp', 'Summer precipitation (mm)', c(100, 700, sm), sum_prcp.calib)
response_curve('pToPET', 'p to PET ratio', c(0, 4, sm), pToPET.calib)
dev.off()


