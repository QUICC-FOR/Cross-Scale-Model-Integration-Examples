library(coda)
load("results/posteriors.rdata")
scaling = readRDS("dat/parameterScaling.rds")
rawDat = readRDS("dat/rawData.rds")

# varnames = c('ddeg', 'an_prcp', 'pToPET')

unscale = function(x, scale.fit)
{
	cen = attr(scale.fit, "scaled:center")
	sc = attr(scale.fit, "scaled:scale")
	(x * sc) + cen
}

predict.mcmc = function(b, newdat)
{
	with(newdat, plogis(b[1] + b[2]*ddeg + b[3]*ddeg^2 + b[4]*ddeg^3 + 
			b[5]*an_prcp + b[6]*an_prcp^2 + b[7]*pToPET + b[8]*pToPET^2))
}

npPr = apply(naivePosterior, 1, predict.mcmc, newdat=newdat)

response_curve = function(v, v.name, lims, calib.range, draw.legend=FALSE)
{
#	color = c('#a6cee3', '#1f78b4', '#b2df8a')
	color = c('#66c2a5', '#fc8d62', '#8da0cb')
	calCols = c('#aaaaaa44', "#444444")
	bcol = paste(color, "88", sep="")
	pcol = paste(color, "66", sep="")
	v.unscaled = seq(lims[1], lims[2], length.out=lims[3])
	v.scaled = scale(v.unscaled, center = attr(scaling[[v]], "scaled:center"), scale = attr(scaling[[v]], "scaled:scale"))
	newdat = data.frame(
		ddeg = rep(0,length(v.scaled)),
		an_prcp = rep(0,length(v.scaled)),
		pToPET = rep(0,length(v.scaled)))
	newdat[,v] = v.scaled
	calib.range.unscaled = unscale(calib.range, scaling[[v]])

	postNames = paste(v, c('1', '2'), sep="")
	postPreds = lapply(list(naivePosterior, intPresPosterior, intFutPosterior), 
			function(x) apply(x, 1, predict.mcmc, newdat=newdat))
	names(postPreds) = c('naive', 'pres', 'fut')
	yy = lapply(postPreds, rowMeans)
	quant = lapply(postPreds, function(x) t(apply(x, 1, quantile, c(0.025, 0.975))))
	
	plot(0,0, type='n', ylim=c(0,1), col=color[1], ylab = "Probability of presence", xlab=v.name, bty='n', xlim=range(v.unscaled))
	polygon(c(calib.range.unscaled, rev(calib.range.unscaled)), c(0,0,1,1), col=calCols[1], border=calCols[2], lwd=0.5)
	lines(v.unscaled, yy[[1]], col=color[1])
	lines(v.unscaled, yy[[2]], col=color[2])
	lines(v.unscaled, yy[[3]], col=color[3])
	polygon(c(v.unscaled, rev(v.unscaled)), c(quant[[1]][,1], rev(quant[[1]][,2])), col=pcol[1], border=bcol[1])
	polygon(c(v.unscaled, rev(v.unscaled)), c(quant[[2]][,1], rev(quant[[2]][,2])), col=pcol[2], border=bcol[2])
	polygon(c(v.unscaled, rev(v.unscaled)), c(quant[[3]][,1], rev(quant[[3]][,2])), col=pcol[3], border=bcol[3])
	if(draw.legend) 
		legend(550, 0.97, bty='o', cex=0.7, legend=c("Naive", "Integrated-Present", "Integrated-Future", "Naive Calibration Range"), fill=c(color, calCols[1]), bg="#FFFFFF")
}

# figure out calibration ranges for the variables
ddeg.calib = with(rawDat$calib[complete.cases(rawDat$calib),], range(ddeg))
an_prcp.calib = with(rawDat$calib[complete.cases(rawDat$calib),], range(an_prcp))
pToPET.calib = with(rawDat$calib[complete.cases(rawDat$calib),], range(pToPET))

pdf(file="ex2_response.pdf", w=7, h=2.75)
par(mfrow=c(1,3), mar=c(4,4,0,0))
sm = 100 # controls how smooth the curves are; higher is smoother but slower
# note: 6000 covers sugar maple all right, but misses the south of NA for the pres and up to 35 degrees lat for the future
# 7000 gets you to 30 degrees. 9000 is needed to get all the way to florida
response_curve('ddeg', 'Degree days', c(0, 7000, sm), ddeg.calib, FALSE)
response_curve('an_prcp', 'Annual precipitation (mm)', c(0, 2000, sm), an_prcp.calib, TRUE)
response_curve('pToPET', 'p to PET ratio', c(0, 4, sm), pToPET.calib)
dev.off()


