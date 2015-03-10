


	
	
	## need to do this, but it will have to happen AFTER the naive model gets fit
# 	priors = data.frame(mean=c(naiveModel$model$coefficients), sd = summary(
# 		naiveModel$model)$coefficients[,2])
# 	# put the priors into the same order as the other pieces
# 	row.names(priors) = sapply(row.names(priors), function(x) {
# 		naiveModel$variables$parameter[naiveModel$variables$coefName == x]})
# 	priors = priors[c(1, order(row.names(priors)))[-(nrow(priors)+1)],]	
# 
# 	inits = data.frame(inits = unlist(naiveModel$starts()))
# 	write.csv(priors, file='dat/integratedPriors.csv', row.names = FALSE)
# 	write.csv(inits, file='dat/integratedInits.csv', row.names = FALSE)


