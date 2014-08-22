setwd("/Users/mtalluto/Documents/git_projects/Cross-Scale-Model-Integration-Examples/example_2/src")
TEST_LEN = 100

inv_logit = function(x) exp(x) / (1 + exp(x))

testDat = data.frame(x1 = rnorm(TEST_LEN))
# starting with only 0 or 1 probabilities for now so we get 100% consistent datasets
testDat$probs = rbinom(TEST_LEN, 1, prob = inv_logit(0.5 + 2.3 * testDat$x1))

# make a model to show us the "right" answers
summary(glm(probs ~ x1, family=binomial, data=testDat))

priors = data.frame(mean=c(0,0), sd = c(100000, 100000))

write.csv(testDat, file='testData.csv', row.names = FALSE)
write.csv(priors, file='testPriors.csv', row.names = FALSE)