# SAA-VAA-2011-2012 Bayesian


library(rjags)
library(R2jags)
library(R2WinBUGS)


d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

#d = d[d$f > 0,]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)

# model
sink('adult-addition.txt')
cat("
model {
	# Priors
	p ~ dunif(0, 1)
	R ~ dunif(0,100)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- R * O[i]
		O[i] ~ dbin(p, F[i])
	}
}
", fill = TRUE)
sink()


# Initial values play with this
inits = function() list(p = 0.5, R = 3.1, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('p', 'R')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12)

# MCMC controls
ni = 50000
nb = 0.1 * ni
thin = 10

# define rjags model
mrj = jags.model('adult-addition.txt', data, inits, n.chains=3, n.adapt=1000)

# burn in
update(mrj, 1000)

# take MCMC samples
orj = coda.samples(mrj, parms, n.iter=ni, thin=10)


# convert confusing mcmc.list object into a nice, happy matrix
orjmat = as.matrix(orj) #iters = TRUE will append the iteration numbers

# look at posterior means and 95% HPDI
apply(orjmat, 2, mean)
HPDinterval(orj)

# look at joint posterior
plot(orjmat[,1] ~ orjmat[,2])

# look at marginal posteriors
hist(orjmat[,1])
hist(orjmat[,2])



summary(orj)
gelman.diag(orj)
autocorr.plot(orj)
crosscorr(orj)
crosscorr.plot(orj)
gelman.plot(orj)

