# SAA-VAA-2011-2012 Bayesian - effect of galls2011 and density dependence

library(rjags)
library(R2jags)
library(R2WinBUGS)


d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

d = d[d$f > 0,]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)

# model
sink('adult-addition-galls-dd.txt')
cat("
model {
	# Priors
	p ~ dunif(0, 1)
	R ~ dunif(-20, 20)
	q ~ dunif(-20, 20)
	b ~ dunif(-20, 20)
	k <- exp(logk)
	logk ~ dnorm(0, 0.0001)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- mu[i] * rho[i]
		mu[i] <- exp(R + q * g11[i] + b * O[i]) * O[i]
		O[i] ~ dbin(p, F[i])
		rho[i] ~ dgamma(k, k)
		#pgr[i] <- lambda[i] / F[i]
	}
}
", fill = TRUE)
sink()


# Initial values play with this
inits = function() list(p = 0.5, R = 2, b=0, q=0, logk=0, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('p', 'R', 'b', 'q', 'k', 'logk', 'pgr')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 50000
nb = 0.1 * ni
thin = 10


####### R2jags #########
nb2fg = jags(data, inits, parms, 'adult-addition-galls-dd.txt', n.chains=3, n.iter=ni)
traceplot(nb2fg)


### diagnostics
# convert rjags object to mcmc.list
nb2fg.mcmc = as.mcmc(nb2fg)
nb2fgmat = as.matrix(nb2fg.mcmc)

# use plotting methods from coda
densityplot(nb2fg.mcmc)

pairs(nb2fg$BUGSoutput$sims.list)

autocorr.plot(nb2fg)
gelman.diag(nb2fg.mcmc)



#### plotting by females
plot(pgr ~ f, data=d)

points(apply(nb2fg$BUGSoutput$sims.list$pgr, 2, mean) ~ d$f, col='green')

newf = seq(0,12,length=80)



# for plants with 0 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * z + nb2fgmat[,'q'] * 0) ))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * z + nb2fgmat[,'q'] * 0) ))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

# for plants with 15 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * z + nb2fgmat[,'q'] * 15) ))
lines(newf, pred.pgr, col='blue')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * z + nb2fgmat[,'q'] * 15) ))
lines(newf, ci.pgr[1,], lty=2, col='blue')
lines(newf, ci.pgr[2,], lty=2, col='blue')

#### plotting by galls
plot(pgr ~ g11, data=d)
newg = seq(0,40,length=80)

# for plants with 0 females
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * 0 + nb2fgmat[,'q'] * z) ))
lines(newg, pred.pgr, col='darkred')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * 0 + nb2fgmat[,'q'] * z) ))
lines(newg, ci.pgr[1,], lty=2, col='darkred')
lines(newg, ci.pgr[2,], lty=2, col='darkred')


# for plants with 12 females
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * 12 + nb2fgmat[,'q'] * z) ))
lines(newg, pred.pgr, col='blue')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgmat[,'R'] + nb2fgmat[,'b'] * 12 + nb2fgmat[,'q'] * z) ))
lines(newg, ci.pgr[1,], lty=2, col='blue')
lines(newg, ci.pgr[2,], lty=2, col='blue')




##################################################
# with rjags
mrj.gdd = jags.model('adult-addition-galls-dd.txt', data, inits, n.chains=3, n.adapt=1000)

# burn in
update(mrj.gdd, 1000)

# take MCMC samples
orj.gdd = coda.samples(mrj.gdd, parms, n.iter=ni, thin=thin)


# convert confusing mcmc.list object into a nice, happy matrix
orjmatgdd = as.matrix(orj.gdd) #iters = TRUE will append the iteration numbers

# look at posterior means and 95% HPDI
apply(orjmatgdd, 2, mean)
HPDI(orj.gdd)

# look at joint posterior
plot(orjmatgdd[,'R'] ~ orjmatgdd[,'p'])
plot(orjmatgdd[,'R'] ~ orjmatgdd[,'b'])

# look at marginal posteriors
hist(orjmatgdd[,'R'])
hist(orjmatgdd[,'p'])
hist(orjmatgdd[,'b'])
hist(orjmatgdd[,'dd'])


plot(pgr ~ f, data=d)
newf = seq(0,12,length=80)
pred.pgr = sapply(newf,
	function(z) mean( exp(orjmatgdd[,'R'] + orjmatgdd[,'b'] * mean(d$g11) + orjmatgdd[,'dd'] * z) * orjmatgdd[,'p'] ))
ci.pgr = sapply(newf,
	function(z) HPDI( exp(orjmatgdd[,'R'] + orjmatgdd[,'b'] * mean(d$g11) + orjmatgdd[,'dd'] * z) * orjmatgdd[,'p'] ))
lines(newf, pred.pgr)
lines(newf, ci.pgr[1,], lty=2)
lines(newf, ci.pgr[2,], lty=2)


summary(orj.gdd)
gelman.diag(orj.gdd)
autocorr.plot(orj.gdd)


######### DIC

dic3 = dic.samples(mrj.gdd, n.iter=ni, thin=thin)

