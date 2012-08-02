# SAA-VAA-2011-2012 Bayesian - effect of galls2011


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
sink('adult-addition-galls.txt')
cat("
model {
	# Priors
	p ~ dunif(0, 1)
	R ~ dunif(0,100)
	b ~ dunif(-10, 10)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- (R + b * g11[i]) * O[i] # do I want a link function or not?
		O[i] ~ dbin(p, F[i])
	}
	
	
}
", fill = TRUE)
sink()


# Initial values play with this
inits = function() list(p = 0.5, R = 9, b=1, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('p', 'R', 'b')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 50000
nb = 0.1 * ni
thin = 10

# define rjags model
mrjg = jags.model('adult-addition-galls.txt', data, inits, n.chains=3, n.adapt=1000)

# burn in
update(mrjg, 1000)

# take MCMC samples
orjg = coda.samples(mrjg, parms, n.iter=ni, thin=10)


# convert confusing mcmc.list object into a nice, happy matrix
orjmat = as.matrix(orj) #iters = TRUE will append the iteration numbers

# look at posterior means and 95% HPDI
apply(orjmat, 2, mean)
HPDinterval(orj)

# look at joint posterior
plot(orjmat[,'R'] ~ orjmat[,'p'])
plot(orjmat[,'R'] ~ orjmat[,'b'])

# look at marginal posteriors
hist(orjmat[,1])
hist(orjmat[,'b'])

plot(pgr ~ g11, data=d)
newg11 = seq(0,40,length=80)
pred.pgr = sapply(newg11,
	function(z) mean( (orjmat[,'R'] + orjmat[,'b'] * z) * orjmat[,'p'] ))
ci.pgr = sapply(newg11,
	function(z) HPDI( (orjmat[,'R'] + orjmat[,'b'] * z) * orjmat[,'p'] ))
lines(newg11, pred.pgr)
lines(newg11, ci.pgr[1,], lty=2)
lines(newg11, ci.pgr[2,], lty=2)

pred.pgr.nodeath = sapply(newg11,
	function(z) mean( (orjmat[,'R'] + orjmat[,'b'] * z) ))
ci.pgr.nodeath = sapply(newg11,
	function(z) HPDI( (orjmat[,'R'] + orjmat[,'b'] * z) ))
lines(newg11, pred.pgr.nodeath)
lines(newg11, ci.pgr.nodeath[1,], lty=2)
lines(newg11, ci.pgr.nodeath[2,], lty=2)


summary(orj)
gelman.diag(orj)
autocorr.plot(orj)
crosscorr(orj)
crosscorr.plot(orj)
gelman.plot(orj)


####### DIC
dic2 = dic.samples(mrjg, n.iter=ni, thin=thin)



#--------------- R2jags --------------------#

# run MCMC
m2 = jags(data, inits=inits, parms, 'adult-addition-galls.txt', n.chains=3, n.iter=ni)

# plot samples
ts.plot(m2$BUGSoutput$sims.list$p)
ts.plot(m2$BUGSoutput$sims.list$R)
ts.plot(m2$BUGSoutput$sims.list$b)
ts.plot(m2$BUGSoutput$sims.list$deviance)

# display output
m2
plot(m2)

# traceplot
traceplot(m2)

# convert rjags object to mcmc.list
m2.mcmc = as.mcmc(m2)

# use plotting methods from coda
xyplot(m2.mcmc)
densityplot(m2.mcmc)

