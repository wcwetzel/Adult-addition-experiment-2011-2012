# SAA-VAA-2011-2012 Bayesian - effect of female density


library(rjags)
library(R2jags)
library(R2WinBUGS)
library(rethinking)

d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

#d = d[d$f > 0,]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)

# model
sink('adult-addition-dd.txt')
cat("
model {
	# Priors
	p ~ dunif(0,1)
	R ~ dunif(-20, 20)
	b ~ dunif(-20, 20)
	
	
	# Likelihood
	for(i in 1:nPlants){
		
		G[i] ~ dpois(lambda[i])
		lambda[i] <- pgr[i] * O[i] # galls per plant
		pgr[i] <- exp(R + b * O[i]) # per capita galling rate of surviving females
		O[i] ~ dbin(p, F[i]) # female mortality
	}
	
	
}
", fill = TRUE)
sink()


# Initial values
inits = function() list(p = 0.5, R = 2, b=-1, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('p', 'R', 'b', 'pgr')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12)

# MCMC controls
ni = 100000
nb = 0.1 * ni
thin = 10


# run MCMC with R2jags
m2 = jags(data, inits=inits, parms, 'adult-addition-dd.txt', n.chains=3, n.iter=ni)


# display output
m2

# traceplot
traceplot(m2exp)

# convert rjags object to mcmc.list
m2.mcmc = as.mcmc(m2)

# use plotting methods from coda
densityplot(m2.mcmc)

pairs(m2$BUGSoutput$sims.list)


m2$BUGSoutput$sims.list$pgr
pgr.mean = apply(m2exp$BUGSoutput$sims.list$pgr,2, mean)

newf = seq(0,12,length=80)

m2mat = as.matrix(m2.mcmc)


plot(pgr ~ f, data=d)
points(pgr.mean ~ d$f, col='blue')

pred.pgr = sapply(newf,
	function(z) mean( exp(m2mat[,'R'] + m2mat[,'b'] * z * m2mat[,'p']) ))
ci.pgr = sapply(newf,
	function(z) HPDI( exp(m2mat[,'R'] + m2mat[,'b'] * z * m2mat[,'p']) ))
lines(newf, pred.pgr)
lines(newf, ci.pgr[1,], lty=2)
lines(newf, ci.pgr[2,], lty=2)

pred.pgr.nodeath = sapply(newg11,
	function(z) mean( (orjmat[,'R'] + orjmat[,'b'] * z) ))
ci.pgr.nodeath = sapply(newg11,
	function(z) HPDI( (orjmat[,'R'] + orjmat[,'b'] * z) ))
lines(newg11, pred.pgr.nodeath)
lines(newg11, ci.pgr.nodeath[1,], lty=2)
lines(newg11, ci.pgr.nodeath[2,], lty=2)



# summary stuff from rjags change orjmat to m2mat
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



















# rjags set up
# define rjags model
mrjg = jags.model('adult-addition-galls.txt', data, inits, n.chains=3, n.adapt=1000)

# burn in
update(mrjg, 1000)

# take MCMC samples
orjg = coda.samples(mrjg, parms, n.iter=ni, thin=10)


# convert confusing mcmc.list object into a nice, happy matrix
orjmat = as.matrix(orj) #iters = TRUE will append the iteration numbers