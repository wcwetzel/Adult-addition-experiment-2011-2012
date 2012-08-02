# SAA-VAA-2011-2012 Bayesian - effect of galls2011 and density dependence

library(rjags)
library(R2jags)
library(R2WinBUGS)


d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

d = d[d$f > 0,]
d = d[-which(d$tag==1417),]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)

# model
sink('adult-addition-galls-dd.txt')
cat("
model {
	# Priors
	R ~ dunif(-20, 20)
	q ~ dunif(-20, 20)
	b ~ dunif(-20, 20)
	k <- exp(logk)
	logk ~ dnorm(0, 0.0001)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- mu[i] * rho[i]
		mu[i] <- exp(R + q * sqrt(g11[i]) + b * F[i]) * F[i]
		rho[i] ~ dgamma(k, k)
		#pgr[i] <- lambda[i] / F[i]
	}
}
", fill = TRUE)
sink()


# Initial values play with this
inits = function() list(R = 2, b=0, q=0, logk=0, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('R', 'b', 'q', 'k', 'logk')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 100000
nb = 0.1 * ni
thin = 100


####### R2jags #########
nb2fgnbsqrt = jags(data, inits, parms, 'adult-addition-galls-dd.txt', n.chains=3, n.iter=ni, n.thin=thin)

traceplot(nb2fgnbsqrt)


### diagnostics
# convert rjags object to mcmc.list
nb2fgnbsqrt.mcmc = as.mcmc(nb2fgnbsqrt)
nb2fgnbsqrtmat = as.matrix(nb2fgnbsqrt.mcmc)

# use plotting methods from coda
densityplot(nb2fgnbsqrt.mcmc)

pairs(nb2fgnbsqrt$BUGSoutput$sims.list)

autocorr.plot(nb2fgnbsqrt)
gelman.diag( nb2fgnb.mcmc)



#### plotting by females
plot(pgr ~ f, data=d)
newf = seq(0,12,length=80)

# for plants with 0 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(0)) ))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(0)) ))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

# for plants with 15 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(15)) ))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(15)) ))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

#### plotting by females not by pgr
plot(g12 ~ f, data=d)
newf = seq(0,12,length=80)

# for plants with 0 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(0)) * z))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(0))  * z))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

# for plants with 15 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(15))  * z))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * z + nb2fgnbsqrtmat[,'q'] * sqrt(15))  * z))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')


#### plotting by galls
plot(pgr ~ g11, data=d)
newg = seq(0,40,length=80)

# for plants with 1 female
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 1 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) ))
lines(newg, pred.pgr, col='darkred')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 1 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) ))
lines(newg, ci.pgr[1,], lty=2, col='darkred')
lines(newg, ci.pgr[2,], lty=2, col='darkred')


# for plants with 12 females
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgnbmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 12 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) ))
lines(newg, pred.pgr, col='blue')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgnbmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 12 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) ))
lines(newg, ci.pgr[1,], lty=2, col='blue')
lines(newg, ci.pgr[2,], lty=2, col='blue')


#### plotting by galls not pgr
plot(g12 ~ g11, data=d)
newg = seq(0,40,length=80)

# for plants with 1 female
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 1 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) * 1))
lines(newg, pred.pgr, col='darkred')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 1 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) * 1))
lines(newg, ci.pgr[1,], lty=2, col='darkred')
lines(newg, ci.pgr[2,], lty=2, col='darkred')

# for plants with 4 females
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 4 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) * 4))
lines(newg, pred.pgr, col='red')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 4 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) * 4))
lines(newg, ci.pgr[1,], lty=2, col='red')
lines(newg, ci.pgr[2,], lty=2, col='red')


# for plants with 12 females
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 12 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) * 12))
lines(newg, pred.pgr, col='blue')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2fgnbsqrtmat[,'R'] + nb2fgnbsqrtmat[,'b'] * 12 + nb2fgnbsqrtmat[,'q'] * sqrt(z)) * 12))
lines(newg, ci.pgr[1,], lty=2, col='blue')
lines(newg, ci.pgr[2,], lty=2, col='blue')


