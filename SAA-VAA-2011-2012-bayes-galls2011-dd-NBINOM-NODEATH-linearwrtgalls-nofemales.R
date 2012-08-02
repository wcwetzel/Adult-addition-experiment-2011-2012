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
	R ~ dunif(-20, 20)
	q ~ dunif(-20, 20)
	k <- exp(logk)
	logk ~ dnorm(0, 0.0001)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- mu[i] * rho[i]
		mu[i] <- (R + q * g11[i]) * F[i]
		rho[i] ~ dgamma(k, k)
	}
}
", fill = TRUE)
sink()


# Initial values play with this
inits = function() list(R = 2, q=0, logk=0, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('R', 'q', 'k', 'logk')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 100000
nb = 0.1 * ni
thin = 100


####### R2jags #########
nb2 = jags(data, inits, parms, 'adult-addition-galls-dd.txt', n.chains=3, n.iter=ni, n.thin=thin)

traceplot(nb2)


### diagnostics
# convert rjags object to mcmc.list
nb2.mcmc = as.mcmc(nb2)
nb2mat = as.matrix(nb2.mcmc)

# use plotting methods from coda
densityplot(nb2.mcmc)

pairs(nb2$BUGSoutput$sims.list)

autocorr.plot(nb2)
gelman.diag( nb2.mcmc)



#### plotting by females
plot(pgr ~ f, data=d)
newf = seq(0,12,length=80)

# for plants with 0 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2mat[,'R'] + nb2mat[,'b'] * z + nb2mat[,'q'] * 0) ))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2mat[,'R'] + nb2mat[,'b'] * z + nb2mat[,'q'] * 0) ))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

# for plants with 15 galls
pred.pgr = sapply(newf,
	function(z) mean( exp(nb2mat[,'R'] + nb2mat[,'b'] * z + nb2mat[,'q'] * 15) ))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( exp(nb2mat[,'R'] + nb2mat[,'b'] * z + nb2mat[,'q'] * 15) ))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

##### by females - not PGR
plot(g12 ~ f, data=d)
newf = seq(0,12,length=80)

# for plants with 0 galls
pred.pgr = sapply(newf,
	function(z) mean( (nb2mat[,'R'] + nb2mat[,'q'] * 0) * z))
lines(newf, pred.pgr, col='darkred')
ci.pgr = sapply(newf,
	function(z) HPDI( (nb2mat[,'R'] + nb2mat[,'q'] * 0) * z))
lines(newf, ci.pgr[1,], lty=2, col='darkred')
lines(newf, ci.pgr[2,], lty=2, col='darkred')

# for plants with 15 galls
pred.pgr = sapply(newf,
	function(z) mean( (nb2mat[,'R'] + nb2mat[,'q'] * 20) * z))
lines(newf, pred.pgr, col='blue')
ci.pgr = sapply(newf,
	function(z) HPDI( (nb2mat[,'R'] + nb2mat[,'q'] * 20) * z))
lines(newf, ci.pgr[1,], lty=2, col='blue')
lines(newf, ci.pgr[2,], lty=2, col='blue')


#### plotting by galls
plot(pgr ~ g11, data=d)
newg = seq(0,40,length=80)

pred.pgr = sapply(newg,
	function(z) mean( (nb2mat[,'R'] + nb2mat[,'q'] * z) ))
lines(newg, pred.pgr, col='darkred')
ci.pgr = sapply(newg,
	function(z) HPDI( (nb2mat[,'R'] + nb2mat[,'q'] * z) ))
lines(newg, ci.pgr[1,], lty=2, col='darkred')
lines(newg, ci.pgr[2,], lty=2, col='darkred')




#### plotting by galls not pgr
plot(g12 ~ g11, data=d)
newg = seq(0,40,length=80)

# for plants with 1 female
pred.pgr = sapply(newg,
	function(z) mean( (nb2mat[,'R'] + nb2mat[,'q'] * z) * 1))
lines(newg, pred.pgr, col='darkred')
ci.pgr = sapply(newg,
	function(z) HPDI( (nb2mat[,'R'] + nb2mat[,'q'] * z) * 1))
lines(newg, ci.pgr[1,], lty=2, col='darkred')
lines(newg, ci.pgr[2,], lty=2, col='darkred')

# for plants with 4 females
pred.pgr = sapply(newg,
	function(z) mean( (nb2mat[,'R'] + nb2mat[,'q'] * z) * 4))
lines(newg, pred.pgr, col='darkred')
ci.pgr = sapply(newg,
	function(z) HPDI( (nb2mat[,'R'] + nb2mat[,'q'] * z) * 4))
lines(newg, ci.pgr[1,], lty=2, col='darkred')
lines(newg, ci.pgr[2,], lty=2, col='darkred')


# for plants with 12 females
pred.pgr = sapply(newg,
	function(z) mean( exp(nb2mat[,'b'] * 12) * (nb2mat[,'R'] + nb2mat[,'q'] * z) * 12))
lines(newg, pred.pgr, col='blue')
ci.pgr = sapply(newg,
	function(z) HPDI( exp(nb2mat[,'b'] * 12) * (nb2mat[,'R'] + nb2mat[,'q'] * z) * 12))
lines(newg, ci.pgr[1,], lty=2, col='blue')
lines(newg, ci.pgr[2,], lty=2, col='blue')


