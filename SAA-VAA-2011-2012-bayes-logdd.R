
# SAA-VAA-2011-2012 Bayesian - effect of female density on log(pgr+0.01)


library(rjags)
library(R2jags)
library(R2WinBUGS)
library(rethinking)

d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

d = d[d$f > 0,]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.5)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)


# model
sink('adult-addition-dd-logPGR.txt')
cat("
model {
	# Priors
	R ~ dunif(-100,100)
	b ~ dunif(-100, 100)
	tau ~ dunif(0, 100)
	
	
	# Likelihood
	for(i in 1:nPlants){
		
		pgrl[i] ~ dnorm(mu[i], tau)
		mu[i] <- R + b * F[i] 
	}
	
	
}
", fill = TRUE)
sink()

# model
sink('adult-addition-dd-logPGR.txt')
cat("
model {
	# Priors
	R ~ dunif(-100,100)
	b ~ dunif(-100, 100)
	tau ~ dunif(0, 100)
	c ~ dunif(-100, 100)
	
	
	# Likelihood
	for(i in 1:nPlants){
		
		pgrl[i] ~ dnorm(mu[i], tau)
		mu[i] <- R + b * F[i] + c * g11[i]
	}
	
	
}
", fill = TRUE)
sink()


# Initial values
inits = function() list(R = 2, b=-1, c=0)

# parms
parms = c('R', 'b', 'c')

# data
data = list(F=d$f, nPlants=nrow(d), pgrl=d$pgrl, g11=d$g11)

# MCMC controls
ni = 50000
nb = 0.1 * ni
thin = 10


# run MCMC with R2jags
mpgrl2 = jags(data, inits=inits, parms, 'adult-addition-dd-logPGR.txt', n.chains=3, n.iter=ni)



traceplot(mpgrl2)
mpgrl.mcmc = as.mcmc(mpgrl2)
densityplot(mpgrl.mcmc)
pairs(mpgrl$BUGSoutput$sims.list)



mpgrlmat = as.matrix(mpgrl.mcmc)


plot(pgrl ~ f, data=d)
newf = seq(0,12,length=80)

pred.pgr = sapply(newf,
	function(z) mean( mpgrlmat[,'R'] + mpgrlmat[,'b'] * z + mpgrlmat[,'c'] * 7))
ci.pgr = sapply(newf,
	function(z) HPDI( mpgrlmat[,'R'] + mpgrlmat[,'b'] * z + mpgrlmat[,'c'] * 7))
lines(newf, pred.pgr)
lines(newf, ci.pgr[1,], lty=2)
lines(newf, ci.pgr[2,], lty=2)



