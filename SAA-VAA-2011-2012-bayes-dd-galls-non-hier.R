# SAA-VAA-2011-2012 Bayesian
# females have binomial survival
# females have poisson eggs
# egg survival depends on egg density and plant quality


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
	R ~ dunif(-20, 20)
	q ~ dunif(-20, 20)
	b ~ dunif(-20, 20)
	
	# Likelihood
	for(i in 1:nPlants){
		
		G[i] ~ dpois(lambda[i])
		lambda[i] <- exp(R + q * g11[i] + b * F[i]) * F[i]
	
	}

}
", fill = TRUE)
sink()


# Initial values
inits = function() list(R = 2, q=0, b=0)

# parms
parms = c('R', 'q', 'b')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 50000
nb = 0.1 * ni
thin = 10


# run MCMC with R2jags
m2fg = jags(data, inits=inits, parms, 'adult-addition-dd.txt', n.chains=3, n.iter=ni)

m2fg = update(m2fg, n.iter=50000)

traceplot(m2fg)

# convert rjags object to mcmc.list
m2fg.mcmc = as.mcmc(m2fg)
m2fgmat = as.matrix(m2fg.mcmc)

# use plotting methods from coda
densityplot(m2fg.mcmc)

pairs(m2fg$BUGSoutput$sims.list)

plot(pgr ~ f, data=d)
newf = seq(0,12,length=80)
pred.pgr = sapply(newf,
	function(z) mean( exp(m2fgmat[,'R'] + m2fgmat[,'b'] * z + m2fgmat[,'q'] * 5) ))
ci.pgr = sapply(newf,
	function(z) HPDI( exp(m2fgmat[,'R'] + m2fgmat[,'b'] * z + m2fgmat[,'q'] * 5) ))
lines(newf, pred.pgr)
lines(newf, ci.pgr[1,], lty=2)
lines(newf, ci.pgr[2,], lty=2)




