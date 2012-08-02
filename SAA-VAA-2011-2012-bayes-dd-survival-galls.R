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
	p ~ dunif(0, 1)
	si ~ dunif(-20, 20)
	R ~ dunif(-20, 20)
	q ~ dunif(-20, 20)
	b ~ dunif(-20, 20)
	
	# Likelihood
	for(i in 1:nPlants){
		
		G[i] ~ dpois(lambda[i])
		lambda[i] <- exp(R) * O[i] * s[i]
		logit(s[i]) <- si + b * O[i] + q * g11[i]
		O[i] ~ dbin(p, F[i]) # female mortality
	}

}
", fill = TRUE)
sink()


# Initial values
inits = function() list(p = 0.5, R = 2, b=-0.5, q=0.5, si=0.5, 
	O = ifelse(d$f>0, 1, 0))

# parms
parms = c('p', 'R', 'b', 'si', 'q')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 10000
nb = 0.1 * ni
thin = 10


# run MCMC with R2jags
m2 = jags(data, inits=inits, parms, 'adult-addition-dd.txt', n.chains=3, n.iter=ni)

m2 = update(m2, n.iter=50000)

traceplot(m2)

# convert rjags object to mcmc.list
m2.mcmc = as.mcmc(m2)

# use plotting methods from coda
densityplot(m2.mcmc)

pairs(m2$BUGSoutput$sims.list)






