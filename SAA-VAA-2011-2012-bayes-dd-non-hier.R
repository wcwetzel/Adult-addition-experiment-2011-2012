# SAA-VAA-2011-2012 Bayesian
# females have binomial survival
# females have poisson eggs
# egg survival depends on egg density and plant quality


library(rjags)
library(R2jags)
library(R2WinBUGS)
library(rethinking)

d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

d = d[d$f > 0,]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)

# model
sink('adult-addition-dd.txt')
cat("
model {
	# Priors
	#si ~ dunif(-20, 20)
	R ~ dunif(-20, 20)
	b ~ dunif(-20, 20)
	
	# Likelihood
	for(i in 1:nPlants){
		
		G[i] ~ dpois(lambda[i])
		lambda[i] <- exp(R + b * F[i]) * F[i] #* s[i]
		#logit(s[i]) <- si + b * F[i]
		#O[i] ~ dbin(p, F[i]) # female mortality
	}

}
", fill = TRUE)
sink()


# Initial values
inits = function() list(R = 2, b=-0.5) #, si=0.5)
	#O = ifelse(d$f>0, 1, 0))

# parms
parms = c('R', 'b')#, 'si')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12)

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
m2mat = as.matrix(m2.mcmc)

# use plotting methods from coda
densityplot(m2.mcmc)

pairs(m2$BUGSoutput$sims.list)

plot(pgr ~ f, data=d)
newf = seq(0,12,length=80)
pred.pgr = sapply(newf,
	function(z) mean( exp(m2mat[,'R'] + m2mat[,'b'] * z) ))
ci.pgr = sapply(newf,
	function(z) HPDI( exp(m2mat[,'R'] + m2mat[,'b'] * z) ))
lines(newf, pred.pgr)
lines(newf, ci.pgr[1,], lty=2)
lines(newf, ci.pgr[2,], lty=2)




