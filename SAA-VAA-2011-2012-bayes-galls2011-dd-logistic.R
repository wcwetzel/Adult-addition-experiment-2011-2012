# SAA-VAA-2011-2012 Bayesian - effect of galls2011 
# and density dependence - logistic

library(rjags)
library(R2jags)
library(R2WinBUGS)
library(rethinking)
library(coda)

d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

#d = d[d$f > 0,]

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

nPlants=nrow(d)

# model
sink('adult-addition-galls-dd-logistic.txt')
cat("
model {
	# Priors
	p ~ dunif(0, 1)
	R ~ dunif(0,100)
	b ~ dunif(-10, 10)
	dd ~ dunif(-10, 10)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- ( exp(R + b * g11[i] + dd * O[i]) / (1 + exp(R + b * g11[i] + dd * O[i])) ) * O[i] # do I want a link function or not?
		O[i] ~ dbin(p, F[i])
	}
}
", fill = TRUE)
sink()


# Initial values play with this
inits = function() list(p = 0.5, R = 9, b=1, dd=0, O = ifelse(d$f>0, 1, 0))

# parms
parms = c('p', 'R', 'b', 'dd')

# data
data = list(F=d$f, nPlants=nrow(d), G=d$g12, g11=d$g11)

# MCMC controls
ni = 50000
nb = 0.1 * ni
thin = 10

# run MCMC
m4 = jags(data, inits, parms, 'adult-addition-galls-dd-logistic.txt', n.chains=3, n.iter=ni)
plot(m4)
traceplot(m4)

m3.mcmc = as.mcmc(m3)
m3.mcmcmat = as.matrix(m3.mcmc)

ts.plot(m3$BUGSoutput$sims.list$p)
ts.plot(m3$BUGSoutput$sims.list$R)
ts.plot(m3$BUGSoutput$sims.list$b)
ts.plot(m3$BUGSoutput$sims.list$dd)
ts.plot(m3$BUGSoutput$sims.list$deviance)

traceplot(m3)

xyplot(m3.mcmc)
densityplot(m3.mcmc)

plot(pgr ~ f, data=d)
newf = seq(0,12,length=80)
pred.pgr = sapply(newf,
	function(z) mean( exp(m4.mcmcmat[,'R'] + m4.mcmcmat[,'b'] * mean(d$g11) + m4.mcmcmat[,'dd'] * z) / (1 + exp(m4.mcmcmat[,'R'] + m4.mcmcmat[,'b'] * mean(d$g11) + m4.mcmcmat[,'dd'] * z)) * m4.mcmcmat[,'p'] ))
ci.pgr = sapply(newf,
	function(z) HPDI( exp(m4.mcmcmat[,'R'] + m4.mcmcmat[,'b'] * mean(d$g11) + m4.mcmcmat[,'dd'] * z) / (1 + exp(m4.mcmcmat[,'R'] + m4.mcmcmat[,'b'] * mean(d$g11) + m4.mcmcmat[,'dd'] * z)) * m4.mcmcmat[,'p'] ))
lines(newf, pred.pgr)
lines(newf, ci.pgr[1,], lty=2)
lines(newf, ci.pgr[2,], lty=2)


summary(orj.gdd)
gelman.diag(orj.gdd)
autocorr.plot(orj.gdd)
