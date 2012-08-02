### SAA-VAA-2011-2012-SIMPLE
# 1 Aug 2012
# simple analysis for ESA2012
# all nbinom
# all except m0 have model per capita galling rate


library(bbmle)
library(lattice)
library(rethinking)

d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)


#---------------- models -----------------#


# model 0: intercept
# if we put females in the cage, then we expect mu galls
m0 = mle2(g12 ~ dnbinom(mu = mu * fpres, size=size), start=list(mu=17, size=0.1), 
	data=d)

# model 2: linear females
# we expect b2 galls per female
m2 = mle2(g12 ~ dnbinom(mu = b2 * f, size=size), 
	start=list(b2=3, size=0.1), 
	data=d)

# model 4: pgr is a linear function of galls
# we expect (b1 + b2 * g11) galls per female
m4 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11) * f, size=size), 
	start=list(b1=1, b2=1, size=0.1), data=d)

# model 5: pgr is a linear function of galls AND females
# we expect (b1 + b2 * g11 + b3 * f) galls per female
m5 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11 + b3 * f) * f, size=size), 
	start=list(b1=1, b2=1, b3=0, size=0.1), data=d)

# model 6: pgr is a linear function of galls AND non-linear function of females
# we expect (b1 + b2 * g11 + exp(b3 * f)) galls per female
m6 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11 + exp(b3 * f)) * f, size=size), 
	start=list(b1=1, b2=1, b3=0, size=0.1), data=d)

# model 7: pgr is a nonlinear function of galls AND linear function of females
# we expect (b1 + b2 * g11 + b3 * f) galls per female
m7 = mle2(g12 ~ dnbinom(mu =  (b1 + (b2 * g11)/(b4 + g11) + b3 * f) * f, size=size), 
	start=list(b1=1, b2=1, b3=0, b4=1, size=0.1), data=d)


# model 8: pgr is a linear function of galls AND females AND interaction
# we expect (b1 + b2 * g11 + b3 * f + b4 * g11 * f) galls per female
m8 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11 + b3 * f + b4 * g11 * f) * f, size=size), 
	start=list(b1=1, b2=1, b3=0, b4=0, size=0.1), data=d)


AICctab(m0, m2, m4, m5, m6, m7, m8, nobs=60, weights=TRUE)



#---------------- sample naive post ---------------------------#


post.m4 = sample.naive.posterior(m4, n=100000)

post.m5 = sample.naive.posterior(m5, n=100000)


#---------------- plotting with x = females --------------------#

newf = seq(0,12, length=40)

##### model 5 #######

### g12 (not PGR) ##### PRESENT THIS WITH TWO LINES
plot(g12 ~ f, data=d)

# 0 galls
mu.m5 = sapply(newf, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * 0 + post.m5[,'b3'] * z) * z )  )
ci.m5 = sapply(newf, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 0 + post.m5[,'b3'] * z) * z )  )
lines(mu.m5 ~ newf, col='blue')
lines(ci.m5[1,] ~ newf, lty=2, col='blue')
lines(ci.m5[2,] ~ newf, lty=2, col='blue')

# 15 galls
mu.m5 = sapply(newf, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * 15 + post.m5[,'b3'] * z) * z )  )
ci.m5 = sapply(newf, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 15 + post.m5[,'b3'] * z) * z )  )
lines(mu.m5 ~ newf)
lines(ci.m5[1,] ~ newf, lty=2)
lines(ci.m5[2,] ~ newf, lty=2)


### PGR 
plot(pgr ~ f, data=d)

# 0 galls
mu.m5 = sapply(newf, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * 0 + post.m5[,'b3'] * z) )  )
ci.m5 = sapply(newf, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 0 + post.m5[,'b3'] * z) )  )
lines(mu.m5 ~ newf)
lines(ci.m5[1,] ~ newf, lty=2)
lines(ci.m5[2,] ~ newf, lty=2)

# 15 galls
mu.m5 = sapply(newf, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * 15 + post.m5[,'b3'] * z) )  )
ci.m5 = sapply(newf, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 15 + post.m5[,'b3'] * z) )  )
lines(mu.m5 ~ newf)
lines(ci.m5[1,] ~ newf, lty=2)
lines(ci.m5[2,] ~ newf, lty=2)



#------------ plotting with x = g11 ---------------------#


newg = seq(0,40, length=80)

##### model 5 #######

### g12 (not PGR)
plot(g12 ~ g11, data=d, subset=d$fpres==1)

# 1 female
mu.m5 = sapply(newg, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 1) * 1 )  )
ci.m5 = sapply(newg, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 1) * 1 )  )
lines(mu.m5 ~ newg)
lines(ci.m5[1,] ~ newg, lty=2)
lines(ci.m5[2,] ~ newg, lty=2)

# 6 females
mu.m5 = sapply(newg, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 6) * 6 )  )
ci.m5 = sapply(newg, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 6) * 6 )  )
lines(mu.m5 ~ newg)
lines(ci.m5[1,] ~ newg, lty=2)
lines(ci.m5[2,] ~ newg, lty=2)


### PGR #### PRESENT THIS WITH JUST ONE LINE
plot(pgr ~ g11, data=d, subset=d$fpres==1)

# 1 female
mu.m5 = sapply(newg, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 1) )  )
ci.m5 = sapply(newg, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 1) )  )
lines(mu.m5 ~ newg)
lines(ci.m5[1,] ~ newg, lty=2)
lines(ci.m5[2,] ~ newg, lty=2)

# 6 females
mu.m5 = sapply(newg, function(z)
	mean( (post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 6) )  )
ci.m5 = sapply(newg, function(z)
	HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 6) )  )
lines(mu.m5 ~ newg)
lines(ci.m5[1,] ~ newg, lty=2)
lines(ci.m5[2,] ~ newg, lty=2)

