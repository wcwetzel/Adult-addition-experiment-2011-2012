# prelim look at SAA
# 4 Jun 2012

library(bbmle)
library(lattice)
library(rgl)
library(rethinking)


d=read.csv("/Users/will/Documents/DATA/2011 DATA/field/SAA_GALLS_2011.csv")

f = d$addition
g11 = d$Galls
g12 = d$galls2012

fpres = as.numeric(f>0)
pcg = g12 / f

hist(pcg)
par(mfrow=c(2,2))
plot(g12 ~ f)
lines(loess.smooth(f, g12))
plot(g12 ~ g11)
lines(loess.smooth(g12, g11))
plot(pcg ~ f)
lines(loess.smooth(f, pcg))
plot(pcg ~ g11)
lines(loess.smooth(g11, pcg))



cloud(g12 ~ f * g11)

plot3d(g12 ~ f + g11, type='s', size=1)
plot3d(g12 ~ f + g11, type='h', add=TRUE)

#------------- models -------------------------#

# model 0: intercept
# if we put females in the cage, then we expect mu galls
m0 = mle2(g12 ~ dnbinom(mu = mu * fpres, size=size), start=list(mu=17, size=0.1), 
	data=data.frame(g12))

# model 1: just galls
# if we put females in the cage, then we expect mu galls + b1 * galls2011
m1 = mle2(g12 ~ dnbinom(mu = mu * fpres + fpres * b1 * g11, size=size), 
	start=list(mu = 17, b1=1, size=0.1), trace=TRUE,
	data=data.frame(g12))

# model 2: linear females
# we expect b2 galls per female
m2 = mle2(g12 ~ dnbinom(mu = b2 * f, size=size), 
	start=list(b2=3, size=0.1), 
	data=data.frame(g12), trace=TRUE)

# model 3: galls and linear females, not sure this one makes sense
# if we put females in the cage, then we expect
m3 = mle2(g12 ~ dnbinom(mu = fpres * b1 * g11 + b2 * f, size=size), 
	start=list(b1=1, b2=1, size=0.1), data=data.frame(g12), trace=TRUE)

# model 4: pgr is a linear function of galls
# we expect (b1 + b2 * g11) galls per female
m4 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11) * f, size=size), 
	start=list(b1=1, b2=1, size=0.1), data=data.frame(g12), trace=TRUE)

# model 5: pgr is a linear function of galls AND females
# we expect (b1 + b2 * g11 + b3 * f) galls per female
m5 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11 + b3 * f) * f, size=size), 
	start=list(b1=1, b2=1, b3=0, size=0.1), data=data.frame(g12), trace=TRUE)

# model 5l: pgr is a log linear function of galls AND females
# we expect log(b1 + b2 * g11 + b3 * f) galls per female
m5l = mle2(g12 ~ dnbinom(mu =  exp(b1 + b2 * g11 + b3 * f) * f, size=size), 
	start=list(b1=0, b2=0, b3=0, size=0.1), data=data.frame(g12), trace=TRUE)


# model 6: pgr is a nonlinear function of females, linear function of galls
# we expect (b1 + b2 * g11) / (1 + b3 * f) galls per female
m6 = mle2(g12 ~ dnbinom(mu =  ( (b1 + b2 * g11) /  (1 + b3 * f) )* f, size=size), 
	start=list(b1=2.388, b2=0.705, b3=0.179, size=1.013), data=data.frame(g12), trace=TRUE)



AICctab(m0, m1, m2, m3, m4, m5, m5l, m6, nobs=39, weights=TRUE)
AICtab(m0, m1, m2, m3, m4, m5, m5l, m6, weights=TRUE)

anova(m0, m1)
anova(m0, m2)
anova(m4, m1)
anova(m5, m4)

d0 = lm(pcg ~ 1)
d1 = lm(pcg ~ f)
d2 = lm(pcg ~ g11)
d3 = lm(pcg ~ f + g11)
d4 = lm(pcg ~ f * g11)


AICtab(d0,d1,d2,d3,d4)
anova(d2, d0)
anova(d1, d0)


#----------------- plot model results ----------------------------#

###### first g12 ~ f

### model 5 ###
post.m5 = sample.naive.posterior(m5)
new.f = seq(0, 15)

# plot g12 ~ f with m5 predictions, separate by galliness
mu.m5.gall0 = sapply(new.f,
	function(z) mean((post.m5[,'b1'] + post.m5[,'b2'] * 0 + post.m5[,'b3'] * z) * z))
ci.m5.gall0 = sapply(new.f,
	function(z) HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 0 + post.m5[,'b3'] * z) * z))

mu.m5.gall5 = sapply(new.f,
	function(z) mean((post.m5[,'b1'] + post.m5[,'b2'] * 5 + post.m5[,'b3'] * z) * z))
ci.m5.gall5 = sapply(new.f,
	function(z) HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 5 + post.m5[,'b3'] * z) * z))

mu.m5.gall15 = sapply(new.f,
	function(z) mean((post.m5[,'b1'] + post.m5[,'b2'] * 15 + post.m5[,'b3'] * z) * z))
ci.m5.gall15 = sapply(new.f,
	function(z) HPDI((post.m5[,'b1'] + post.m5[,'b2'] * 15 + post.m5[,'b3'] * z) * z))


## plot
plot(g12 ~ f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1,
	main='Predictions from density-dependent model')

# lines for predictions by galliness
lines(new.f, mu.m5.gall0, col='blue')
lines(new.f, ci.m5.gall0[1,], lty=2, col='blue')
lines(new.f, ci.m5.gall0[2,], lty=2, col='blue')

lines(new.f, mu.m5.gall5, col='orange')
lines(new.f, ci.m5.gall5[1,], lty=2, col='orange')
lines(new.f, ci.m5.gall5[2,], lty=2, col='orange')

lines(new.f, mu.m5.gall15, col='red')
lines(new.f, ci.m5.gall15[1,], lty=2, col='red')
lines(new.f, ci.m5.gall15[2,], lty=2, col='red')

legend(-0.65, 101, legend = c('0 galls in 2011', '5 galls in 2011', '15 galls in 2011'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')

### model 4 ###
post.m4 = sample.naive.posterior(m4)
new.f = seq(0, 15)

# plot g12 ~ f with m4 predictions, separate by galliness
mu.m4.gall0 = sapply(new.f,
	function(z) mean((post.m4[,'b1'] + post.m4[,'b2'] * 0) * z))
ci.m4.gall0 = sapply(new.f,
	function(z) HPDI((post.m4[,'b1'] + post.m4[,'b2'] * 0) * z))

mu.m4.gall5 = sapply(new.f,
	function(z) mean((post.m4[,'b1'] + post.m4[,'b2'] * 5) * z))
ci.m4.gall5 = sapply(new.f,
	function(z) HPDI((post.m4[,'b1'] + post.m4[,'b2'] * 5) * z))

mu.m4.gall15 = sapply(new.f,
	function(z) mean((post.m4[,'b1'] + post.m4[,'b2'] * 15) * z))
ci.m4.gall15 = sapply(new.f,
	function(z) HPDI((post.m4[,'b1'] + post.m4[,'b2'] * 15) * z))


## plot
plot(g12 ~ f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1,
	main='Predictions from model w/o density-dependence')

# lines for predictions by galliness
lines(new.f, mu.m4.gall0, col='blue')
lines(new.f, ci.m4.gall0[1,], lty=2, col='blue')
lines(new.f, ci.m4.gall0[2,], lty=2, col='blue')

lines(new.f, mu.m4.gall5, col='orange')
lines(new.f, ci.m4.gall5[1,], lty=2, col='orange')
lines(new.f, ci.m4.gall5[2,], lty=2, col='orange')

lines(new.f, mu.m4.gall15, col='red')
lines(new.f, ci.m4.gall15[1,], lty=2, col='red')
lines(new.f, ci.m4.gall15[2,], lty=2, col='red')

legend(-0.65, 101, legend = c('0 galls in 2011', '5 galls in 2011', '15 galls in 2011'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')


##### second g12 ~ g11

## model 5 with density dependence
new.g = seq(0, 40)

# plot g12 ~ f with m5 predictions, separate by galliness
mu.m5.f1 = sapply(new.g,
	function(z) mean((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 1) * 1))
ci.m5.f1 = sapply(new.g,
	function(z) HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 1) * 1))

mu.m5.f5 = sapply(new.g,
	function(z) mean((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 5) * 5))
ci.m5.f5 = sapply(new.g,
	function(z) HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 5) * 5))

mu.m5.f12 = sapply(new.g,
	function(z) mean((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 12) * 12))
ci.m5.f12 = sapply(new.g,
	function(z) HPDI((post.m5[,'b1'] + post.m5[,'b2'] * z + post.m5[,'b3'] * 12) * 12))


## plot
plot(g12 ~ g11, ylab='Galls in 2012', xlab='Galls in 2011', pch=20, las=1,
	main='Predictions from density-dependent model')

# lines for predictions by galliness
lines(new.g, mu.m5.f1, col='blue')
lines(new.g, ci.m5.f1[1,], lty=2, col='blue')
lines(new.g, ci.m5.f1[2,], lty=2, col='blue')

lines(new.g, mu.m5.f5, col='orange')
lines(new.g, ci.m5.f5[1,], lty=2, col='orange')
lines(new.g, ci.m5.f5[2,], lty=2, col='orange')

lines(new.g, mu.m5.f12, col='red')
lines(new.g, ci.m5.f12[1,], lty=2, col='red')
lines(new.g, ci.m5.f12[2,], lty=2, col='red')

legend(-2, 101, legend = c('1 female', '5 females', '12 females'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')

## model 4 with density dependence
new.g = seq(0, 40)

# plot g12 ~ f with m4 predictions, separate by galliness
mu.m4.f1 = sapply(new.g,
	function(z) mean((post.m4[,'b1'] + post.m4[,'b2'] * z) * 1))
ci.m4.f1 = sapply(new.g,
	function(z) HPDI((post.m4[,'b1'] + post.m4[,'b2'] * z) * 1))

mu.m4.f5 = sapply(new.g,
	function(z) mean((post.m4[,'b1'] + post.m4[,'b2'] * z) * 5))
ci.m4.f5 = sapply(new.g,
	function(z) HPDI((post.m4[,'b1'] + post.m4[,'b2'] * z) * 5))

mu.m4.f12 = sapply(new.g,
	function(z) mean((post.m4[,'b1'] + post.m4[,'b2'] * z) * 12))
ci.m4.f12 = sapply(new.g,
	function(z) HPDI((post.m4[,'b1'] + post.m4[,'b2'] * z) * 12))


## plot
plot(g12 ~ g11, ylab='Galls in 2012', xlab='Galls in 2011', pch=20, las=1,
	main='Predictions from model w/o density-dependence')

# lines for predictions by galliness
lines(new.g, mu.m4.f1, col='blue')
lines(new.g, ci.m4.f1[1,], lty=2, col='blue')
lines(new.g, ci.m4.f1[2,], lty=2, col='blue')

lines(new.g, mu.m4.f5, col='orange')
lines(new.g, ci.m4.f5[1,], lty=2, col='orange')
lines(new.g, ci.m4.f5[2,], lty=2, col='orange')

lines(new.g, mu.m4.f12, col='red')
lines(new.g, ci.m4.f12[1,], lty=2, col='red')
lines(new.g, ci.m4.f12[2,], lty=2, col='red')

legend(-2, 101, legend = c('1 female', '5 females', '12 females'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')


