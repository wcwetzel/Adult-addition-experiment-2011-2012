# Adult addition experiment 2011-2012
# 23 Jul 2012

library(bbmle)
library(lattice)
library(rethinking)

d=read.csv("/Users/will/Documents/DATA/2012 DATA/FIELD/SAA-VAA-2011-2012.csv")

d$pgr = d$g12 / d$f
d$pgrl = log(d$pgr+0.01)

d$fpres = as.numeric(d$f>0)

dSAA = d[d$site=='snarl',]
dVAA = d[d$site=='val',]


# plot overview for SAA and VAA combined
par(mfrow=c(2,2))
plot(g12 ~ f, data=d)
lines(loess.smooth(d$f, d$g12))
plot(g12 ~ g11, data=d)
lines(loess.smooth(d$g11, d$g12))
plot(pgr ~ f, data=d)
lines(loess.smooth(d$pgr, d$f))
plot(pgr ~ g11, data=d)
lines(loess.smooth(d$g11, d$g12))

# SAA alone
par(mfrow=c(2,2))
plot(g12 ~ f, data=dSAA)
lines(loess.smooth(dSAA$f, dSAA$g12))
plot(g12 ~ g11, data=dSAA)
lines(loess.smooth(dSAA$g11, dSAA$g12))
plot(pgr ~ f, data=dSAA)
lines(loess.smooth(dSAA$pgr, dSAA$f))
plot(pgr ~ g11, data=dSAA)
lines(loess.smooth(dSAA$g11, dSAA$g12))

# VAA alone
par(mfrow=c(2,2))
plot(g12 ~ f, data=dVAA)
lines(loess.smooth(dVAA$f, dVAA$g12))
plot(g12 ~ g11, data=dVAA)
lines(loess.smooth(dVAA$g11, dVAA$g12))
plot(pgr ~ f, data=dVAA)
lines(loess.smooth(dVAA$pgr, dVAA$f))
plot(pgr ~ g11, data=dVAA)
lines(loess.smooth(dVAA$g11, dVAA$g12))



#library(rgl)
#plot3d(d$g12 ~ d$f + d$g11, type='s', size=1)
#plot3d(d$g12 ~ d$f + d$g11, type='h', add=TRUE)


#------------- models -------------------------#

# model 0: intercept
# if we put females in the cage, then we expect mu galls
m0 = mle2(g12 ~ dnbinom(mu = mu * fpres, size=size), start=list(mu=17, size=0.1), 
	data=d)

# model 1: just galls
# if we put females in the cage, then we expect mu + b1 * g11 galls
m1 = mle2(g12 ~ dnbinom(mu = mu * fpres + fpres * b1 * g11, size=size), 
	start=list(mu = 17, b1=1, size=0.1), trace=TRUE,
	data=d)
	
# model 1l: just log(galls)
# if we put females in the cage, then we expect mu + b1 * g11 galls
# m1l = mle2(g12 ~ dnbinom(mu = mu * fpres + fpres * b1 * log(g11+0.01), size=size), 
	# start=list(mu = 17, b1=1, size=0.1), trace=TRUE,
	# data=d)

# model 2: linear females
# we expect b2 galls per female
m2 = mle2(g12 ~ dnbinom(mu = b2 * f, size=size), 
	start=list(b2=3, size=0.1), 
	data=d, trace=TRUE)

# model 3: galls and linear females, not sure this one makes sense
# if we put females in the cage, then we expect
m3 = mle2(g12 ~ dnbinom(mu = fpres * b1 * g11 + b2 * f, size=size), 
	start=list(b1=1, b2=1, size=0.1), data=d, trace=TRUE)

# model 4: pgr is a linear function of galls
# we expect (b1 + b2 * g11) galls per female
m4 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11) * f, size=size), 
	start=list(b1=1, b2=1, size=0.1), data=d, trace=TRUE)

# model 5: pgr is a linear function of galls AND females
# we expect (b1 + b2 * g11 + b3 * f) galls per female
m5 = mle2(g12 ~ dnbinom(mu =  (b1 + b2 * g11 + b3 * f) * f, size=size), 
	start=list(b1=1, b2=1, b3=0, size=0.1), data=d, trace=TRUE)

# model 5l: pgr is a log linear function of galls AND females
# we expect log(b1 + b2 * g11 + b3 * f) galls per female
m5l = mle2(g12 ~ dnbinom(mu =  exp(b1 + b2 * g11 + b3 * f) * f, size=size), 
	start=list(b1=0, b2=0, b3=0, size=0.1), data=d, trace=TRUE)

# model 6: pgr is a rational function of females, linear function of galls
# we expect (b1 + b2 * g11) / (1 + b3 * f) galls per female
m6 = mle2(g12 ~ dnbinom(mu =  ( (b1 + b2 * g11) /  (1 + b3 * f) )* f, size=size), 
	start=list(b1=4.277, b2=0.91, b3=0.471, size=0.88), data=d, trace=TRUE,
	method='Nelder-Mead')

# model 6.5: pgr is a rational function of females, linear function of galls
# we expect (b1 + b2 * g11) / (1 + b3 * f) galls per female
m6.5 = mle2(g12 ~ dnbinom(mu =  ( (b1 + b2 * g11) /  (b4 + b3 * f) )* f, size=size), 
	start=list(b1=2.388, b2=0.705, b3=0.179, b4=1, size=1.013), data=d, trace=TRUE)

# model 6.7
m6.7 = mle2(g12 ~ dnbinom(mu =  ( (a + R * g11) /  (c + g11 + b * f) ) * f, size=size), 
	start=list(a = 2, R=12, b=1, c=1, size=1.013), data=d, trace=TRUE)


# model 7: pgr is a rational function of females, rational function of galls
# we expect (b1 + b2 * g11) / (1 + b3 * f + b4 * g11) galls per female
m7 = mle2(g12 ~ dnbinom(mu =  ( (b1 + b2 * g11) /  (1 + b3 * f + b4 * g11) )* f, 
	size=size), start=list(b1=14.785, b2=16.803, b3=3.825, b4=1.89, size=0.9849), 
	data=d, trace=TRUE, method='Nelder-Mead')

# model 8: pgr is a neg exp function of females, linear function of galls
# we expect exp(r * (1 - f / k) + b1 * g11) galls per female
m8 = mle2(g12 ~ dnbinom(mu =  exp( r * (1 - f / k) + b1 * g11 ) * f, size=size), 
	start=list(r=1, k=15, b1=0, size=1.013), data=d, trace=TRUE)

# model 9: pgr is a michaelis-menten function of galls
# we expect ( a * (g11+1) ) / ( b + (g11+1) ) galls per female
m9 = mle2(g12 ~ dnbinom(mu =  (( b1 * (g11+1) ) / ( b2 + (g11+1) )) * f, size=size), 
	start=list(b1=1, b2=1, size=1.013), data=d, trace=TRUE)

# model 10: as proposed in QE presentation
m10 = mle2(g12 ~ dnbinom(mu =  f * exp(r * (1 - f / k) + b1 * g11), size=size), 
	start=list(r=6, k=25, b1=0, size=1.013), data=d, trace=TRUE)

# model 10ng: as proposed in QE presentation - no galls
m10ng = mle2(g12 ~ dnbinom(mu =  f * exp(r * (1 - f / k)), size=size), 
	start=list(r=6, k=25, size=1.013), data=d, trace=TRUE)

# model 69: neg exp females, rational galls
m69 = mle2(g12 ~ dnbinom(mu =  ( exp(r + b * f) * (g11+1) /  (c + g11) ) * f, size=size), 
	start=list(r = 2, b=-0.5, c=1, size=1.013), data=d, trace=TRUE)


AICctab(m0, m1, m2, m3, m4, m5, m5l, m6, m6.5, m6.7, m8, m9, m10, m10ng, m69, nobs=60, weights=TRUE)
AICtab(m0, m1, m2, m3, m4, m5, m5l, m6, weights=TRUE)

anova(m0, m1)
anova(m0, m2)
anova(m4, m1)
anova(m5, m4)
anova(m1, m6)
anova(m6, m6.7)

d0 = lm(pgrl ~ 1, data=d)
d1 = lm(pgrl ~ f, data=d)
d2 = lm(pgrl ~ g11, data=d)
d3 = lm(pgrl ~ f + g11, data=d)
d4 = lm(pgrl ~ f * g11, data=d)


AICtab(d0,d1,d2,d3,d4)
anova(d2, d0)
anova(d1, d0)


#----------------- plot model results ----------------------------#

###### first g12 ~ f
new.f = seq(0, 15)


post.m69 = sample.naive.posterior(m69, n=100000)
mu.m69 = sapply(new.f,
	function(z) mean( (exp(post.m69[,'r'] + post.m69[,'b'] * z) * (1+1) /  (post.m69[,'c'] + 1) ) * z) ) 
ci.m69 = sapply(new.f,
	function(z) HPDI( (exp(post.m69[,'r'] + post.m69[,'b'] * z) * (1+1) /  (post.m69[,'c'] + 1) ) * z) ) 
plot(d$g12 ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1)
lines(new.f, mu.m69, col='blue')
lines(new.f, ci.m69[1,], lty=2, col='blue')
lines(new.f, ci.m69[2,], lty=2, col='blue')

mu.m69pgr = sapply(new.f,
	function(z) mean( (exp(post.m69[,'r'] + post.m69[,'b'] * z) * (1+1) /  (post.m69[,'c'] + 1) ) ) ) 
ci.m69pgr = sapply(new.f,
	function(z) HPDI( (exp(post.m69[,'r'] + post.m69[,'b'] * z) * (1+1) /  (post.m69[,'c'] + 1) ) ) ) 
plot(d$pgr ~ d$f, ylab='PGR', xlab='Females in 2011', pch=20, las=1)
lines(new.f, mu.m69pgr, col='blue')
lines(new.f, ci.m69pgr[1,], lty=2, col='blue')
lines(new.f, ci.m69pgr[2,], lty=2, col='blue')

mu.m69pgr = sapply(new.f,
	function(z) mean( (exp(post.m69[,'r'] + post.m69[,'b'] * z) * (1+10) /  (post.m69[,'c'] + 10) ) ) ) 
ci.m69pgr = sapply(new.f,
	function(z) HPDI( (exp(post.m69[,'r'] + post.m69[,'b'] * z) * (1+10) /  (post.m69[,'c'] + 10) ) ) ) 
lines(new.f, mu.m69pgr, col='orange')
lines(new.f, ci.m69pgr[1,], lty=2, col='orange')
lines(new.f, ci.m69pgr[2,], lty=2, col='orange')


### model m10ng ###
post.m10ng = sample.naive.posterior(m10ng, n=100000)
mu.m10ng = sapply(new.f,
	function(z) median( z * exp(post.m10ng[,'r'] * (1 - z / post.m10ng[,'k'])) ) )

ci.m10ng = sapply(new.f,
	function(z) HPDI( z * exp(post.m10ng[,'r'] * (1 - z / post.m10ng[,'k'])) ) )

plot(d$g12 ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1)
lines(new.f, mu.m10ng, col='blue')
lines(new.f, ci.m10ng[1,], lty=2, col='blue')
lines(new.f, ci.m10ng[2,], lty=2, col='blue')

### model m10ng PGR ###
mu.m10ngpgr = sapply(new.f,
	function(z) median( exp(post.m10ng[,'r'] * (1 - z / post.m10ng[,'k'])) ) )

ci.m10ngpgr = sapply(new.f,
	function(z) HPDI( exp(post.m10ng[,'r'] * (1 - z / post.m10ng[,'k'])) ) )

plot(d$pgr ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1)
lines(new.f, mu.m10ngpgr, col='green')
lines(new.f, ci.m10ngpgr[1,], lty=2, col='green')
lines(new.f, ci.m10ngpgr[2,], lty=2, col='green')



### model m10 ###
post.m10 = sample.naive.posterior(m10, n=100000)
mu.m10 = sapply(new.f,
	function(z)( z * exp(post.m10[,'r'] * (1 - z / post.m10[,'k'])) ) )

ci.m10ng = sapply(new.f,
	function(z) HPDI( z * exp(post.m10ng[,'r'] * (1 - z / post.m10ng[,'k'])) ) )

plot(d$g12 ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1)
lines(new.f, mu.m10ng, col='blue')
lines(new.f, ci.m10ng[1,], lty=2, col='blue')
lines(new.f, ci.m10ng[2,], lty=2, col='blue')



### model 6.7 ###
post.m67 = sample.naive.posterior(m6.7, n=100000)
mu.m67.gall0 = sapply(new.f,
	function(z) mean(((post.m67[,'a'] + post.m67[,'R'] * 0 ) / (post.m67[,'c'] + 0 + post.m67[,'b'] * z)) * z))
ci.m67.gall0 = sapply(new.f,
	function(z) HPDI(((post.m67[,'a'] + post.m67[,'R'] * 0 ) / (post.m67[,'c'] + 0 + post.m67[,'b'] * z)) * z))

mu.m67.gall5 = sapply(new.f,
	function(z) mean(((post.m67[,'a'] + post.m67[,'R'] * 5 ) / (post.m67[,'c'] + 5 + post.m67[,'b'] * z)) * z))
ci.m67.gall5 = sapply(new.f,
	function(z) HPDI(((post.m67[,'a'] + post.m67[,'R'] * 5 ) / (post.m67[,'c'] + 5 + post.m67[,'b'] * z)) * z))

mu.m67.gall10 = sapply(new.f,
	function(z) mean(((post.m67[,'a'] + post.m67[,'R'] * 10 ) / (post.m67[,'c'] + 10 + post.m67[,'b'] * z)) * z))
ci.m67.gall10 = sapply(new.f,
	function(z) HPDI(((post.m67[,'a'] + post.m67[,'R'] * 10 ) / (post.m67[,'c'] + 10 + post.m67[,'b'] * z)) * z))


plot(d$g12 ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1,
	main='Predictions from model with rational females, linear males')

# lines for predictions by galliness
lines(new.f, mu.m67.gall0, col='blue')
lines(new.f, ci.m67.gall0[1,], lty=2, col='blue')
lines(new.f, ci.m67.gall0[2,], lty=2, col='blue')

lines(new.f, mu.m67.gall5, col='orange')
lines(new.f, ci.m67.gall5[1,], lty=2, col='orange')
lines(new.f, ci.m67.gall5[2,], lty=2, col='orange')

lines(new.f, mu.m67.gall10, col='red')
lines(new.f, ci.m67.gall10[1,], lty=2, col='red')
lines(new.f, ci.m67.gall10[2,], lty=2, col='red')

legend(-0.65, 101, legend = c('0 galls in 2011', '5 galls in 2011', '15 galls in 2011'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')





### model 6 ###
post.m6 = sample.naive.posterior(m6, n=100000)
mu.m6.gall0 = sapply(new.f,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * 0 ) / (1 + post.m6[,'b3'] * z)) * z))
ci.m6.gall0 = sapply(new.f,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * 0 ) / (1 + post.m6[,'b3'] * z)) * z))

mu.m6.gall5 = sapply(new.f,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * 5 ) / (1 + post.m6[,'b3'] * z)) * z))
ci.m6.gall5 = sapply(new.f,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * 5 ) / (1 + post.m6[,'b3'] * z)) * z))

mu.m6.gall10 = sapply(new.f,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * 10 ) / (1 + post.m6[,'b3'] * z)) * z))
ci.m6.gall10 = sapply(new.f,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * 10 ) / (1 + post.m6[,'b3'] * z)) * z))


plot(d$g12 ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1,
	main='Predictions from model with rational females, linear males')

# lines for predictions by galliness
lines(new.f, mu.m6.gall0, col='blue')
lines(new.f, ci.m6.gall0[1,], lty=2, col='blue')
lines(new.f, ci.m6.gall0[2,], lty=2, col='blue')

lines(new.f, mu.m6.gall5, col='orange')
lines(new.f, ci.m6.gall5[1,], lty=2, col='orange')
lines(new.f, ci.m6.gall5[2,], lty=2, col='orange')

lines(new.f, mu.m6.gall10, col='red')
lines(new.f, ci.m6.gall10[1,], lty=2, col='red')
lines(new.f, ci.m6.gall10[2,], lty=2, col='red')

legend(-0.65, 101, legend = c('0 galls in 2011', '5 galls in 2011', '15 galls in 2011'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')



### model 5 ###
post.m5 = sample.naive.posterior(m5, n=100000)

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
plot(d$g12 ~ d$f, ylab='Galls in 2012', xlab='Females in 2011', pch=20, las=1,
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
plot(g12 ~ f, data=d, ylab='Galls in 2012', xlab='Females in 2011', 
	pch=20, las=1, main='Predictions from model w/o density-dependence')

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



##### second g12 ~ g11 #####
new.g = seq(0, 40)

### model 1 no PGR
post.m1 = sample.naive.posterior(m1, n=1000000)


mu.m1 = sapply(new.g,
	function(z) mean( post.m1[,'mu'] + post.m1[,'b1'] * z ))
ci.m1 = sapply(new.g,
	function(z) HPDI( post.m1[,'mu'] + post.m1[,'b1'] * z ))
plot(d$g12 ~ d$g11, ylab='Galls in 2012', xlab='Galls in 2011', pch=20, las=1, subset=d$fpres==1)
lines(new.g, mu.m1)
lines(new.g, ci.m1[1,], lty=2)
lines(new.g, ci.m1[2,], lty=2)

### model 1 PGR

mu.m1 = sapply(new.g,
	function(z) mean( post.m1[,'mu'] + post.m1[,'b1'] * z ))
ci.m1 = sapply(new.g,
	function(z) HPDI( post.m1[,'mu'] + post.m1[,'b1'] * z ))
plot(d$g12 ~ d$g11, ylab='Galls in 2012', xlab='Galls in 2011', pch=20, las=1, subset=d$fpres==1)
lines(new.g, mu.m1)
lines(new.g, ci.m1[1,], lty=2)
lines(new.g, ci.m1[2,], lty=2)


## model 6
mu.m6.f1 = sapply(new.g,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 1)) * 1))
ci.m6.f1 = sapply(new.g,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 1)) * 1))

mu.m6.f5 = sapply(new.g,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 5)) * 5))
ci.m6.f5 = sapply(new.g,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 5)) * 5))

mu.m6.f12 = sapply(new.g,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 12)) * 12))
ci.m6.f12 = sapply(new.g,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 12)) * 12))


plot(d$g12 ~ d$g11, ylab='Galls in 2012', xlab='Galls in 2011', pch=20, las=1,
	main='Predictions from model with rational females, linear males', data=d)

# lines for predictions by galliness
lines(new.g, mu.m6.f1, col='blue')
lines(new.g, ci.m6.f1[1,], lty=2, col='blue')
lines(new.g, ci.m6.f1[2,], lty=2, col='blue')

lines(new.g, mu.m6.f5, col='orange')
lines(new.g, ci.m6.f5[1,], lty=2, col='orange')
lines(new.g, ci.m6.f5[2,], lty=2, col='orange')

lines(new.g, mu.m6.f12, col='red')
lines(new.g, ci.m6.f12[1,], lty=2, col='red')
lines(new.g, ci.m6.f12[2,], lty=2, col='red')

legend(-0.65, 101, legend = c('1 female', '5 females', '12 females'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')




## model 5 with density dependence

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
	main='Predictions from density-dependent model', data=d)

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
	main='Predictions from model w/o density-dependence', data=d)

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


##### plotting pgr #####
### model 6.7 ###
mu.m67.f1pgr = sapply(new.g,
	function(z) mean(((post.m67[,'a'] + post.m67[,'R'] * z ) / (post.m67[,'c'] + z + post.m67[,'b'] * 1)) ))
ci.m67.f1pgr = sapply(new.g,
	function(z) HPDI(((post.m67[,'a'] + post.m67[,'R'] * z ) / (post.m67[,'c'] + z + post.m67[,'b'] * 1)) ))

mu.m67.f5pgr = sapply(new.g,
	function(z) mean(((post.m67[,'a'] + post.m67[,'R'] * z ) / (post.m67[,'c'] + z + post.m67[,'b'] * 5)) ))
ci.m67.f5pgr = sapply(new.g,
	function(z) HPDI(((post.m67[,'a'] + post.m67[,'R'] * z ) / (post.m67[,'c'] + z + post.m67[,'b'] * 5)) ))

mu.m67.f12pgr = sapply(new.g,
	function(z) mean(((post.m67[,'a'] + post.m67[,'R'] * z ) / (post.m67[,'c'] + z + post.m67[,'b'] * 12)) ))
ci.m67.f12pgr = sapply(new.g,
	function(z) HPDI(((post.m67[,'a'] + post.m67[,'R'] * z ) / (post.m67[,'c'] + z + post.m67[,'b'] * 12)) ))


plot(d$pgr ~ d$g11, ylab='PGR', xlab='Galls in 2011', pch=20, las=1,
	main='')

# lines for predictions by galliness
lines(new.g, mu.m67.f1pgr, col='blue')
lines(new.g, ci.m67.f1pgr[1,], lty=2, col='blue')
lines(new.g, ci.m67.f1pgr[2,], lty=2, col='blue')

lines(new.g, mu.m67.f5pgr, col='orange')
lines(new.g, ci.m67.f5pgr[1,], lty=2, col='orange')
lines(new.g, ci.m67.f5pgr[2,], lty=2, col='orange')

lines(new.g, mu.m67.f12pgr, col='red')
lines(new.g, ci.m67.f12pgr[1,], lty=2, col='red')
lines(new.g, ci.m67.f12pgr[2,], lty=2, col='red')

legend(-0.65, 101, legend = c('0 galls in 2011', '5 galls in 2011', '15 galls in 2011'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')



## model 6
mu.m6.f1.pgr = sapply(new.g,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 1)) ))
ci.m6.f1.pgr = sapply(new.g,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 1)) ))

mu.m6.f5.pgr = sapply(new.g,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 5)) ))
ci.m6.f5.pgr = sapply(new.g,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 5)) ))

mu.m6.f12.pgr = sapply(new.g,
	function(z) mean(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 12)) ))
ci.m6.f12.pgr = sapply(new.g,
	function(z) HPDI(((post.m6[,'b1'] + post.m6[,'b2'] * z ) / (1 + post.m6[,'b3'] * 12)) ))


plot(d$pgr ~ d$g11, ylab='Per capita galling rate', xlab='Galls in 2011', pch=20, las=1,
	main='Predictions from model with rational females, linear males', data=d)

# lines for predictions by galliness
lines(new.g, mu.m6.f1.pgr, col='blue')
lines(new.g, ci.m6.f1.pgr[1,], lty=2, col='blue')
lines(new.g, ci.m6.f1.pgr[2,], lty=2, col='blue')

lines(new.g, mu.m6.f5.pgr, col='orange')
lines(new.g, ci.m6.f5.pgr[1,], lty=2, col='orange')
lines(new.g, ci.m6.f5.pgr[2,], lty=2, col='orange')

lines(new.g, mu.m6.f12.pgr, col='red')
lines(new.g, ci.m6.f12.pgr[1,], lty=2, col='red')
lines(new.g, ci.m6.f12.pgr[2,], lty=2, col='red')

legend(-0.65, 20, legend = c('1 female', '5 females', '12 females'),
	col=c('blue', 'orange', 'red'), lty=1, cex=0.8, bty='n')


