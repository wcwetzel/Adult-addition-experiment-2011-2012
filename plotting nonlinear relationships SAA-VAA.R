#### 3-D plot of SAA-VAA predictions

library(lattice)

g11 = rep(seq(0, 40, length=41), 12)
f = sort(rep(seq(1,12, length=12), 41))


R = 12
c = 1
a = 1
b = 10
b2 =0.5
r = 3


pgr = a + R * g11 / (c + g11 + b * f)

wireframe(pgr ~ g11 + f, scales=list(arrows = FALSE, draw=TRUE, tick.number=8))

(1.4687 + 6.779 * 10) / (-0.857 + 10 + 1.037 * 1)


m10pgr = function() exp(coef(m10)['r'] * (1 - f / coef(m10)['k']) + coef(m10)['b1'] * g11)

wireframe(m10pgr() ~ g11 + f)

pgr = (a + exp(r - b2 * f) * g11) / (c + g11)

wireframe(pgr ~ g11 + f, scales=list(arrows = FALSE, draw=TRUE, tick.number=8))



