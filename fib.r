
library(changepoint)
set.seed(10)
m.data=c(rnorm(100,0,1),rnorm(100,1,1),rnorm(100,0,1),rnorm(100,0.2,1))
ts.plot(m.data,xlab="Index")
m.pelt=cpt.mean(m.data,method="PELT")
plot(m.pelt,type="l",cpt.col="blue",xlab="Index",cpt.width=4)
cpts(m.pelt)
m.binseg=cpt.mean(m.data,method="BinSeg")
plot(m.binseg,type="l",xlab="Index",cpt.width=4)
cpts(m.binseg)
m.pm=cpt.mean(m.data,penalty="Manual",pen.value="1.5*log(n)",method="PELT")
plot(m.pm,type="l",cpt.col="blue",xlab="Index",cpt.width=4)
cpts(m.pm)

#bayesian change point
library(bcp)
plot(bcp(m.data))

#structural change
library(strucchange)
ocus.n <- efp(m.data ~ 1, type = "OLS-CUSUM")
plot(ocus.n)

fs.n <- Fstats(m.data ~ 1)
plot(fs.n)

bp.n <- breakpoints(m.data ~ 1)
plot(bp.n)
fm0.n <- lm(m.data ~ 1)
coef(fm0.n)
c<-2

n.fac <- breakfactor(bp.n)
fm1.n <- lm(m.data ~ n.fac - 1)
coef(fm1.n)

ts.plot(m.data,xlab="Index")
lines(fitted(fm0.n))
lines(fitted(fm1.n))
