data=read.table("test.txt") 

lcases= log(data[,2])
data=cbind(data,lcases)
colnames(data)=c("income","cases","CrCards","lcases")
data

## Fit the data with poisson linear model with offset=lcases

log.fit=glm(CrCards~income+offset(lcases),family=poisson,data=data)
summary(log.fit)

## Fitted Values under poisson linear model

fitted(log.fit)


library(changepoint)
#set.seed(10)
#m.data=c(rnorm(100,0,1),rnorm(100,1,1),rnorm(100,0,1),rnorm(100,0.2,1))
ts.plot(data[,3],xlab="Index")
m.pelt=cpt.mean(data[,3],method="PELT")
plot(m.pelt,type="l",cpt.col="blue",xlab="Index",cpt.width=4)
cpts(m.pelt)
m.binseg=cpt.mean(data[,3],method="BinSeg")
plot(m.binseg,type="l",xlab="Index",cpt.width=4)
cpts(m.binseg)
m.pm=cpt.mean(data[,3],penalty="Manual",pen.value="1.5*log(n)",method="PELT")
plot(m.pm,type="l",cpt.col="blue",xlab="Index",cpt.width=4)
cpts(m.pm)

#bayesian change point
library(bcp)
plot(bcp(data[,3]))

#structural change
library(strucchange)
ocus.n <- efp(data[,3] ~ 1, type = "OLS-CUSUM")
plot(ocus.n)

fs.n <- Fstats(data[,3] ~ 1)
plot(fs.n)

bp.n <- breakpoints(data[,3] ~ 1)
plot(bp.n)
fm0.n <- lm(data[,3] ~ 1)
coef(fm0.n)


n.fac <- breakfactor(bp.n)
fm1.n <- lm(data[,3] ~ n.fac - 1)
coef(fm1.n)

ts.plot(data[,3],xlab="Index")
lines(fitted(fm0.n))
lines(fitted(fm1.n))
