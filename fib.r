N <- 20;
f <- 0:2;
for (i in 3:N) f <- c(f, sum(f[i-(0:2)]))
f

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