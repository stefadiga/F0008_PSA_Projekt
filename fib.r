N <- 20;
f <- 0:2;
for (i in 3:N) f <- c(f, sum(f[i-(0:2)]))
f
