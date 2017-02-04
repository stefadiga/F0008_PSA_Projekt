N <- 30;
f <- 0:1;
for (i in 2:N) f <- c(f, sum(f[i-(0:1)]))
f
