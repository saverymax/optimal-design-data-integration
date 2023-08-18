
r <- seq(0, 10, 0.1)
r
x <- exp(-10*((r-7)/5)^2)
min(x)
max(x)
min(0.98*x)
max(0.98*x)
min(0.98*x+.01)
max(0.98*x+.01)
plot(r, x)
plot(r, 0.98*x)
plot(r, 0.98*x +0.01)
plot(r, qnorm(x))
plot(r, qnorm(0.98*x))
plot(r, qnorm(0.98*x + 0.01))

