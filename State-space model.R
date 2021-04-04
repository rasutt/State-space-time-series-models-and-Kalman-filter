# State-space model

# White noise
n_t <- 20
sigma <- 1
e <- rnorm(n_t, sd = sigma)

# Random walk
tau <- 1
a <- rnorm(n_t, sd = tau)
mu <- cumsum(a)

# Random walk plus white noise
y <- mu + e
plot(y, type = 'l', xlab = "n_time", ylab = "Value",
     main = "Random walk plus white noise")

# Quarterly random walks
mu_q <- as.vector(apply(matrix(a, nrow = 4), 2, cumsum))

# Quarterly random walks plus noise
y_q <- mu_q + e
plot(y_q, type = 'l', xlab = "n_time", ylab = "Value",
     main = "Quarterly random walks plus white noise")
grid(nx = n_t / 4)

# White noise with evolving variance
h <- cumsum(a)
y_v <- exp(h / 2) * e
plot(y_v, type = 'l', xlab = "n_time", ylab = "Value",
     main = "White noise with evolving variance")

