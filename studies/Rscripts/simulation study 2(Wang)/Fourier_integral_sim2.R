set.seed(1)
library(matrixStats)
library(CholWishart)
library(mvtnorm)
library(ggplot2)
library(patchwork)

y = as.matrix(read.csv("observations.csv")[, 2:3])
m = 200
R = 20
S2 = (m - 1) * cov(y)
sample.mean = apply(y, 2, mean)
## hyper-parameters
mu0 = c(0, 0)
k0 = 0.01
nu0 = 3
lambda0 = matrix(data = c(1, 0.7, 0.7, 1), nrow = 2)
## parameters
mu = c(0, 0)
sigma = matrix(data = c(1, 0.7, 0.7, 1), nrow = 2)
## analytical solution computation
d = 2
p = 5
sample.size = 10^4
vm = nu0 + m
km = k0 + m
lambdam = lambda0 + S2 + (k0 * m / (k0 + m)) * (mu0 - sample.mean) %*% t(mu0 - sample.mean)
log.true.c = -0.5 * m * d * log(pi) + lmvgamma(x = 0.5 * vm, p = d) - lmvgamma(x = 0.5 * nu0, p = d) + 0.5 * nu0 * log(det(lambda0)) - 0.5 * vm * log(det(lambdam)) + 0.5 * d * log(k0) - 0.5 * d * log(km)
mum = (k0 * mu0 + m * sample.mean)/(k0 + m) # See Bayesian Data Analysis by Andrew Gelman, pg 73

# joint density of data and parameter (posterior kernel)
correlation.transform = function(cov12, sigmas1, sigmas2){
  rho = cov12 / sqrt(sigmas1 * sigmas2)
  alpha = log((rho + 1) / (1 - rho))
  return(alpha)
}

lprior = function(mu, mu0, k0, nu0, sigma, lambda0){
  lprior.sigma = dInvWishart(x = sigma, df = nu0, Sigma = lambda0, log = T)
  lprior.mu = dmvnorm(x = mu, mean = mu0, sigma = (1 / k0) * sigma, log = T)
  lp = lprior.sigma + lprior.mu
  return(lp)
}

llike = function(y, mu, sigma){
  ll = sum(dmvnorm(x = y, mean = mu, sigma = sigma, log = T))
  return(ll)
}

lpost = function(y, mu, sigma, mu0, k0, nu0, lambda0){
  lpst = lprior(mu, mu0, k0, nu0, sigma, lambda0) + llike(y = y, mu = mu, sigma = sigma)
  return(lpst)
}
# computation of posterior kernel value at mu and sigma
ljacobian <- function(alpha, beta1, beta2) {
  1.5 * (beta1 + beta2) + log(2) + alpha - 2 * log(1 + exp(alpha))
}
lq = lpost(y, mu, sigma, mu0, k0, nu0, lambda0) + ljacobian(correlation.transform(0.7, 1, 1), 0, 0)

# we choose to evaluate the posterior density at the point specified by mu and sigma (1, 0.7, 1, 0, 0)
x = c(0, correlation.transform(0.7, 1, 1), 0, mu)
x.mat = matrix(rep(x, sample.size), nrow = sample.size, byrow = T)

simulation.results = rep(NA, 1000)
ptm <- proc.time()
for (i in 1:1000){
  sigma.post = rInvWishart(n = sample.size, df = vm, Sigma = lambdam)
  mupost.dispersion = sigma.post / km
  mu.post = apply(mupost.dispersion, 3, rmvnorm, n = 1, mean = mum)
  transformed.sigma.post.vec = cbind(log(sigma.post[1,1,]),
                                     correlation.transform(sigma.post[1,2,], sigma.post[1,1,], sigma.post[2,2,]),
                                     log(sigma.post[2,2,]))
  transformed.post.sample.vec = cbind(transformed.sigma.post.vec, t(mu.post))
  sample.eval.diff = transformed.post.sample.vec - x.mat
  a = rowProds(sin(R * sample.eval.diff) / sample.eval.diff)
  logpostdens = log(sum(a)) - (log(sample.size) + p * log(pi))
  simulation.results[i] <- lq - logpostdens
  if (is.nan(simulation.results[i]) == T){
    stop("NaN produced in the last iteration")
  }
  print(paste0("iteration", i))
}
proc.time() - ptm
(sim.mean = mean(simulation.results))
(sim.se = sd(simulation.results) / sqrt(1000))
abs(sim.mean - log.true.c) / sim.se
square.diff = (simulation.results - rep(log.true.c, 1000))^2
(1 / 1000) * sum(square.diff)

fi.result.density = ggplot(data = as.data.frame(simulation.results), aes(x = simulation.results)) +
  geom_density(alpha=.2, fill="lightblue") +
  geom_vline(aes(xintercept = log.true.c, color = "true value")) +
  xlim(-505.65, -503.8) +
  labs(title = "The Fourier Integral Estimates") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )

fi.result.box = ggplot(data = as.data.frame(simulation.results), aes(x = simulation.results)) +
  geom_boxplot() +
  geom_vline(xintercept = log.true.c, color = "red") +
  xlim(-505.65, -503.8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pdf("pics2/FI.png")
fi.result.density / fi.result.box +
  plot_layout(heights = c(3, 1)) +
  labs(x = "estimates of marginal likelihood",
       color = "legend")
dev.off()

write.csv(simulation.results, "transfi_sim2.csv")

plot(density(simulation.results), xlab = "estimates of marginal likelihood",
     main = "The Fourier Integral Estimates")
abline(v = log.true.c, col = "red")
legend("topright", legend = c("density of estimates", "true value"), col = c("black", "red"), lwd = 2, lty = 1)

