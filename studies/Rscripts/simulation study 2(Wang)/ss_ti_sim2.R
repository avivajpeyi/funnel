set.seed(1)
# To install powModSel R-package
devtools::install_github('pmat747/powModSel');
library(powModSel);
library(CholWishart)
library(mvtnorm)
library(ggplot2)
library(patchwork)

#####################################
### parameters of case background ###
#####################################
y = as.matrix(read.csv("simulation_study_2/observations.csv")[, 2:3])
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

m = 200
d = 2
num = nu0 + m
km = k0 + m
lambdam = lambda0 + S2 + (k0 * m / (k0 + m)) * (mu0 - sample.mean) %*% t(mu0 - sample.mean)
log.true.c = -0.5 * m * d * log(pi) + lmvgamma(x = 0.5 * num, p = d) - lmvgamma(x = 0.5 * nu0, p = d) + 0.5 * nu0 * log(det(lambda0)) - 0.5 * num * log(det(lambdam)) + 0.5 * d * log(k0) - 0.5 * d * log(km)

##########################
### functions for SSTI ###
##########################
# power likelihood function
quadratic.form = function(yj, mu, sigma){ # this is the quadratic form inside the exponential for the likelihood
  # not in use
  t(yj - mu) %*% solve(sigma) %*% (yj - mu)
}

llike_pow = function(sigma, mu, y){
  # not in use
  m = nrow(y) # sample size
  quadratic = apply(y, 2, quadratic.form, mu = mu, sigma = sigma)
  result = -0.5*sum(quadratic) - (m / 2) * determinant(sigma, logarithm = T)
  return(result)
}

llike_pow_evaluation = function(sample.vector, y){
  # this function is prepared to evaluate the likelihood density value of a point where the input is a vector (ccccmm)
  sample.covariance = matrix(sample.vector[1:4], nrow = 2, ncol = 2, byrow = T)
  sample.mean = sample.vector[5:6]
  ll = sum(dmvnorm(x = y, mean = sample.mean, sigma = sample.covariance, log = T))
  return(ll)
}

##############
# power posterior function
pow_sampling = function(beta, lambda0, kappa0, nu0, mu0, m, ybar, S, n){
  # this function samples from power posteriors of power beta
  mu_bm = (kappa0 * mu0 + beta * m * ybar) / (beta * m + kappa0)
  kappa_bm = kappa0 + beta * m
  lambda_bm = lambda0 + beta*S + (beta*m*kappa0)/(beta*m+kappa0)*(ybar-mu0) %*% t(ybar-mu0)
  sigma.sample = rInvWishart(n = n,
                             df = nu0 + beta * m,
                             Sigma = lambda_bm)
  mu.dispersion = sigma.sample / kappa_bm
  mu.sample = apply(mu.dispersion, 3, rmvnorm, n = 1, mean = mu_bm)
  
  sigma.result = matrix(c(sigma.sample), nrow = 4, byrow = F)
  mu.result = matrix(c(mu.sample), nrow = 2, byrow = F)
  result = rbind(sigma.result, mu.result)# the output will be a matrix where each row is an observation (ccccmm)
  return(result)
}

######################
# function to prepare for the matrix of log likelihoods and temperatures
f = function(k, n, y, lambda0, kappa0, nu0, mu0){ # using sufficient statistic to summarize y
  
  # f returns a vector with the log-likelihoods and temperatures 
  # of n power posterior samples
  
  # k: number of temperatures
  # n: number of samples per temperature
  # lambda0: prior parameter of sigma
  # kappa0: prior sample size
  # nu0: prior degrees of freedom
  # mu0: prior location parameter of mu
  
  m = nrow(y)
  ybar = apply(y, 2, mean) #observation mean
  S = (m - 1) * cov(y) #sample covariance multiplied by sample size. Note that cov() uses m-1 as denom.
  
  beta = qbeta(seq(0,1,length = k), 0.3, 1)
  
  pow_samples.list = lapply(beta, pow_sampling, lambda0 = lambda0, kappa0 = kappa0,
                            nu0 = nu0, mu0 = mu0, m = m, ybar = ybar, S = S, n = n)
  pow_samples = matrix(unlist(pow_samples.list), ncol = 6, byrow = T)
  ls = apply(pow_samples, 1, llike_pow_evaluation, y = y)
  ls = as.data.frame(cbind(ls, rep(beta, each = n)))
  colnames(ls) = c("logL", "invTemp")
  
  return(ls)
  
}

#r = f(k=200, n=50, y=y, lambda0, kappa0 = k0, nu0 = nu0, mu0); # log-likelihoods and temperatures

#ss(r); # Stepping-stone sampling estimate
#ti(r) # Thermodynamic integration estimate

###################
### simulations ###
###################
ssti.results = matrix(nrow = 1000, ncol = 2)
ptm <- proc.time()
for (i in 1:1000){
  r = f(k=1000, n=10, y=y, lambda0, kappa0 = k0, nu0 = nu0, mu0)
  ssti.results[i, 1] = ss(r)
  ssti.results[i, 2] = ti(r)
  print(paste0("iteration", i))
}
(proc.time() - ptm) / 1000

#################################################
### post processing of the simulation results ###
#################################################
(sim.mean = apply(ssti.results, 2, mean))
(sim.se = apply(ssti.results, 2, sd) / sqrt(1000))
abs(sim.mean - rep(log.true.c, 2)) / sim.se
square.diff1 = (ssti.results[, 1] - rep(log.true.c, 1000))^2
(1 / 1000) * sum(square.diff1)
square.diff2 = (ssti.results[, 2] - rep(log.true.c, 1000))^2
(1 / 1000) * sum(square.diff2)
write.csv(x = as.data.frame(ssti.results), file = "ssti.csv")

colours <- c("true value" = "red")
ss.result.density = ggplot(data = as.data.frame(ssti.results[, 1]), aes(x = ssti.results[, 1])) +
  geom_density(alpha=.2, fill="lightblue") +
  xlim(-505.65, -504.5) + 
  geom_vline(aes(xintercept = log.true.c, color = "true value")) +
  labs(title = "The Stepping-Stone Estimates") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ss.result.box = ggplot(data = as.data.frame(ssti.results[, 1]), aes(x = ssti.results[, 1])) +
  geom_boxplot() +
  xlim(-505.65, -504.5) + 
  geom_vline(xintercept = log.true.c, color = "red") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ss.result.density / ss.result.box +
  plot_layout(heights = c(3, 1)) +
  labs(x = "estimates of marginal likelihood",
       color = "legend")

ti.result.density = ggplot(data = as.data.frame(ssti.results[, 2]), aes(x = ssti.results[, 2])) +
  geom_density(alpha=.2, fill="lightblue") +
  xlim(-505.65, -504.5) + 
  geom_vline(aes(xintercept = log.true.c, color = "true value")) +
  labs(title = "The Thermodynamic Integration Estimates") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ti.result.box = ggplot(data = as.data.frame(ssti.results[, 2]), aes(x = ssti.results[, 2])) +
  geom_boxplot() +
  xlim(-505.65, -504.5) + 
  geom_vline(xintercept = log.true.c, color = "red") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ti.result.density / ti.result.box +
  plot_layout(heights = c(3, 1)) +
  labs(x = "estimates of marginal likelihood",
       color = "legend")

#plot(density(ssti.results[, 1]), xlab = "estimates of marginal likelihood",
#     main = "The SSTI Estimates")
#lines(density(ssti.results[, 2]), col = "blue")
#abline(v = log.true.c, col = "red")
#legend("topright", legend = c("SS estimates", "TI estimates", "true value"),
#       col = c("black", "blue", "red"), lwd = 2, lty = 1)