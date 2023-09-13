set.seed(1)
library("Rcpp")
library(CholWishart)
library(mvtnorm)
logminus = function(x, y){
  if(x >= y){
    return(x + log(1.0 - exp(y - x)));
  }else{
    return(NA);
  }
}
sourceCpp("logsum.cpp")
single.radius.indicator <- T
optimal.radius = 0.1 # this should be the radius that is believed to be optimal based on some pilot runs
radius.seq = c(0.09, 0.1, 0.25, 0.5, 0.75, 1)

#####################################
### parameters of case background ###
#####################################
y = as.matrix(read.csv("simulation_study_2/observations.csv")[, 2:3])
sample.mean = apply(y, 2, mean)
S2 = (m - 1) * cov(y)
## hyper-parameters
mu0 = c(0, 0)
k0 = 0.01
nu0 = 3
lambda0 = matrix(data = c(1, 0.7, 0.7, 1), nrow = 2)
## parameters
mu = c(0, 0)
sigma = matrix(data = c(1, 0.7, 0.7, 1), nrow = 2)
## analytical solution computation
m = 200
d = 2
p = 5 # dimension of parameter space
sample.size = 10^4
num = nu0 + m
km = k0 + m
lambdam = lambda0 + S2 + (k0 * m / (k0 + m)) * (mu0 - sample.mean) %*% t(mu0 - sample.mean)
log.true.c = -0.5 * m * d * log(pi) + lmvgamma(x = 0.5 * num, p = d) - lmvgamma(x = 0.5 * nu0, p = d) + 0.5 * nu0 * log(det(lambda0)) - 0.5 * num * log(det(lambdam)) + 0.5 * d * log(k0) - 0.5 * d * log(km)
mum = (k0 * mu0 + m * sample.mean)/(k0 + m) # See Bayesian Data Analysis by Andrew Gelman, pg 73

#########################
### functions for IDR ###
#########################
## IDR-structural functions
norm.vec = function(x){
  sqrt(sum(x^2))
}

log.ball.volume = function(p, radius){ # see wikipedia "volume of an n-ball"
  if (p == 0) {volume = 0}
  if (p == 1){volume = log(2) + log(radius)}
  if (p > 1 &  p <= 10){
    volume = log(2) + log(pi) - log(p) + 2*log(radius) + log.ball.volume(p-2, radius = radius)
  } else {
    volume = (p / 2) * log(pi) + p * log(radius) - lgamma(p / 2 + 1)
  }
  return(volume)
}

h = function(xi, xi.0, p, r){
  # note that if the distance is too small, then the point is inside the circle and h won't be calculated for this point.
  distance = norm.vec(xi - xi.0)
  output = (1 - (r^p)/(distance^p))^(1/p) * (xi - xi.0)
  return(output)
}

estimated.rmse = function(log.k, log.c.hat, n, linflated.densities, original.densities){
  # this function uses (3.21) in "IDR for marginal likelihood in Bayesian phylogenetics"
  density.ratio = exp(linflated.densities - original.densities)
  rmse = (1 / sqrt(n)) * (exp(log.c.hat - log.k)) * sd(density.ratio)
  return(rmse)
}

#################
## case-specific functions (log posterior density function)
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
  lpst = lprior(mu, mu0, k0, nu0, sigma, lambda0) + llike(y, mu, sigma)
  return(lpst)
}

###################
## density functions
qg = function(y, xi, mu0, k0, nu0, lambda0){ # using the notation on page 246 of Wang (2020).
  mu = xi[4:5]
  alpha = xi[2]
  beta1 = xi[1]
  beta2 = xi[3]
  
  ljacobian = 1.5 * (beta1 + beta2) + log(2) + alpha - 2 * log(1 + exp(alpha))
  rhosigma1sigma2 = ((exp(alpha)-1)/(1+exp(alpha))) * exp(0.5 * (beta1 + beta2))
  conumat = matrix(c(exp(beta1), rhosigma1sigma2, rhosigma1sigma2, exp(beta2)), nrow = 2, byrow = T)
  lpst = lpost(y, mu, conumat, mu0, k0, nu0, lambda0) + ljacobian
  return(lpst)
}

qr = function(y, xi, xi.0, p, radius, mu0, k0, nu0, lambda0, qg.top){ # p=5 as the dimension of parameter space is 5
  # determine whether the particle is inside the unit circle from the centre (sample mean)
  #! note that qr(xi=, xi.0=, ..., lambda0=) need to be specified when used, otherwise there will be errors
  if (radius <= 0) {stop("Inflation radius must be positive")}
  distance = norm.vec(xi - xi.0)
  if (distance <= radius){ # the density value will be the density at which the inflation occurs
    est.dens = qg.top
  } 
  if (distance > radius){ 
    h.xi.diff = h(xi = xi, xi.0 = xi.0, p = p, r = radius)
    est.dens = qg(y, xi = (xi.0 + h.xi.diff), mu0, k0, nu0, lambda0)
  }
  return(est.dens)
}

#####################
## estimation functions
c.hat.estimation.single.radius = function(y, p, n = sample.size,
                                          radius, 
                                          original.densities, 
                                          target.sample, 
                                          mu0, k0, nu0, lambda0, qg.top, xi.0, single.radius.indicator){ # = the wrap.linflatedpost function
  linflated.densities = apply(X = target.sample, 1,
                              FUN = qr, y = y, p = p, mu0 = mu0, k0 = k0, nu0 = nu0, lambda0 = lambda0,
                              qg.top = qg.top, radius = radius,
                              xi.0 = xi.0)
  # marginal likelihood
  log.k = log.ball.volume(p = p, radius = radius) + qg.top # or + qg(y, xi = xi.0, mu0, k0, nu0, lambda0)
  log.dens.ratio = logplusvec(linflated.densities - original.densities)
  log.denominator = logminus(log.dens.ratio - log(n), 0)
  log.c.hat = log.k - log.denominator
  
  if (single.radius.indicator == T){ # only returning the estimation result is sufficient
    return(log.c.hat)
  } else {
    # RMSE
    est.rmse = estimated.rmse(log.k = log.k, 
                              log.c.hat = log.c.hat, 
                              n = n, 
                              linflated.densities = linflated.densities, 
                              original.densities = original.densities)
    result.vec = c(radius, log.c.hat, est.rmse)
    names(result.vec) <- c("radius", "log c hat", "estimated RMSE")
    return(result.vec)
  }
}

c.hat.estimation.auto = function(y, p, radius.seq, 
                                 target.sample, original.densities,
                                 mu0, k0, nu0, lambda0, qg.top,
                                 n = sample.size, xi.0){
  group.results = t(sapply(radius.seq,
                           c.hat.estimation.single.radius,
                           y = y, p = p, qg.top = qg.top,
                           target.sample = target.sample, 
                           original.densities = original.densities, 
                           n = n, xi.0 = xi.0,
                           mu0 = mu0, k0 = k0, nu0 = nu0, lambda0 = lambda0,
                           single.radius.estimator = F))
  min.rmse.index = which.min(group.results[, 3])
  final.estimate = group.results[min.rmse.index, ]
  return(final.estimate)
}

########################
## other functions
correlation.transform = function(cov12, sigmas1, sigmas2){
  rho = cov12 / sqrt(sigmas1 * sigmas2)
  alpha = log((rho + 1) / (1 - rho))
  return(alpha)
}

###################
### simulations ###
###################
ptm <- proc.time()
if (single.radius.indicator == T){ # only using the optimal radius for estimations
  simulation.results = rep(NA, 1000)
  for (i in 1:1000){
    sigma.post = rInvWishart(n = sample.size, df = num, Sigma = lambdam)
    mupost.dispersion = sigma.post / km
    mu.post = apply(mupost.dispersion, 3, rmvnorm, n = 1, mean = mum)
    transformed.sigma.post.vec = cbind(log(sigma.post[1,1,]),
                                       correlation.transform(sigma.post[1,2,], sigma.post[1,1,], sigma.post[2,2,]),
                                       log(sigma.post[2,2,]))
    transformed.post.sample.vec = cbind(transformed.sigma.post.vec, t(mu.post))
    transformed.post.sample.mean = apply(transformed.post.sample.vec, 2, mean)
    original.densities = apply(transformed.post.sample.vec, 1, qg, y = y, mu0 = mu0, k0 = k0, nu0 = nu0, lambda0 = lambda0)
    qg.top = qg(y, xi = transformed.post.sample.mean, mu0, k0, nu0, lambda0)
    simulation.results = c.hat.estimation.single.radius(y = y, p = p, radius = optimal.radius,
                                                        original.densities = original.densities, 
                                                        target.sample = transformed.post.sample.vec,
                                                        mu0 = mu0, k0 = k0, nu0 = nu0, lambda0 = lambda0, qg.top = qg.top,
                                                        xi.0 = transformed.post.sample.mean, single.radius.indicator = T)
    print(paste0("iteration", i))
    #print(simulation.results[i, ])
  }
} else { # using multiple radius for estimations and then choosing the result using estimated RMSE
  simulation.results = matrix(nrow = 1000, ncol = 3)
  for (i in 1:1000){
    sigma.post = rInvWishart(n = sample.size, df = num, Sigma = lambdam)
    mupost.dispersion = sigma.post / km
    mu.post = apply(mupost.dispersion, 3, rmvnorm, n = 1, mean = mum)
    transformed.sigma.post.vec = cbind(log(sigma.post[1,1,]),
                                       correlation.transform(sigma.post[1,2,], sigma.post[1,1,], sigma.post[2,2,]),
                                       log(sigma.post[2,2,]))
    transformed.post.sample.vec = cbind(transformed.sigma.post.vec, t(mu.post))
    transformed.post.sample.mean = apply(transformed.post.sample.vec, 2, mean)
    original.densities = apply(transformed.post.sample.vec, 1, qg, y = y, mu0 = mu0, k0 = k0, nu0 = nu0, lambda0 = lambda0)
    qg.top = qg(y, xi = transformed.post.sample.mean, mu0, k0, nu0, lambda0)
    simulation.results[i, ] = c.hat.estimation.auto(y = y, p = p, radius.seq = radius.seq,
                                                    target.sample = transformed.post.sample.vec, original.densities = original.densities,
                                                    mu0 = mu0, k0 = k0, nu0 = nu0, lambda0 = lambda0, qg.top = qg.top,
                                                    xi.0 = transformed.post.sample.mean)
    
    print(paste0("iteration", i))
    #print(simulation.results[i, ])
  }
}
(proc.time() - ptm) / 1000

#################################################
### post processing of the simulation results ###
#################################################
if (single.radius.indicator == T){
  mean(simulation.results)
  sd(simulation.results) / sqrt(1000)
  square.diff = (simulation.results - rep(log.true.c, 1000))^2
  (1 / 1000) * sum(square.diff)
  
  plot(density(simulation.results), xlab = "estimates of marginal likelihood",
       main = "The Inflated Density Ratio Estimates, radius = 0.1")
  abline(v = log.true.c, col = "red")
  legend("topright", legend = c("density of estimates", "true value"), col = c("black", "red"), lwd = 2, lty = 1)
} else {
  apply(simulation.results[, 2:3], 2, mean)
  apply(simulation.results[, 2:3], 2, sd) / sqrt(1000)
  simulation.results[, 1]
  square.diff = (simulation.results[, 2] - rep(log.true.c, 1000))^2
  (1 / 1000) * sum(square.diff)
  
  
  mean(simulation.results[, 1] == 0.1)
  which(simulation.results[, 1] == 0.1)
  plot(density(simulation.results[, 2]), xlab = "estimates of marginal likelihood",
       main = "The Inflated Density Ratio Estimates")
  abline(v = log.true.c, col = "red")
  legend("topright", legend = c("density of estimates", "true value"), col = c("black", "red"), lwd = 2, lty = 1)
  lines(density(simulation.results[which(simulation.results[, 1] == 0.1), 2]), col = "blue")
  lines(density(simulation.results[which(simulation.results[, 1] != 0.1), 2]), col = "purple")
  legend("topright", legend = c("density of all estimates", "true value", "radius 0.1 estimates", "radius non-0.1 estimates"), 
         col = c("black", "red", "blue", "purple"), lwd = 2, lty = 1)
}
