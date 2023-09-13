library(mvtnorm)
library("Rcpp")
library(CholWishart)
library(ncdf4)
#library(haven)
#library(LaplacesDemon)
# Wang example
## hyper-parameters
mu0 = c(0, 0)
k0 = 0.01
nu0 = 3
lambda0 = matrix(data = c(1, 0.7, 0.7, 1), nrow = 2)

## parameters
mu = c(0, 0)
sigma = matrix(data = c(1, 0.7, 0.7, 1), nrow = 2)

## analytical solution computation (1)
m = 200
d = 2
num = nu0 + m
km = k0 + m

## sampling of observed data
set.seed(1) # this is the set of data we actually "observe"
y = rmvnorm(n = 200, mean = mu, sigma = sigma)
write.csv(y, file = "observations.csv")
sample.mean = apply(y, 2, mean)
(rho = cor(y)[1, 2])
S2 = (m - 1) * cov(y)

## analytical solution computation (2)
lambdam = lambda0 + S2 + (k0 * m / (k0 + m)) * (mu0 - sample.mean) %*% t(mu0 - sample.mean)
log.true.c = -0.5 * m * d * log(pi) + lmvgamma(x = 0.5 * num, p = d) - lmvgamma(x = 0.5 * nu0, p = d) + 0.5 * nu0 * log(det(lambda0)) - 0.5 * num * log(det(lambdam)) + 0.5 * d * log(k0) - 0.5 * d * log(km)
mum = (k0 * mu0 + m * sample.mean)/(k0 + m) # See Bayesian Data Analysis by Andrew Gelman, pg 73


###################### NO NEED TO RUN THE FOLLOWING CODES
## sampling from the posterior distribution
sigma.post = rInvWishart(n = 10^7, df = num, Sigma = lambdam)
mupost.dispersion = sigma.post / km
### estimating the amount of time needed for sampling
mu.post0 = matrix(nrow = 2, ncol = 10^5)
a = Sys.time()
for (j in 1:10^5){
  mu.post0[, j] = rmvnorm(1, mum, mupost.dispersion[,,j])
  if (j %% 1000 == 0){
    print(paste0("iteration", j))
  }
}
b = Sys.time()    
paste0(round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 3), " Seconds")
### 40.102 seconds, so the amount of time needed using apply() should be no longer than 100*40.102 seconds = 4010.2 seconds = 66.8 minutes = 1.11 hours
mu.post = apply(mupost.dispersion, 3, rmvnorm, n = 1, mean = mum) # in fact: 1 hour 19 minutes

# posterior sample storage
sigma.post.mat = matrix(sigma.post, ncol = 2, byrow = T)
mu_post_vec = c(mu.post)
post.sample.mat = cbind(sigma.post.mat, mu_post_vec)
post.sample.frame = as.data.frame(post.sample.mat)
write.csv(post.sample.frame, file = "posterior_sample.csv") # 6 min

# netCDF (not used)
sigma.post.vec = cbind(sigma.post[1,1,], sigma.post[1,2,], sigma.post[2,2,])
post.sample.vec = cbind(sigma.post.vec, mu.post)
d1 = ncdim_def("covdim1", "unit1", as.double(-1:1), unlim = T)
varsigma21 = ncvar_def("sigma1", "", d1, prec = "double")
netfile = nc_create("posterior_sample0.nc", list(varsigma21))


## MH sampling (not used)
burnin <- 2000
Niter <- 12000
thin <- 1
mu.sample = cbind(rep(NA, Niter), rep(NA, Niter)) # one column for mu1, one column for mu2
sigma.sample = matrix(nrow = 2 * Niter, ncol = 2)
keep = seq(burnin, Niter, thin)
initial.mu = c(1, 1)
initial.sigma = diag(nrow = 2, ncol = 2)

# the first iteration
mu.star = c(rnorm(1), rnorm(1))
sigma21.star = exp(rnorm(1))
sigma22.star = exp(rnorm(1))
alpha.star = rnorm(1)
rho.star = correlation.backtransform(alpha.star)
covariance.star = rho.star * sqrt(sigma21.star) * sqrt(sigma22.star)
sigma.sample= matrix(c(sigma21.star, covariance.star, covariance.star, sigma22.star), nrow = 2, ncol = 2)
log.postdensratio = lpost(y, mu.star, sigma.star, mu0, k0, lambda0) - lpost(y, initial.mu, initial.sigma, mu0, k0, lambda0)
if (log(runif(1)) < log.postdensratio) { # ACCEPT
  mu.sample[1, ] <- mu.star
  sigma.sample[1:2, ] <- sigma.star
} else { # REJECT
  mu.sample[1, ] <- initial.mu
  sigma.sample[1:2, ] <- initial.sigma
}

## iterations
for (i in 2:Niter) {
  
  
  # step 1: propose a new value:
  mu.star = c(rnorm(1), rnorm(1))
  sigma21.star = exp(rnorm(1))
  sigma22.star = exp(rnorm(1))
  alpha.star = rnorm(1)
  rho.star = correlation.backtransform(alpha.star)
  covariance.star = rho.star * sqrt(sigma21.star) * sqrt(sigma22.star)
  sigma.star = matrix(c(sigma21.star, covariance.star, covariance.star, sigma22.star), nrow = 2, ncol = 2)
  
  # step 2: evaluate the posterior kernels
  log.postdensratio = lpost(y, mu.star, sigma.star, mu0, k0, lambda0) - lpost(y, mu.sample[(i-1), ], sigma.sample[(2*i-3):(2*i-2), ], mu0, k0, lambda0)
  
  # step 3: accept or reject
  if (log(runif(1)) < log.postdensratio) { # ACCEPT
    mu.sample[i, ] <- mu.star
    sigma.sample[(2*i-1):(2*i), ] <- sigma.star
  } else { # REJECT
    mu.sample[i, ] <- mu.sample[(i - 1), ]
    sigma.sample[(2*i-1):(2*i), ] <- sigma.sample[(2*i-3):(2*i-2), ]
  }
}

## try block Gibbs - sample theta using the conditional distribution given in the DrIDR and sample xi using adaptive MH (but could be conjugate) or using the conjugate that the posterior of sigma is again IW.
## or try JAGS or use the likelihood function in the TISS file
## try to derive the analytical solution of the power posteriors in TI and SS.




## outcome
mcmc.outcome = y[keep]
mcmc.size = length(mcmc.outcome)

plot(density(mcmc.outcome))
