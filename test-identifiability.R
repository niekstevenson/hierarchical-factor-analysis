#Script for testing the recovery of factor models with custom 

rm(list = ls())
source("pmwg/variants/factor.R")

ll <- function(pars, data){
  n_pars <- length(pars)
  return(sum(dmvnorm(as.matrix(data[,2:(1 + n_pars)]), pars, diag(1, nrow = n_pars), log = T)))
}

n_pars <- 15 # number of parameters
n_subjects <-20 # number of subjects
n_trials <- 200 #number trials per subject per condition
n_factors <- 2 # number of factors

parNames <- paste0("V", 1:n_pars)

lambda <- matrix(rep(c(.5,0,.8,.4,-.3), n_factors*3), nrow = n_pars, ncol = n_factors, byrow = TRUE)

sigma <- diag(rep(c(0.15, 0.25, 0.3, 0.2, 0.1), 3))
eta <- rmvnorm(n_subjects, sigma = diag(rep(1, n_factors)))
epsilon <- rmvnorm(n_subjects, sigma = sigma)
y <- sweep(eta %*% t(lambda) + epsilon, 2, -seq(-2,2, length.out = n_pars))

all_data <- data.frame()
for(j in 1:n_subjects){
  sub_data <- data.frame(subject = rep(j, n_trials))
  sub_data[,2:(n_pars + 1)] <- data.frame(rmvnorm(n_trials, y[j,], diag(1, nrow = n_pars)))
  all_data <- rbind(all_data, sub_data)
}

##Factor 
# Create the Particle Metropolis within Gibbs sampler object ------------------
sampler <- pmwgs(
  data = all_data,
  pars = parNames,
  ll_func = ll,
  n_factors = n_factors
)

n_cores <- 8

# start the sampler ---------------------------------------------------------
# sampler <- init(sampler, n_cores = n_cores) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 1500, n_cores = n_cores)
adapted <- run_stage(burned, stage = "adapt",iter = 5000,n_cores = n_cores)
sampled <- run_stage(adapted, stage = "sample",iter = 1500,n_cores = n_cores)

