rm(list = ls())
# We source the appropriate group level here
source("pmwg/variants/factor.R")
source("plots.R")

# Define the joint likelihood
ll_joint <- function(x,data,sample=FALSE) {
  x_cog <- x[grep("cog", names(x))]
  x_MRI <- x[grep("MRI", names(x))]
  names(x_cog) <- gsub("cog_", "", names(x_cog))
  names(x_MRI) <- gsub("MRI_", "", names(x_MRI))
  ll_cog <- ll_cog(x=x_cog,data=data[[2]][[1]])
  ll_MRI <- ll_MRI(x=x_MRI,data=data[[3]][[1]])
  if(ll_cog > -1e+10 && ll_MRI > -1e+10){
    out = sum(ll_cog,ll_MRI)
  } else {
    out = -1e10
  }
  return(out)
}

# The likelihood of the cognitive model
ll_cog <- function(x,data){
  # Simple DDM
  out1 <- rtdists::ddiffusion(data[data$cond == 1, ], a = exp(x["a"]), t0 = exp(x["t0"]),
                              v = x[1])
  out2 <- rtdists::ddiffusion(data[data$cond == 2, ], a = exp(x["a"]), t0 = exp(x["t0"]),
                              v = x[3])
  out <- c(out1, out2)
  bad <- (out<1e-10)|(!is.finite(out))
  out[bad] <- 1e-10
  return(sum(log(pmax(out, 1e-10))))
}

# The likelihood of the neural activity
ll_MRI <- function(x, data){
  # The first column is the total BOLD signal
  y <- data[,1]
  # The next set of columns are the individual regions activity
  X <- as.matrix(data[,2:(ncol(data)-1)])
  # The first set of parameters are the regressors one for each region
  betas <- x[1:(length(x) -1)]
  # The last parameter is the standard deviation
  sigma <- exp(x[length(x)])
  return(max(sum(dnorm(y, mean=X %*% betas, sd = sigma, log = T)), -10000))
}

# Some constants
n_subjects <-50
n_trials_cog <- 200   
n_timepoints_mri <- 200
n_factors <- 2

pars <- c("cog_v1", "cog_a", "cog_v2", "cog_t0", "MRI_c1", "MRI_c2", "MRI_sd")
group_means <- c(.5, log(1.2), 1, log(0.25), -.3, .5, log(.3)) #a, t0 and neural sd are estimated on log-scale to have support on the natural scale
group_factors <- matrix(c(.4, 0,
                          -.4, .4,
                          .35, -0.05,
                          .05, 0.1,
                          .45, -.05,
                          .4, -.2,
                          .1, .25), ncol = n_factors, byrow = T)

rownames(group_factors) <- pars

group_resid_variance <- c(.05, .1, .15, .2, .15, .1, .05)

group_cov <- group_factors %*% t(group_factors) + diag(group_resid_variance)
# Inspect our simulated factors
corrplot::corrplot(group_factors, is.corr = F, col.lim = c(-1,1))

# Inspect our simulated correlations matrix
corrplot::corrplot(cov2cor(group_cov), is.corr = T)

# seems reasonable, now make some random effects
random_effects <- mvtnorm::rmvnorm(n_subjects, group_means, group_cov)
colnames(random_effects) <- pars


behav_data <- data.frame()
neural_data <- data.frame()

for(sub in 1:n_subjects){
  pars_cog <- random_effects[sub,1:4]
  pars_mri <- random_effects[sub,5:7]
  # Generate 2 difficulty conditions
  behav1 <- rtdists::rdiffusion(n_trials_cog/2, a = exp(pars_cog[2]), t0 = exp(pars_cog[4]),
                                v = pars_cog[1])
  behav1$cond <- 1
  behav2 <- rtdists::rdiffusion(n_trials_cog/2, a = exp(pars_cog[2]), t0 = exp(pars_cog[4]),
                                v = pars_cog[3])
  behav2$cond <- 2
  behav <- rbind(behav1, behav2)
  behav$subject <- sub
  behav_data <- rbind(behav_data, behav)
  # Just generate some random neural time signal
  time_signal <-  mvtnorm::rmvnorm(n_timepoints_mri, mean = rep(0, 2))
  BOLD <- rnorm(n_timepoints_mri, mean = 
                  time_signal %*% pars_mri[1:2], sd = exp(pars_mri[3]))
  neural <- data.frame(y = BOLD, X = time_signal, subject = sub)
  neural_data <- rbind(neural_data, neural)
}


# This makes one big data frame of them, which can be used in the sampler
all_data <- data.frame(subject = 1:n_subjects)
all_data['cog'] <- I(list(split(behav_data, f = behav_data$subject)))
all_data['MRI'] <- I(list(split(neural_data, f = neural_data$subject)))

##Factor analysis
# Create the Particle Metropolis within Gibbs sampler object ------------------
# This is with standard diagonal constraint
# Here we don't have to supply a constraint matrix, since diagonal constraint is the default. But if we would it would look like this:
constraintMat <- matrix(1, nrow = length(pars), ncol = n_factors)
constraintMat[upper.tri(constraintMat, diag = T)] <- 0 

sampler <- pmwgs(
  data = all_data,
  pars = pars,
  ll_func = ll_joint,
  n_factors = n_factors
) 

# This is how we set up the sampler with minimal constraint
constraintMat <- matrix(1, nrow = length(pars), ncol = n_factors)
constraintMat[upper.tri(constraintMat, diag = F)] <- 0

sampler <- pmwgs(
  data = all_data,
  pars = pars,
  ll_func = ll_joint,
  n_factors = n_factors,
  constraintMat = constraintMat
)

# Lastly we can also apply the additional constraint that for example t0 is not loading on factor 1.
# We again use minimal constraint
constraintMat <- matrix(1, nrow = length(pars), ncol = n_factors)
constraintMat[upper.tri(constraintMat, diag = F)] <- 0
constraintMat[pars == "cog_t0",1] <- 0

sampler <- pmwgs(
  data = all_data,
  pars = pars,
  ll_func = ll_joint,
  n_factors = n_factors,
  constraintMat = constraintMat
)

#We'll use this constraint form to actually sample from this time.


n_cores <- 12

# start the sampler ---------------------------------------------------------
sampler <- init(sampler, n_cores = n_cores) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 1000,  n_cores = n_cores)
adapted <- run_stage(burned, stage = "adapt",iter = 1500, n_cores = n_cores)
sampled <- run_stage(adapted, stage = "sample",iter = 2000,  n_cores = n_cores)

# since we did minimal constraint we have to check the bimodality of the samples
# Here I only take samples from the chain in the sample stage
hist(sampled$samples$theta_lambda[1,1,sampled$samples$stage == "sample"])
hist(sampled$samples$theta_lambda[2,2,sampled$samples$stage == "sample"])
# Potentially some small bimodality in the second factor, since the loadings are flipped once they're negative, and 0 seems to be cut-off
# Safer to run with standard constraint

# Check chains
idx <- sampled$samples$stage == "sample"
matplot(t(sampled$samples$theta_mu[,idx]), type = "l")
# Factor 1, we want theta_lambda (lambda_true in the paper) and not lambda_untransf
matplot(t(sampled$samples$theta_lambda[,1,idx]), type = "l")
# Factor 2, these look autocorrelated, could thin them
matplot(t(sampled$samples$theta_lambda[,2,idx]), type = "l")

# Mean recovery, this is more dependent on the DDM being sloppy rather than the sampling method, but looks very good
plot(group_means, rowMeans(sampled$samples$theta_mu[,idx]), xlab = "simulated", ylab = "recovered")
abline(coef = c(0,1))

# Loadings recovery
mean_loadings_recovered <- apply(sampled$samples$theta_lambda[,,idx], 1:2, mean)

# Looks like good recovery
lambda_recovery_plot(sampled, group_factors)

# This is how to plot the loadings with corrplot, depending on your monitor set up you will need to play with the x and y of the legend in the function
plot_relations(sampled, cred = T) # will also plot 95% credible intervals

# We can do the same for the implied correlation matrix
plot_relations(sampled, cred = T, do_corr = T)

