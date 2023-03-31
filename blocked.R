#Script for testing the recovery of the blocked covariance matrix

rm(list = ls())
source("pmwg/variants/blocked.R")

library(rtdists)
library(MCMCpack)

# Define the joint ll
joint_ll <- function(x, data){
  parPreFixs <- gsub("[|].*", "", names(x))
  totalSum <- 0
  for(mod in unique(parPreFixs)){
    currentPars <- x[which(parPreFixs == mod)]
    names(currentPars) <- gsub(".*[|]", "", names(currentPars))
    modelData <- data[[mod]][[1]]
    totalSum <- totalSum + log_likelihood(currentPars, modelData, sample = F)
  }
  return(totalSum)
}

# Define the individual LBA log_likelihood, this could be any likelihood
log_likelihood=function(x,data, sample=TRUE) {
  x <- exp(x)
  bPars <- grep("b", names(x))
  bs <- x["A"]+x[bPars][data$condition]
  if (sample) { #for sampling
    out=rLBA(n=nrow(data),A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
  } else { #for calculating density
    out=dLBA(rt=data$rt,response=data$resp,A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
    bad=(out<1e-10)|(!is.finite(out))
    out[bad]=1e-10
    out=sum(log(out))
  }
  out
}


n.trials <- 80      #number trials per subject per conditions
n.subj <- 20 #number of subjects
n.cond <- 3
n.exp <- 2

allparameters <- numeric()
alldata <- list()

# Simulate a joint model with different experiments that each have different conditions that effect the threshold, but the same amount per experiment. 
# The data is fit with an lba for each experiment

for(i in 1:n.exp){
  names=c("subject","rt","resp","condition") #names of columns
  data <- data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
  names(data) <- names
  data$condition <- rep(1:n.cond,times = n.trials) #filling in condition
  data$subject <- rep(1:n.subj, each = n.trials*n.cond) #filling in subjects
  alldata[[i]] <- data
  parameter.names <- c(c("A","v1","v2","t0"), c(paste0("b", 1:n.cond)))
  parameters <- c((0.1 + .02*i)*(-1)^i + c(0.4, 1.2, .7, -2.4), seq(0.1, 0.5, length.out = n.cond))
  names(parameters) <- paste0("Mod", i, "|", parameter.names)
  allparameters <- c(allparameters, parameters)
}

covMat <- riwish(2.5*length(allparameters),diag(length(allparameters)))


# Here I set one block per experiment, but other blocks could also be conceived of course
Cov_blocked <- matrix(0, nrow = nrow(covMat), ncol = ncol(covMat))
for(i in 1:n.exp){
  idx <-((i-1)*(4+n.cond) + 1):(i*(4+n.cond))
  Cov_blocked[idx, idx] <- covMat[idx, idx]
}
corrplot::corrplot(cov2cor(Cov_blocked))


subj_random_effects <- t(mvtnorm::rmvnorm(n.subj,mean=allparameters,Cov_blocked))
exp(subj_random_effects)
parPreFixs <- gsub("[|].*", "", rownames(subj_random_effects))
df <- data.frame(subject = 1:n.subj)
for(i in 1:n.exp){
  modIdx <- which(parPreFixs == paste0("Mod", i))
  task_random_effects <- subj_random_effects[modIdx,]
  rownames(task_random_effects) <- gsub(".*[|]", "", rownames(task_random_effects))
  data <- alldata[[i]]
  for (j in 1:n.subj){
    tmp<- log_likelihood(task_random_effects[,j],sample=TRUE,data=data[data$subject==j,])
    data$rt[data$subject==j]=tmp$rt
    data$resp[data$subject==j]=tmp$response
  }
  df[paste0("Mod", i)] <- I(list(split(data, f = data$subject)))
}


pars <- rownames(subj_random_effects)
# Create the Particle Metropolis within Gibbs sampler object ------------------


sampler <- pmwgs(
  data = df,
  pars = pars,
  ll_func = joint_ll,
  par_groups = rep(1:n.exp, each = 4 + n.cond)
)
# start the sampler ---------------------------------------------------------
sampler <- init(sampler, n_cores = 12) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 500, n_cores =12)
adapted <- run_stage(burned, stage = "adapt", iter = 5000, n_cores = 12) # set up high, will drop-out early anyways
sampled <- run_stage(sampled, stage = "sample", iter = 1000,n_cores = 12)

