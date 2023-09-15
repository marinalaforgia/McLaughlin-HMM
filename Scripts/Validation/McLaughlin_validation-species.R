# simulate mclaughin
rm(list=ls())

library(foreach)
library(doParallel)
library(tidyverse)

source("Scripts/Validation/simulation.R")
source("Scripts/HMM_Functions.R")

mcl.df.q <- readRDS("Data/20221212_mcl_HMM_quadrat.RDS") # quadrat level HMM
mcl.df.q <- filter(mcl.df.q, iter < 150, n.plots >= 50, FunGroup != "Native Grass")

#### simulation ###
row.names(mcl.df.q) <- mcl.df.q$Species_Name
param <- mcl.df.q[,2:6]

trueParam = list() #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.8      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

tmp <- expand.grid(Species_Name = rownames(mcl.df.q), sim = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA)

#### Run HMM ####

sim.df <- foreach(j = rownames(mcl.df.q), .combine = rbind) %dopar% { # for each species
 
  sim <- list()
  
  for(k in 1:100){ # simulate a dataset (400 patches, 19 years)
    sim[[k]] <- simulation(param[j,])
    
    X = as.matrix(sim[[k]]) # use their time series PA data
    n = 5 
    p0Results = rep(0,n)
    gResults = rep(0,n)
    cResults = rep(0,n)
    sResults = rep(0,n)
    rResults = rep(0,n)
    llResults = rep(0,n)
    lltrueParam = rep(0,n)
    
    for (i in 1:n) {
      #for (k in 1:3) print(i) 
      #print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"
      EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
      #print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"; why is this on here twice?
      # fill in results with new params from the best log likelihood model
      p0Results[i] = EMresult$param$p0 # initial seed bank prob
      gResults[i] = EMresult$param$g # germ
      cResults[i] = EMresult$param$c #col
      sResults[i] = EMresult$param$s #surv
      rResults[i] = EMresult$param$r #prod
      llResults[i] = EMresult$ll # log-likelihood
      lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
    }
    
  # individual output for each sim of each species
  tmp[tmp$Species_Name == j & tmp$sim == k,]$p0 <- EMresult$param$p0 
  tmp[tmp$Species_Name == j & tmp$sim == k,]$g <- EMresult$param$g
  tmp[tmp$Species_Name == j & tmp$sim == k,]$c <- EMresult$param$c
  tmp[tmp$Species_Name == j & tmp$sim == k,]$s <- EMresult$param$s
  tmp[tmp$Species_Name == j & tmp$sim == k,]$r <- EMresult$param$r
  tmp[tmp$Species_Name == j & tmp$sim == k,]$iter <- which.max(EMresult$lllist)

  }
  
   # mean sim output per species
  
      cbind(j, 
            mean(tmp[tmp$Species_Name == j,]$p0), 
            mean(tmp[tmp$Species_Name == j,]$g), 
            mean(tmp[tmp$Species_Name == j,]$c), 
            mean(tmp[tmp$Species_Name == j,]$s),
            mean(tmp[tmp$Species_Name == j,]$r),
            mean(tmp[tmp$Species_Name == j,]$iter),
            sd(tmp[tmp$Species_Name == j,]$p0), 
            sd(tmp[tmp$Species_Name == j,]$g), 
            sd(tmp[tmp$Species_Name == j,]$c), 
            sd(tmp[tmp$Species_Name == j,]$s),
            sd(tmp[tmp$Species_Name == j,]$r),
            sd(tmp[tmp$Species_Name == j,]$iter)
      )
}
  


#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("Species_Name", "p0", "g", "c", "s", "r", "iter","p0.sd", "g.sd", "c.sd", "s.sd", "r.sd", "iter.sd")

saveRDS(sim.df, "McL-sim-species-validation.RDS")