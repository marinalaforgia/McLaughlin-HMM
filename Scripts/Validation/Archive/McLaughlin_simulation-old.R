# simulate mclaughin
rm(list=ls())

source("Scripts/HMMs/simulation/simulation.R")
source("Scripts/HMMs/simulation/HMM_Functions.R")

mcl.df.q <- readRDS("Scripts/HMMs/simulation/20221212_mcl_HMM_quadrat.RDS") # quadrat level HMM

#### simulation ###
param <- data.frame(matrix(ncol = 5, nrow = 5000))

colnames(param) <- c('p0', 'g', 'c', 's', 'r')

param$p0 <- sample(mcl.df.q$p0, size = 5000, replace = T)
param$g <- sample(mcl.df.q$g, size = 5000, replace = T)
param$c <- sample(mcl.df.q$c, size = 5000, replace = T)
param$s <- sample(mcl.df.q$s, size = 5000, replace = T)
param$r <- 1

sim <- list()

for(i in 1:nrow(param)){
  sim[[i]] <- simulation(param[i,])
}

#### Run HMM ####
trueParam = list() #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.8      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)


sim.df <- expand.grid(Species_Name = 1:5000, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates

system.time(for(j in 1:3) { # for each species
  X = as.matrix(sim[[j]]) # use their time series PA data
  n = 5 
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  lltrueParam = rep(0,n)
 
  print(j)  
  
  for (i in 1:n) {
    for (k in 1:5) print(i) 
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"; why is this on here twice?
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
    lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
  }

  sim.df[sim.df$Species_Name == j,]$p0 <- EMresult$param$p0
  sim.df[sim.df$Species_Name == j,]$g <- EMresult$param$g
  sim.df[sim.df$Species_Name == j,]$c <- EMresult$param$c
  sim.df[sim.df$Species_Name == j,]$s <- EMresult$param$s
  sim.df[sim.df$Species_Name == j,]$r <- EMresult$param$r
  sim.df[sim.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
})
# little over 7 min to do THREE 'species' using 1 core
#    user  system elapsed 
# 411.186   4.806 426.316 

saveRDS(sim.df, "McL-sim.RDS")