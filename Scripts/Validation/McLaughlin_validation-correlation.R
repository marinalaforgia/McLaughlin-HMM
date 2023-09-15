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
param <- data.frame(matrix(ncol = 5, nrow = 5000))

colnames(param) <- c('p0', 'g', 'c', 's', 'r')

# independently sample from existing species parameters
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

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

sim.df <- foreach(j = 1:5000, .combine = rbind) %dopar% { # for each species
  X = as.matrix(sim[[j]]) # use their time series PA data
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
  
  cbind(j, 
        EMresult$param$p0, 
        EMresult$param$g, 
        EMresult$param$c, 
        EMresult$param$s,
        EMresult$param$r,
        which.max(EMresult$lllist))
}
  
#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("Species_Name", "p0", "g", "c", "s", "r", "iter")

#saveRDS(sim.df, "McL-sim.RDS")

sim.df <- filter(sim.df, iter < 150)

sim.cor <- data.frame(cor = NA, p.value = NA)

for(i in 1:50000){
  tmp <- slice_sample(sim.df, n = 59, replace = F)
  tmp.test <- cor.test(tmp$s, tmp$c)
  sim.cor[i,]$cor <- tmp.test$estimate
  sim.cor[i,]$p.value <- tmp.test$p.value
}


ggplot(sim.cor, aes(x = cor)) +
  geom_histogram(bins = 20, col = "black", fill = "white") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  geom_vline(xintercept = -0.553, col = "red") + # correlation from real data is -0.553
  labs(x = "simulated correlation coefficient")

mean(sim.cor$cor)
nrow(sim.cor[sim.cor$cor <= -0.553,])/nrow(sim.cor) #p value < 0.0001...?