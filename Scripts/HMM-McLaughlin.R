#### McLaughlin HMM ####
rm(list = ls())
source("Scripts/HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(tidyverse)

# Prep Data ####
trait <- read.csv("Data/McL_80SitesSpeciesTraits_012615.csv")
mcl <- read.csv("Data/Core_Community_Data2019.csv")
serp <- read.csv("Data//McL_Abiotic-Data.csv")

rm <- c("Logfia spp./Micropus californicus", "Aegilops triuncialis", "Trifolium bifidum/gracilentum", "Torilis arvensis/nodosa", "Hesperolinon spp.", "Lactuca serriola", "Cryptantha hispidula", "Cuscuta californica", "Silene gallica", "Erodium botrys", "Erodium brachycarpum", "Trifolium dubium", "Trifolium microdon", "Trifolium microcephalum", "Logfia gallica") # rm commonly mis-ID'd species, species that hadnt always been censused (cuscuta) and goatgrass which is treated at the reserve and occasionally pulled from transects

mcl <- filter(mcl, Species_Name %in% trait[trait$Annual.Perennial == "Annual",]$Species_Name, !(Species_Name %in% rm)) 
mcl <- merge(serp[,1:2], mcl, by = "Site")

mcl[mcl$Species_Name == "Geranium molle",]$Species_Name <- "Geranium dissectum"
mcl <- ddply(mcl, .(Year, Site, Quadrat, Serpentine, Species_Name), summarize, Cover = sum(Cover))


##### Transect Level ####
# mcl <- ddply(mcl, .(Site, Year, Species_Name), summarize, PA = 1)
# mcl <- pivot_wider(mcl, names_from = "Year", values_from = "PA")
# mcl <- as.data.frame(mcl)
# mcl[is.na(mcl)] <- 0
# mcl <- merge(serp[,1:2], mcl, by = "Site")
# 
# mcl <- filter(mcl, Serpentine == "S")
# 
# # full <- data.frame()
# # for(s in names(sites)){
# #  mcl <- filter(mcl.tmp, Site %in% sites[[s]])
# 
# # looking at species averages
# mcl$avg.P.site <- apply(mcl[,3:22], 1, mean)
# mcl.sp.avg <- ddply(mcl, .(Species_Name), summarize, avg.P.sp = mean(avg.P.site))

##### Quadrat Level ####
#mcl[which(mcl$Cover == 0),] # checks out

mcl$PA <- 1
mcl$ID <- paste(mcl$Site, mcl$Quadrat, sep = ".")

mcl <- pivot_wider(mcl[,c(1,5,7,8)], names_from = "Year", values_from = "PA")
mcl[is.na(mcl)] <- 0

# looking at species averages
mcl$avg.P.site <- apply(mcl[,3:ncol(mcl)], 1, mean)
mcl.sp.avg <- ddply(mcl, .(Species_Name), summarize, avg.P.sp = mean(avg.P.site))

mcl <- mcl[,-23]

#### Create Species List ####
length(unique(mcl$ID))

#Full 400
# serp 190
# nonserp: 210
species.list <- list()

for(i in unique(mcl$Species_Name)) {
  tmp <- filter(mcl, Species_Name == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  #if(nrow(tmp) >= length(unique(mcl$ID))*.2) species.list[[i]] <- tmp 
  if(nrow(tmp) >= 20) species.list[[i]] <- tmp 
}


df <- data.frame(names = NA, obs = NA)

for(n in names(species.list)) {
  df <- data.frame(rbind(df, c(names = n, obs = nrow(species.list[[n]]))))
}


df$obs <- as.numeric(as.character(df$obs))


#### Run HMM ####
trueParam = list() #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.8      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)


mcl.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates


for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
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
    for (k in 1:5) print(i) # I don't understand what's going on here, why print it 5 times?
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

  mcl.df[mcl.df$Species_Name == j,]$p0 <- EMresult$param$p0
  mcl.df[mcl.df$Species_Name == j,]$g <- EMresult$param$g
  mcl.df[mcl.df$Species_Name == j,]$c <- EMresult$param$c
  mcl.df[mcl.df$Species_Name == j,]$s <- EMresult$param$s
  mcl.df[mcl.df$Species_Name == j,]$r <- EMresult$param$r
  mcl.df[mcl.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  mcl.df[mcl.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  mcl.df[mcl.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  mcl.df[mcl.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
  mcl.df[mcl.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
}


#### Explore Output ####
trait <- read.csv("Data/McL_80SitesSpeciesTraits_012615.csv")
trait$FunGroup <- paste(trait$Native.Exotic, trait$Grass.Forb.Shrub)
mcl.df <- merge(mcl.df, trait[,c(1,3,24)], by = "Species_Name", all.y = F)

mcl <- merge(mcl, trait[,c(1,3,24)], by = "Species_Name", all.y = F)
mcl$site.id <- paste(mcl$Site, mcl$Quadrat, sep = ".")

saveRDS(mcl.df, "Scripts/HMMs/McLaughlin/20230704_mcl_HMM_quadrat.RDS")

