#! /usr/bin/Rscript

#Estimation of the negative binomial parameters in bash mode with RScript
#setwd() #if required

library(scModels)
library(dplyr)
library(numDeriv)
source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required


#Import data
x <- read.csv("E28_RT_test_processing1_counts.txt",sep="\t",header=TRUE) #Example with E28 count file
#Adapt the previous file name accordingly

samples <- unique(x$sample) %>% sort()
aL <- split(x, x$sample)
aL <- aL[samples]


#Parameter inference
optX <- lapply(aL, function(y) FUN_NBmle(y, method="BFGS", init=c(1,1), control=list(maxit=10000)))


#Hessian matrix estimation
hesX <- lapply(samples, function(sa){
  y <- aL[[sa]]
  optN <- optX[[sa]]
  hesN <- FUN_HesM(y=y, optN=optN)
  return(hesN)
})
names(hesX) <- samples

#Perform KS test
pvalX <- lapply(samples, function(sa){
  y <- aL[[sa]]
  optN <- optX[[sa]]
  pval <- FUN_KStest(y=y, optN=optN)
  return(pval)
})

names(pvalX) <- samples


#Arrange results into a data frame
dataX <- lapply(optX, function(y){
  df <- sapply(y, function(y2) y2$par)
  
  df <- t(df) %>% as.data.frame()
  colnames(df) <- c("exp", "kon")
  df$r <- df$bfreq <- df$kon
  df$bsize <- df$exp/df$kon
  rownames(df) <- names(y)
  return(df)
})


#Then get the standard deviation values from the hessian matrices
aL2 <- lapply(hesX, function(y){
  std <- sapply(y, FUN_hes2sd)
  return(std)
})

for (i in 1:length(dataX)){
  df <- dataX[[i]]
  std <- aL2[[i]]
  df$r.delta <- df$bfreq.delta <- std[2,]
  df$bsize.delta <- sqrt( ( std[1,]/df$bfreq )**2 + ( std[2,]*df$exp/df$bfreq**2 )**2 ) #NA values will give NA on this expression
  df$mu.delta <- std[1,]
  
  #Order the data frame
  df <- df[,c("bfreq", "bsize", "bfreq.delta", "bsize.delta", "exp", "mu.delta", "r", "r.delta")]
  
  dataX[[i]] <- df
}


#Save files
aL <- list(optLN = optX, hesN = hesX, pvalL=pvalX, bDataN=dataX)
saveRDS(aL, "E28_NB_Optim_results.rds") #This is just an example, adapt file name accordingly
