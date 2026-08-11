#! /usr/bin/Rscript

#Estimation of the negative binomial parameters in bash mode with RScript


#library(scModels)
library(dplyr)
library(numDeriv)
source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required


#Import data
x <- read.csv("E28_RT_test_processing1_counts.txt",sep="\t",header=TRUE)

samples <- unique(x$sample) %>% sort()
aL <- split(x, x$sample)
aL <- aL[samples]

#Add an extra slot with the pooled sample
f <- grepl("RT2_r(2|3)", x$sample)
aL$pooled <- x[f,]

samples <- names(aL)

#Parameter inference
optX <- lapply(aL, function(y) FUN_NBmle(y, method="BFGS", init=c(1,1), control=list(maxit=10000)))


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
  colnames(df) <- c("exp", "bfreq")
  df$bsize <- df$exp/df$bfreq
  rownames(df) <- names(y)
  #Arrange the columns of the data frame
  df <- df[,c("bfreq", "bsize", "exp")]
  
  return(df)
})

#Add the p-values of the KS test
for (i in 1:length(dataX)){
  dataX[[i]]$KS_pval <- pvalX[[i]]
}


#Export the files
library(WriteXLS)
#WriteXLS(dataX, ExcelFileName = "E28_burst_parameters_NB.xls", SheetNames =names(dataX),row.names=TRUE)
