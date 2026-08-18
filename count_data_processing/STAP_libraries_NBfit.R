#! /usr/bin/Rscript

#Analysis of promoter libraries
#Estimation of the 0-truncated negative binomial parameters in bash mode with RScript

library(dplyr)
library(numDeriv)
source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required


#Import datasets of the promoter libraries
dir_files <- "./data" #adapt the path with the directory containing the files with raw UMI counts of the promoter libraries (STAP-seq experiments)
count_files <- list.files(path=dir_files)
count_files <- grep("STAP", count_files, value=TRUE) #subsetting to file with only STAP-seq experiments

#Import datasets into a list object
data_list <- lapply(count_files, function(y){
  file <- paste(dir_files, y, sep="/")
  x <- read.csv(file, sep="\t",header=TRUE)
  return(x)
})

#Pool samples of the same library into new slots on the data_list object
data_list$libFL001_STAP_zfh1_pool <- Reduce(rbind, data_list[1:2])
data_list$libFL002_STAP_CG4822_pool <- Reduce(rbind, data_list[3:4])
data_list$libFL003_STAP_Pld_pool <- Reduce(rbind, data_list[5:6])

#Arrange slots in the data_list object in alphabetical order:
samples <- names(data_list)
samples <- sort(samples)
data_list <- data_list[samples]


#Parameter inference
optX <- lapply(data_list, function(y) FUN_NBmle(y, method="BFGS", init=c(1,1), control=list(maxit=10000)))


#Perform KS test
pvalX <- lapply(samples, function(sa){
  y <- data_list[[sa]]
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
WriteXLS(dataX, ExcelFileName = "Promoter_burst_parameters_NB.xls", SheetNames =names(dataX),row.names=TRUE)
