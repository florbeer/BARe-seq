#! /usr/bin/Rscript

#Analysis of enhancer libraries with pasted motifs (library libFL007_STARR_DSCP)
#Estimation of the 0-truncated negative binomial parameters in bash mode with RScript

library(dplyr)
library(numDeriv)
source("R01_Opitmization_Functions.R") #adapt the path if required


#Import datasets of the enhancer libraries
dir_files <- "./data" #adapt the path with the directory containing the files with raw UMI counts of the enhancer libraries (STARR-seq experiments)
count_files <- list.files(path=dir_files)
count_files <- grep("007_STARR", count_files, value=TRUE) #subsetting to file with the libFL007_STARR_DSCP library

#count_files <- c("libFL007_STARR_DSCP_replicate1_batch1_counts.txt", "libFL007_STARR_DSCP_replicate1_batch2_counts.txt", "libFL007_STARR_DSCP_replicate2_batch3_counts.txt", "libFL007_STARR_DSCP_replicate2_batch4_counts.txt")

#Import datasets into a list object
data_list0 <- lapply(count_files, function(y){
  file <- paste(dir_files, y, sep="/")
  x <- read.csv(file, sep="\t",header=TRUE)
  return(x)
})

names(data_list0) <- gsub("_batch._counts.txt", "", count_files)

#Adjust the data, by combining different batches of the data
data_list <- c( Reduce(rbind, data_list0[1:2]), #batch 1 & 2
                Reduce(rbind, data_list0[3:4]) #batch 3 & 4
                )

names(data_list) <- paste("libFL007_STARR_DSCP_replicate", c(1,2), sep="")

#Pool samples of the same library into new slots within the same library
data_list$libFL007_STARR_DSCP_pool <- Reduce(rbind, data_list[1:2])

#Arrange slots in the list in alphabetical order:
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
WriteXLS(dataX, ExcelFileName = "Enhancer_libFL007_burst_parameters_NB.xls", SheetNames =names(dataX),row.names=TRUE)
