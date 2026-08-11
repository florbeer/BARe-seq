#! /usr/bin/Rscript

#Computing the profile likelihood and confidence intervals in bash mode with RScript
#setwd() #if required

#library(scModels)
library(dplyr)
library(numDeriv)
library(rootSolve)
source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required


#Import data
x <- read.csv("E28_RT_test_processing1_counts.txt",sep="\t",header=TRUE) #Example with E28 count file
#Adapt the previous file name accordingly
samples <- unique(x$sample) %>% sort()

############

#Importing the inferred parameters
aL <- readRDS("E28_NB_Optim_results.rds") #Example with E28 count file, adapt file name accordingly

names(aL)
optLN <- aL$optLN
hesN <- aL$hesN
pvalL <- aL$pvalL
bDataN <- aL$bDataN

#Run
system.time(CIN <- FUN_CI95(x, optLN, hesN, samples))

#Save results
saveRDS(CIN, "E28_NB_ConfidenceIntervals.rds") #This is just an example, adapt file name accordingly
