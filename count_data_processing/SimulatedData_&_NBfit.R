#! /usr/bin/Rscript

#Generating a simulated data for testing accuracy of the model
#Script run with Rscript

library(dplyr)
library(numDeriv)

library(cowplot)
library(LearnBayes)
#library(pheatmap)

source("R01_Opitmization_Functions.R") #adapt the path if required

#Set a range of Negative Binomial parameters mu and r (dispersion parameter)
#Selected range of values:
#mu: values between: -2.85 and 1.15 in log10 scale
#r: values between: -2.85 and 2.25 in log10 scale
#Fixed number of drawn observations: 1700
#Number of simulated sequences: 5000

#mu values
a <- seq(from=-2.85, to=1.15, length.out=10000)
set.seed(15)
sm <- sample(a, 5000, replace=TRUE)
sm <- 10^sm
summary(sm)

#r values
a <- seq(from=-2.85, to=2.25, length.out=10000)
set.seed(90)
sr <- sample(a, 5000, replace=TRUE)
sr <- 10^sr
summary(sr)

#Number of observations
sn <- rep(1700, 5000)

#Making a data frame with ground truth values
df_GT <- data.frame(sim_mu=sm, sim_bfreq=sr, drawn_obs=sn)
df_GT$sim_bsize <- sm/sr
rownames(df_GT) <- df_GT$seq_ID <- paste("Seq", 1:nrow(df_gt), sep="")
head(df_GT)
write.table(df_GT, "SimData_NB_GTpars.txt", quote = FALSE)


#Custom function for computing the probability of 0-truncated NB distribution and sampling counts from it
#It allows to generate non-0-counts
r0T_NB <- function(x, n, mu, size){
  d <- dnbinom(x, mu=mu, size=size)/(1-dnbinom(0, mu=mu, size=size))
  sa <- sample(x, size=n, replace=TRUE, prob=d)
  return(sa)
}

#Generating counts
WD  <- apply(cbind(sm,sr, sn),1,function(y){
  r0T_NB(1:1000, n=y[3], mu=y[1], size=y[2])
})
dim(WD)

#Arranging counts into a data frame, similar to the count matrices with experimental data
sim_counts <- data.frame(seq_ID = paste0("Seq", rep(1:ncol(WD), each=nrow(WD) )),
                 plasmid_bc = paste0("barcode", rep(1:nrow(WD), times=ncol(WD) )),
                 N=expand.grid(WD))
colnames(sim_counts) <- c("seq_ID", "plasmid_bc", "N")
dim(sim_counts)
head(sim_counts)

#Saving the count matrix
write.table(sim_counts, "SimData_counts.txt", quote = FALSE)


#Parameter inference
#dfw$sample <- "simDat"
optW <- FUN_NBmle(sim_counts, method="BFGS", init=c(1,1), control=list(maxit=10000))


#Perform KS test
pvalW <- FUN_KStest(sim_counts, optW)


#Arrange inferred parameters into a data frame
dataW <- sapply(optW, function(y2) y2$par)
dim(dataW)
dataW <- t(dataW) %>% as.data.frame()
colnames(dataW) <- c("exp", "bfreq")
dataW$r <- dataW$bfreq
dataW$bsize <- dataW$exp/dataW$bfreq
rownames(dataW) <- names(optW)

#Arrange the columns of the data frame
dataW <- dataW[,c("bfreq", "bsize", "exp")]

#Add the p-values of the KS test
dataW$KS_pval <- pvalW


#Join data frame with ground truth values (df_GT) and inferred parameters (dataW)
dataW$seq_ID <- rownames(dataW)

dataW2 <- full_join(df_GT, dataW, join_by(seq_ID == seq_ID))
dim(dataW2)
head(dataW2)


#Save table with parameters
saveRDS(dataW2, "SimData_NBFit.rds")