#Script for initial data analysis for a given BARe-seq experiment
#It can be either promoter or enhancer libraries
#It contains:
#Quality control checks (detected sequences per samples, barcodes per samples, etc..)
#Data arrangement into data frames. It includes parameter values, filtering criteria, standard deviation values, precision score values, etc.
#Plots for data exploration, checking different filtering criteria, etc.

#setwd() #if required

#library(scModels)
library(dplyr)
library(numDeriv)

library(ggplot2)
library(ggpubr)
library(patchwork)

source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required

#Import data
x <- read.csv("E28_RT_test_processing1_counts.txt",sep="\t",header=TRUE) #Example with E28 count file
#Adapt the previous file name accordingly

head(x)
dim(x)
unique(x$sample)

#Initial quality measurements
#Plasmid barcodes per regulatory element
aL <- lapply(samples, function(y){
  y2 <- x[x$sample == y,]
  y2 <- aggregate( rep(1,nrow(y2)), by=list(y2$seq_ID), sum)
  y3 <- y2[,2]
  names(y3) <- y2[,1]
  return(y3)
})
sapply(aL, length); sapply(aL, summary)
sapply(aL, head)
names(aL) <- samples
#Color code, when having 3 libraries and 3 samples per library
#Adjust accordingly
a <- scales::hue_pal()(3)
cols <- rep(a, each=3)

boxplot(aL, cex=0.2, main="Plasmid barcodes per regulatory element", log="y", border=cols)
abline(h = median(aL[[1]]), lty=2)


#Estimate UMI counts per (sequence) regulatory element
aL2 <- lapply(samples, function(y){
  y2 <- x[x$sample == y,]
  y2 <- aggregate(y2$N, by=list(y2$seq_ID), sum)
  y3 <- y2[,2]
  names(y3) <- y2[,1]
  return(y3)
})
sapply(aL2, length); sapply(aL2, summary)
sapply(aL2, head)
names(aL2) <- samples
boxplot(aL2, cex=0.2, main="UMI counts per (sequence) regulatory element", log="y", border=cols)
abline(h = median(aL2[[1]]), lty=2)

#Estimate normalized counts per sequence, i.e. UMI counts per sequence, divided by the number of sequences per sample
aL3 <- lapply(samples, function(y){
  y2 <- x[x$sample == y,]
  y3 <- aggregate( y2$N, by=list(y2$seq_ID), sum)
  
  y4 <- aggregate( rep(1,nrow(y2)), by=list(y2$seq_ID), sum)
  y4 <- y3[,2] / y4[,2]
  names(y4) <- y3[,1]
  return(y4)
})

sapply(aL3, length); sapply(aL3, summary)
sapply(aL3, head)
names(aL3) <- samples
boxplot(aL3, cex=0.2, main="Normalized counts per (sequence) regulatory element", log="y", border=cols)
abline(h = median(aL3[[1]]), lty=2)

#Library size:
df <- aggregate(x$N, by=list(x$sample), sum)
colnames(df) <- c("sample", "libsize")
head(x)
#Add library / condition information
df$condition <- gsub(".+_E28_", "", df$sample) %>% gsub("_r.", "", .) #adjust this code accordingly

ggplot(df, aes(x = sample, y = libsize)) +
  geom_bar(aes(color = condition, fill = condition), stat = "identity", position = position_dodge(0.8),
           width = 0.7) +
  geom_text(aes(label = libsize, group = condition), position = position_dodge(0.8), vjust = -0.3, size = 3.5)


#Detected sequences (regulatory elements) per sample
aL4 <- sapply(samples, function(y){
  y2 <- x[x$sample == y,]
  y3 <- unique(y2$seq_ID) %>% length()
  return(y3)
})

df$detected_seqs <- aL4

ggplot(df, aes(x = sample, y = detected_seqs)) +
  geom_bar(aes(color = condition, fill = condition),
           stat = "identity", position = position_dodge(0.8),
           width = 0.7) +
  labs(title="Detected sequences (regulatory sequences) per samples", x ="Sample", y = "number sequences") +
  geom_text(aes(label = detected_seqs, group = condition), position = position_dodge(0.8), vjust = -0.3, size = 3.5)

#Check also number of barcodes detected in each sample
a$detected_bc <- sapply(samples, function(y) sum(x$sample == y) )
ggplot(a, aes(x = samples, y = detected_bc)) +
  geom_bar(aes(color = lib, fill = lib),
           stat = "identity", position = position_dodge(0.8),
           width = 0.7) +
  labs(title="Detected plasmid barcodes per samples", x ="Sample", y = "number plasmid barcodes") +
  geom_text(aes(label = detected_bc, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5)


#Apply filtering criteria:
#Option 1: keep sequences at least 3 UMI counts in at least 3 barcodes
aL <- split(x, x$sample)
aL2 <- lapply(aL, function(y){
  y2 <- split(y$N, y$seq_ID)
  f <- sapply(y2, function(y3){
    sum(y3 == 3) > 2
  })
  y4 <- names(y2)[f]
  return(y4)
})
sapply(aL2, length)
head(aL2[[1]])

#Option 2: keep sequences with >2 UMI counts in at least 3 barcodes
aL3 <- lapply(aL, function(y){
  y2 <- split(y$N, y$seq_ID)
  f <- sapply(y2, function(y3){
    sum(y3 > 2) > 2
  })
  y4 <- names(y2)[f]
  return(y4)
})
sapply(aL3, length)
head(aL3[[1]])



##################
#Parameter inference is performed in bash mode with R03_NBFit_Exe.R file
#It includes: inferred parameters mu and r for the Negative binomial model
#hessian matrices, used to estimate standard deviation of the parameters
#p adjusted values when performing Kolmogorov-Smirnov test
#Data frames with different quantities summarized from the previous scores. There is a data frame per sample, arranged into a big list


#Import results
aL <- readRDS("E28_NB_Optim_results.rds") #Example with E28 experiment. Adapt accordingly
names(aL)
optLN <- aL$optLN #inferred parameter values
hesN <- aL$hesN #hessian matrices
pvalL <- aL$pvalL #p values of the KS test
bDataN <- aL$bDataN #list of data frames


#Add additional information to the data frames
#Compute the sd thresholds: sd < 5*estimate
bDataN2 <- lapply(bDataN, as.data.frame)
bDataN2 <- lapply(bDataN2, function(y) {
  #y$r <- y$bfreq
  #y$r.delta <- y$bfreq.delta
  y$bfreq.thres <- y$bfreq.delta < 5*y$bfreq & !is.na(y$bfreq) & !is.na(y$bfreq.delta)
  y$bsize.thres = y$bsize.delta < 5*y$bsize & !is.na(y$bsize) & !is.na(y$bsize.delta)
  y$mu_sd.thres <- y$mu.delta < 5*y$exp & !is.na(y$exp) & !is.na(y$mu.delta)
  y$r_sd.thres <- y$r.delta < 5*y$r & !is.na(y$r) & !is.na(y$r.delta)
  y$sd_thres2 <- y$mu_sd.thres & y$r_sd.thres
  return(y)
})

#Add p value from Kolmogorov-Smirnov test
sapply(pvalL, summary)
head(pvalL[[1]])
sapply(pvalL, length)
sapply(bDataN2, nrow)
#Since the p values are in the same order as optLN, we can just add them directly to the data frame
for (i in 1:length(bDataN2)){
  bDataN2[[i]]$KS_pval <- pvalL[[i]]
}
head(bDataN2[[1]])

#Add the two filtering criteria for passing the quality control:
#Option 1: keep sequences at least 3 UMI counts in at least 3 barcodes
aL <- split(x, x$sample)
aL2 <- lapply(aL, function(y){
  y2 <- split(y$N, y$seq_ID)
  f <- sapply(y2, function(y3){
    sum(y3 == 3) > 2
  })
  y4 <- names(y2)[f]
  return(y4)
})
sapply(aL2, length)
head(aL2[[1]])

#Option 2: keep sequences with >2 UMI counts in at least 3 barcodes
aL3 <- lapply(aL, function(y){
  y2 <- split(y$N, y$seq_ID)
  f <- sapply(y2, function(y3){
    sum(y3 > 2) > 2
  })
  y4 <- names(y2)[f]
  return(y4)
})
sapply(aL3, length)
head(aL3[[1]])

for (i in 1:length(bDataN2)){
  y <- bDataN2[[i]]
  y$filter1 <- rownames(y) %in% aL2[[i]]
  y$filter2 <- rownames(y) %in% aL3[[i]]
  bDataN2[[i]] <- y
}

#Add the convergence value from optLN
length(optLN)
length(optLN[[1]])
length(optLN[[1]][[1]])
head(bDataN2[[1]])
for (i in 1:length(optLN)){
  identical(names(optLN[[i]]), rownames(bDataN2[[i]])) %>% print()
}

for (i in 1:length(bDataN2)){
  y <- bDataN2[[i]]
  y$convergence <- sapply(optLN[[i]], function(y2) y2$convergence)
  bDataN2[[i]] <- y
}

#Comparison between technical replicates, to get reproducible values

#For experiments with up to 3 samples, make pair-wise comparisons of one sample versus another one at a time
#Use self-custom function FUN_repVal2 from ../main_functions/R01_Opitmization_Functions.R script.
names(bDataN2)
sapply(bDataN2, dim) #Checking which are the samples with higher number of observations to use as a query

at1 <- 2
at2 <- 0.5
#Just some examples:
bDataN2[[1]] <- FUN_repVal2(listdf=bDataN2, query=1, bkg=2, thr1=at1, thr2=at2) #Compares samples 1 versus sample 2 and returns columns of reproducible values for sample 1
bDataN2[[2]] <- FUN_repVal2(listdf=bDataN2, query=2, bkg=3, thr1=at1, thr2=at2)
bDataN2[[3]] <- FUN_repVal2(listdf=bDataN2, query=3, bkg=2, thr1=at1, thr2=at2)

bDataN2[[4]] <- FUN_repVal2(listdf=bDataN2, query=4, bkg=5, thr1=at1, thr2=at2)
bDataN2[[5]] <- FUN_repVal2(listdf=bDataN2, query=5, bkg=4, thr1=at1, thr2=at2)
bDataN2[[6]] <- FUN_repVal2(listdf=bDataN2, query=6, bkg=5, thr1=at1, thr2=at2)

bDataN2[[7]] <- FUN_repVal2(listdf=bDataN2, query=7, bkg=9, thr1=at1, thr2=at2)
bDataN2[[8]] <- FUN_repVal2(listdf=bDataN2, query=8, bkg=9, thr1=at1, thr2=at2)
bDataN2[[9]] <- FUN_repVal2(listdf=bDataN2, query=9, bkg=7, thr1=at1, thr2=at2)


head(bDataN2[[1]])
sapply(bDataN2, function(y) sum(y$pair_rep))

#Confidence intervals
#This is computed in bash mode with R04_ProfLik_Exe.R script

#Import results
CIN <- readRDS("E28_NB_Optim_results.rds") #This is just an example, adapt file name accordingly
class(CIN)
length(CIN)
sapply(CIN, length)
length(bDataN2)
length(CIN[[1]][[1]])
sapply(CIN[[1]][[1]], length)

for (i in 1:length(bDataN2)){
  identical(rownames(bDataN2[[i]]), names(CIN[[i]])) %>% print()
}
#Check how many threshold values are in the samples
sapply(CIN, function(y) {
  #sapply(y, function(y2) length(y2$mu_roots)) %>% table()
  sapply(y, function(y2) length(y2$r_roots)) %>% table()
})
#There are some examples with 1 and 0 interval values, instead of 2
#Maybe these are sequences with low quality. Lets make a quick check
head(bDataN2[[1]])
length(bDataN2)
dataX <- bDataN2
for (i in 1:length(dataX)){
  df <- dataX[[i]]
  df$nthres_mu <- sapply(CIN[[i]], function(y) length(y$mu_roots)) %>% paste0("n", .)
  df$nthres_r <- sapply(CIN[[i]], function(y) length(y$r_roots)) %>% paste0("n", .)
  dataX[[i]] <- df
}
yD <- lapply(dataX, function(y) y[y$filter2,])
sapply(yD, nrow)
sapply(yD, function(y) table(y$nthres_mu)) #all of them with 2 values
sapply(yD, function(y) table(y$nthres_r)) #some examples still with 1 value

rm(dataX, yD)

#Proceed with annotating the threshold values into the data frame
for (i in 1:length(bDataN2)){
  df <- bDataN2[[i]]
  df$thres_r2 <- df$thres_r1 <- df$thres_mu2 <- df$thres_mu1 <- NA
  df$thres_mu1 <- sapply(CIN[[i]], function(y) y$mu_roots[1])
  df$thres_mu2 <- sapply(CIN[[i]], function(y) y$mu_roots[2])
  df$thres_r1 <- sapply(CIN[[i]], function(y) y$r_roots[1])
  df$thres_r2 <- sapply(CIN[[i]], function(y) y$r_roots[2])
  bDataN2[[i]] <- df
}
head(bDataN2[[1]])

#Estimate the precision score, defined as: estimate / (CI_high - CI_low)
sapply(bDataN2, function(y){ sum(is.na(y$thres_r2)) })
bDataN2 <- lapply(bDataN2, function(y){
  #Precision for mu
  y$mu_prec <- apply(y[,c("exp", "thres_mu1", "thres_mu2")], 1,function(ay) {
    if(anyNA(ay)) { a <- NA}
    else {
      #a <- (ay[3] - ay[2]) / ay[1]
      a <- ay[1] / (ay[3] - ay[2])
    }
    return(a)
  })
  #Precision for r
  y$r_prec <- apply(y[,c("r", "thres_r1", "thres_r2")], 1,function(ay) {
    if(anyNA(ay)) { a <- NA}
    else {
      #a <- (ay[3] - ay[2]) / ay[1]
      a <- ay[1] / (ay[3] - ay[2])
    }
    return(a)
  })
  #Geometric average
  y$full_prec <- apply(y[,c("mu_prec", "r_prec")], 1, function(ay){
    if(anyNA(ay)) { a <- NA}
    else { exp(mean( log(ay) )) }
  })
  
  return(y)
})

sapply(bDataN2, function(y) summary(y$mu_prec) )
sapply(bDataN2, function(y) summary(y$r_prec) )
sapply(bDataN2, function(y) summary(y$full_prec) )


#Set a threshold for the precision score
#Make the precision score to allow for a range of 0.5 of relative error with respect to their parameter values mu and r
#Therefore: expected precision to be > 1

bDataN2 <- lapply(bDataN2, function(y){
  y$mu_prec2 <- y$mu_prec > 1 & !is.na(y$mu_prec)
  y$r_prec2 <- y$r_prec > 1 & !is.na(y$r_prec)
  y$full_prec2 <- y$mu_prec2 & y$r_prec2
  
  return(y)
})

sapply(bDataN2, function(y) sum(y$mu_prec2))
sapply(bDataN2, function(y) sum(y$r_prec2))
sapply(bDataN2, function(y) sum(y$full_prec2))



#Determine mu threshold based on the bimodal distribution of mu values

#Determining arbitrary values by visual inspection. Once the thresholds were chosen, they were applied to all datasets
#Choose just thresholds for mu and r, individually

#Arrange all the information into a big data frame
#Use sequences after filtering (filtering criteria 2), reproducible and with precision score > 1
yD <- lapply(bDataN2, function(y) y[y$filter2 & y$pair_rep & y$full_prec2,])
#Optional use all the sequences without any filter
#yD <- bDataN2
#Another option: apply only filtering criteria 2 and reproducible values
#yD <- lapply(bDataN2, function(y) y[y$filter2 & y$pair_rep,])

samples <- names(bDataN2)

for (i in 1:length(yD)){
  y <- yD[[i]]
  if(nrow(y) == 0){
    y <- y
  }
  else{
    y$seq_ID <- rownames(y)
    y$sample <- samples[i]
    
    y$log2_mu <- log2(y$exp)
    y$log2_r <- log2(y$r)
    y$log2_bfreq <- log2(y$bfreq)
    y$log2_bsize <- log2(y$bsize)
    y$var <- y$exp + (y$exp^2/y$r)
    y$log2_var <- log2(y$var)
  }
  yD[[i]] <- y
}
yD <- Reduce(rbind, yD)
head(yD)

#Make a scatter plot for testing different thresholds
a <- sapply(bDataN2, function(y) log2(y$exp)) %>% unlist() %>% range(., na.rm=TRUE)
a2 <- sapply(bDataN2, function(y) log2(y$r)) %>% unlist() %>% range(., na.rm=TRUE)
plot(yD$log2_mu, yD$log2_r, pch=16, col="gray", xlab="log2 mean", ylab="log2 r", main="NB parameters, BFGS", xlim=a, ylim=a2)
#abline(h=c(-4, -3.5, 2, -4.5, -5), col="orange", lty=2)
#abline(v=c(-1.8, -2, -4), col="orange", lty=2)
abline(h=c(-5, 2), col="orange", lty=2)
abline(v=c(-4), col="orange", lty=2)
f <- yD$log2_mu < -4 | yD$log2_r < -5 | yD$log2_r > 2
points(yD$log2_mu[f], yD$log2_r[f], pch=16, col="red")
legend("topleft", c("sequences", "accepted", "rejected"), col=c("white", "gray", "red"), pch=16, bty="n", y.intersp = 0.5)


#Write the cloud annotation into the data frames contained in the main list
head(bDataN2[[1]])
sapply(bDataN2, function(y) anyNA(y$exp))
sapply(bDataN2, function(y) anyNA(y$r))
anyNA(yD$log2_mu)
anyNA(yD$log2_r)

bDataN2 <- lapply(bDataN2, function(y){
  f <- log2(y$exp) < -4 | log2(y$r) < -5 | log2(y$r) > 2
  y$cloud <- ifelse(f, "cl_low", "cl_high")
  return(y)
})
sapply(bDataN2, function(y) table(y$cloud))


#Determining mu threshold based on the distribution of mean expression values.
#Selecting a threshold that splits the two modes of the distribution
#This is run for each individual dataset, libraries or conditions

#Filter values by precision score
yD <- lapply(bDataN2, function(y) y[y$full_prec2,])
sapply(yD, dim)
samples <- names(bDataN2)

#The next vector "condition" helps to aggregate samples into libraries / conditions / experiments
#Change accordingly based on the dataset
condition <- gsub("_E.+", "", samples) #Example for E28 experiment

#Annotate information into the data frames
for (i in 1:length(yD)){
  y <- yD[[i]]
  
  if(nrow(y) == 0){
    next
  }
  y$seq_ID <- rownames(y)
  y$sample <- samples[i]
  #y$sample2 <- samples2[i]
  #y$library <- libs[i]
  #y$lib_condition <- libs_condition[i]
  y$condition <- condition[i]
  y$log2_mu <- log2(y$exp)
  y$log2_r <- log2(y$r)
  y$log2_bfreq <- log2(y$bfreq)
  y$log2_bsize <- log2(y$bsize)
  y$var <- y$exp + (y$exp^2/y$r)
  y$log2_var <- log2(y$var)
  yD[[i]] <- y
}
yD <- Reduce(rbind, yD)
head(yD)

#Estimate the lower value of the mean
unique(yD$condition)
condition_ref <- unique(condition)
fd <- yD$condition == condition_ref[8]
log2mu <- yD$log2_mu[fd]
length(log2mu)
mu_d <- density(log2mu)
plot(mu_d)
#Add vertical lines, that can be used as the range of values to look for the minimum
abline(v=c(-5, -1.5), col="red")


#Test different bandwidths and sample sizes
a <- seq(from=0.1, to=0.8, by=0.05) #bandwidth values
a2 <- 2^(10:13) # sample size
log2mu <- yD$log2_mu[fd]

aL <- sapply(a, function(y){
  y3 <- sapply(a2, function(y2){
    a3 <- density(log2mu, bw=y, n=y2)
    f <- a3$x > -5 & a3$x < -1.5
    f2 <- which(a3$y == min(a3$y[f]))
    a4 <- a3$x[f2] %>% round(., 3)
    return(a4)
  })
  return(y3)
})
colnames(aL) <-paste("bw", a, sep="_")
rownames(aL) <-paste("n", a2, sep="_")
aL

plot(density(log2mu, bw=0.45), main="density mu") #test different band widths. Choose the one that looks smooth, but without over smoothing it
#lines(density(ab, bw=0.5), col="red")
abline(v=-3, col="red")
abline(v=-3.5, col="orange", lty=2)

#Condition, threshold (in log2 scale) and bandwidth (example for E57 v2 experiment)
#"libFL010_ctrl": -3, 0.45
#"libFL010_Ecd": -3.4, 0.35
#"libFL030": -2.9, 0.3
#"libFL031": -3.4, 0.35
#"libFL033": -4, arbitrary thres. There are two mu modes, but their values are high (around -1 and 0.25) compared
#to drosophila data sets (remember this is the human library)
#"libFL034": -3.7, 0.4
#"libFL035": -2.1, 0.6
#"libFL036": -2.4, 0.45

#Write it down in the data frames
head(bDataN2[[1]])
colnames(bDataN2[[1]])
ab <- c(-3, -3.4, -2.9, -3.4, -4, -3.7, -2.1, -2.4) %>% rep.int(., c(2,2,4,4,2,2,4,4)) #Adjust this information based on the experiment, number of samples per condition, selected threshold values
length(ab)
for (i in 1:length(bDataN2)){
  bDataN2[[i]]$mu_bimodal_thres <- ab[i]
}


#########################################################
#Export data frame as .xls file.
df <- bDataN2
sapply(df, class)
head(df[[1]])

#Export
library(WriteXLS)
#WriteXLS(df, ExcelFileName = "E28_burst_parameters_NB.xls", SheetNames =names(df),row.names=TRUE) #Example for E28 experiment, adjust name accordingly

#Another option: export also into the data frames into the compressed R format (.rds)
#saveRDS(df, "E28_burst_parameters_NB.rds")



#########################################################
#########################################################
#Diagnostic plots
#Make another data frame that summarizes number of sequences and so on

#Reproducible sequences (based on pair-wise comparisons):
df <- sapply(bDataN2, function(y) sum(y$pair_rep)) %>% as.data.frame()
colnames(df)[1] <- c("repro")
df$repro_fr <- sapply(bDataN2, function(y) sum(y$pair_rep)/nrow(y)) %>% round(.,3)

samples <- names(bDataN2)
df$sample <- names(samples)
df$lib <- gsub("_E38.+", "", samples) #Example with E38 experiment samples, change accordingly
df$lib_condition <- gsub("_E38.+", "", samples) #Extra option, mainly used for pooled samples and complex data set designs. Adjust accordingly.
#In this case it's equal to the library information

#Vector of colors. Adjust accordingly
a <- scales::hue_pal()(4)
cols <- rep(a, each=3)


ggplot(df, aes(x = sample, y = repro)) +
  labs(title="Reproducible values, absolute numbers", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = repro, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = repro_fr)) +
  labs(title="Fraction of reproducible values", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = repro_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

#Check reproducible sequences after filtering criterion 1 (not used anymore):
#yD <- lapply(bDataN2, function(y) y[y$filter1,])
#Check reproducible sequences after filtering criterion 2:
yD <- lapply(bDataN2, function(y) y[y$filter2,])
#Number of sequences passing QC filtering
df$filter <- sapply(yD, nrow)

#Estimate the total sequences to compute also the fraction
aL4 <- sapply(samples, function(y){
  y2 <- x[x$sample == y,]
  y3 <- unique(y2$seq_ID) %>% length()
  return(y3)
})
sapply(yD, nrow) / aL4
df$filter_fr <- (sapply(yD, nrow) / aL4) %>% round(.,3)

#Now check how many of them are reproducible
df$repro2 <- sapply(yD, function(y) sum(y$pair_rep))
df$repro2_fr <- sapply(yD, function(y) sum(y$pair_rep)/nrow(y)) %>% round(.,3)


#Number of sequences with an estimate
df$n_seqs <- sapply(bDataN2, nrow)
df$n_seqs_fr <- (sapply(bDataN2, nrow) / aL4) %>% round(.,3) #aL4 was computed in the previous step, it has de number of sequences detected in the experiment


#Also check how many sequences pass the Kolmogorov-Smirnov test
head(bDataN2[[1]])
sapply(bDataN2, function(y) sum(y$KS_pval > 0.001))
df$KS_05 <- sapply(bDataN2, function(y) sum(y$KS_pval > 0.05))
df$KS_01 <- sapply(bDataN2, function(y) sum(y$KS_pval > 0.01))
df$KS_001 <- sapply(bDataN2, function(y) sum(y$KS_pval > 0.001))
#df$KS_fr <- sapply(bDataN2, function(y) sum(y$KS_pval > 0.01) / nrow(y)) %>% round(.,3)
#Add also what happen when filtering sequences first
yD <- lapply(bDataN2, function(y) y[y$filter2,]) #filtering criteria 2
df$KS_01_2 <- sapply(yD, function(y) sum(y$KS_pval > 0.01))
#df$KS_fr2 <- sapply(yD, function(y) sum(y$KS_pval > 0.001) / nrow(y)) %>% round(.,3)

#Plots
#Sequences with estimates
ggplot(df, aes(x = sample, y = n_seqs)) +
  labs(title="Sequences with estimated parameters", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = n_seqs, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = n_seqs_fr)) +
  labs(title="Fraction of sequences with estimated parameters", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = n_seqs_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

#Sequences passing filtering criteria
ggplot(df, aes(x = sample, y = filter)) +
  labs(title="Sequences passing filtering criteria, absolute numbers", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_fill_brewer(palette="Dark2") + #or "Spectral"
  scale_colour_brewer(palette="Dark2") +
  geom_text(aes(label = filter, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5)

ggplot(df, aes(x = sample, y = filter_fr)) +
  labs(title="Fraction of sequences passing filtering criteria", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_fill_brewer(palette="Dark2") + #or "Spectral"
  scale_colour_brewer(palette="Dark2") +
  geom_text(aes(label = filter_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5)

#Reproducible values
ggplot(df, aes(x = sample, y = repro2)) +
  labs(title="Reproducible values, absolute numbers", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = repro2, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = repro2_fr)) +
  labs(title="Fraction of reproducible values", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = repro2_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

#Sequences passing the KS test
aL <- list()
aL[[1]] <- ggplot(df, aes(x = sample, y = KS_05)) +
  labs(title="Well fit sequences, KS test, padj > 0.05", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7, show.legend=FALSE) +
  geom_text(aes(label = KS_05, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 2.5) +
  ylim(0, 1000) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

aL[[2]] <- ggplot(df, aes(x = sample, y = KS_01)) +
  labs(title="Well fit sequences, KS test, padj > 0.01", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7, show.legend=FALSE) +
  geom_text(aes(label = KS_01, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 2.5) +
  ylim(0, 1000) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

aL[[3]] <- ggplot(df, aes(x = sample, y = KS_001)) +
  labs(title="Well fit sequences, KS test, padj > 0.001", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = KS_001, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 2.5) +
  ylim(0, 1000) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

wrap_plots(aL)

#Fraction of sequences passing KS test
#ggplot(df, aes(x = sample, y = KS_fr)) +
#  labs(title="Fraction well fit sequences, KS test, padj > 0.001", x ="Sample", y = "Fraction of sequences") +
#  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
#  geom_text(aes(label = KS_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
#  scale_fill_manual(values = unique(cols)) +
#  scale_color_manual(values = unique(cols))

#KS after filtering
ggplot(df, aes(x = sample, y = KS_01_2)) +
  labs(title="Well fit sequences, KS test, padj > 0.01", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = KS_01_2, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))


#Also check the convergence values (not required indeed)
head(bDataN2[[1]])
df$convergence <- sapply(bDataN2, function(y) sum(y$convergence == 1))
ggplot(df, aes(x = sample, y = convergence)) +
  labs(title="Sequences with convergence = 1", x ="Sample", y = "Number sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = convergence, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5)
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))


#Add to df, how many sequences have a precision score > 1
yD <- bDataN2 #values without filtering
#yD <- lapply(bDataN2, function(y) y[y$filter2,]) #Using filtering criteria 2
#yD <- lapply(bDataN2, function(y) y[y$filter2 & y$pair_rep,]) #Using filtering criteria 2 and reproducibe values
sapply(yD, nrow)

head(yD[[1]])
df$mu_prec <- sapply(yD, function(y) sum(y$mu_prec2))
df$mu_prec_fr <- (sapply(yD, function(y) sum(y$mu_prec2))/sapply(yD, nrow)) %>% round(., 3)
df$r_prec <- sapply(yD, function(y) sum(y$r_prec2))
df$r_prec_fr <- (sapply(yD, function(y) sum(y$r_prec2))/sapply(yD, nrow)) %>% round(., 3)
df$full_prec <- sapply(yD, function(y) sum(y$full_prec2))
df$full_prec_fr <- (sapply(yD, function(y) sum(y$full_prec2))/sapply(yD, nrow)) %>% round(., 3)


#Plots
ggplot(df, aes(x = sample, y = mu_prec)) +
  labs(title="mu precision score, sequences accepted", x ="Sample", y = "Number of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = mu_prec, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = mu_prec_fr)) +
  labs(title="mu precision score, fraccion of sequences accepted", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = mu_prec_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = r_prec)) +
  labs(title="r precision score, sequences accepted", x ="Sample", y = "Number of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = r_prec, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = r_prec_fr)) +
  labs(title="r precision score, fraccion of sequences accepted", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = r_prec_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = full_prec)) +
  labs(title="Intersect mu and r precision score, sequences accepted", x ="Sample", y = "Number of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = full_prec, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = full_prec_fr)) +
  labs(title="Intersect mu and r precision score, fraccion of sequences accepted", x ="Sample", y = "Fraction of sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = full_prec_fr, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))


#Add how many sequences pass the mu threshold, based on the bimodal distribution of mean expression values
#With and without filtering by precision score
df$mu_bimod1 <- sapply(bDataN2, function(y) sum(log2(y$exp) > y$mu_bimodal_thres))
df$mu_bimod2 <- sapply(bDataN2, function(y) sum(log2(y$exp) > y$mu_bimodal_thres & y$full_prec2))

ggplot(df, aes(x = sample, y = mu_bimod1)) +
  labs(title="Mu > thres., individual threshold for condition", x ="Sample", y = "Number of accepted sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = mu_bimod1, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = mu_bimod2)) +
  labs(title="Intercept Mu > thres. & prec. score > 1", x ="Sample", y = "Number of accepted sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = mu_bimod2, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))


#Add how many sequences pass the mu threshold, arbitrary threshold
#With and without filtering by precision score
df$mu_cloud1 <- sapply(bDataN2, function(y) sum(y$cloud == "cl_high"))
df$mu_cloud2 <- sapply(bDataN2, function(y) sum(y$cloud == "cl_high" & y$full_prec2))

ggplot(df, aes(x = sample, y = mu_cloud1)) +
  labs(title="Mu > 2^-4", x ="Sample", y = "Number of accepted sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = mu_cloud1, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))

ggplot(df, aes(x = sample, y = mu_bimod2)) +
  labs(title="Intercept Mu > 2^-4 & prec. score > 1", x ="Sample", y = "Number of accepted sequences") +
  geom_bar(aes(color = lib_condition, fill = lib_condition), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = mu_cloud2, group = lib), position = position_dodge(0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = unique(cols)) +
  scale_color_manual(values = unique(cols))


####################
#Scatterplots with the bursting parameters and NB parameters
head(bDataN2[[1]])

#Copy data frames as they are or apply different filters (optional)
yD <- bDataN2 #no filter
#yD <- lapply(bDataN2, function(y) y[y$filter2 & y$pair_rep,])
#yD <- lapply(bDataN2, function(y) y[y$full_prec2,])

samples <- names(bDataN2)
samples2 <- gsub("_E38.+_rep", "_rep", samples) #Example for E38 experiment, change accordingly
libs <- gsub("_E38.+", "", samples) #Example for E38 experiment, change accordingly
contition <- gsub("_E38.+", "", samples) #Extra information, mainly used for complex datasets or so. Here it's the same as liraries. Adjust accordingly
lib_condition <- gsub("_E38.+", "", samples) #Extra information, mainly used for complex datasets or so. Here it's the same as liraries. Adjust accordingly
reps <- gsub(".+_rep", "rep", samples) #Extra information, mainly used for complex datasets or so. Here it's the same as liraries. Adjust accordingly

for (i in 1:length(yD)){
  y <- yD[[i]]
  y$seq_ID <- rownames(y)
  y$sample <- samples[i]
  y$sample2 <- samples2[i]
  y$library <- libs[i]
  y$lib_condition <- libs_condition[i]
  y$condition <- condition[i]
  y$rep <- reps[i]
  y$log10_mu <- log10(y$exp)
  y$log10_r <- log10(y$r)
  y$log10_bfreq <- log10(y$bfreq)
  y$log10_bsize <- log10(y$bsize)
  y$var <- y$exp + (y$exp^2/y$r)
  y$log10_var <- log10(y$var)
  yD[[i]] <- y
}
yD <- Reduce(rbind, yD)
head(yD)

#Color values, adjust accordingly to the number of samples and conditions
a <- scales::hue_pal()(4)
names(a) <- unique(libs)
yD$cols <- a[yD$library]
cols <- rep(a, each=3)

#General scatter plot with mu and r
head(yD)
ggplot(yD, aes(x=log10_mu, y=log10_r, color=library)) +
  labs(x="log10 mu", y = "log10 r", color="Library", title="NB parameters") +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(. ~ library)


#Plots based on precision score acceptance
#All points
ggplot(yD, aes(x=log10_mu, y=log10_r, color=full_prec2)) +
  labs(x="log10 mu", y = "log10 r", title="NB parameters") +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(rep ~ lib_condition) +
  scale_color_manual(values = c("gray", "black"))

#Remaining points
f <- yD$full_prec2
ax <- range(yD$log10_mu, na.rm=TRUE)
ay <- range(yD$log10_r, na.rm=TRUE)
ggplot(yD[f,], aes(x=log10_mu, y=log10_r)) +
  labs(x="log10 mu", y = "log10 r", title="NB parameters") +
  xlim(ax) +
  ylim(ay) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(rep ~ lib_condition) +
  scale_color_manual(values = c("black"))


#Plots with the mu cloud
#All points
ggplot(yD, aes(x=log10_mu, y=log10_r, color=log10_mu > mu_bimodal_thres)) +
  labs(x="log10 mu", y = "log10 r", title="NB parameters", color="mu > thres.") +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(rep ~ lib_condition) +
  scale_color_manual(values = c("gray", "black"))

#Remaining points
f <- yD$full_prec2
ax <- range(yD$log10_mu, na.rm=TRUE)
ay <- range(yD$log10_r, na.rm=TRUE)
ggplot(yD[f,], aes(x=log10_mu, y=log10_r, color=log10_mu > mu_bimodal_thres)) +
  labs(x="log10 mu", y = "log10 r", title="Estimates with prec. score > 1 \n & mu > threshold", color="mu > thres.") +
  xlim(ax) +
  ylim(ay) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(rep ~ lib_condition) +
  scale_color_manual(values = c("gray", "black"))

