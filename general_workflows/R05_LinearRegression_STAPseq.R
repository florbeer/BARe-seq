#Script for comparing different approaches for motif analysis for BARe-seq experiments prepared with the promoter libraries
#It contains analysis of all the datasets at a time, selecting either pooled or individual sequences of the main experiments E28 and E38
#The linear regression is performed with all the selected motifs included in one model, and not one regression per motif.

#Motif analysis strategies
#Part 1: fold change estimation, comparison of group-containing and group-non containing sequences
#Part 2: conventional linear regression with only linear terms (without crossed terms)
#Part 3: Lasso regression with linear and interaction terms
#Part 4: Lasso regression with only linear terms (without crossed terms)

#There are complementary analysis not included in this script, but in other scripts:
#Conventional linear regression with crossed terms -> to be prepared
#Comparison between elastic, lasso and Ridge regression -> to be prepared

#setwd() #if required

#library(scModels)
library(dplyr)
library(numDeriv)

library(ggplot2)
library(ggpubr)
library(patchwork)

library("glmnet")
#library("caret")

source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required
source("../main_functions/R02_LinearReg_Functions.R") #adapt the path if required


######
#Import data of interest
data_E38v3 <- readRDS("../E38_v3_STAPseq/E38_NB_dRT2_pools.rds") #E38 experiment v3, with the replicates of 6, lib008
data_E28 <- readRDS("../E28_STAPseq/E28_NB_dataframes.rds") #use E28 experiment RT2, lib008, instead of E38. It just doesn't have pooled samples
data_E28_p2 <- readRDS("../E28_STAPseq/E28_STAPseq_RT2_pools.rds") #these are the pooled samples of E28
data_E38v1_16 <- readRDS("../E38_STAPseq/E38_STAPseq_lib016_pools.rds")
data_E38v1_17 <- readRDS("../E38_STAPseq/E38_STAPseq_lib017_pools.rds")

#Make the analyses of all samples together
names(data_E38v3) #Let's use the pairs 1&2, and the pool4
names(data_E28) #Use RT2, replicates 2&3
names(data_E28_p2) #Use "RT2_pair"
names(data_E38v1_16) #Let's use rep 1&2, and pair
names(data_E38v1_17) #let's use rep 5&6 and pair


#Selecting the indicated samples (4 sets of samples)
dataP <- c(data_E28[5:6], data_E28_p2[6], data_E38v3[c(7,8,12)], data_E38v1_16[c(1,2,6)], data_E38v1_17[c(2,3,6)])
names(dataP)[1:3] <- c("E28_lib008_rep2", "E28_lib008_rep3", "E28_lib008_pair")
names(dataP)[4:6] <- c("E38_lib008_pair1", "E38_lib008_pair2", "E38_lib008_pool4")
names(dataP)[9] <- c("E38_lib016_pair")
names(dataP)[12] <- c("E38_lib017_pair")

#Option selecting only samples from experiment E38 (excluding E28)
#dataP <- c(data_E38v3[c(7,8,12)], data_E38v1_16[c(1,2,6)], data_E38v1_17[c(2,3,6)])
#names(dataP)[1:3] <- c("lib008_pair1", "lib008_pair2", "lib008_pool4")
#names(dataP)[6] <- c("lib016_pair")
#names(dataP)[9] <- c("lib017_pair")

rm(data_E38v3, data_E38v1_16, data_E38v1_17, data_E28, data_E28_p2)

#Arranging the data for the plots
names(dataP)
sapply(dataP, dim)
yP <- dataP
length(yP)
head(yP[[1]])

samples <- names(dataP)
libs <- rep.int(c("lib008", "lib016", "lib017"), c(6, 3, 3)) #option for the 4 sets of samples (E28 and E38)
#libs <- c("lib008", "lib016", "lib017") %>% rep(., each=3) #option for the 3 sets of samples (only E38)

#Add individual mu thresholds (in log2 scale, STAPseq data, E38 library): lib_008: -2.7; lib_016: -3.94, lib_017: -3.46
#Threshold for E28 experiments (lib008): -3.3
thres <- c(-3.3, -2.7, -3.94, -3.46) %>% rep(., each=3) #option for the 4 sets of samples (E28 and E38)
#thres <- c(-2.7, -3.94, -3.46) %>% rep(., each=3) #option for the 3 sets of samples (only E38)


for (i in 1:length(yP)){
  y <- yP[[i]]
  y$seq_ID <- rownames(y)
  y$sample <- samples[i]
  y$library <- libs[i]
  #y$condition <- conditions[i]
  y$log2_mean <- log2(y$exp)
  y$log2_r <- log2(y$r)
  y$log2_bf <- log2(y$bfreq)
  y$log2_bs <- log2(y$bsize)
  y$var <- y$exp + (y$exp^2/y$r)
  y$log2_var <- log2(y$var)
  #Quantities in log10
  y$log10_mu <- log10(y$exp)
  y$log10_r <- log10(y$r)
  y$log10_bfreq <- log10(y$bfreq)
  y$log10_bsize <- log10(y$bsize)

  y$threshold <- thres[i]
  
  yP[[i]] <- y
}
head(yP[[1]])
sapply(yP, dim)

a <- colnames(yP[[1]]) #to get common columns among the datasets. E28 has 2-3 columns less.
#sapply(yP, function(y) sum(a %in% colnames(y)))

#Filter samples, if required
f <- 1:12 #to keep everything
#f <- c(1:3,7:12) #to remove E38_lib008 from this analysis

#Select samples of interest (vector f) and columns found in E28 (vector a)
yD <- lapply(yP[f], function(y) y[,a]) %>% Reduce(rbind, .) %>% as.data.frame()
anyDuplicated(rownames(yD))
anyDuplicated(yD$seq_ID)
dim(yD)

##
#Filtering sequences.
#There are 4 options, run only one of them

#Option 1: filter by precision score and mu cloud
yD <- yD[yD$full_prec2 & yD$log2_mean > yD$threshold,]

#Option 2: filter only by precision score
#yD <- yD[yD$full_prec2,]

table(yD$full_prec2)
summary(yD$log2_mean)


#Option 3: filter by precision score and mu cloud and reproducible values across replicates
#Make a function that applies this for the individual samples.
#And for the pooled samples, get the sequences reported as reproducible in individual samples
a <- colnames(yP[[1]])
sapply(yP, function(y) table(y$pair_rep))
sapply(yP, function(y) anyNA(y$pair_rep))
yD <- lapply(yP, function(y) {
  if (anyNA(y$pair_rep)){
    y2 <- y[y$full_prec2 & y$log2_mean > y$threshold,a]
    #y2 <- y[,a]
  }
  else{
    y2 <- y[y$full_prec2 & y$log2_mean > y$threshold & y$pair_rep,a]
  }
  return(y2)
})
#Select sequences for the pooled samples
sapply(yD, dim)
f <- seq(from=1, to=12, by=3)
f2 <- seq(from=3, to=12, by=3)

aL <- lapply(yP[f], function(y) rownames(y[y$pair_rep,]))
for (i in 1:length(f2)){
  f3 <- f2[i]
  f4 <- aL[[i]]
  yD[[f3]] <- yD[[f3]][rownames(yD[[f3]]) %in% f4,]
}
sapply(yD, dim)
sapply(aL, length)

yD <- Reduce(rbind, yD) %>% as.data.frame()

#Option 4: filter by KS and mu threshold
yD <- yD[yD$KS_pval > 0.01 & yD$log2_mean > yD$threshold,]
dim(yD)


#############
#Import motif data
motif <- read.csv("libFL008_anno.txt",sep="\t",header=TRUE) #used for lib008
dim(motif)

#Define motif names
#MO <- colnames(motif)[25:43] #motif names (all of them)
#Optional: select a subset, corresponding to motifs of interest (used in my thesis)
#MO <- c("TATAbox", "INR", "Ebox", "DPE", "MTE", "DPE_extended", "DRE", "Ohler1", "Ohler7", "Ohler8")
MO <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT") #motifs in the figures

motif2 <- left_join(yD[,c("seq_ID", "sample", "library","bfreq", "bsize","exp", "log2_mean", "log2_bf", "log2_bs")],
                    motif[,c("oligo_id", "class2",MO)], by = join_by("seq_ID" == "oligo_id") )
head(motif2)
dim(motif2)
table(motif2$class2, useNA ="ifany")

##
#Quick exploration of the parameters
plot(motif2$exp, main="mean expression", ylab="mean")
plot(motif2$log2_mean, main="log2 mean expression", ylab="log2 mean")
plot(density(motif2$log2_mean), main="density log2 mean expression")

plot(motif2$bfreq, main="r disp. param", ylab="r")
plot(density(motif2$bfreq), main="density r disp. param")
plot(motif2$log2_bf, main="log2 r disp. param", ylab="log2 r")
plot(density(motif2$log2_bf), main="density log2 r disp. param")

plot(motif2$bsize, main="bursting size", ylab="bursting size")
plot(density(motif2$bsize), main="density bursting size")
plot(motif2$log2_bs, main="log2 bursting size", ylab="log2 bursting size")
plot(density(motif2$log2_bs), main="density bursting size")
#There is no batch effect, so we can proceed with the normal analysis with the values as they are


#Make a list of the motif2 data frame, with different combinations to test:
#Individual samples without putting them together
#2 individual samples, putting them together (not recommended)
#Pooled samples

a <- unique(motif2$sample)
unique(motif2$library)

f <- c(seq(from=1, to=length(a), by=3), seq(from=2, to=length(a), by=3)) %>% sort()
f2 <- a[f]
f3 <- matrix(f2, ncol=2, byrow=TRUE)

#Extracting individual samples without joining them together
aL_mo <- lapply(f2, function(y) motif2[motif2$sample == y,])
names(aL_mo) <- f2
length(aL_mo)
sapply(aL_mo, dim)

#Extracting individual samples, and joining the ones that belong to the same experiment
aL_mo <- c(aL_mo, apply(f3, 1, function(y) {motif2[motif2$sample %in% y,]} ) )
#names(aL_mo)[9:12] <- paste(c("E28_lib008", "E38_lib008", "E38_lib016", "E38_lib017"), "joined", sep="_")
names(aL_mo)[7:9] <- paste(c("E28_lib008", "E38_lib016", "E38_lib017"), "joined", sep="_")

#Now extracting the pooled samples
f4 <- setdiff(a, f2)
aL_mo <- c(aL_mo, lapply(f4, function(y) motif2[motif2$sample == y,]) )
#names(aL_mo)[13:16] <- f4
names(aL_mo)[10:12] <- f4
sapply(aL_mo, dim)


############################################
#Part 1
#Compute the fold changes, comparing group-containing and group-non containing sequences

dim(motif)
#MO <- colnames(motif)[25:43] #all motif names
#MO <- c("TATAbox", "INR", "Ebox", "DPE", "MTE", "DPE_extended", "DRE", "Ohler1", "Ohler7", "Ohler8") #motifs in my thesis
#MO <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT") #motifs in current figures

#Info to add to the data frames in yP
#Individual mu thresholds (in log2 scale): lib_008 (E38): -2.7; lib_016: -3.94, lib_017: -3.46
#thres <- c(-2.7, -3.94, -3.46) %>% rep(., each=3)
#When having 4 libraries: lib_008 (E28): -3.3; lib_008 (E38): -2.7; lib_016: -3.94, lib_017: -3.46
#thres <- c(-3.3, -2.7, -3.94, -3.46) %>% rep(., each=3)

#Samples, shorter name, and condition
samples <- names(dataP)
samples2 <- gsub("libFL", "E38_lib",samples) %>% gsub("_E.+_r", "_r", .)
#condition <- c("lib008", "lib016", "lib017") %>% rep(., each=3)
condition <- c("lib008", "lib008", "lib016", "lib017") %>% rep(., each=3)

yP2 <- list()

for (i in 1:length(yP)){
  y <- yP[[i]]
  #y$threshold <- thres[i]
  y$sample2 <- samples2[i]
  y$condition <- condition[i]
  
  y2 <- left_join(y, motif[,c("oligo_id", MO)], join_by(seq_ID == oligo_id))
  rownames(y2) <- y2$seq_ID
  
  for (j in MO){
    f <- y2[,j] > 0
    y2[,j] <- ifelse(f, j, paste0("no", j))
    a <- c(j, paste0("no", j)) #motif names
    y2[,j] <- factor(y2[,j], levels=a )
  }
  
  yP2[[i]] <- y2
}
names(yP2) <- names(yP)
head(yP2[[1]])

#Select samples of interest
#f <- c(3,6,9) #Select pooled samples: all E38 samples
f <- c(3,9,12) #Select pooled samples: E28 for lib008, and E38 for lib016, 017
#f <- c(3,6,9,12) #Select pooled samples: E28, and E38 all libraries
#f <- c(1,2,4,5,7,8) #Select individual samples

#Filtering. There are several options, select only one
#Option 1 Filter by precision score, and mu threshold
yD2 <- lapply(yP2[f], function(y) y[y$full_prec2 & y$log2_mean > y$threshold,]) %>% Reduce(rbind, .)
dim(yD2)
unique(yD2$sample)
unique(yD2$sample2)


#Option 2 Filter by precision score and mu cloud and reproducible values across replicates
#Make a function that applies this for the individual samples.
#And for the pooled samples, get the sequences reported as reproducible in individual samples
a <- colnames(yP2[[1]])
sapply(yP2, function(y) table(y$pair_rep))
sapply(yP2, function(y) anyNA(y$pair_rep))
yD2 <- lapply(yP2, function(y) {
  if (anyNA(y$pair_rep)){
    y2 <- y[y$full_prec2 & y$log2_mean > y$threshold,a]
    #y2 <- y[,a]
  }
  else{
    y2 <- y[y$full_prec2 & y$log2_mean > y$threshold & y$pair_rep,a]
  }
  return(y2)
})
#Select sequences for the pooled samples
sapply(yD2, dim)
f <- seq(from=1, to=12, by=3)
f2 <- seq(from=3, to=12, by=3)

aL <- lapply(yP2[f], function(y) rownames(y[y$pair_rep,]))
for (i in 1:length(f2)){
  f3 <- f2[i]
  f4 <- aL[[i]]
  yD2[[f3]] <- yD2[[f3]][rownames(yD2[[f3]]) %in% f4,]
}
sapply(yD2, dim)
sapply(aL, length)

f <- c(3,6,9,12) #Select pooled samples: E28, and E38 all libraries
yD2 <- Reduce(rbind, yD2[f]) %>% as.data.frame()

#Optiion 3 Filter by KS and mu threshold
#f <- c(3,6,9) #Select pooled samples: all E38 samples
f <- c(3,9,12) #Select pooled samples: E28 for lib008, and E38 for lib016, 017
#f <- c(3,6,9,12) #Select pooled samples: E28, and E38 all libraries
#f <- c(1,2,4,5,7,8) #Select individual samples
yD2 <- lapply(yP2[f], function(y) y[y$KS_pval > 0.01 & y$log2_mean > y$threshold,]) %>% Reduce(rbind, .)
dim(yD2)
unique(yD2$sample)
unique(yD2$sample2)



#Make a data frame for the heatmaps. It should contain FC of the medians of presence / absence of the motif groups
#Function to analyze one parameter for all motifs: analyze_parameter, located in ../main_functions/R01_Opitmization_Functions.R
a <- c("sample", "library", "condition")
am <- analyze_parameter(cp_data=yD2, param_col="exp", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="mu")
af <- analyze_parameter(cp_data=yD2, param_col="bfreq", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="b_freq")
as <- analyze_parameter(cp_data=yD2, param_col="bsize", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="b_size")

#Quick adjustment for matching with the plotting code
colnames(am)[4:5] <- colnames(af)[4:5] <- colnames(as)[4:5] <- c("samples2", "samples")

#Combine all data frames to make only one plot
ab2 <- rbind(af, as, am)
#ab2$samples2 <- factor(ab2$samples2, levels=samples2) #Check if required
#Some extra parameter of lines:
#a <- seq(from=2.5, to=6, by=2)

#Optional: selecting motifs of interest
#Motifs in the current figures:
#a <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT")
#head(ab2)
#ab2 <- ab2[ab2$motif %in% a,]
#ab2$motif <- factor(ab2$motif, levels=a)

#Make a list of heatmap plots for future comparisons between models
pL2 <- list()


pL2$motifFC <- ggplot(ab2, aes(samples2, motif, fill=log2_fc, label=wilcox_pval)) +
  #ggplot(ab2, aes(motif, samples2, fill=log2_fc, label=wilcox_pval)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "log2 FC", title= "Motif analysis, fold change") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  #geom_text(size=10) +
  geom_text(size=7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  #facet_grid(param ~.) #+
  facet_grid(. ~ param)
#geom_hline(yintercept=a, color = "black", size=0.25)
#geom_hline(yintercept=20, linetype="dashed", color = "red", size=2)



############################################
#Part 2
#Conventional linear regression
#Only linear terms, without any interaction (no crossed terms)

#Generic function for linear regression
#lreg <- lm(log2_mean ~ Ohler1 + DRE + TATAbox + INR + Ebox + Ohler6 + Ohler7 + Ohler8 + DPE + MTE + TCT + DPE_extended + DPE_Kadonaga + BREu + BREd + DCE1 + DCE2 + DCE3 + TC_17_Zabidi, data = motif2)
#summary(lreg)
#aL3 <- summary(lreg)[[4]] %>% as.data.frame() #extract coefficients
#head(aL3)
#dim(aL3)

#Running the regression for each data frame contained in aL_mo
sapply(aL_mo, dim)
#MO <- colnames(motif)[25:43] #motif names
#MO <- c("TATAbox", "INR", "Ebox", "DPE", "MTE", "DPE_extended", "DRE", "Ohler1", "Ohler7", "Ohler8") #motifs selected in my my thesis
#MO <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT") #Motifs selected in the figures

aL3 <- lapply(aL_mo, function(y){
  #Make a list that contains the regression for each parameter at a time: mean, bf and bs
  lreg <- list()
  #All motifs
  #lreg$mu <- lm(log2_mean ~ Ohler1 + DRE + TATAbox + INR + Ebox + Ohler6 + Ohler7 + Ohler8 + DPE + MTE + TCT + DPE_extended + DPE_Kadonaga + BREu + BREd + DCE1 + DCE2 + DCE3 + TC_17_Zabidi, data = y)
  #lreg$bf <- lm(log2_bf ~ Ohler1 + DRE + TATAbox + INR + Ebox + Ohler6 + Ohler7 + Ohler8 + DPE + MTE + TCT + DPE_extended + DPE_Kadonaga + BREu + BREd + DCE1 + DCE2 + DCE3 + TC_17_Zabidi, data = y)
  #lreg$bs <- lm(log2_bs ~ Ohler1 + DRE + TATAbox + INR + Ebox + Ohler6 + Ohler7 + Ohler8 + DPE + MTE + TCT + DPE_extended + DPE_Kadonaga + BREu + BREd + DCE1 + DCE2 + DCE3 + TC_17_Zabidi, data = y)
  #Only a selection of motifs (motifs selected in my my thesis)
  #lreg$mu <- lm(log2_mean ~ TATAbox + INR + Ebox + DPE + MTE + DPE_extended + DRE + Ohler1 + Ohler7 + Ohler8, data = y)
  #lreg$bf <- lm(log2_bf ~ TATAbox + INR + Ebox + DPE + MTE + DPE_extended + DRE + Ohler1 + Ohler7 + Ohler8, data = y)
  #lreg$bs <- lm(log2_bs ~ TATAbox + INR + Ebox + DPE + MTE + DPE_extended + DRE + Ohler1 + Ohler7 + Ohler8, data = y)
  #Motifs selected in the figures
  lreg$mu <- lm(log2_mean ~ DPE + INR + Ebox + Ohler6 + DRE + MTE + Ohler7 + Ohler1 + TATAbox + TCT, data = y)
  lreg$bf <- lm(log2_bf ~ DPE + INR + Ebox + Ohler6 + DRE + MTE + Ohler7 + Ohler1 + TATAbox + TCT, data = y)
  lreg$bs <- lm(log2_bs ~ DPE + INR + Ebox + Ohler6 + DRE + MTE + Ohler7 + Ohler1 + TATAbox + TCT, data = y)
  #Extract coefficients and arrange the data frames in the lists
  #a <- c("mu", "bf", "bs")
  #for (i in 1:length(lreg)){
  #  
  #}
  
  lreg2 <- lapply(lreg, function(ay){
    a <- summary(ay)[[4]] %>% as.data.frame() 
    a <- a[MO,]
    a$motif <- MO
    a$significance <- ifelse(a$'Pr(>|t|)' < 0.05, "*", NA)
    return(a)
  })
  
  #Add a column with the parameter
  a <- c("mu", "bf", "bs")
  for (i in 1:length(lreg2)){
    lreg2[[i]]$parameter <- a[i]
  }
  lreg3 <- Reduce(rbind, lreg2)
  return(lreg3)
})
head(aL3[[1]])

a3 <- names(aL3)
for(i in 1:length(aL3)){
  aL3[[i]]$sample <- a3[i]
}


#Select samples of interest
names(aL3)
#f <- 9:12 #select the joined samples
f <- 10:12 #select the pooled samples (analysis with E28, E38 libs 16 and 17)
#f <- 13:16 #select the pooled samples
#f <- 9:16 #Or select both pooled and joined


aL4 <- Reduce(rbind, aL3[f])
dim(aL4)
head(aL4)
unique(aL4$motif)
aL4$motif <- factor(aL4$motif, levels=MO)


#Save the plot into the pL2 list
#a2 <- seq(from=2.5, to=7, by=2) #extra parameter for plotting vertical lines, not used here
pL2$Lin_reg_ind <- ggplot(aL4, aes(sample, motif, fill=Estimate, label=significance)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Coefficients", title="Linear regression, individual motifs (no interaction)") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  #geom_vline(xintercept = c(a2)) +
  geom_text(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  facet_grid(. ~ parameter)



####################################################
#Part 3
#Lasso regression, using functions from ../main_functions/R02_LinearReg_Functions.R
#Regression for mean
#For one data frame of the list aL_mo:
a <- FUN_netr_modsel(df=aL_mo[[1]], resp_var="log2_mean", motifs=MO, alpha_test=1, square_err=FALSE)
names(a)
head(a[[1]])
head(a[[2]])

#Regression for each data frame in aL_mo
aL_tests <- lapply(aL_mo, function(y) {
  aL <- list()
  aL$df <- y #Copy the data frame to have it in the next permutation function
  aL$index <- FUN_netr_modsel(df=y, resp_var="log2_mean", motifs=MO, alpha_test=1, square_err=FALSE)$index #extract the index once
  
  #Extract alpha and lambda values into a different list
  aL2 <- list()
  aL2$mean <- FUN_netr_modsel(df=y, resp_var="log2_mean", motifs=MO, alpha_test=1, square_err=FALSE)$parameters
  aL2$bfreq <- FUN_netr_modsel(df=y, resp_var="log2_bf", motifs=MO, alpha_test=1, square_err=FALSE)$parameters
  aL2$bsize <- FUN_netr_modsel(df=y, resp_var="log2_bs", motifs=MO, alpha_test=1, square_err=FALSE)$parameters
  
  #Add the parameter names to the data frames
  a <- c("log2_mean", "log2_bf", "log2_bs")
  for(i in 1:length(aL2)){
    aL2[[i]]$par <- a[i]
  }
  #Simplify aL2 and add it to the main list aL
  aL$params <- Reduce(rbind, aL2)
  
  return(aL)
})
#So we have as an output a list with different slots
#First level the 16 tested data frames from aL_mo
#Second level: 3 slots
#df: Original data frame where the model is applied. It is important to keep it, otherwise, the indexes returned by FUN_netr_modsel does not make sense later
#index: indexes used for the validation model
#params: another data frame with the chosen alpha and lambda parameter values

length(aL_tests)
sapply(aL_tests, length)
sapply(aL_tests, names)
sapply(aL_tests, function(y) dim(y[[1]])) #copied data frames from aL_mo
sapply(aL_tests, function(y) length(y[[2]])) #indexes used as validation
sapply(aL_tests, function(y) dim(y[[3]])) #data frame with the chosen alpha and lambda parameter values
aL_tests[[1]][[3]]

#Permutation with FUN_per version 3, which has the data frame motif2 and the response variable as input, instead of the model
aL_perm <- lapply(aL_tests, function(aL_data){
  df_data <- aL_data[[1]] #data frame with dependent and predictive variables, copied from aL_mo
  f <- aL_data[[2]] #indexes of the training dataset
  pars <- aL_data[[3]] #data frame with parameters alpha and lambda
  
  aL_perm <- apply(pars, 1, function(y){
    ap <- y[2]
    lb <- y[3]
    resp <- y[4]
    y2 <- FUN_per(df=df_data, resp_var=resp, motifs=MO, n_per=1000, lambda=lb, alpha=ap, index=f)
    return(y2)
  })
  names(aL_perm) <- pars[,4]
  return(aL_perm)
})

length(aL_perm)
sapply(aL_perm, length)
#sapply(aL_perm, names)
sapply(aL_perm, function(y) sapply(y, class))
sapply(aL_perm, function(y) sapply(y, dim))
head(aL_perm[[1]][[1]])


#Check how many coefficients are zero or different from 0
#Write indexes in the order of each experiment, just for better comparison
f <- grep("E28", names(aL_perm))
#f2 <- grep("E38_lib008", names(aL_perm))
f3 <- grep("016", names(aL_perm))
f4 <- grep("017", names(aL_perm))
#f5 <- c(f, f2, f3, f4)
f5 <- c(f, f3, f4)

sapply(aL_perm[f5], function(y){
  #sapply(y, function(y2) sum(y2$coefficients == 0))
  sapply(y, function(y2) sum(y2$coefficients != 0))
})
#Check how many coefficients are significant
sapply(aL_perm[f5], function(y){
  sapply(y, function(y2) sum(y2$padj < 0.05))
})
#Also checking the number of sequences in each data frame
sapply(aL_mo[f5], dim)


#Since there are many samples where only one coefficient is significant, check if it is the same across
#arranged samples of the same experiment, in other words, check which coefficients are consistent
aL <- lapply(aL_perm, function(y){
  sapply(y, function(y2) rownames(y2[y2$padj < 0.05,]))
})
sapply(aL, function(y) sapply(y, length))

f <- grep("E28", names(aL))
f2 <- grep("E38_lib008", names(aL))
f3 <- grep("016", names(aL))
f4 <- grep("017", names(aL))

aL[f]
aL[f2]
aL[f3]
aL[f4]


###
#Heatmap with the significant coefficients, one sample at a time
#For the heatmaps with several samples, scroll down a little bit more
head(aL_perm[[1]][[1]], 30)
sapply(aL_perm, names)
identical(rownames(aL_perm[[1]][[3]]), rownames(aL_perm[[2]][[3]]) )
identical(rownames(aL_perm[[1]][[3]]), rownames(aL_perm[[3]][[3]]) )

#Select one sample
names(aL_perm)
ab <- 12 #Joined samples are indexes: 9:12, or 7:9 (when removing E38_lib008 from the analysis)
aL <- aL_perm[[ab]]
length(aL)
names(aL)
head(aL[[1]])

a <- c("mu", "bf", "bs")
for (i in 1:length(aL)){
  aL[[i]]$parameter <- a[i]
}

aL2 <- lapply(aL, function(y){
  a <- rownames(y) %>% strsplit(., ":")
  a <- sapply(a, function(y2){
    if (length(y2) == 1){
      y2 <- c(y2, "1_motif")
    }
    else {y2}
    return(y2)
  }) %>% t()
  #return(a)
  y$motif1 <- a[,1]
  y$motif2 <- a[,2]
  return(y)
})
sapply(aL2, dim)
head(aL2[[1]], 30)
aL2 <- Reduce(rbind, aL2)
dim(aL2)
head(aL2, 30)
tail(aL2, 10)

#Some arrangements for plotting
f <- aL2$padj > 0.05
aL2$coeff2 <- aL2$coefficients
aL2$coeff2[f] <- NA #For showing only the coefficients with significant value

aL2$motif1_par <- paste(aL2$motif1, aL2$parameter, sep="_") #For separating observations based on the parameters
#MO
a <- c("mu", "bf", "bs")
f2 <- paste(rep(c("(Intercept)", MO), each=length(a)), rep(a, times=length(MO)+1), sep="_")
aL2$motif1_par <- factor(aL2$motif1_par, levels=f2)
aL2$motif2 <- factor(aL2$motif2, levels=c("1_motif", MO))

#List for saving plots, in case of selecting different samples in the above line: `ab <- 12`
#pL <- list()
#pL$E28_lib008 <- ggplot(data = aL2, aes(x=motif2, y=motif1_par, fill=coeff2)) + #Change the order of x and y, so that the values go into the lower half of the matrix
#pL$E38_lib016 <- ggplot(data = aL2, aes(x=motif2, y=motif1_par, fill=coeff2)) +
pL$E38_lib017 <- ggplot(data = aL2, aes(x=motif2, y=motif1_par, fill=coeff2)) +
  geom_tile(color = "white") +
  labs(fill="coefficient") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", na.value="#d8d8d8") +
  geom_hline(yintercept = seq(from=3.5, to=56, by=3), colour='black') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust=0.8), axis.title.x = element_blank(), axis.title.y = element_blank())
wrap_plots(pL)


####
#Get the interactions that are relevant and make boxplots
#Select interactions that have a significant change, at least in one of the parameters
#Go back to the previous code just to select the model of interest (Lasso)
head(aL2)
a <- aL2[aL2$padj < 0.05, c("motif1", "motif2")]
dim(a)
a$motif_inter <- paste(a$motif1, a$motif2, sep=":")
a <- a[!grepl("1_motif",a$motif_inter),] #Remove coefficients with one motif
a <- a[!duplicated(a$motif_inter),] #Remove repeated interactions
a <- a[,1:2] #Remove the interaction column, to match with the next piece of code

#Use aL_mo list, which has subsets of the motif2 data frame
names(aL_mo)
#ab <- 12 #select a slot in aL_mo if required
head(aL_mo[[ab]])

aL5 <- aL_mo[[ab]]
aL3 <- apply(a, 1, function(y){
  #interaction(motif2[,y])
  interaction(aL5[,y])
  
})
class(aL3)
dim(aL3)
colnames(aL3) <- paste(a$motif1, a$motif2, sep=":")
aL3[1:5,1:13]
#aL3 <- data.frame(motif2[,c("log2_mean", "log2_bf", "log2_bs")], aL3)
aL3 <- data.frame(aL5[,c("log2_mean", "log2_bf", "log2_bs")], aL3)


for (i in 4:ncol(aL3)){
  y <- aL3[,i]
  aL3[,i] <- factor(y, levels=c("0.0", "1.0", "0.1", "1.1")) 
}

head(aL3)
class(aL3$Ohler1.TATAbox)
class(aL3$log2_mean)

#One plot at a time
ggplot(aL3, aes(x=aL3[,8], y=log2_mean)) +
  #geom_boxplot(fill=c('#0D748C', '#F9FFAF', '#0D748C', '#F9FFAF')) +
  geom_boxplot(outlier.shape = NA, color=c("#fd950d","#39ab00","#8c0d74","#0f87a3")) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_jitter(width = 0.25, height = 0, size=.5, color="#4c4c4c") +
  labs(title=colnames(aL3)[8],x="Interaction groups", y = "log2 mean") +
  theme_classic() +
  theme(legend.position="none")

#List of plots
a <- colnames(aL3)[4:ncol(aL3)]
a2 <- colnames(aL3)[1:3]
a3 <- c("log2 mu", "log2 bfreq", "log2 bsize")

aL4 <- lapply(a, function(ay){
  pL <- list()
  
  for (i in 1:length(a2)){
    pL[[i]] <- ggplot(aL3, aes(x=.data[[ay]], y=.data[[a2[i]]]) ) +
      #geom_boxplot(fill=c('#0D748C', '#F9FFAF', '#0D748C', '#F9FFAF')) +
      geom_boxplot(outlier.shape = NA, color=c("#fd950d","#39ab00","#8c0d74","#0f87a3")) +
      #geom_jitter(shape=16, position=position_jitter(0.2)) +
      geom_jitter(width = 0.25, height = 0, size=.5, color="#4c4c4c") +
      labs(title=paste(ay, a3[i], sep=",\n"), x="Interaction groups", y = a3[i]) +
      theme_minimal() +
      theme(legend.position="none")
  }
  return(pL)
})
length(aL4)
names(aL4) <- a
sapply(aL4[[1]], class)
#wrap_plots(aL4[1:12]) + plot_layout(ncol = 4)

#Select plots of interest, choose one option at a time:
f <- grep("TATA", a) #interactions containing TATAbox
f <- c(grep("Ohler1.DCE3", a), grep("INR.Ohler8", a), grep("INR.Ebox", a), grep("MTE", a)) #individual interactions
f <- grep("Ohler8", a) #interactions containing Ohler8
a[f]

#Reduce(c, aL4[f]) %>% wrap_plots(.) + plot_layout(ncol = 6)
Reduce(c, aL4[f]) %>% wrap_plots(.) + plot_layout(ncol = 3)

f2 <- f[c(1:3,5)]
Reduce(c, aL4[f2]) %>% wrap_plots(.) + plot_layout(ncol = 3)


#####
#Heatmaps with the significant coefficients, taking several samples

#Save the plot into the pL2 list, which compares different strategies for motif analysis

#For this, select samples of interest
names(aL_perm)
head(aL_perm[[1]][[1]], 30)
sapply(aL_perm, names)
identical(rownames(aL_perm[[1]][[3]]), rownames(aL_perm[[2]][[3]]) )
identical(rownames(aL_perm[[1]][[3]]), rownames(aL_perm[[3]][[3]]) )

#Select the pooled samples
f <- 10:12
#f <- 13:16
#Select the joined samples
#f <- 9:12
#Or select all of them together
#f <- 9:16
aL <- aL_perm[f]

#Extract coefficients of individual motifs
#MO <- colnames(motif)[25:43]
#MO <- c("TATAbox", "INR", "Ebox", "DPE", "MTE", "DPE_extended", "DRE", "Ohler1", "Ohler7", "Ohler8") #motifs in my thesis
MO <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT") #motifs in the figures
aL <- lapply(aL, function(y){
  a <- c("mu", "bf", "bs")
  for (i in 1:length(y)){
    y[[i]]$parameter <- a[i]
  }
  
  y2 <- lapply(y, function(ab) {
    a2 <- ab[MO,]
    a2$motif <- MO
    return(a2)
  }) %>% Reduce(rbind,.)
  return(y2)
})
sapply(aL, dim)
head(aL[[1]])

a3 <- names(aL)

for(i in 1:length(aL)){
  aL[[i]]$sample <- a3[i]
}
aL2 <- Reduce(rbind, aL)
aL2$motif <- factor(aL2$motif, levels=MO)
aL2$significance <- ifelse(aL2$padj < 0.05, "*", NA)
head(aL2)

#Plot
#Additional parameter for adding vertical lines, not used here
#a4 <- seq(from=2.5, to=7, by=2)
#a4 <- seq(from=0.5, to=9, by=2)


pL2$Lasso_inter <- ggplot(aL2, aes(sample, motif, fill=coefficients, label=significance)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Coeff.", title="Lasso regression, interaction terms") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  #geom_vline(xintercept = a4) +
  geom_text(size=7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  facet_grid(. ~ parameter)




####################################################
#Part 4
#Running Lasso regression for each data frame contained in aL_mo, and only the linear terms, without interaction terms
sapply(aL_mo, dim)

#Using different functions of the ../main_functions/R02_LinearReg_Functions.R script

#Regression for each data frame in aL_mo

#Selecting one set of motifs
#MO <- colnames(motif)[25:43] #all motif names
#MO <- c("TATAbox", "INR", "Ebox", "DPE", "MTE", "DPE_extended", "DRE", "Ohler1", "Ohler7", "Ohler8") #my thesis' motifs
#MO <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT") #motifs in the figures

aL_tests2 <- lapply(aL_mo, function(y) {
  aL <- list()
  aL$df <- y #Copy the data frame to have it in the next permutation function
  aL$index <- FUN_netr_modsel(df=y, resp_var="log2_mean", motifs=MO, alpha_test=1, crossed_terms=FALSE, square_err=FALSE)$index #extract the index vector once
  
  #Extract alpha and lambda values into a different list
  aL2 <- list()
  aL2$mean <- FUN_netr_modsel(df=y, resp_var="log2_mean", motifs=MO, alpha_test=1, crossed_terms=FALSE, square_err=FALSE)$parameters
  aL2$bfreq <- FUN_netr_modsel(df=y, resp_var="log2_bf", motifs=MO, alpha_test=1, crossed_terms=FALSE, square_err=FALSE)$parameters
  aL2$bsize <- FUN_netr_modsel(df=y, resp_var="log2_bs", motifs=MO, alpha_test=1, crossed_terms=FALSE, square_err=FALSE)$parameters
  
  #Add the parameter names to the data frames in aL2
  a <- c("log2_mean", "log2_bf", "log2_bs")
  for(i in 1:length(aL2)){
    aL2[[i]]$par <- a[i]
  }
  #Simplify aL2 and add it to the main list aL
  aL$params <- Reduce(rbind, aL2)
  
  return(aL)
})
#So we have as an output a list with different slots
#First level the 16 tested data frames from aL_mo
#Second level: 3 slots
#df: Original data frame where the model is applied. It is important to keep it, otherwise, the indexes returned by FUN_netr_modsel does not make sense later
#index: indexes used for the validation model
#params: another data frame with the chosen alpha and lambda parameter values
length(aL_tests2)
sapply(aL_tests2, length)
sapply(aL_tests2, names)
sapply(aL_tests, names) #comparison with the first version
sapply(aL_tests2, function(y) dim(y[[1]])) #copied data frames from aL_mo
sapply(aL_tests2, function(y) length(y[[2]])) #indexes used as validation
sapply(aL_tests2, function(y) dim(y[[3]])) #data frame with the chosen alpha and lambda parameter values
aL_tests2[[1]][[3]]
aL_tests[[1]][[3]]
#Compared to the full model (with crossed terms, there are not so much differences)

#Permutation with FUN_per version 3, which has the data frame with the motifs and the response variable as input, instead of the model
aL_perm2 <- lapply(aL_tests2, function(aL_data){
  df_data <- aL_data[[1]] #data frame with dependent and predictive variables, copied from aL_mo
  f <- aL_data[[2]] #indexes of the training dataset
  pars <- aL_data[[3]] #data frame with parameters alpha and lambda
  
  aL_perm <- apply(pars, 1, function(y){
    ap <- y[2]
    lb <- y[3]
    resp <- y[4]
    y2 <- FUN_per(df=df_data, resp_var=resp, motifs=MO, n_per=1000, lambda=lb, alpha=ap, index=f, crossed_terms=FALSE)
    return(y2)
  })
  names(aL_perm) <- pars[,4]
  return(aL_perm)
})

length(aL_perm2)
sapply(aL_perm2, length)
#sapply(aL_perm, names)
sapply(aL_perm2, function(y) sapply(y, class))
sapply(aL_perm2, function(y) sapply(y, dim))
head(aL_perm2[[1]][[1]])
head(aL_perm[[1]][[1]]) #comparison with the crossed-terms' model


#Check how many coefficients are zero or different from 0
#Write indexes in the order of each experiment
f <- grep("E28", names(aL_perm2))
#f2 <- grep("E38_lib008", names(aL_perm2))
f3 <- grep("016", names(aL_perm2))
f4 <- grep("017", names(aL_perm2))
#f5 <- c(f, f2, f3, f4)
f5 <- c(f, f3, f4)

sapply(aL_perm2[f5], function(y){
  #sapply(y, function(y2) sum(y2$coefficients == 0))
  sapply(y, function(y2) sum(y2$coefficients != 0))
})
#Check how many coefficients are significant
sapply(aL_perm2[f5], function(y){
  sapply(y, function(y2) sum(y2$padj < 0.05))
})
#Check the number of sequences in each data frame
sapply(aL_mo[f5], dim)


#Since there are many samples where only one coefficient is significant, check if it is the same across
#arranged samples of the same experiment, in other words, check which coefficients are consistent
aL3 <- lapply(aL_perm2, function(y){
  sapply(y, function(y2) rownames(y2[y2$padj < 0.05,]))
})
sapply(aL3, function(y) sapply(y, length))

f <- grep("E28", names(aL3))
f2 <- grep("E38_lib008", names(aL3))
f3 <- grep("016", names(aL3))
f4 <- grep("017", names(aL3))

aL3[f]
aL3[f2]
aL3[f3]
aL3[f4]

#Make a heatmaps
#Select samples of interest
names(aL_perm2)
head(aL_perm2[[1]][[1]], 30)
sapply(aL_perm2, names)
identical(rownames(aL_perm2[[1]][[3]]), rownames(aL_perm2[[2]][[3]]) )
identical(rownames(aL_perm2[[1]][[3]]), rownames(aL_perm2[[3]][[3]]) )

#Select the pooled samples
f <- 10:12
#f <- 13:16
#Select the joined samples
#f <- 9:12
#Or select all of them together
#f <- 9:16
#Select the pooled samples, but only E28_lib008, and not E38_lib008
#f <- c(13,15:16)
aL3 <- aL_perm2[f]
head(aL3[[1]])
#Extract coefficients of individual motifs
#MO <- colnames(motif)[25:43]
#MO <- c("TATAbox", "INR", "Ebox", "DPE", "MTE", "DPE_extended", "DRE", "Ohler1", "Ohler7", "Ohler8") #my thesis' motifs
MO <- c("DPE", "INR", "Ebox", "Ohler6", "DRE", "MTE", "Ohler7", "Ohler1", "TATAbox", "TCT") #motifs in the figures

head(aL3[[1]][[1]])
aL3[[1]][[1]][MO,]

aL3 <- lapply(aL3, function(y){
  a <- c("mu", "bf", "bs")
  for (i in 1:length(y)){
    y[[i]]$parameter <- a[i]
  }
  
  y2 <- lapply(y, function(ab) {
    a2 <- ab[MO,]
    a2$motif <- MO
    return(a2)
  }) %>% Reduce(rbind,.)
  return(y2)
})
sapply(aL3, dim)
head(aL3[[1]])

a3 <- names(aL3)

for(i in 1:length(aL3)){
  aL3[[i]]$sample <- a3[i]
}
aL4 <- Reduce(rbind, aL3)
aL4$motif <- factor(aL4$motif, levels=MO)
aL4$significance <- ifelse(aL4$padj < 0.05, "*", NA)
head(aL4)

#Plot
#Extra parameter for vertical lines, if desired
#a4 <- seq(from=2.5, to=7, by=2)
#a4 <- seq(from=0.5, to=9, by=2)
pL2$Lasso_ind <- ggplot(aL4, aes(sample, motif, fill=coefficients, label=significance)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Coeff.", title="Lasso regression, individual motifs (no interaction)") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  #geom_vline(xintercept = a4) +
  geom_text(size=7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  facet_grid(. ~ parameter)


############
#Plotting all the heatmaps from different models at once
names(pL2)
f <- c(1,2,4,3) #Arrange slots if desired
wrap_plots(pL2[f])
