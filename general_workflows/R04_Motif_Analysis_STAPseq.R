#Script for motif analysis for BARe-seq experiment prepared with the promoter libraries
#It contains the fold change estimation between motif-containing and motif-non-containing groups
#First part: analysis for one dataset
#Second part: analysis for multiple datasets


#setwd() #if required

#library(scModels)
library(dplyr)
library(numDeriv)

library(ggplot2)
library(ggpubr)
library(patchwork)

source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required

#######################################
#Analysis for one dataset. In the second part there is the code for analyzing several datasets, with individual and pooled samples

#Import data
x <- read.csv("E28_RT_test_processing1_counts.txt",sep="\t",header=TRUE) #Example with E28 count file
#Adapt the previous file name accordingly

head(x)
dim(x)
unique(x$sample)

#Import also processed data
dataP <- readRDS("E28_NB_dataframes.rds") #Example with E28 count file, adapt accordingly


#######################################
#Motif analyses with fold changes
#Import motif annotation
motifP <- read.csv("libFL008_anno.txt",sep="\t",header=TRUE) #used for lib008
dim(motifP)
MO <- colnames(motifP)[25:43] #motif names


#Arrange the data
#Samples, shorter names (samples2 vector), and condition
samples <- names(dataP)
samples2 <- gsub("FL", "",samples) %>% gsub("_E.+_r", "_r", .) 
condition <- paste0("RT", 1:3) %>% rep(., each=3)

yP2 <- list()

for (i in 1:length(yP)){
  y <- yP[[i]]
  #y$threshold <- thres[i]
  y$sample2 <- samples2[i]
  y$condition <- condition[i]
  
  y2 <- left_join(y, motifP[,c("oligo_id", MO)], join_by(seq_ID == oligo_id))
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

#Select samples and filtering sequences
f <- 1:9 #Select all samples
f <- c(4:6) #Select only samples with RT2

#Filter sequences: precision score and mu > specific threshold
yD2 <- lapply(yP2[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
#Option with no filter
#yD2 <- Reduce(rbind, yP2[f])

#Make a data frame for the heatmaps. It should contain FC of the medians of presence / absence of the motif groups
#Function to analyze one parameter for all motifs. Use the analyze_parameter function from ../main_functions/R01_Opitmization_Functions.R
a <- c("sample", "lib", "condition")
am <- analyze_parameter(cp_data=yD2, param_col="exp", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="mu")
af <- analyze_parameter(cp_data=yD2, param_col="bfreq", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="b_freq")
as <- analyze_parameter(cp_data=yD2, param_col="bsize", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="b_size")


#Combine all data frames to make only one plot
ab2 <- rbind(af, as, am)


#Optional: select motifs of interest, in case of willing to remove non-informative motifs
#Example:
#a <- c("INR", "TATAbox", "DPE", "MTE", "Ebox", "DRE", "Ohler1", "Ohler6", "Ohler7", "TCT", "DCE3", "DCE2", "DCE1", "BREd", "BREu") %>% rev()
#head(ab2)
#ab2 <- ab2[ab2$motif %in% a,]
#ab2$motif <- factor(ab2$motif, levels=a)

ggplot(ab2, aes(samples2, motif, fill=log2_fc, label=wilcox_pval)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "log2 FC", title= "Motif analyses") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  geom_text(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  facet_grid(. ~ param)



###########################
#Boxplots for individual motifs and parameters

#Repeated code: select samples and filtering sequences if required
#Select samples and filtering sequences
f <- 1:9 #Select all samples
#f <- c(4:6) #Select only samples with RT2

#Filter sequences: precision score and mu > especific threshold
yD2 <- lapply(yP2[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
#Option with no filter
#yD2 <- Reduce(rbind, yP2[f])

aL <- lapply(MO, function(y) {
  a <- c(y, paste0("no", y))
  a2 <- unique(yD2$sample)
  a3 <- paste(rep(a2, each=length(a)), rep(a, times=length(a2)), sep=":")
  
  yD2$group <- paste(yD2$sample, yD2[,y], sep=":")
  yD2$group <- factor(yD2$group, levels=a3)
  
  df <- data.frame(bfreq=yD2[,"bfreq"], log2_bf=yD2[,"log2_bfreq"], bsize=yD2[,"bsize"],
                   log2_bs=yD2[,"log2_bsize"], exp=yD2[,"exp"], log2_mu=yD2[,"log2_mu"], group=yD2$group,
                   motif=y)
  return(df)
})

names(aL) <- MO
length(aL)
sapply(aL, dim)
head(aL[[1]])


#Setting up colors
#Colors for 3 libraries with 3 samples each
cols <- c(RColorBrewer::brewer.pal(n = 3, name = "Blues")[2:3] %>% rep(., times=3),
         RColorBrewer::brewer.pal(n = 3, name = "Oranges")[2:3] %>% rep(., times=3),
         RColorBrewer::brewer.pal(n = 5, name = "PuRd")[4:5] %>% rep(., times=3) )
#Colors for 1 library with 3 samples
#cols <- c(RColorBrewer::brewer.pal(n = 6, name = "Blues")[c(4,6)] %>% rep(., times=3) )
#Colors for 3 pooled libraries:
#cols <- scales::hue_pal()(3) %>% rep(., each=2)


pFL <- lapply(aL, function(yd){
  n <- unique(yd$motif)
  
  #p <- ggplot(data = yd, aes(x = group, y = bfreq, color=group)) +
  p <- ggplot(data = yd, aes(x = group, y = log10_bf, color=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0, size=.5, color="gray") +
    scale_color_manual(values = cols) +
    xlab("") + ylab("log10 bursting frequency") + ggtitle(paste0(n, ", bursting frequency")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) #+
  #tryCatch({ stat_compare_means( comparisons = a2 , size = 3, method="wilcox.test",tip.length = 0.01,
  #                               bracket.size = 0.1) },
  #         error = function(err) {NA} )
  return(p)
})


pSL <- lapply(aL, function(yd){
  n <- unique(yd$motif)
  
  #p <- ggplot(data = yd, aes(x = group, y = bfreq, color=group)) +
  p <- ggplot(data = yd, aes(x = group, y = log10_bs, color=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0, size=.5, color="gray") +
    scale_color_manual(values = cols) +
    xlab("") + ylab("log10 bursting size") + ggtitle(paste0(n, ", bursting size")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) #+
  #tryCatch({ stat_compare_means( comparisons = a2 , size = 3, method="wilcox.test",tip.length = 0.01,
  #                               bracket.size = 0.1) },
  #         error = function(err) {NA} )
  return(p)
})

pEL <- lapply(aL, function(yd){
  n <- unique(yd$motif)
  
  #p <- ggplot(data = yd, aes(x = group, y = bfreq, color=group)) +
  p <- ggplot(data = yd, aes(x = group, y = log10_mu, color=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0, size=.5, color="gray") +
    scale_color_manual(values = cols) +
    xlab("") + ylab("log10 mu") + ggtitle(paste0(n, ", mean expression")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) #+
  #tryCatch({ stat_compare_means( comparisons = a2 , size = 3, method="wilcox.test",tip.length = 0.01,
  #                               bracket.size = 0.1) },
  #         error = function(err) {NA} )
  return(p)
})

names(pFL) <- names(pSL) <- names(pEL) <- MO

length(pFL)
sapply(pFL, class)
#Plot bursting frequency values for the first 2 motifs
wrap_plots(pFL[1:2])

#Or plot all parameter values for a given motif
names(pFL)
i <- 2 #index for the motif name, in this case TATAbox
wrap_plots(pFL[[i]] + theme(legend.position="none"), pSL[[i]] + theme(legend.position="none"), pEL[[i]])



#######################################
#Second part: analysis of several datasets, with the option to focus on individual, pooled samples or all of them

#Import data of interest
#STAPseq - promoter library
#Import all the samples we need:
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
yP <- dataP
length(yP)
head(yP[[1]])

samples <- names(dataP)
libs <- rep.int(c("lib008", "lib016", "lib017"), c(6, 3, 3)) #option for the 4 sets of samples (E28 and E38)
#libs <- c("lib008", "lib016", "lib017") %>% rep(., each=3) #option for the 3 sets of samples (only E38)

#Add individual mu thresholds (in log2 scale, STAPseq data, E38 library): lib_008: -3.5; lib_016: -4.1, lib_017: -3.7
#Threshold for E28 experiments (lib008): -3.5
thres <- c(-3.5, -3.5, -4.1, -3.7) %>% rep(., each=3) #option for the 4 sets of samples (E28 and E38)
#thres <- c(-3.5, -4.1, -3.7) %>% rep(., each=3) #option for the 3 sets of samples (only E38)

for (i in 1:length(yP)){
  y <- yP[[i]]
  y$seq_ID <- rownames(y)
  y$sample <- samples[i]
  y$library <- libs[i]
  #y$condition <- conditions[i]
  y$log2_mu <- log2(y$exp)
  y$log2_r <- log2(y$r)
  y$log2_bfreq <- log2(y$bfreq)
  y$log2_bsize <- log2(y$bsize)
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

###
#Motif analyses with fold changes
#Repeated code for importing motif annotation
motifP <- read.csv("libFL008_anno.txt",sep="\t",header=TRUE) #used for lib008
dim(motifP)
MO <- colnames(motifP)[25:43] #motif names


#Samples, shorter name (samples2 vector), and condition
samples <- names(dataP)
samples2 <- gsub("FL", "",samples) %>% gsub("_E.+_r", "_r", .) 
condition <- c("lib008", "lib008", "lib016", "lib017") %>% rep(., each=3) #option for the 4 sets of samples (E28 and E38)
#condition <- c("lib008", "lib016", "lib017") %>% rep(., each=3) #option for the 3 sets of samples (only E38)

yP2 <- list()

for (i in 1:length(yP)){
  y <- yP[[i]]
  #y$threshold <- thres[i]
  y$sample2 <- samples2[i]
  y$condition <- condition[i]
  
  y2 <- left_join(y, motifP[,c("oligo_id", MO)], join_by(seq_ID == oligo_id))
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

#Select samples and filtering sequences
f <- 1:length(yP2) #Select everything
#f <- c(3,6,9) #Select pooled samples: all E38 samples
#f <- c(3,9,12) #Select pooled samples: E28 for lib008, and E38 for lib016, 017
#f <- c(1,2,4,5,7,8) #Select individual samples

#Filtering
yD2 <- lapply(yP2[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
dim(yD2)
unique(yD2$sample)
unique(yD2$sample2)

#Make a data frame for the heatmaps. It should contain FC of the medians of presence / absence of the motif groups
#Function to analyze one parameter for all motifs. Use the analyze_parameter function from ../main_functions/R01_Opitmization_Functions.R
a <- c("sample", "lib", "condition")
am <- analyze_parameter(cp_data=yD2, param_col="exp", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="mu")
af <- analyze_parameter(cp_data=yD2, param_col="bfreq", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="b_freq")
as <- analyze_parameter(cp_data=yD2, param_col="bsize", qr_motifs=MO, qr_sample="sample2", add_sample=a, param_name="b_size")

#Function for one parameter, i.e. bursting frequency
#ggplot(af, aes(motif, library, fill=log2_fc, label=wilcox_pval)) +
ggplot(af, aes(samples2, motif, fill=log2_fc, label=wilcox_pval)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "log2 FC", title="Bursting frequency") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  geom_text(size=10) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8))

#Combine all data frames to make only one plot
ab2 <- rbind(af, as, am)
ab2$samples2 <- factor(ab2$samples2, levels=samples2)

#Optional: selecting motifs of interest
a <- c("INR", "TATAbox", "DPE", "MTE", "Ebox", "DRE", "Ohler1", "Ohler6", "Ohler7", "TCT", "DCE3", "DCE2", "DCE1", "BREd", "BREu") %>% rev()
head(ab2)
ab2 <- ab2[ab2$motif %in% a,]
ab2$motif <- factor(ab2$motif, levels=a)

#Some extra parameter for lines (to see better the split between samples):
#a <- seq(from=2.5, to=6, by=2)


#Different options to plot for:
#ggplot(ab2, aes(motif, library, fill=log2_fc, label=wilcox_pval)) +
#ggplot(ab2, aes(library, motif, fill=log2_fc, label=wilcox_pval)) +
ggplot(ab2, aes(samples2, motif, fill=log2_fc, label=wilcox_pval)) +
#ggplot(ab2, aes(motif, samples2, fill=log2_fc, label=wilcox_pval)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "log2 FC", title= "Motif analyses") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  #geom_text(size=10) +
  geom_text(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  #facet_grid(param ~.) #+
  facet_grid(. ~ param)
#geom_hline(yintercept=a, color = "black", size=0.25)
#geom_hline(yintercept=20, linetype="dashed", color = "red", size=2)


##
#Boxplots for individual motifs and parameters

#Select samples and filtering sequences (repeated code)
f <- 1:length(yP2) #Select everything
#f <- c(3,6,9) #Select pooled samples: all E38 samples
#f <- c(3,9,12) #Select pooled samples: E28 for lib008, and E38 for lib016, 017
#f <- c(1,2,4,5,7,8) #Select individual samples

#Filtering
yD2 <- lapply(yP2[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
dim(yD2)
unique(yD2$sample)
unique(yD2$sample2)
head(yD2)

aL <- lapply(MO, function(y) {
  a <- c(y, paste0("no", y))
  a2 <- unique(yD2$sample)
  a3 <- paste(rep(a2, each=length(a)), rep(a, times=length(a2)), sep=":")
  
  yD2$group <- paste(yD2$sample, yD2[,y], sep=":")
  yD2$group <- factor(yD2$group, levels=a3)
  
  df <- data.frame(bfreq=yD2[,"bfreq"], log2_bf=yD2[,"log2_bfreq"], bsize=yD2[,"bsize"],
                   log2_bs=yD2[,"log2_bsize"], exp=yD2[,"exp"], log2_mu=yD2[,"log2_mu"], group=yD2$group,
                   motif=y)
  return(df)
})

names(aL) <- MO
length(aL)
sapply(aL, dim)
head(aL[[1]])

#Arranging colors
#Colors for 3 pooled libraries:
cols <- scales::hue_pal()(3) %>% rep(., each=2)
#For 3 libraries with 3 samples each
#cols <- c(RColorBrewer::brewer.pal(n = 3, name = "Blues")[2:3] %>% rep(., times=3),
#         RColorBrewer::brewer.pal(n = 3, name = "Oranges")[2:3] %>% rep(., times=3),
#         #RColorBrewer::brewer.pal(n = 5, name = "Purples")[c(2,4)] %>% rep(., times=3) )
#         RColorBrewer::brewer.pal(n = 5, name = "PuRd")[4:5] %>% rep(., times=3) )
#For 1 library with 3 samples
#cols <- c(RColorBrewer::brewer.pal(n = 6, name = "Blues")[c(4,6)] %>% rep(., times=3) )


pFL <- lapply(aL, function(yd){
  n <- unique(yd$motif)
  
  #p <- ggplot(data = yd, aes(x = group, y = bfreq, color=group)) +
  p <- ggplot(data = yd, aes(x = group, y = log2_bf, color=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0, size=.5, color="gray") +
    scale_color_manual(values = cols) +
    xlab("") + ylab("log2 bursting frequency") + ggtitle(paste0(n, ", bursting frequency")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) #+
  #tryCatch({ stat_compare_means( comparisons = a2 , size = 3, method="wilcox.test",tip.length = 0.01,
  #                               bracket.size = 0.1) },
  #         error = function(err) {NA} )
  return(p)
})


pSL <- lapply(aL, function(yd){
  n <- unique(yd$motif)
  
  #p <- ggplot(data = yd, aes(x = group, y = bfreq, color=group)) +
  p <- ggplot(data = yd, aes(x = group, y = log2_bs, color=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0, size=.5, color="gray") +
    scale_color_manual(values = cols) +
    xlab("") + ylab("log2 bursting size") + ggtitle(paste0(n, ", bursting size")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) #+
  #tryCatch({ stat_compare_means( comparisons = a2 , size = 3, method="wilcox.test",tip.length = 0.01,
  #                               bracket.size = 0.1) },
  #         error = function(err) {NA} )
  return(p)
})

pEL <- lapply(aL, function(yd){
  n <- unique(yd$motif)
  
  #p <- ggplot(data = yd, aes(x = group, y = bfreq, color=group)) +
  p <- ggplot(data = yd, aes(x = group, y = log2_mu, color=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0, size=.5, color="gray") +
    scale_color_manual(values = cols) +
    xlab("") + ylab("log2 mu") + ggtitle(paste0(n, ", mean expression")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) #+
  #tryCatch({ stat_compare_means( comparisons = a2 , size = 3, method="wilcox.test",tip.length = 0.01,
  #                               bracket.size = 0.1) },
  #         error = function(err) {NA} )
  return(p)
})

names(pFL) <- names(pSL) <- names(pEL) <- MO

length(pFL)
sapply(pFL, class)
wrap_plots(pFL[1:2]) #plot changes in bursting frequency of the first two motifs

#Combine changes in parameters for one motif
names(pFL)
i <- 3
wrap_plots(pFL[[i]] + theme(legend.position="none"), pSL[[i]] + theme(legend.position="none"), pEL[[i]])



