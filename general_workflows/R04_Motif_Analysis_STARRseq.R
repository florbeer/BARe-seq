#Script for motif analysis for BARe-seq experiment prepared with the promoter libraries
#It contains:
#Estimation of fold changes

#setwd() #if required

#library(scModels)
library(dplyr)
library(numDeriv)

library(ggplot2)
library(ggpubr)
library(patchwork)

source("../main_functions/R01_Opitmization_Functions.R") #adapt the path if required

#Import data
x <- read.csv("E39_STARR_v4_processing1_counts.txt",sep="\t",header=TRUE) #Example with E39 v4 count file, adapt file name accordingly
head(x)
dim(x)
unique(x$sample)

#Import also processed data
dataX <- readRDS("E39_v4_NB_dataframes.rds") #Example with E39 v4 count file, adapt accordingly


#######################################
#Motif analyses with fold changes
#Option 1: import the original STARR_TWIST_metadata_v2_corrected.txt file from Franzi, and correct motif names
motif <- read.csv("STARR_TWIST_metadata_v2_corrected.txt",sep="\t",header=TRUE) #updated in April 2025

#Change column names from the motif annotation into motif names
colnames(motif) #indexes: 42:50
a <- colnames(motif)[42:50]
#Transformation given by Franzi:
TF_motifs <- list(GATA=data.frame(Motif="GATA", ID="flyfactorsurvey__srp_SANGER_5_FBgn0003507"),
                  AP1=data.frame(Motif="AP1", ID="flyfactorsurvey__kay_Jra_SANGER_5_FBgn0001291"),
                  twist=data.frame(Motif="twist", ID="flyfactorsurvey__twi_da_SANGER_5_FBgn0000413"),
                  Trl=data.frame(Motif="Trl", ID="flyfactorsurvey__Trl_FlyReg_FBgn0013263"),
                  ETS=data.frame(Motif="ETS", ID="flyfactorsurvey__Ets97D_SANGER_10_FBgn0004510"),
                  SREBP=data.frame(Motif="SREBP", ID="flyfactorsurvey__HLH106_SANGER_10_FBgn0015234"),
                  # hk
                  Dref=data.frame(Motif="Dref", ID="homer__AVYTATCGATAD_DREF"),
                  Ohler1=data.frame(Motif="Ohler1", ID="homer__MYGGTCACACTG_Unknown1"),
                  Ohler6=data.frame(Motif="Ohler6", ID="homer__AAAAATACCRMA_Unknown4"))
#Simplifying the list:
TF_motifs <- Reduce(rbind, TF_motifs)
#Check if colnames in motif are in order, if so, we can just directly replace them
identical(colnames(motif)[42:50], TF_motifs[,2])

f <- match(colnames(motif)[42:50], TF_motifs[,2]) #matching elements from the firs vector into the second. So it gives: element 1 in v1 matches with 5 in v2, and so on
colnames(motif)[42:50] <- TF_motifs[f,1]
#Order the motif data frame in the order I have been using it
MO <- c("GATA", "AP1", "twist", "Trl", "ETS", "SREBP", "Dref", "Ohler1", "Ohler6")
sum(colnames(motif) %in% MO)
motif <- cbind(motif[,1:41], motif[,MO])


#Option 2: use the corrected file in rds format. Not need to modify anything else.
#motif <- readRDS("STARRseq_motif_annot_corrected.rds") #This is the corrected file Franzi sent, with the motif names in the right order. The original file comes from: STARR_TWIST_metadata_v2_corrected.txt, and still I need to edit colnames and so on
MO <- colnames(motif)[42:50]


#Arrange data
#Make a list with the motif information
yR <- dataX
#samples <- unique(x$sample)
samples <- names(dataX)
samples2 <- gsub("E39_polyA_", "", samples)

for (i in 1:length(yR)){
  y <- yR[[i]]
  y$sample2 <- samples2[i]
  if (is.null(y$condition)) {
    y$condition <- y$sample2
  }
  
  y2 <- left_join(y, motif[,c("ID_seq", MO)], join_by(seq_ID == ID_seq))
  #rownames(y2) <- y2$seq_ID
  
  for (j in MO){
    f <- y2[,j] > 0
    y2[,j] <- ifelse(f, j, paste0("no", j))
    a <- c(j, paste0("no", j)) #motif names
    y2[,j] <- factor(y2[,j], levels=a )
  }
  
  yR[[i]] <- y2
}

#Select samples of interest
f <- 1:length(yR) #option to select everything
#f <- 1:4 #option to select only samples from lib010

#Filter sequences: precision score and mu > especific threshold
yD2 <- lapply(yR[f], function(y) y[y$full_prec2 & log2(y$exp) > y$mu_bimodal_thres,]) %>% Reduce(rbind, .)
#Option with no filter
#yD2 <- Reduce(rbind, yR[f])


#Making data frames with a self-made function for making the heatmaps. Function found in ../R01_Optimization_Functions script
#The function selects one parameter at a time.
#It splits the data into groups with and without the indicated motifs across samples or groups of samples

af <- analyze_parameter(cp_data=yD2, param_col="bfreq", qr_motifs=MO, qr_sample="sample2", param_name="b_freq")
as <- analyze_parameter(cp_data=yD2, param_col="bsize", qr_motifs=MO, qr_sample="sample2", param_name="b_size")
am <- analyze_parameter(cp_data=yD2, param_col="exp", qr_motifs=MO, qr_sample="sample2", param_name="mu")

#Joining all the data frames
ab2 <- rbind(af, as, am)

#Heat map with ggplot2
a <- 4.5

ggplot(ab2, aes(motif, samples2, fill=log2_fc, label=wilcox_pval)) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "log2 FC", title= "Motif analyses") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446") +
  #geom_text(size=10) +
  geom_text(size=5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8)) +
  facet_grid(param ~.) +
  geom_hline(yintercept=a, color = "black", size=0.25)


###########################
#Boxplots for individual motifs and parameters

#Repeated code: select samples and filtering sequences if required
#Select samples and filtering sequences
f <- 1:8 #Select all samples

#Filter sequences: precision score and mu > especific threshold
yD2 <- lapply(yR[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
#Option with no filter
#yD2 <- Reduce(rbind, yR[f])

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
#Colors for 2 libraries with 4 samples each
cols <- c(RColorBrewer::brewer.pal(n = 3, name = "Blues")[2:3] %>% rep(., times=4),
         RColorBrewer::brewer.pal(n = 3, name = "Oranges")[2:3] %>% rep(., times=4) )
#Colors for 1 library with 4 samples
#cols <- c(RColorBrewer::brewer.pal(n = 6, name = "Blues")[c(4,6)] %>% rep(., times=4) )
#Colors for 2 pooled libraries:
#cols <- scales::hue_pal()(2) %>% rep(., each=2)


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
i <- 1 #index for the motif name, in this case GATA
wrap_plots(pFL[[i]] + theme(legend.position="none"), pSL[[i]] + theme(legend.position="none"), pEL[[i]])




#######################################
#Second part: analysis of several datasets, with the option to focus on individual, pooled samples or all of them

#Import data of interest
#STARRseq - enhancer libraries
#Import all the samples we need:
data_E24_lib010 <- readRDS("../E24_STARRseqE24_STARRseq_lib010_pools.rds")
data_E39_lib013 <- readRDS("../E39_STARRseq/E39_STARRseq_lib013_pools.rds")
data_E39_lib014 <- readRDS("../E39_STARRseq/E39_STARRseq_lib014_pools.rds")

#Make the analyses of all samples together
names(data_E24_lib010) #use rep1 & rep 2, lib010, and the pair
names(data_E39_lib013) #Pair1 (r1&2) & pair2 (r5&6) -> pool of 4
names(data_E39_lib014) #Triplet1 (reps 1,2,3) & triplet2 (reps 4,5,6)-> pool of 6

#Selecting the indicated samples
dataR <- c(data_E24_lib010[c(4,5,7)], data_E39_lib013[c(6,7, 9)], data_E39_lib014[c(10:12)])
names(dataR)[1:3] <- c("lib010_rep1", "lib010_rep2", "lib010_pair")
names(dataR)[4:6] <- paste("lib013", c("pair1", "pair2", "pool4"), sep="_")
names(dataR)[7:9] <- paste("lib014", c("triplet1", "triplet2", "pool6"), sep="_")


rm(data_E24_lib010, data_E39_lib013, data_E39_lib014)


#Arranging the data for the plots
names(dataR)
sapply(dataR, dim)
f <- colnames(dataR[[1]])
colnames(dataR[[4]])
yR <- lapply(dataR, function(y) y[,f])
length(yR)
head(yR[[1]])
sapply(yR, dim)

samples <- names(dataR)
libs <- c("lib010", "lib013", "lib014") %>% rep(., each=3)
#Add individual mu thresholds (in log2 scale, STARRseq data): lib010: -2.8; lib013: -2.6, lib014: -2.4
thres <- c(-2.8, -2.6, -2.4) %>% rep(., each=3)

for (i in 1:length(yR)){
  y <- yR[[i]]
  y$seq_ID <- rownames(y)
  y$sample <- samples[i]
  y$lib <- libs[i]
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
  
  #Size of the confidence intervals ?? maybe later. We can just plot the precision score
  y$threshold <- thres[i]
  
  yR[[i]] <- y
}
head(yR[[1]])

###
#Motif analyses with fold changes
#Repeated code for importing motif annotation. Importing just the corrected file in rds format
motif <- readRDS("STARRseq_motif_annot_corrected.rds") #This is the corrected file Franzi sent, with the motif names in the right order. The original file comes from: STARR_TWIST_metadata_v2_corrected.txt, and still I need to edit colnames and so on
MO <- colnames(motif)[42:50]


#Samples, shorter name (samples2 vector), and condition
samples <- names(dataR)
samples2 <- gsub("FL", "",samples) %>% gsub("_E.+_r", "_r", .) 
condition <- c("lib010", "lib013", "lib014") %>% rep(., each=3)

yR2 <- list()

for (i in 1:length(yR)){
  y <- yR[[i]]
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
  
  yR2[[i]] <- y2
}
names(yR2) <- names(yR)
head(yR2[[1]])

#Select samples and filtering sequences
f <- 1:length(yR2) #Select everything
#f <- c(3,6,9) #Select pooled samples
#f <- c(1,2,4,5,7,8) #Select individual samples

#Filtering
yD2 <- lapply(yR2[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
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


########
#Boxplots for individual motifs and parameters

#Select samples and filtering sequences (repeated code)
f <- 1:length(yR2) #Select everything
#f <- c(3,6,9) #Select pooled samples
#f <- c(1,2,4,5,7,8) #Select individual samples

#Filtering
yD2 <- lapply(yR2[f], function(y) y[y$full_prec2 & y$log2_mu > y$threshold,]) %>% Reduce(rbind, .)
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



