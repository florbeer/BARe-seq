#Script for motif analysis for BARe-seq experiment prepared with the promoter libraries
#It contains:
#Estimation of fold changes

#setwd() #if required

library(scModels)
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



