library(scModels)
library(dplyr)
library(numDeriv)
library(cowplot)
library(LearnBayes)

library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggpmisc)

source("E12_FunctionsContourplot_v2.R")


#Making contour plots with Simulated data 3.2
#And the axis are in log10 scale, instead of log2

#Option 1
#Contour plots by using different data frames provided
#Import list of data frames (simulated data):

aL2 <- readRDS("SimData_03_2_DFs_Contour_v2_NBPar_FC_14.rds")
names(aL2)

df_nb <- aL2$df_nb #data frame with the quantities arranged in bins and ready for plotting with ggplot functions
#df_mu_r <- aL2$df_mu_r # data frame with individual quantities for each simulated sequence. This was used to create the previous df_np data frame. This is not used for this section of code
l_breaks_FC_14 <- aL2$l_breaks_FC_14 #list with two data frames: one with the breaks for the contour plot, the other with labels of the intervals and the colors to use
#l_breaks_FC_14 contains the breaks estimated for the avg_fold_change_14 score


#Plot function (only the contours):
geom_contour(data=dfw_nb, aes(x = x_mid, y = y_mid, z=avg_fold_change_14), breaks=l_breaks_FC_14$df$bk_log2) +
  geom_contour_filled(data=dfw_nb, aes(x = x_mid, y = y_mid, z=avg_fold_change_14), breaks=l_breaks_FC_14$df$bk_log2) +
  labs(title=n2, x ="log10 mu", y = "log10 r") +
  #geom_point(data=y2, aes(x=log10_mu, y=log10_r), color=alpha("black", 0.4), size=1.5, shape=16, show.legend=FALSE) +
  scale_fill_manual(values = l_breaks_FC_14$df2$cols, labels=l_breaks_FC_14$df2$labs_log2) +
  theme_minimal()

#To overlay the real data, use the line starting with "geom_point(...)"


###################
#Option 2: adding the fraction of points into each of the regions
#It requires the data frames of the simulated data contained in "SimData_03_2_DFs_Contour_v2_NBPar_FC_14.rds" file

#Import list of data frames (repeated code):
aL2 <- readRDS("SimData_03_2_DFs_Contour_v2_NBPar_FC_14.rds")
names(aL2)
df_nb <- aL2$df_nb #data frame with the quantities arranged in bins and ready for plotting with ggplot functions
df_mu_r <- aL2$df_mu_r # data frame with individual quantities for each simulated sequence. This was used to create the previous df_np data frame.
l_breaks_FC_14 <- aL2$l_breaks_FC_14 #list with two data frames: one with the breaks for the contour plot, the other with labels of the intervals and the colors to use
#l_breaks_FC_14 contains the breaks estimated for the avg_fold_change_14 score
df_breaks <- l_breaks_FC_14$df

#Estimate the fraction of real data points contained into in each region of the contour
#You require yD2, a data frame with the real values. And values should be previously selected (i.e. passing different filtering criteria, coming from specific samples, libraries, etc.)
#The arguments for the real data are in the last line: yD2: data frame, param_ed="log10_mu", param_ed_y="log10_r" are the colnames in yD2
#Use the fraction_real_data function from the E12_FunctionsContourplot.R script
df_fr <- fraction_real_data(df_metrics=df_mu_r, param_x="log10_smu", param_y="log10_sr", bins = 20,
                            df_avg_metric=df_nb, z="avg_fold_change_14", df_breaks=df_breaks, scale="log2",
                            df_expdat=yD2, param_ed="log10_mu", param_ed_y="log10_r")

df_fr <- df_fr[,c(1,3)]
colnames(df_fr) <- c("bin", "percent")
df_fr[,2] <- strsplit(df_fr[,2], " ") %>% sapply(., function(y) y[1])
df_fr_tb <- tibble(x = -7.5, y = 10, df_fr = list(df_fr))

#Plot function:
ggplot() +
  geom_contour(data=df_nb, aes(x = x_mid, y = y_mid, z=avg_fold_change_14), breaks=l_breaks_FC_14$df$bk_log2) +
  geom_contour_filled(data=df_nb, aes(x = x_mid, y = y_mid, z=avg_fold_change_14), breaks=l_breaks_FC_14$df$bk_log2) +
  labs(title="avg_fold_change_14", x ="log10 mu", y = "log10 r") +
  geom_point(data=yD2, aes(x=x=log10_mu, y=log10_r), color=alpha("black", 0.4), size=1.5, shape=16, show.legend=FALSE) +
  scale_fill_manual(values = l_breaks_FC_14$df2$cols, labels=l_breaks_FC_14$df2$labs_log2) +
  theme_minimal() +
  annotate(geom="table", x=7.5, y=5, label=df_fr)
