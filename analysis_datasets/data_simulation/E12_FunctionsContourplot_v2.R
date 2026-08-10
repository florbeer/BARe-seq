#Bursting dynamics. Alex Stark Collab.
#Functions for making the contour plot and handling the data quickly
#setwd("/home/gruengroup/reyna/varID02/21Stark_Collab/E28_STAPseq")
list.files()

library(scModels)
library(dplyr)
library(numDeriv)
#library(tidyverse)
library(cowplot)
library(LearnBayes)
library(pheatmap)

library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggpmisc)


#Quick function to run the loess with na values and remove them
#We can also use the argument na.action="na.exclude" in the loess function, but the point of this is to compare the predicted values
#If they change or not
#The default is na.action="na.omit", and I would prefer to use it, since I used it for the first versions of the plots
loess_fNA <- function(formula, data, predicted_var, ...){
  predicted <- loess(formula, data=data, ...) %>% predict()
  pred_correction <- rep(NA, nrow(data))
  if (!is.null(predicted_var) ){
    f <- !is.na(predicted_var)
  } else {
    f <- !is.na(data[[predicted_var]])
  }
  pred_correction[f] <- predicted
  return(pred_correction)
}


##########################################
#Function for computing different distance metrics
#Arguments:
#data: data frame with columns of parameter values.
#param1_gt and pram2_gt are colnames in the data frame containing the ground truth values of the parameters to evaluate
#param1_est and param2_est, are colnames in the data frame containing the estimated values of the parameters to evaluate
#param1_loess and param2_loess are colnames in the data frame of the parameters to use in the loess regression. This gives
#flexibility to estimate the predictor variable in one side with linear or log2 scale, and to regress out in linear
#or log2 scales. This is just for testing which option would better to use
#metric: choose between: 'fold_change_111', 'fold_change_112', 'fold_change_12', 'fold_change_13', 'fold_change_14', 'percentage_error', or 'precision_score'
#The differences between the foldchange methods are in the order of arranging the data. The most preferred one can be fold_change_112 or fold_change_14
#If choosing metric == "precision_score"
#comb_pred_scr indicates which column in the data frame has the combined precision score
#distance: it indicates another column in the data frame with a different metric to be used for loess regression
#... used for additional arguments of the loess function. I recommend to test different values for span argument (default: 0.7)

#Output:
#Vector with the predicted values of the loess regression
compute_distance_loess <- function(data, param1_gt=NULL, param2_gt=NULL, param1_est=NULL, param2_est=NULL, param1_loess, param2_loess, metric, comb_prec_scr=NULL, distance=NULL, ...) {
  require(dplyr)
  
  #Apply one of the different distance metrics and local prediction with loess function
  if (is.null(distance)){
    #if (metric == "fold_change_111") {
    if (grepl("fold_change", metric)) {
      #Estimate fold changes
      data[["param1_FC"]] <- data[[param1_est]] / data[[param1_gt]]
      data[["param2_FC"]] <- data[[param2_est]] / data[[param2_gt]]
      
      if (metric == "fold_change_111"){
        #Fold change 1.1.1:
        #Fold change per parameter -> combine (geometric average) -> loess -> prediction
        comb_FC <- data[,c("param1_FC", "param2_FC")] %>% apply(., 1, function(x) exp(mean(log(x), na.rm=TRUE) ))
        predicted <- loess_fNA(comb_FC ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var=comb_FC, ...)
      }
      else if (metric == "fold_change_112") {
        #Option 1.1.2 (the favorite one)
        #Fold change per parameter -> combine (geometric average) -> log2 -> loess -> prediction
        log2FC <- data[,c("param1_FC", "param2_FC")] %>% apply(., 1, function(x) exp(mean(log(x), na.rm=TRUE) )) %>% log2()
        predicted <- loess_fNA(log2FC ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= log2FC, ...)
      }
      else if (metric == "fold_change_12") {
        #Option 1.2
        #Fold change per parameter -> combine (geometric average) -> log2 -> abs -> loess -> prediction
        abslog2FC <- data[,c("param1_FC", "param2_FC")] %>% apply(., 1, function(x) exp(mean(log(x), na.rm=TRUE) )) %>% log2() %>% abs()
        predicted <- loess_fNA(abslog2FC ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= abslog2FC, ...)
      }
      else if (metric == "fold_change_13") {
        #Option 1.3
        #Fold change per parameter -> log2 -> abs -> combine (geometric average) -> loess -> prediction
        abslog2FC <- data[,c("param1_FC", "param2_FC")] %>% log2() %>% abs() %>% apply(., 1, function(x) exp(mean(log(x), na.rm=TRUE) ))
        predicted <- loess_fNA(abslog2FC ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= abslog2FC , ...)
      }
      else if (metric == "fold_change_14") {
        #Option 1.4
        #This is the second alternative for avoiding regions of over- and under estimation
        #Fold change per parameter -> combine (geometric average) -> log2 -> loess -> prediction -> abs
        log2FC <- data[,c("param1_FC", "param2_FC")] %>% apply(., 1, function(x) exp(mean(log(x), na.rm=TRUE) )) %>% log2()
        predicted <- loess_fNA(log2FC ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= log2FC , ...) %>% abs()
      }
      else {
        stop("Invalid metric version for 'fold_change'. Choose from '111', '112', '12', '13', '14'.")
      }
    }
    
    else if (metric == "percentage_error") {
      #percentage error per parameter -> abs -> combine (geometric average) -> log2 -> loess -> prediction
      data[["param1_perr"]] <- abs((data[[param1_gt]] - data[[param1_est]]) / data[[param1_gt]]) * 100
      data[["param2_perr"]] <- abs((data[[param2_gt]] - data[[param2_est]]) / data[[param2_gt]]) * 100
      
      comb_perr <- data[,c("param1_perr", "param2_perr")] %>% apply(., 1, function(x) exp(mean(log(x), na.rm=TRUE) )) %>% log2()
      predicted <- loess_fNA(comb_perr ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= comb_perr, ...)
    }
    else if (metric == "precision_score") {
      #precision score per parameter -> combine (geometric average) -> log2 -> loess -> prediction
      #Usually I have the precision score already estimated, so use it directly:
      log2ps <- log2(data[[comb_prec_scr]])
      predicted <- loess_fNA(log2ps ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= log2ps, ...)
    } 
    else {
      stop("Invalid metric. Choose from 'fold_change', 'percentage_error', or 'precision_score'.")
    }
  }
  else {
    #Use a combined metric already provided:
    combined_metric <- data[[distance]]
    predicted <- loess_fNA(combined_metric ~ data[[param1_loess]] + data[[param2_loess]], data=data, predicted_var= combined_metric, ...)
  }
  
  # Return the predicted values from the loess fit
  return(predicted)
}

####################################
#Function for arranging the data, in order to make the contour plot
#Arguments
#data: data frame with the parameters to plot
#param_x and param_y characters, which are colnames in the data frame that will be binned and plot in the x and y axis, respectively
#metric_cols: character vector with the colnames in the data frame, that will be averaged based on the binning of parameters x and y
#bins: number of bins to split the parameters into
prepare_contour_data <- function(data, param_x, param_y, metric_cols, bins = 20) {
  require(dplyr)
  
  # Bin the x parameter
  x_range <- range(data[[param_x]], na.rm = TRUE)
  x_range <- c(floor(min(x_range)), ceiling(max(x_range)))
  x_bks <- seq(x_range[1], x_range[2], length.out = bins + 1)
  #x_range <- range(data[[param_x]], na.rm = TRUE) + c(-0.5, 0.5)
  #x_range <- round(x_range, 0)
  x_bins <- cut(data[[param_x]], breaks = x_bks, labels = FALSE)
  #Correct values in case of having NA's.
  #We already tried to correct for them with the + or - 0.5, but still not good for certain values
  #It's a very long code, trying to cover both scenarios: when having only one NA value, or more than one
  if (anyNA(x_bins)){
    f <- which(is.na(x_bins))
    checks <- sapply(data[[param_x]][f], function(y) diff(c(y, x_range[1])) ) %>% abs()
    checks <- rbind(checks, abs(sapply(data[[param_x]][f], function(y) diff(c(y, x_range[2])) )) )
    f2 <- apply(checks, 2, function(y) which(y == min(y)))
    x_bins[f] <- ifelse(f2 == 1, 1, bins)
  }
  
  x_midpoints <- x_bks
  x_midpoints <- head(x_midpoints, -1) + diff(x_midpoints) / 2
  x_midpoints <- round(x_midpoints, 2)
  #Adjust the first and last elements in the midpoints vector. It allows to have all real values within the plotting area
  x_midpoints[1] <- x_range[1]
  x_midpoints[length(x_midpoints)] <- x_range[2]
  
  # Bin the y parameter
  #y_range <- range(data[[param_y]], na.rm = TRUE) + c(-0.5, 0.5)
  #y_range <- round(y_range, 0)
  y_range <- range(data[[param_y]], na.rm = TRUE)
  y_range <- c(floor(min(y_range)), ceiling(max(y_range)))
  y_bks <- seq(y_range[1], y_range[2], length.out = bins + 1)
  y_bins <- cut(data[[param_y]], breaks = y_bks, labels = FALSE)
  
  #Correct values in case of having NA's.
  if (anyNA(y_bins)){
    f <- which(is.na(y_bins))
    checks <- sapply(data[[param_y]][f], function(y) diff(c(y, y_range[1])) ) %>% abs()
    checks <- rbind(checks, abs(sapply(data[[param_y]][f], function(y) diff(c(y, y_range[2])) )) )
    f2 <- apply(checks, 2, function(y) which(y == min(y)))
    y_bins[f] <- ifelse(f2 == 1, 1, bins)
  }
  
  y_midpoints <- y_bks
  y_midpoints <- head(y_midpoints, -1) + diff(y_midpoints) / 2
  y_midpoints <- round(y_midpoints, 2)
  #Adjust the first and last elements in the midpoints vector
  y_midpoints[1] <- y_range[1]
  y_midpoints[length(y_midpoints)] <- y_range[2]
  
  # Add bins and bin combinations to the data frame
  data$x_bin <- x_bins
  data$y_bin <- y_bins
  #data$bin_combination <- paste(data$x_bin, data$y_bin, sep = ".")
  
  # Compute the average of the selected metrics for each bin combination
  aggregated_data <- data %>%
    group_by(x_bin, y_bin) %>%
    summarise(across(all_of(metric_cols), ~ mean(.x, na.rm = TRUE), .names = "avg_{col}"), .groups = "drop")
  
  # Add midpoints for plotting
  aggregated_data <- aggregated_data %>%
    mutate(
      x_mid = x_midpoints[x_bin],
      y_mid = y_midpoints[y_bin]
    )
  
  # Add the bin combination
  aggregated_data$bin_combination <- paste(aggregated_data$x_bin, aggregated_data$y_bin, sep = ".")
  #df_bks <- data.frame(x_breaks=x_bks, y_breaks=y_bks)
  #aL <- list(aggregated_data=aggregated_data, df_bks=df_bks)
  return(aggregated_data)
  #return(aL)
}


####################################################
#Function for defining the breaks in the contour plot
#Maybe it is good to define a criteria based on the selected metric
#Start checking the fold change options:

#Option 1.1.1
#Fold change per parameter -> combine (geometric average) -> loess -> prediction
#Best value around 1, values go from 0 to infinite. Linear
#Diverging palette

#Option 1.1.2
#this is the most used version from previous plots
#Fold change per parameter -> combine (geometric average) -> log2 -> loess -> prediction
#Best value around 0, values go from negative to positive. Log2
#Diverging palette

#Option 1.2
#Fold change per parameter -> combine (geometric average) -> log2 -> abs -> loess -> prediction
#Best value around 0, values go from 0 to infinite (ideally, but the loess gives some negative values). Log2
#Diverging palette

#Option 1.3
#Fold change per parameter -> log2 -> abs -> combine (geometric average) -> loess -> prediction
#Best value around 0, values go from 0 to infinite (ideally, but the loess gives some negative values). Log2
#Maybe gradual or sequential palette, when we have only positive values. Diverging palette when having negative values

#Option 1.4
#This is the second alternative for avoiding regions of over- and under estimation
#Fold change per parameter -> combine (geometric average) -> log2 -> loess -> prediction -> abs
#Best value around 0, values go from 0 to infinite. Log2
#Gradual or sequential palette

#Now let;s go to the function
#Arguments:
#z: vector to plot. It is usually one of the metrics used to determine accuracy, and averaged across bins
#scale: indicates if the z vector is in "linear", "log2" scale. Default is log2
#col_scale: choose between "diverging" or "gradual". default is "gradual"
#remove.one: logical argument, it indicates whether to remove 1 from the breaks or not. Useful for removing it from
#the fold changes-derived metrics, and to have a bin around 1, instead of 2 bins with values higher and
#lower than 1. Default is FALSE

prepare_breaks <- function(z, scale="log2", col_scale="gradual", remove.one=FALSE){
  require(dplyr)
  
  if (scale == "log2") {
    range_z <- range(z, na.rm=TRUE)
  } else if (scale == "linear") {
    range_z <- range(z, na.rm=TRUE) %>% log2()
  } else {
    stop("Invalid scale. Choose between 'log2' or 'linear'.")
  }
  
  #range_z <- range_z + c(-0.4, 0.4)
  #range_z <- round(range_z, 0)
  range_z <- c(floor(min(range_z)), ceiling(max(range_z)))
  bk_linear <- seq(from=range_z[1], to=range_z[2]) %>% 2^.
  #Conditional to remove the 1 value for the fold change metrics, but only if it is not the lower value
  if (bk_linear[1] != 1 & remove.one){
    f <- bk_linear != 1
    bk_linear <- bk_linear[f]
  }
  #Condition to add one or two points around 1, depending if the scale is negative or positive
  if (range_z[1] < log2(1/1.2)) {
    bk_linear <- c(bk_linear, 1/1.2, 1.2) %>% sort() #maybe consider to add 1.5 as well 
  }
  else {
    bk_linear <- c(bk_linear, 1.2) %>% sort()
  }
  
  
  #Checking if there is a bin without any observation. This cause some problems in the labeling when making the plot
  if (scale == "log2") {
    z_binned <- cut(z, breaks=log2(bk_linear), labels=FALSE)
  } else if (scale == "linear") {
    z_binned <- cut(z, breaks=bk_linear, labels=FALSE)
  } 
  an <- length(bk_linear)-1 #number of labels to be expected
  #Remove one of the breaks if not populated
  #Old code, only useful for one bin
  #if (all(1:an %in% unique(z_binned)) == FALSE) {
  #  a <- setdiff(1:an, unique(z_binned))
  #  if (a == 1) { #this is when the lower bin is not populated
  #    bk_linear <- bk_linear[-1] 
  #  } else if (a == length(bk_linear)) { #this is when the upper bin is not populated
  #    bk_linear <- bk_linear[-length(bk_linear)]
  #  } else { #When another bin in between is not populated, just mix it with the lower bin, let's see how it looks
  #    bk_linear <- bk_linear[-a]
  #  }
  #}
  if (all(1:an %in% unique(z_binned)) == FALSE) {
    a <- setdiff(1:an, unique(z_binned)) #vector with empty bins
    a2 <- setdiff(a, an) #this removes the last bin
    if (an %in% a) { #this is when the upper bin is not populated
      bk_linear <- bk_linear[-length(bk_linear)]
    }
    if (length(a2) > 0) { #all the other cases: when the first, or other bins in between are not populated, just remove the lower threshold
      bk_linear <- bk_linear[-a2]
    }
  }
  
  bk_log2 <- log2(bk_linear)
  labs_linear <- cut(bk_linear, breaks=bk_linear) %>% levels()
  labs_log2 <- cut(bk_log2, breaks=bk_log2) %>% levels()
  
  
  #Color palette
  if (col_scale == "diverging"){
    #Code considering an asymmetric number of bins, where the center is not exactly the medium value
    n <- length(labs_linear)
    f2 <- grep("1.2", labs_linear) #Define where is the center
    f2 <- f2[1] #use the first bin containing 1.2
    #bins in the upper direction
    n2 <- n-f2
    #bins in the lower direction
    n3 <- n-n2-1
    #number of bins for the color palette
    n4 <- max(c(n3, n2))
    n5 <- n4*2+1
    #Indexes for the color palette
    n6 <- n4+1 #center of the palette
    n7 <- n6-n3
    n8 <- n6+n2
    
    cols <- colorRampPalette(c("blue", "#FAF1FF","red"))(n5)[n7:n8]
  }
  else if (col_scale == "gradual") {
    n <- length(labs_linear)
    cols <- colorRampPalette(c("#FAF1FF", "blue"))(n)
  }
  else {
    stop("Invalid col_scale. Choose between 'divergent' or 'gradual'.")
  }
  
  aL <- list()
  aL$df <- data.frame(bk_linear=bk_linear, bk_log2=bk_log2)
  aL$df2 <- data.frame(labs_linear=labs_linear, labs_log2=labs_log2, cols=cols)
  return(aL)
}

######################
#Function to estimate the fraction of real data points into the regions in the contour plot

#df_metrics: data frame with individual values (no average), and all the info: parameter values and metric values
#param_x, param_y: colnames in df_metrics, it is a character indicating the parameters to be analysed, we need to bin them again
#bins: number of bins for param_x and param_y. Default is 20
#df_avg_metric: data frame or tibble with the metric values averaged across the different bins

#z: character (only one value), indicating the averaged metrics to bin, it should be included in the colnames in the df_avg_metric data frame
#metrics: character (or character vector) indicating the averaged metrics to bin, arranged in columns in the df_avg_metric data frame #Thinking this would not be possible, to run it for several metrics at a time. It depends on the last function prepare_breaks() which requires a data frame design to be run for each metric
#df_breaks: data frame with the breaks to use for the averaged metrics. It is the first slot from the list given as output by the prepare_breaks() function
#scale: character indicating whether the averaged metric is in logarithmic or linear scale. Options: "log2" (the default), and "linear"
#df_expdat: data frame with the experimental data to be binned
#param_ed_x, param_ed_y: parameters for plotting in x and y, contained in the df_expdat data frame as colnames

fraction_real_data <- function(df_metrics, param_x, param_y, bins = 20, df_avg_metric, z, df_breaks, scale="log2", df_expdat, param_ed_x, param_ed_y){
  require (dplyr)
  
  #Breaks for the parameters param_x and param_y contained in df_metrics data frame
  x_range <- range(df_metrics[[param_x]], na.rm = TRUE)
  x_range <- c(floor(min(x_range)), ceiling(max(x_range)))
  x_bks <- seq(x_range[1], x_range[2], length.out = bins + 1)
  
  y_range <- range(df_metrics[[param_y]], na.rm = TRUE)
  y_range <- c(floor(min(y_range)), ceiling(max(y_range)))
  y_bks <- seq(y_range[1], y_range[2], length.out = bins + 1)
  
  #Then bin the experimental data with the breaks of the previous step
  df_expdat$param_x_bin <- cut(df_expdat[[param_ed_x]], breaks=x_bks, labels=FALSE)
  df_expdat$param_y_bin <- cut(df_expdat[[param_ed_y]], breaks=y_bks, labels=FALSE)
  df_expdat$bin_combination <- paste(df_expdat$param_x_bin, df_expdat$param_y_bin, sep=".")
  
  
  #Bin the average metric vector z from df_avg_metric tibble
  bk_linear<- df_breaks$bk_linear
  
  if (scale == "log2") {
    df_avg_metric$z_binned <- cut(df_avg_metric[[z]], breaks=log2(bk_linear), labels=FALSE)
  } else if (scale == "linear") {
    df_avg_metric$z_binned <- cut(df_avg_metric[[z]], breaks=bk_linear, labels=FALSE)
  }
  
  #Combine the bin of the z metric with the real data
  #The ideas is that z metric has two binning codes: the combined bin according to x and y parameters. And the bins based on the range of values of itself (z)
  #So one bin correspond to the others and this information is used for annotating z-metric bins into the real data
  df_expdat <- left_join(df_expdat, df_avg_metric[,c("bin_combination", "z_binned")], by="bin_combination")
  
  #Make another data frame with the relative frequencies of each region
  a <- unique(df_expdat[["z_binned"]]) %>% sort()
  a2 <- sapply(a, function(x) sum(df_expdat[["z_binned"]] == x & !is.na(df_expdat[["z_binned"]]))/nrow(df_expdat) )
  a2 <- c(a2, sum(is.na(df_expdat[["z_binned"]]))/nrow(df_expdat))
  df3 <- data.frame(z_binned = c(a,NA), freq=round(a2, 4))
  df3$percent <- paste(df3$freq*100, "%")
  
  #Nest task would be to add the labels, such as:
  #Let's see later how to make it
  #df3$mulab <- c(0.7,2.2,2.4,1.6,-0.2,-2.8, -4.5, -5.9, -7, NA)
  #df3$rlab <- c(8,7,4.5,-6,-7.2,-7.2, -8.3, -8.3, -8.3, NA)
  #Maybe with ggrepel::geom_text_repel or geom_labels_repel
  #Or some functions of the ggpp package, such as annotate. But I haven't explored it that much
  
  return(df3)
}
