#Setting up a threshold for mu values for each library
#STARR-seq experiment with the library with inserted motifs (library libFL007_STARR_DSCP)

library(dplyr)
library(readxl)

#Import parameter values of the different enhancer samples contained in Enhancer_libFL007_burst_parameters_NB.xls
samples <- excel_sheets("Enhancer_libFL007_burst_parameters_NB.xls") #adjust path if required
data_list <- lapply(samples, function(y) {
  y2 <- read_excel("Enhancer_libFL007_burst_parameters_NB.xls", sheet = y) %>%
    as.data.frame()
  rownames(y2) <- y2[,1]
  y2 <- y2[,2:ncol(y2)]
  return(y2)
})
names(data_list) <- samples


#Estimate threshold for each library
#Take the replicate samples and exclude the pooled samples
samples2 <- grep("rep", samples, value=TRUE)
data_list2 <- data_list[samples2]

#Take the mean expression values per library and transform them into log2 scale
libs <- gsub("_rep.+", "", samples2) %>% unique()
log2mu_list <- lapply(libs, function(y){
  f <- grep(y, names(data_list2))
  mu <- lapply(data_list2[f], function(y2) y2$exp) %>% Reduce(c, .)
  log2mu <- log2(mu)
  return(log2mu)
})

names(log2mu_list) <- libs


#Function to estimate the lower value of mean expression
#The intended function will take mean expression values and estimate kernel density with "density" function
#It will test different smoothing bandwidths (bw argument) and give a vector as output for the different bandwidth tested values
#Argument:
#log2_mean_expr: numerical vector of mean expression value sin log2 scale
#low and high: numerical value, indicating the range of values to screen for the getting the mean expression value with minimal density

finding_mean_thres <- function(log2_mean_expr, low=-4, high=-1.5){
  bw_test <- seq(from=0.1, to=0.8, by=0.05) #bandwidth values to test
  
  mu_thres <- sapply(bw_test, function(y){
    dens <- density(log2_mean_expr, bw=y)
    f <- dens$x > low & dens$x < high
    f2 <- which(dens$y == min(dens$y[f]))
    mu_thres <- dens$x[f2] %>% round(., 3)
    return(mu_thres)
  })
  names(mu_thres) <- paste("bw", bw_test, sep="_")
  return(mu_thres)
}


#Plot the density of mean expression values for each library, in order to have an estimate of the range of values to screen for the threshold
par(mfrow=c(2,2))
for (i in 1:length(log2mu_list)){
  plot(density(log2mu_list[[i]]), main=names(log2mu_list)[i])
}
dev.off()

#Apply the finding_mean_thres function to each library
finding_mean_thres(log2mu_list[[1]], low=-4, high=-1.5)

#Plot density with different bandwidths. This is just for visualization and helping to find a good bandwidth
#The chosen bandwidth depends on how smooth / spiky the plot looks like. It's suggested to choose a smooth line, but without over smoothing
i <- 1 #change to plot different libraries
plot(density(log2mu_list[[i]], bw=0.4), main="density mu") #test different bandwidths. 
abline(v=-3.7, col="red") #adjust value according to the output of finding_mean_thres

#Chosen mean expression thresholds (rounded values) and bandwidths per library:
#libFL007_STARR_DSCP			-3.2			0.4

#Go back to the corresponding datasets, and annotate which mena expression values are higher than the thresholds
#data_list <- lapply()
mean_thres <- c(-3.2) %>% rep(., each=3)
for (i in 1:length(data_list)){
  df <- data_list[[i]]
  f <- log2(df$exp) > mean_thres[i]
  df$high_mu <- f
  data_list[[i]] <- df
}

#Export the dataframes
library(WriteXLS)
WriteXLS(data_list, ExcelFileName = "Enhancer_libFL007_burst_parameters_NB_&_muthres.xls", SheetNames =names(data_list),row.names=TRUE)

