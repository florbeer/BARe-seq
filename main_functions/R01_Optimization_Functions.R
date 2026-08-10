#STAPseq and STARRseq data analysis
#Side functions for the analysis

#Function for running MLE with one of the optimization methods
#Arguments:
#x: a data frame with at least three columns: "seq_ID", "plasmid_UMI_main" (or "plasmid_bc"), and "N"
#For real data from STAPseq and so on, split the x information by individual samples
#method: choose one option: "Nelder-Mead", "L-BFGS-B", "BFGS"
#If selecting L-BFGS-B, it requires also lower = c(3e-03,3e-03), upper = c(50,8e+05) arguments (values are only examples, used for the simulated data)
#init: vector with initial values, usually I used 1,1 in most cases

FUN_NBmle <- function(x, method=c("Nelder-Mead", "L-BFGS-B", "BFGS"), init, ...){
  optL <- zo <- list()
  z  <- aggregate(x$N,by=list(x$seq_ID),mean)
  zo <- z[order(z$x,decreasing=TRUE),]
  init <- init
  
  fn <- function(x2){
    mu <- x2[1]
    size <- x2[2]
    -sum( log(dnbinom(k, size=size, mu=mu )/( 1 - dnbinom(0, size=size, mu=mu ) ) ) )
  }
  
  for ( gr in zo$Group.1 ){ #sequences
    cat(gr,"\n")
    f <- x$seq_ID == gr
    
    k <- x[f,"N"]
    
    suppressWarnings( w <- tryCatch( {
      opt <- optim( init, fn, method = method, ...)
      TRUE
    }, error = function(err){ FALSE } ))
    
    if ( w ){
      optL[[gr]] <- opt
    }
  }
  return(optL)
}


#######
#Function to estimate hessian matrices
#Version 3: adapted function to run it with either 1 sample (i.e. optLN is the direct output of FUN_NBmle function
#and not a list with several FUN_NBmle outputs). And if running it with several samples, the adapt it with
#a lapply function
#y: a data frame with at least three columns: "seq_ID", "plasmid_UMI_main" (or "plasmid_bc"), and "N", and corresponding to the data of a unique sample (not several ones)
#optN: the list of optimization values, it is the output of the FUN_NBmle function
#Internally it computes the hessian matrices with a for loop, instead of parallelization with lapply (not sure why it's so different)

FUN_HesM <- function(y, optN){
  require(numDeriv)
  #Internal function for computing MLE of the 0-truncated negative binomial model
  fn <- function(z){
    mu <- z[1]
    size <- z[2]
    -sum( log(dnbinom(k, size=size, mu=mu )/( 1 - dnbinom(0, size=size, mu=mu ) ) ) )
  }
  
  hesN2 <- list()
  for ( gr in names(optN) ){
    #cat(gr,"\n")
    f  <- y$seq_ID == gr
    k <- y[f,"N"]
    z <- optN[[gr]]$par
    
    suppressWarnings( w <- tryCatch( {
      h <- hessian(fn,z) #estimate the hessian matrix based on the function fn, estimated at the points in z
      TRUE
    }, error = function(err){ FALSE } ))
    
    if ( w ){ hesN2[[gr]] <- h }
  }
  
  return(hesN2)
}


#Function to extract the standard deviation of mu and r from the corresponding hessian matrix
#Input:
#y: a 2x2 numerical matrix, it is the hessian matrices coming from the FUN_HesM function, which are arranged into a list of values
#Output: a 2 numerical vector with the standard deviation of mu and r, respectively
#To apply this function to each matrix individually, it is required to use lapply function
FUN_hes2sd <- function(y){
  if (length(y) > 0){
    h <- y
    suppressWarnings( v <- tryCatch( {
      std <- sqrt(diag(solve(h))) 
      TRUE
    }, error = function(err){ FALSE } ))
    if ( v ){
      mu.delta <- std[1]
      r.delta <- std[2]
    }else{
      mu.delta <- NA
      r.delta <- NA
    }
  }else{
    mu.delta <- NA
    r.delta <- NA
  }
  return(c(mu.delta=mu.delta, r.delta=r.delta))
}


#Function for computing the Kolmogorov-Smirnov test
#Version 2 with inputs arranged for one sample
#Input:
#y: a data frame with at least three columns: "seq_ID", "plasmid_UMI_main" (or "plasmid_bc"), and "N", and corresponding to the data of a unique sample (not several ones)
#optN: the list of optimization values, it is the output of the FUN_NBmle function
#Internally it computes the Kolmogorov Smirnov test by comparing the empirical counts versus hypothettcal counts given the inferred parameters in optN

FUN_KStest <- function(y,optN){
  fitF <- function(k,opt){
    mu <- opt$par[1]
    size<- opt$par[2]
    dnbinom(k, size=size, mu=mu )/( 1 - dnbinom(0, size=size, mu=mu ) )
  }
  
  pvl <- c()
  
  for ( gr in names(optN) ){
    f  <- y$seq_ID == gr
    
    if (max(y[f,"N"]) < 1000) h <- hist(y[f,"N"],breaks=0:1000,plot=FALSE)
    else h <- hist(y[f,"N"],breaks=0:max(y[f,"N"]), plot=FALSE)
    k <- y[f,"N"]
    
    opt <- optN[[gr]] 
    
    x1 <- h$counts/sum(h$counts)
    x2 <- fitF(h$breaks + 1, opt)
    b <- h$breaks + 1
    ra <- 1:15
    #plot(b[ra],x1[ra],type="l")
    #lines(b[ra],x2[ra],col="red")
    ks <- ks.test(x1[ra],x2[ra])
    pvl <- c(pvl,ks$p.value)
  }
  padj <- p.adjust(pvl,method="BH")
  return(padj)
}



###################
#Profile likelihood functions for estimation of confidence intervals
#All these functions were tested on "21Stark_Collab/E39_STARRseq/R07E39_ProfileLikelihood.R" script

#Function to evaluate likelihood function of the 0-truncated negative binomial model. It only computes the likelihood value
#by providing specific parameter values, without performing optimization anymore
#x1: mu parameter
#x2: r or size parameter
#obs: observations or UMI counts
fn2 <- function(x1, x2, obs){
  mu <- x1
  size <- x2
  obs <- obs
  -sum( log(dnbinom(obs, size=size, mu=mu )/( 1 - dnbinom(0, size=size, mu=mu ) ) ) )
}

##
#Function that generates the ratio profile likelihood when fixing mu
#Mu is fixed over a range of values, and the likelihood expressions are evaluated (with fn2 function), without optimization
#mu_exp: range of possible mu values to test, given into log2 scale to cover a better spectrum of the parameter space
#obs: observations or UMI counts in these case coming from the original data
#mle_pars: vector with mu and r parameter values inferred by maximum likelihood
#Output: ratio likelihood of mu
ratio_mu3 <- Vectorize( function(mu_exp, obs, mle_pars){
  obs <- obs
  lk_mle <- fn2(mle_pars[1], mle_pars[2], obs)
  
  mu <- 2^mu_exp
  r_mle <- mle_pars[2]
  lk_mu <- fn2(mu, r_mle, obs)
  
  r_mu <- exp(lk_mle - lk_mu)
  return(r_mu)
}, vectorize.args = "mu_exp")

#Function that search for the thresholds of the 95% confidence interval
#0.147 is the critical point of a chi squared distribution with 1 degree of freedom and significance = 0.05. These conditions describe the distribution of the profile likelihood.
#Same parameters as before:
#mu_exp: range of possible mu values to test, given into log2 scale to cover a better spectrum of the parameter space
#obs: observations or UMI counts in these case coming from the original data
#mle_pars: vector with mu and r parameter values inferred by maximum likelihood
#output: critical values when the profile likelihood is equal to 0.147. They are in log2 scale, so they need to be transformed into linear scale
fn_muC3 <- Vectorize(function(mu_exp, obs, mle_pars) {
  #root_val <- ratio_mu(mu = mu, init=init, obs=obs, mle_pars=mle_pars, method=method, ... ) - 0.147
  a <- ratio_mu3(mu_exp = mu_exp, obs=obs, mle_pars=mle_pars) #%>% unlist()
  root_val <- a - 0.147
  return(root_val)
}, vectorize.args = "mu_exp")


###
#Function that generates the ratio profile likelihood when fixing r (size or dispersion parameter)
#r is fixed over a range of values, and the likelihood expressions are evaluated (with fn2 function), without optimization
#size_exp: range of possible r (size) values to test, given into log2 scale to cover a better spectrum of the parameter space
#obs: observations or UMI counts in these case coming from the original data
#mle_pars: vector with mu and r parameter values inferred by maximum likelihood
#Output: ratio likelihood of r (size) parameter
ratio_r3 <- Vectorize( function(size_exp, obs, mle_pars){
  obs <- obs
  lk_mle <- fn2(mle_pars[1], mle_pars[2], obs)
  
  mu_mle <- mle_pars[1]
  size <- 2^size_exp
  lk_r <- fn2(mu_mle, size, obs)
  
  rat_r <- exp(lk_mle - lk_r)
  return(rat_r)
}, vectorize.args = "size_exp")

#Function that search for the thresholds of the 95% confidence interval
#0.147 is the critical point of a chi squared distribution with 1 degree of freedom and significance = 0.05. These conditions describe the distribution of the profile likelihood.
#Same parameters as before:
#size_exp: range of possible r (size) values to test, given into log2 scale to cover a better spectrum of the parameter space
#obs: observations or UMI counts in these case coming from the original data
#mle_pars: vector with mu and r parameter values inferred by maximum likelihood
#output: critical values when the profile likelihood is equal to 0.147. They are in log2 scale, so they need to be transformed into linear scale
fn_rC3 <- Vectorize(function(size_exp, obs, mle_pars) {
  #root_val <- ratio_r(size=size, init=init, obs=obs, mle_pars=mle_pars, method=method, ... ) - 0.147
  a <- ratio_r3(size_exp=size_exp, obs=obs, mle_pars=mle_pars)
  root_val <- a - 0.147
  return(root_val)
}, vectorize.args = "size_exp")


##
#Function for computing the confidence intervals across diferent samples in the datasets
#x: table with original UMI counts across different regulatory elements and samples, coming from the massive parallel reporter assay
#optLN: list with the parameters inferred by maximum likelihood per regulatory element (or sequence) and per sample
#hesN: list with the hessian matrices estimated from optLN
#samples: character vector with the names of the samples to evaluate
#output: a list with 3 levels: samples, regulatory element (sequence) and roots. Roost contain two slots: mu_roots and r_roots for mu and r, respectively
FUN_CI95 <- function(x, optLN, hesN, samples){
  require(dplyr)
  require(rootSolve)
  CIN <- list()
  err <- 2.5
  
  for ( sa in samples ){
    CIN[[sa]] <- list() 
    
    for ( gr in names(optLN[[sa]]) ){
      f  <- x$sample == sa
      y  <- x[f,]
      f  <- y$seq_ID == gr
      sum(f)
      k <- y[f,"N"]
      
      z <- optLN[[sa]][[gr]]$par
      hesm <- hesN[[sa]][[gr]]
      prop_sigma <- sqrt(diag(solve(hesm)) )
      
      #Apply low and upper values for mu
      if(!is.na(prop_sigma[1]) ){
        a <- z[1] - err*prop_sigma[1]
        if(a > 0) low_mu <- log2(a)
        else low_mu <- -16
        
        up_mu <- log2(z[1] + err*prop_sigma[1])
      }
      else{
        low_mu <- -16
        up_mu <- log2(z[1] + err)
      }
      
      
      suppressWarnings(w <- tryCatch( {
        ab <- uniroot.all(fn_muC3, c(low_mu,up_mu), obs=k, mle_pars=z) %>% 2^.
        TRUE
      }, error = function(err){ FALSE }
      ) )
      
      if (w & length(ab) < 2) {
        if (z[1] < 0.1) {
          low_mu <- -16
          up_mu <- 1
        }
        if (z[1] >= 0.1 & z[1] < 1) {
          low_mu <- -13
          up_mu <- 4
        }
        if (z[1] >= 1) {
          low_mu <- -2.5
          up_mu <- 6
        }
        
        suppressWarnings(w <- tryCatch( {
          ab <- uniroot.all(fn_muC3, c(low_mu,up_mu), obs=k, mle_pars=z)  %>% 2^.
          TRUE
        }, error = function(err){ FALSE }
        ) )
      }
      
      #Apply low and upper values for r (size) parameter
      if(!is.na(prop_sigma[2]) ){
        a2 <- z[2] - err*prop_sigma[2]
        if(a2 > 0) {low_r <- log2(a2)}
        else {low_r <- -16}
        
        up_r <- log2(z[2] + err*prop_sigma[2])
      }
      else{
        low_r <- -16
        up_r <- log2(z[2] + err)
      }
      
      suppressWarnings(w2 <-tryCatch( {
        ab2 <- uniroot.all(fn_rC3, c(low_r,up_r), obs=k, mle_pars=z) %>% 2^.
        TRUE
      }, error = function(err){ FALSE }
      ) )
      
      if(w2 & length(ab2) < 2){
        if (z[2] < 0.01) {
          low_r <- -16
          up_r <- 0
        }
        if (z[2] >= 0.01 & z[2] < 2) {
          low_r <- -13
          up_r <- 4
        }
        if (z[2] >= 2) {
          low_r <- -10
          up_r <- 20
        }
        
        suppressWarnings(w2 <-tryCatch( {
          ab2 <- uniroot.all(fn_rC3, c(low_r,up_r), obs=k, mle_pars=z) %>% 2^.
          TRUE
        }, error = function(err){ FALSE }
        ) )
      }
      
      #Third adjustment for r parameter
      if(w2 & length(ab2) < 2){
        low_r <- -13
        up_r <- 20
        
        suppressWarnings(w2 <-tryCatch( {
          ab2 <- uniroot.all(fn_rC3, c(low_r,up_r), obs=k, mle_pars=z) %>% 2^.
          TRUE
        }, error = function(err){ FALSE }
        ) )
      }
      
      CIN[[sa]][[gr]] <- list()
      if (w) { CIN[[sa]][[gr]]$mu_roots <- ab }
      else { CIN[[sa]][[gr]]$mu_roots <- NA }
      
      if (w2) { CIN[[sa]][[gr]]$r_roots <- ab2 }
      else { CIN[[sa]][[gr]]$r_roots <- NA }
    }
  }
  return(CIN)
}


#Functions for motif analyses
#Wilcoxon test, copied and pasted from ../MultiomicsPBMCs/R26SideFunctions_2.R script

#Having a list x, with numeric values in each element, make paired comparisons between element 1 and 2, 3 and 4, and so on...
pairedWilcox <- function(x, alternative="two.sided", ...){
  j <- 2
  wt <- c()
  for (i in seq(from=1, to=length(x), by=2)){
    a <- x[[i]]
    a2 <- x[[j]]
    wt <- c(wt, tryCatch({
      wilcox.test(a, a2, alternative="two.sided", ...)$p.value },
      error = function(err) {NA}
    ))
    j <- j+2
  }
  return(wt)
}

#Transforming the numeric p values (contained into a numeric vector) into asterisk, for easier interpretation
num2star <- function(wt){
  wt2 <- wt
  wt2[wt < 0.001] <- "***"
  wt2[wt < 0.01 & wt > 0.001] <- "**"
  wt2[wt < 0.05 & wt > 0.01] <- "*"
  wt2[wt > 0.05] <- "ns"
  return(wt2)
}

# Function to analyze one parameter for all motifs
#Inputs:
#cp_data: "compiled data" data frame, with all the parameters and motif content from different samples of interest
#param_col: character vector, corresponding to the parameter to be analysed, and equal to one of the colnames in the list of data frames
#qr_motifs: "query motifs", character vector that includes all the motif names to be analyzed
#qr_sample: "query sample", character with the colnames of the data frame to split into groups, it can be sample, sample2 (short names), condition, library, etc.
#add_sample: character vector indicating additional colnames to include in the output data frame, together with "qr_sample".Default is NULL
#param_name: character, indicating the parameter name to use for plotting afterwards

analyze_parameter <- function(cp_data, param_col, qr_motifs, qr_sample, add_sample=NULL, param_name) {
  MO <- qr_motifs
  
  aL <- lapply(MO, function(motif) {
    a <- c(motif, paste0("no", motif))
    a2 <- unique(cp_data[[qr_sample]])
    a3 <- paste(rep(a2, each=length(a)), rep(a, times=length(a2)), sep=":")
    
    cp_data$group <- paste(cp_data[[qr_sample]], cp_data[,motif], sep=":")
    cp_data$group <- factor(cp_data$group, levels=a3)
    
    y2 <- split(cp_data[[param_col]], cp_data$group)
    return(y2)
  })
  names(aL) <- MO
  
  # Create reference data frame for samples
  #f <- !duplicated(cp_data$sample)
  #df2 <- cp_data[f,c("sample", "sample2", "lib", "condition")]
  if(is.null(add_sample)){
    ac <- qr_sample
  }
  else{ ac <- c(qr_sample, add_sample) }
  
  ac2 <- split(cp_data[,ac], cp_data[[qr_sample]])
  df2 <- sapply(ac2, function(y) apply(y,2,unique)) %>% t() %>% as.data.frame()
  
  # Statistical analysis for each motif
  ab <- lapply(aL, function(y) {
    # Wilcoxon test
    y2 <- pairedWilcox(y)
    y3 <- ifelse(y2 < 0.05, "*", NA)
    
    # Fold change based on the median values
    j <- 2
    fc <- c()
    for (i in seq(from=1, to=length(y), by=2)){
      a <- y[[i]]
      a2 <- y[[j]]
      fc <- c(fc, median(a)/median(a2))
      j <- j+2
    }
    log2_fc <- log2(fc)
    
    df <- data.frame(fc=fc, log2_fc=log2_fc, wilcox_pval=y3)
    df <- cbind(df, df2)
    return(df)
  })
  
  # Add motif names to each result
  for (i in 1:length(ab)){
    ab[[i]]$motif <- MO[i]
  }
  
  # Combine all motifs for this parameter
  result <- Reduce(rbind, ab)
  result$motif <- factor(result$motif, levels=MO)
  result$param <- param_name
  
  return(result)
}