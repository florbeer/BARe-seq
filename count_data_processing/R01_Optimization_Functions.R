#BARe-seq analysis of the count data
#Complementary functions for fitting the 0-truncated NB model and inferring bursting parameters

#Function for running MLE with one of the optimization methods
#Arguments:
#x: a data frame with at least three columns: "seq_ID", "plasmid_UMI_main" (or "plasmid_bc"), and "N"
#All the counts in the data frame are processed as if the were one sample. So arranging the data frames into samples prior to running this function is required
#method: choose one option: "Nelder-Mead", "L-BFGS-B", "BFGS"
#If selecting L-BFGS-B, it requires also lower = c(3e-03,3e-03), upper = c(50,8e+05) arguments (values are only examples)
#init: vector with initial values, in most cases we used 1,1

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


#Function for computing the Kolmogorov-Smirnov test
#Input:
#y: a data frame with at least three columns: "seq_ID", "plasmid_UMI_main" (or "plasmid_bc"), and "N", and corresponding to the data of a unique sample, similar to the previous FUN_NBmle function (and not several ones, unless willing to pool the data)
#optN: the list of optimization values, it is the output of the FUN_NBmle function
#Internally it computes the Kolmogorov-Smirnov test by comparing the empirical counts versus hypothettcal counts given the inferred parameters in optN

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

    ks <- ks.test(x1[ra],x2[ra])
    pvl <- c(pvl,ks$p.value)
  }
  padj <- p.adjust(pvl,method="BH")
  return(padj)
}
