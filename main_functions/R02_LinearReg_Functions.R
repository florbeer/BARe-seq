#Functions for linear regression models


#Functions for elastic net regression


#Function to get the best alpha and lambda values for performing elastic net regression
#Current version: 1.3
#Main difference with previous versions 1.1 and 1.2 is that version 1.3 gives the option whether to run a full model with interaction terms, or just predictor variables without any interaction
#The function will compute a model for different alpha values.
#The best alpha value will be selected based on the minimum value of the mean square errors (real - predicted)
#After this, we select lambda values for 3 alpha values: best alpha, alpha=1 ( for lasso reg.), alpha=0 (for ridge reg.)
#Lambda values are selected by cross validation and getting the values that deviates 1 standard error from the minimum lambda.
#Inputs:
#df: data frame with the values to regress (response variable) and the explanatory variables (motif presence/absence)
#resp_var: character, column name with the response variable to be predicted by the model
#motifs: character vector with all the motifs that should be included in the model
#seed1: for random sampling and making the training and validation partitions. Default is NULL and then we use 578
#seed2: for computing the cross validation for each value of alpha Default is NULL and then we use 253
#alpha_test: alpha values to test for model selection. Default is NULL and then we use seq(from=0, to=1, by=0.05)
#crossed_terms: logical value indicating whether to compute linear regression with crossed terms (interaction terms), or only the predictor variables. Default is TRUE
#square_err: logical value indicating whether to compute mean squared errors to determine the best alpha value or not. Default is TRUE
#Output: list with a data frame and a vector.
#The data frame contains the values for alpha and lambda to test. If square_err=TRUE, data frame has 3 alpha values: ridge, best alpha and lasso regression
#If square_err=FALSE, data frame has values for all the alpha's tested in alpha_test
#The slot with the vector contains the index of the values selected for the training data. So we can use it later for getting the validation partition and compute coefficients.

FUN_netr_modsel <- function(df, resp_var, motifs, seed1=NULL, seed2=NULL, alpha_test=NULL, crossed_terms=TRUE, square_err=TRUE){
  require("dplyr")
  require("glmnet")
  FUN_inter <- as.formula( ~ .*.)
  
  if(!is.null(seed1))
    seedn <- seed1
  else seedn <- 578
  
  if(!is.null(seed2))
    seedn2 <- seed2
  else seedn2 <- 253
  
  if(!is.null(alpha_test))
    alpha_test <- alpha_test
  else alpha_test <- seq(from=0, to=1, by=0.05)
  
  ay <- df[,resp_var]
  x_var <- df[,motifs]
  
  if(crossed_terms==TRUE){
    x_var <- model.matrix(FUN_inter, x_var)[, -1]
  } else{
    x_var <- as.matrix(x_var)
  }
  
  #Make a subset of the data, 70% for training, and 30% for testing
  #For convenience, use the same training and validation data partitions for each model assessment
  set.seed(seedn)
  f <- sample(1:nrow(x_var), .7*nrow(x_var), replace = FALSE)
  
  aL <- lapply(alpha_test, function(y){
    set.seed(seedn2)
    y2 <- cv.glmnet(x_var[f,], ay[f], type.measure="mse", alpha=y, family="gaussian")
    return(y2)
  })
  names(aL) <- paste("alpha", alpha_test, sep="_")
  
  #Compute mean square errors:
  if(square_err==TRUE){
    aL2 <- lapply(aL, function(y){
      y2 <- predict(y, s=y$lambda.1se, newx=x_var[f,])
      y3 <- mean((ay[f] - y2)^2)
      return(y3)
    }) %>% unlist()
    
    f2 <- which(aL2 == min(aL2)) %>% names()
    ab <- gsub(".+_", "", f2) %>% as.numeric()
    
    #Vector with different lambda values, corresponding to the selected alpha values
    lb <- c(aL[["alpha_0"]]$lambda.1se, aL[[f2]]$lambda.1se, aL[["alpha_1"]]$lambda.1se)
    df2 <- data.frame(alpha=c("alpha_0", f2, "alpha_1"), alpha_val=c(0, ab, 1), lambda_se=lb)
    rownames(df2) <- c("ridge", "alpha_mse", "lasso")
  }
  else {
    f2 <- names(aL)
    ab <- gsub(".+_", "", f2) %>% as.numeric()
    
    lb <- sapply(aL, function(y) y$lambda.1se)
    df2 <- data.frame(alpha=f2, alpha_val=ab, lambda_se=lb)
    rownames(df2) <- c(f2)
  }
  
  aL3 <- list(parameters=df2, index=f)
  return(aL3)
}



########
#Estimation of p value by permutation method

#Arguments:
#df: data frame with predictor and predicted variables
#resp_var: character value indicating the column name in df which represents the dependent or response variable that is aim to be predicted by the model
#motifs: character vector indicating the columns' names in df which represents the predictor variables to be used by the model
#alpha: numeric value, alpha value for the glmnet function, computing elastic net regression
#lambda: numeric value, lambda value for the glmnet function, computing elastic net regression
#index: numeric vector, indicating the samples (rows or regulatory elements in df) that are selected for the
#training data. If provided, the coefficients of model computed in the function are fit with the validation
#data, instead of the training data. Default is NULL, and then the model would be computed across all observations
#n_per=1000, number of permutations to perform, default = 1000
#crossed_terms: logical value indicating whether to compute linear regression with crossed terms (interaction terms), or only the predictor variables. Default is TRUE
#The function computes internally the design matrix with predictor variables and the vector with predicted
#variable
#The function also computes internally the model with the provided alpha and lambda values
#Output: a data frame containing the coefficient values, p value, and adjusted p value
FUN_per <- function(df, resp_var, motifs, n_per=1000, lambda, alpha, index=NULL, crossed_terms=TRUE){
  require("dplyr")
  require("glmnet")
  
  if(!is.null(index))
    index <- index
  else index <- 1:nrow(df)
  
  ab <- alpha
  lb <- lambda
  FUN_inter <- as.formula( ~ .*.)
  
  ay <- df[,resp_var]
  x_var <- df[,motifs]
  #x_var <- model.matrix(FUN_inter, x_var)[, -1]
  if(crossed_terms==TRUE){
    x_var <- model.matrix(FUN_inter, x_var)[, -1]
  } else{
    x_var <- as.matrix(x_var)
  }
  
  #Get the default model
  f <- index
  model <- glmnet(x_var[-f,], ay[-f], alpha = ab, lambda = lb)
  coefs <- coef(model)
  coefs_n <- coef(model) %>% rownames()
  
  #coefs <- coef(model)
  #if(!is.null(lambda))
  #  lb <- lambda
  #else lb <- model$lambda.1se
  
  #Get null distribution for each coefficient
  null_distribution <- matrix(NA, nrow=length(coefs), ncol=n_per)
  
  for (i in 1:n_per){
    #Permute the response variable
    y_perm <- ay[sample(1:length(ay), replace = TRUE)]
    # Fit elastic net model
    model_perm <- glmnet(x_var, y_perm, alpha = ab, lambda = lb)
    null_distribution[,i] <- coef(model_perm)[,1]
  }
  
  #Compare the null distribution against the model coefficients and estimate p value and p adjusted value
  a <- cbind(coefs, null_distribution)
  a2 <- apply(a, 1, function(w){
    ab <- w[1]
    ab2 <- w[2:length(w)]
    pval <- mean(abs(ab2) >= abs(ab))
    padj <- p.adjust(pval) #p adjusted value
    w2 <- cbind(pval=pval, padj=padj)
    return(w2)
  }) %>% t() %>% as.data.frame()
  colnames(a2) <- c("pval", "padj")
  
  #Make a data frame with coefficients and p values
  netreg <- cbind(data.frame(coefficients=coefs[,1]), a2)
  return(netreg)
}



