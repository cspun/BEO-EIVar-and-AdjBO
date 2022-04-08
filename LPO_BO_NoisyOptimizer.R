############################################################
# This file contains functions to find optimal lambda using Bayes Opt
# (with EIVar as acquisition function and time series adjustment)
# for LPO with sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate
############################################################


library('DiceKriging')
source('noisy_optimizer_withStopping.R')


# Loss function in choosing lambda: Minimize negative sharpe ratio
# negsharperatio<-function(r,eta,gamma=1){
#   # if eta is almost zero or null, set loss to 100000
#   if(sum(eta^2)<1e-06 | is.null(eta)){nsr=100000} else {w=eta/sum(eta);dailyr=r%*%w;nsr=-mean(dailyr)/sd(dailyr)}
#   return (nsr)
# }

negCEQ<-function(r,eta,gamma=1){
  if(sum(eta^2)<1e-06 | is.null(eta)){nceq=100000} else {w=eta/sum(eta);dailyr=r%*%w;nceq=-(mean(dailyr)-gamma/2*var(dailyr))}
  nceq
}

# Get indexes of partition to pseudo train and pseudo test 
cv.part<-function (n, k) 
{
  ntest <- floor(n/k)
  ntrain <- n - ntest
  #    ind <- sample(n)
  ind<-1:n
  trainMat <- matrix(NA, nrow = ntrain, ncol = k)
  testMat <- matrix(NA, nrow = ntest, ncol = k)
  nn <- 1:n
  for (j in 1:k) {
    sel <- ((j - 1) * ntest + 1):(j * ntest) #index for test
    testMat[, j] <- ind[sel] #save test indices 
    sel2 <- nn[!(nn %in% sel)] # index for train  (complement of index of test)
    trainMat[, j] <- ind[sel2] # save train indices
  }
  return(list(trainMat=trainMat, testMat=testMat))
}

# Get indexes of partition for rolling windows evaluation
rollwindows.part <- function(n.train, n.test)
{
  train = matrix(NA, nrow = n.train, ncol = n.test)
  test = matrix(NA, nrow = 1, ncol= n.test)
  
  for (i in 1:n.test)
  {
    trainInd = (i):(i + n.train - 1)
    testInd = i + n.train
    train[,i] = trainInd
    test[,i] = testInd
  }
  return(list(train=train, test=test))
}



# Bayes Opt for LPO with Sample Estimates with EIVar as acquisition function and Time Series adjustments
BO_noisy_optimizer_l1fast <- function(r, use.prevlambda = TRUE, prevlambda= NULL,
                                      coef.trend = NULL, coef.cov = NULL, coef.var = NULL,noise.var = NULL,
                                      random.start=FALSE, delta,step,fold,normalized = FALSE,
                                      perturb=FALSE,lambda=NULL,logspaced=TRUE,coef=0.8,
                                      nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  p<-ncol(r);n<-nrow(r)
  
  if (is.null(lambda)) { 
    coef=ifelse(coef>1,1,coef)
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    lambda.min = lambda.max/50
  }
  
  # Define bounds of the search
  bounds = list(lambda = c(lambda.min,lambda.max))
  
  # Number of initial observations
  n.initlambda = 15
  
  evaluate.initlambda = FALSE
  add.flag=FALSE
  
  # Generate few lambda values as initial observations
  if (use.prevlambda) # Use previous time steps' results
  {
    init_observations = as.data.frame(prevlambda)
    #cat("Using previous lambda with length",nrow(init_observations),"\n")
    if (nrow(init_observations) == 1)
    {
      init.lambda = seq(lambda.min, lambda.max, length.out =  n.initlambda)
      evaluate.initlambda = TRUE
      add.flag = TRUE
      #cat("Evaluating more lambdas with nlambda ", n.initlambda, "\n")
    }
  }
  else # Pick new start points  
  {
    if (random.start==TRUE)
    {
      init.lambda = runif(n= n.initlambda, min = lambda.min, max = lambda.max)
      evaluate.initlambda = TRUE
      cat('Using random start points \n')
    }
    else
    {
      init.lambda = seq(lambda.min, lambda.max, length.out = n.initlambda)
      evaluate.initlambda = TRUE
      cat('Using uniformly spaced start points \n')
    }
  }
  
  if(evaluate.initlambda)
  {
    part.list<-cv.part(n,fold) # Get indexes of partition to pseudo train and pseudo test 
    SR.re.initlambda<-matrix(0,nrow=fold,ncol=n.initlambda)
    
    for (j in 1:fold) {
      indtrain<-part.list$trainMat[,j] #get index of training data for the current fold
      r.train<-r[indtrain,] # training data
      l1fast.init<-l1fast(r.train,delta=delta,step=step,normalized=normalized,lambda= init.lambda, shrink=TRUE)
      indtest<-part.list$testMat[,j] #get index of test data for the current fold
      r.test <- r[indtest,]
      
      # betalist[[jl]] is corresponding eta for current lambda used in the iteration
      for (jl in 1:n.initlambda) {
        SR.re.initlambda[j,jl]<- negsharperatio(r.test, l1fast.init$betalist[[jl]]) } ###
    }
    
    # Create data frame for storing initial observations (i.e. lambda and corresponding f(lambda))
    # Use robust SR as objective function
    f.initlambda = apply(SR.re.initlambda,2,max) 

    init_observations = data.frame(init.lambda, f.initlambda)
  }
    
  if (add.flag) # bind prev and new lambda
  {
    names(init_observations) <- c("lambda","Value")
    prevdf = as.data.frame(prevlambda)
    if (nrow(prevdf)==1)
    {
      init_observations = init_observations
    }
    else
    {
      names(prevdf) <- c("lambda","Value")
      init_observations = rbind(prevdf, init_observations)
    }
  }
  
  # Rename columns
  names(init_observations) <- c("lambda","Value")
  
  # Function to be optimized by BayesOpt   
  l1fast_sr <- function(lambda)
  {
    part.list.lambda<-cv.part(n,fold) # Get indexes of partition to pseudo train and pseudo test 
    vec.re.lambda = rep(0, fold)
    
    for (f in 1:fold) {
      indtrain.lambda<-part.list.lambda$trainMat[,f] #get index of training data for the current fold
      r.train<-r[indtrain.lambda,] # training data
      l1fast.lambda<-l1fast(r.train,delta=delta,step=step,normalized=normalized,lambda= lambda,shrink=TRUE)
      indtest.lambda<-part.list.lambda$testMat[,f] #get index of test data for the current fold
      r.test <- r[indtest.lambda,]
      vec.re.lambda[f] = negsharperatio(r.test, l1fast.lambda$betalist[[1]])
    }
    
    # Calculate objective function
    SR.test = max(vec.re.lambda)

    return (SR.test)
  }
  
  # Fit initial GP model with initial observations. Later will be updated inside the optimizer
  model = km( formula=~1, design= as.data.frame(init_observations[,1]), response= as.data.frame(init_observations[,2]),
              covtype="matern5_2", coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
              noise.var = rep(noise.var,nrow(init_observations)),
              estim.method="MLE",optim.method = "BFGS", control = list(trace = FALSE))
  
  # Run Bayes Opt
  cat("BO starts at ",  as.character(Sys.time()),"\n")
  sink('null')
  optim.param = list()
  #optim.param$plugin.type = "min_km_mean" #Note: plugin.type is needed for EI criterion only
  
  optim.result = noisy.optimizer_withStopping(optim.crit = "EIwithVar", optim.param = optim.param,
                                 model = model,n.ite=15,
                                 funnoise = l1fast_sr,noise.var = noise.var,
                                 lower = lambda.min, upper = lambda.max,
                                 CovReEstimate = TRUE, NoiseReEstimate = TRUE
                                 )
  sink()
  cat("BO ends at ",  as.character(Sys.time()),"\n")
  
  cat("OPTIM LAMBDA IS ", optim.result$best.x,"\n")

  # Get final GP hyperparameter settings
  km.coef = coef(optim.result$model)
  noise.var = optim.result$history.noise.var[length(optim.result$history.noise.var)]

  return (list(best.x= optim.result$best.x, ybestvec= optim.result$ybestvec, kmmeanvec = optim.result$kmmeanvec, 
               history.x = optim.result$history.x, history.y = optim.result$history.y,
               coef.trend = km.coef$trend, coef.cov = km.coef$range, coef.var = km.coef$sd2, noise.var =noise.var )) 
}



# Bayes Opt for LPO with Bayes Stein Estimates with EIVar as acquisition function and Time Series adjustments
BO_noisy_optimizer_l1fastbs <- function(r,mu_0, tau, use.prevlambda = TRUE, prevlambda= NULL,
                                        coef.trend = NULL, coef.cov = NULL, coef.var = NULL,noise.var = NULL,
                                        random.start = FALSE, delta,step,fold, shrinkage= FALSE,
                                        normalized = FALSE,
                                        lambda=NULL,logspaced=TRUE,coef=0.8,
                                        nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  p<-ncol(r);n<-nrow(r)
  
  if (is.null(lambda)) { 
    coef=ifelse(coef>1,1,coef)
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    lambda.min = lambda.max/50
  }
  
  # Define bounds of the search
  bounds = list(lambda = c(lambda.min,lambda.max))
  
  # Number of initial observations
  evaluate.initlambda = FALSE
  add.flag = FALSE
  n.initlambda = 15
  
  # Generate few lambda values as initial observations
  if (use.prevlambda) # use previous day's results
  {
    init_observations = as.data.frame(prevlambda)
    cat("Using previous lambda with length",nrow(init_observations),"\n")
    if (nrow(init_observations) == 1)
    {
      init.lambda = seq(lambda.min, lambda.max, length.out = n.initlambda)
      evaluate.initlambda = TRUE
      add.flag = TRUE
      #cat("Evaluating more lambdas with nlambda ", n.initlambda, "\n")
    }
  }
  else # Pick new start points
  {
    if (random.start==TRUE)
    {
      init.lambda = runif(n= n.initlambda, min = lambda.min, max = lambda.max)
      evaluate.initlambda = TRUE
      cat('Using random start points \n')
    }
    else
    {
      init.lambda = seq(lambda.min, lambda.max, length.out = n.initlambda)
      evaluate.initlambda = TRUE
      cat('Using uniformly spaced start points \n')
    }
  }
  
  if (evaluate.initlambda)
  {
    part.list<-cv.part(n,fold) # Get indexes of partition to pseudo train and pseudo test 
    SR.re.initlambda<-matrix(0,nrow=fold,ncol=n.initlambda)
    
    for (j in 1:fold) {
      indtrain<-part.list$trainMat[,j] #get index of training data for the current fold
      r.train<-r[indtrain,] # training data
      l1fast.init<- l1fast.bayesstein(r.train,mu_0=mu_0,tau=tau,delta=delta,
                                      step=step,shrinkage = shrinkage,
                                      normalized=normalized,lambda=init.lambda)
      indtest<-part.list$testMat[,j] #get index of test data for the current fold
      r.test <- r[indtest,]
      
      # betalist[[jl]] is corresponding eta for current lambda used in the iteration
      for (jl in 1:n.initlambda) {
        SR.re.initlambda[j,jl]<- negCEQ(r.test, l1fast.init$betalist[[jl]]) } ###
    }
    
    # Create data frame for storing initial observations (i.e. lambda and corresponding f(lambda))
    # Use robust SR as objective function
    f.initlambda = apply(SR.re.initlambda,2,mean) 
    
    init_observations = data.frame(init.lambda, f.initlambda)
  }
  
  if (add.flag) # bind prev and new lambda
  {
    names(init_observations) <- c("lambda","Value")
    prevdf = as.data.frame(prevlambda)
    if (nrow(prevdf)==1)
    {
      init_observations = init_observations
    }
    else
    {
      names(prevdf) <- c("lambda","Value")
      init_observations = rbind(prevdf, init_observations)
    }
  }
  
  # Rename columns
  names(init_observations) <- c("lambda","Value")
  
  # Function to be optimized by BayesOpt 
  l1fast_sr <- function(lambda)
  {
    part.list.lambda<-cv.part(n,fold) # Get indexes of partition to pseudo train and pseudo test 
    
    vec.re.lambda = rep(0, fold)
    
    for (f in 1:fold) {
      indtrain.lambda<-part.list.lambda$trainMat[,f] #get index of training data for the current fold
      r.train<-r[indtrain.lambda,] # training data
      l1fast.lambda<- l1fast.bayesstein(r.train,mu_0=mu_0,tau=tau,delta=delta,
                                        step=step,shrinkage = shrinkage,
                                        normalized=normalized,lambda= lambda)
      indtest.lambda<-part.list.lambda$testMat[,f] #get index of test data for the current fold
      r.test <- r[indtest.lambda,]
      vec.re.lambda[f] = negCEQ(r.test, l1fast.lambda$betalist[[1]])
    }
    
    # Calculate objective function
    SR.test = mean(vec.re.lambda)

    return (SR.test)
  }
  
  # Fit initial GP model with initial observations. Later will be updated inside the optimizer
  model = km( formula=~1, design= as.data.frame(init_observations[,1]), response= as.data.frame(init_observations[,2]),
              covtype="matern5_2",coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
              noise.var=rep(noise.var,nrow(init_observations)),
              estim.method="MLE",optim.method = "BFGS", control = list(trace = FALSE))
  
  # Run Bayes Opt
  cat("BO starts at ",  as.character(Sys.time()),"\n")
  sink("null")
  optim.param = list()
  #optim.param$plugin.type = "min_km_mean" #Note: plugin.type is needed for EI criterion only
  optim.result = noisy.optimizer_withStopping(optim.crit = "EIwithVar", optim.param = optim.param,
                                              model = model,n.ite=15, noise.var= noise.var,
                                              funnoise = l1fast_sr,
                                              lower = lambda.min, upper = lambda.max,
                                              CovReEstimate = TRUE, NoiseReEstimate = TRUE,
                                              control = list(trace = FALSE))
  sink()
  cat("BO ends at ",  as.character(Sys.time()),"\n")
  
  cat("OPTIM LAMBDA IS ", optim.result$best.x,"\n")

  # Get final GP hyperparameter settings
  km.coef = coef(optim.result$model)
  noise.var = optim.result$history.noise.var[length(optim.result$history.noise.var)]
  
  return (list(best.x= optim.result$best.x, ybestvec= optim.result$ybestvec, kmmeanvec = optim.result$kmmeanvec,
               history.x = optim.result$history.x, history.y = optim.result$history.y,
               coef.trend = km.coef$trend, coef.cov = km.coef$range, coef.var = km.coef$sd2, noise.var =noise.var ))
  
}



# Bayes Opt for LPO with Bayes Asset Pricing Estimates with EIVar as acquisition function and Time Series adjustments
BO_noisy_optimizer_l1fastassetp <- function(r,f, sigma_alpha,nu, use.prevlambda = TRUE, prevlambda= NULL,
                                            coef.trend = NULL, coef.cov = NULL, coef.var = NULL,noise.var = NULL,
                                            random.start=FALSE,
                                            delta,step,fold,
                                            normalized = FALSE,
                                            lambda=NULL,logspaced=TRUE,coef=0.8,
                                            nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  p<-ncol(r);n<-nrow(r)
  
  if (is.null(lambda)) { 
    coef=ifelse(coef>1,1,coef)
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    lambda.min = lambda.max/50
  }
  
  # Define bounds of the search
  bounds = list(lambda = c(lambda.min,lambda.max))
  
  # Number of initial observations
  evaluate.initlambda = FALSE
  add.flag = FALSE
  n.initlambda = 15
  
  # Generate few lambda values as initial observations
  if (use.prevlambda) # use previous day's results
  {
    init_observations = as.data.frame(prevlambda)
    cat("Using previous lambda with length",nrow(init_observations),"\n")
    if (nrow(init_observations) == 1)
    {
      init.lambda = seq(lambda.min, lambda.max, length.out =  n.initlambda)
      evaluate.initlambda = TRUE
      add.flag = TRUE
      #cat("Evaluating more lambdas with nlambda ", n.initlambda, "\n")
    }
  }
  else # Pick new start points
  {
    if (random.start==TRUE)
    {
      init.lambda = runif(n= n.initlambda, min = lambda.min, max = lambda.max)
      evaluate.initlambda = TRUE
      cat('Using random start points \n')
    }
    else
    {
      init.lambda = seq(lambda.min, lambda.max, length.out = n.initlambda)
      evaluate.initlambda = TRUE
      cat('Using uniformly spaced start points \n')
    }
  }
  
  if (evaluate.initlambda)
  {
    part.list<-cv.part(n,fold) # Get indexes of partition to pseudo train and pseudo test 
    SR.re.initlambda<-matrix(0,nrow=fold,ncol=n.initlambda)
    
    for (j in 1:fold) {
      indtrain<-part.list$trainMat[,j] #get index of training data for the current fold
      r.train<-r[indtrain,] # training data
      f.train = f[indtrain,] #training data for factor returns
      l1fast.init<- l1fast.bayesassetpricing(r.train,f.train, sigma_alpha = sigma_alpha,
                                             nu=nu, delta=delta,
                                             step=step,normalized=normalized,lambda=init.lambda)
      indtest<-part.list$testMat[,j] #get index of test data for the current fold
      r.test <- r[indtest,]
      
      # betalist[[jl]] is corresponding eta for current lambda used in the iteration
      for (jl in 1:n.initlambda) {
        SR.re.initlambda[j,jl]<- negCEQ(r.test, l1fast.init$betalist[[jl]]) } ###
    }
    
    # Create data frame for storing initial observations (i.e. lambda and corresponding f(lambda))
    # Use robust SR as objective function
    f.initlambda = apply(SR.re.initlambda,2,mean)
    init_observations = data.frame(init.lambda, f.initlambda)
  }
  
   
  if (add.flag) # bind prev and new lambda
  {
    names(init_observations) <- c("lambda","Value")
    prevdf = as.data.frame(prevlambda)
    if (nrow(prevdf)==1)
    {
      init_observations = init_observations
    }
    else
    {
      names(prevdf) <- c("lambda","Value")
      init_observations = rbind(prevdf, init_observations)
    }
  }
 
  # Rename columns
  names(init_observations) <- c("lambda","Value")
  
  # Function to be optimized by BayesOpt 
  l1fast_sr <- function(lambda)
  {
    part.list.lambda<-cv.part(n,fold) # Get indexes of partition to pseudo train and pseudo test 
    
    vec.re.lambda = rep(0, fold)
    
    for (fo in 1:fold) {
      indtrain.lambda<-part.list.lambda$trainMat[,fo] #get index of training data for the current fold
      r.train<-r[indtrain.lambda,] # training data
      f.train = f[indtrain.lambda, ] #training data for factor returns
      l1fast.lambda<-l1fast.bayesassetpricing(r.train,f.train, sigma_alpha = sigma_alpha,
                                              nu=nu, delta=delta,
                                              step=step,normalized=normalized,lambda=lambda)
      
      indtest.lambda<-part.list.lambda$testMat[,fo] #get index of test data for the current fold
      r.test <- r[indtest.lambda,]
      vec.re.lambda[fo] = negCEQ(r.test, l1fast.lambda$betalist[[1]])
    }
  
    # Calculate objective function
    SR.test = mean(vec.re.lambda)
    
    return (SR.test)
  }
  
  # Fit initial GP model with initial observations. Later will be updated inside the optimizer
  model = km( formula=~1, design= as.data.frame(init_observations[,1]), response= as.data.frame(init_observations[,2]),
              covtype="matern5_2", coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
              noise.var = rep(noise.var,nrow(init_observations)),
              estim.method="MLE",optim.method = "BFGS", control = list(trace = FALSE))
  
  # Run Bayes Opt
  cat("BO starts at ",  as.character(Sys.time()),"\n")
  sink('null')
  optim.param = list()
  #optim.param$plugin.type = "min_km_mean" #Note: plugin.type is needed for EI criterion only
  optim.result = noisy.optimizer_withStopping(optim.crit = "EIwithVar", optim.param = optim.param,
                                 model = model,n.ite=15,
                                 funnoise = l1fast_sr,noise.var = noise.var,
                                 lower = lambda.min, upper = lambda.max,
                                 CovReEstimate = TRUE, NoiseReEstimate = TRUE)
  sink()
  cat("BO ends at ",  as.character(Sys.time()),"\n")
  
  cat("OPTIM LAMBDA IS ", optim.result$best.x,"\n")
  
  # Get final GP hyperparameter settings
  km.coef = coef(optim.result$model)
  noise.var = optim.result$history.noise.var[length(optim.result$history.noise.var)]
  
  return (list(best.x= optim.result$best.x, ybestvec= optim.result$ybestvec, kmmeanvec = optim.result$kmmeanvec,
               history.x = optim.result$history.x, history.y = optim.result$history.y,
               coef.trend = km.coef$trend, coef.cov = km.coef$range, coef.var = km.coef$sd2, noise.var =noise.var))
}



# Bayes Opt for LPO with Bayes Averaging Estimates with EIVar as acquisition function and Time Series adjustments
BO_noisy_optimizer_l1fastbayesave <- function(r,prob_old, mu_old, sigma_old, k_old, delta_old,
                                              random.subset = FALSE, 
                                              use.prevlambda = TRUE, prevlambda= NULL,
                                              coef.trend = NULL, coef.cov = NULL, coef.var = NULL,noise.var = NULL,
                                              random.start=FALSE,
                                              delta,step,fold,normalized = FALSE,
                                              lambda=NULL,logspaced=TRUE,coef=0.8)
{
  p<-ncol(r)
  t<-nrow(r)
  sample_cov = cov(r)
  
  if (is.null(lambda)) { 
    coef=ifelse(coef>1,1,coef)
    lambda.max=coef*max(abs(delta))
    lambda.min = lambda.max/50
  }
  
  # Define bounds of the search
  bounds = list(lambda = c(lambda.min,lambda.max))
  
  # Number of initial observations
  n.initlambda = 15
  evaluate.initlambda = FALSE
  add.flag=FALSE
  
  # Generate few lambda values as initial observations
  if (use.prevlambda) # use previous day's results
  {
    init_observations = as.data.frame(prevlambda)
    cat("Using previous lambda with length ", nrow(init_observations), "\n")
    if (nrow(init_observations) == 1)
    {
      init.lambda = seq(lambda.min, lambda.max, length.out =  n.initlambda)
      evaluate.initlambda = TRUE
      add.flag = TRUE
      cat("Evaluating more lambdas with nlambda ", n.initlambda, "\n")
    }
  }
  else # Pick new start points
  {
    if (random.start==TRUE)
    {
      init.lambda = runif(n= n.initlambda, min = lambda.min, max = lambda.max)
      evaluate.initlambda = TRUE
      cat('Using random start points \n')
    }
    else
    {
      init.lambda = seq(lambda.min, lambda.max, length.out = n.initlambda)
      evaluate.initlambda = TRUE
      cat('Using uniformly spaced start points \n')
    }
  }
  
  # Note: For Bayes Averaging, only 1 pseudo-test set is used (taken to be last f points). Since must preserve time dependency.
  if (evaluate.initlambda)
  {
    returns.mat <- matrix(0,nrow=fold,ncol= n.initlambda)
    not.fea = FALSE
    
    for (f in fold:1)
    {
      # Use data from time 1 to t-f as training
      indtrain = 1:(t-f)
      n.train = length(indtrain)
      r.train = r[indtrain,]
      
      # Use data at time t-f + 1 as test point in this fold
      indtest = (t-f) + 1
      r.test = r[indtest,]
      
      # Compute estimates for test point: E(r_t|F_t-1) and V(r_t|F_t-1)
      # Estimate for E(r_t|F_t-1)
      mu_hat_t = rep(0,p)
      for (m in 1:n.train)
      {
        mu_hat_t = mu_hat_t + mu_old[[m]]*prob_old[m]
      }
      # Estimate for V(r_t|F_t-1)
      sigma_hat_t = matrix(0,p,p) 
      for (m in 1:n.train)
      {
        sigma_hat_t = sigma_hat_t + ( ((1+ k_old[m])/k_old[m])*sigma_old[[m]] 
                                      + mu_old[[m]]%*%t(mu_old[[m]]) )*prob_old[m]
      }
      sigma_hat_t = sigma_hat_t - mu_hat_t%*%t(mu_hat_t)
      
      # Compute LPO weights for test point
      sigma_hat_t = sigma_hat_t / step 
      if(normalized){
        d=diag(sigma_hat_t)^(-1/2);D=diag(d)
        # normalize delta
        delta=d*delta;maxdelta=max(abs(delta));delta=delta/maxdelta 
        # normalize sigma: becomes correlation matrix
        sigma_hat_t=D%*%sigma_hat_t%*%D
      }

      betalist<-vector("list", n.initlambda) 
      # betalist is initialized with NULL and will be updated to store optimal eta associated with certain lambda.
      # if a lambda yields no feasible sol, the entry in betalist will still be the default NULL
      nzlist<-vector("list",  n.initlambda)
      card<-rep(p, n.initlambda)
      fea<-rep(NA, n.initlambda)
      status<-rep(0, n.initlambda)
      betaeta<-rep(0, n.initlambda)
      diff<-rep(0, n.initlambda)
      gross<-rep(0, n.initlambda)
      if(t<p){
        A<-cbind(sigma_hat_t,-sigma_hat_t);A<-rbind(A,-A);b=c(delta,-delta);c=rep(-1,2*p);bbar=rep(1,2*p)
        for (j in  n.initlambda:1) { # start with maximum lambda and decrease 
          temp<-fastlpnocat(c,A,b+bbar*init.lambda[j]) 
          status[j]=temp$status
          if(status[j]==1){init.lambda=init.lambda[-(1:j)];break} # if no feasible sol found, stop and output only lambda from previous iterations
          betalist[[j]]<-temp$sol[1:p]-temp$sol[-(1:p)] #optimal eta for this lambda
          fea[j]=max(abs(sigma_hat_t%*%betalist[[j]]-delta))-init.lambda[j] # check constraint
          if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}  # revert back from previous normalizations
          if(fea[j]>1e-03){init.lambda=init.lambda[-(1:j)];break} # if violated by > 1e-3 stop and output only lambda from previous iterations
          nzlist[[j]]=which(betalist[[j]]!=0) #indices where value nonzero for this eta
          card[j]=length(nzlist[[j]]) #count number of nonzero = sparsity of this eta (cardinality of this portfolio)
          betaeta[j]=t(delta)%*%betalist[[j]]
          diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/t)/step)%*%betalist[[j]]
          gross[j]=sum(abs(betalist[[j]])) #l1 norm of eta
          # if(card[j]>p*0.15){init.lambda=init.lambda[-(1:j)];break}
        }
      } else {
        for (j in  n.initlambda:1) {
          temp<-lpsimplex(sigma_hat_t,delta,init.lambda[j])
          status[j]=temp$status
          if(status[j]==2){init.lambda=init.lambda[-(1:j)];break} 
          betalist[[j]]<-temp$sol #optimal eta?
          fea[j]=max(abs(sigma_hat_t%*%betalist[[j]]-delta))-init.lambda[j]
          if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}
          if(fea[j]>1e-03){init.lambda=init.lambda[-(1:j)];break}
          nzlist[[j]]=which(betalist[[j]]!=0)
          card[j]=length(nzlist[[j]])
          betaeta[j]=t(delta)%*%betalist[[j]]
          diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/t)/step)%*%betalist[[j]]
          gross[j]=sum(abs(betalist[[j]])) 
          # if(card[j]>p*0.15){init.lambda=init.lambda[-(1:j)];break}
        }
      }
      
      #cat("nlambda before fold",f,"is", n.initlambda,"\n")
      #cat("length lambda after fold",f,"is",length(init.lambda),"\n")
      if (length(init.lambda) ==0)
      {
        cat("no feasible lambda \n")
        add.flag = FALSE
        not.fea = TRUE
        break
      }
      
      if ((n.initlambda-length(init.lambda)) > 0)
      {
        returns.mat = as.matrix(returns.mat[,-(1:(n.initlambda-length(init.lambda)))])
        betalist = betalist[-(1:(n.initlambda-length(init.lambda)))]
      }
      #cat("lossre ncol", ncol(returns.mat),"\n")
      #cat("betalist length",length(betalist),"\n")
      n.initlambda = length(init.lambda)
      #cat("nlambda after fold",f,"is",n.initlambda,"\n")
      for (jl in 1:n.initlambda) 
      {
        returns.mat[f,jl]<- r.test%*%(betalist[[jl]]/ sum(betalist[[jl]]))
      }
      
    } ##### END FOR LOOP FOLD
    
    # Compute SR
    returns.mean<-apply(returns.mat, 2, mean)
    #returns.sd <- apply(returns.mat, 2, sd)
    #f.initlambda <- -(returns.mean/returns.sd)
    returns.var = apply(returns.mat,2,var)
    f.initlambda = -(returns.mean - 0.5*returns.var)
    
    # Create data frame for storing initial observations (i.e. lambda and corresponding f(lambda))
    if (!not.fea)
    {
      init_observations = data.frame(init.lambda, f.initlambda)
    }
    
  } ##### END if evaluate lambda
  
  
  if (add.flag) # bind prev and new lambda
  {
    names(init_observations) <- c("lambda","Value")
    prevdf = as.data.frame(prevlambda)
    if (nrow(prevdf)==1)
    {
      init_observations = init_observations
    }
    else
    {
      names(prevdf) <- c("lambda","Value")
      init_observations = rbind(prevdf, init_observations)
    }
  }
  
  # Rename columns
  names(init_observations) <- c("lambda","Value")
  
  # Function to be optimized by BayesOpt 
  l1fast_sr <- function(lambda)
  {
    n.initlambda = length(lambda)
    init.lambda = lambda
    
    vec.re.lambda = rep(0, fold)
    
    for (f in fold:1)
    {
      # Use data from time 1 to t-f as training
      indtrain = 1:(t-f)
      n.train = length(indtrain)
      r.train = r[indtrain,]
      
      # Use data at time t-f + 1 as test point in this fold
      indtest = (t-f) + 1
      r.test = r[indtest,]
      
      # Compute estimates for test point: E(r_t|F_t-1) and V(r_t|F_t-1)
      # Estimate for E(r_t|F_t-1)
      mu_hat_t = rep(0,p)
      for (m in 1:n.train)
      {
        mu_hat_t = mu_hat_t + mu_old[[m]]*prob_old[m]
      }
      # Estimate for V(r_t|F_t-1)
      sigma_hat_t = matrix(0,p,p) 
      for (m in 1:n.train)
      {
        sigma_hat_t = sigma_hat_t + ( ((1+ k_old[m])/k_old[m])*sigma_old[[m]] 
                                      + mu_old[[m]]%*%t(mu_old[[m]]) )*prob_old[m]
      }
      sigma_hat_t = sigma_hat_t - mu_hat_t%*%t(mu_hat_t)
      
      # Compute LPO weights for test point
      sigma_hat_t = sigma_hat_t / step 
      if(normalized){
        d=diag(sigma_hat_t)^(-1/2);D=diag(d)
        # normalize delta
        delta=d*delta;maxdelta=max(abs(delta));delta=delta/maxdelta 
        # normalize sigma: becomes correlation matrix
        sigma_hat_t=D%*%sigma_hat_t%*%D
      }
      # vector("list", n) creates vector of length n with null values (a list of lists) 
      betalist<-vector("list", n.initlambda) 
      # betalist is initialized with NULL and will be updated to store optimal eta associated with certain lambda.
      # if a lambda yields no feasible sol, the entry in betalist will still be the default NULL
      nzlist<-vector("list",  n.initlambda)
      card<-rep(p, n.initlambda)
      fea<-rep(NA, n.initlambda)
      status<-rep(0, n.initlambda)
      betaeta<-rep(0, n.initlambda)
      diff<-rep(0, n.initlambda)
      gross<-rep(0, n.initlambda)
      if(t<p){
        A<-cbind(sigma_hat_t,-sigma_hat_t);A<-rbind(A,-A);b=c(delta,-delta);c=rep(-1,2*p);bbar=rep(1,2*p)
        for (j in  n.initlambda:1) { # start with maximum lambda and decrease 
          temp<-fastlpnocat(c,A,b+bbar*init.lambda[j]) 
          status[j]=temp$status
          if(status[j]==1){init.lambda=init.lambda[-(1:j)];break} # if no feasible sol found, stop and output only lambda from previous iterations
          betalist[[j]]<-temp$sol[1:p]-temp$sol[-(1:p)] #optimal eta for this lambda
          fea[j]=max(abs(sigma_hat_t%*%betalist[[j]]-delta))-init.lambda[j] # check constraint
          if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}  # revert back from previous normalizations
          if(fea[j]>1e-03){init.lambda=init.lambda[-(1:j)];break} # if violated by > 1e-3 stop and output only lambda from previous iterations
          nzlist[[j]]=which(betalist[[j]]!=0) #indices where value nonzero for this eta
          card[j]=length(nzlist[[j]]) #count number of nonzero = sparsity of this eta (cardinality of this portfolio)
          betaeta[j]=t(delta)%*%betalist[[j]]
          diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/t)/step)%*%betalist[[j]]
          gross[j]=sum(abs(betalist[[j]])) #l1 norm of eta
          # if(card[j]>p*0.15){init.lambda=init.lambda[-(1:j)];break}
        }
      } else {
        for (j in  n.initlambda:1) {
          temp<-lpsimplex(sigma_hat_t,delta,init.lambda[j])
          status[j]=temp$status
          if(status[j]==2){init.lambda=init.lambda[-(1:j)];break} 
          betalist[[j]]<-temp$sol #optimal eta?
          fea[j]=max(abs(sigma_hat_t%*%betalist[[j]]-delta))-init.lambda[j]
          if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}
          if(fea[j]>1e-03){init.lambda=init.lambda[-(1:j)];break}
          nzlist[[j]]=which(betalist[[j]]!=0)
          card[j]=length(nzlist[[j]])
          betaeta[j]=t(delta)%*%betalist[[j]]
          diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/t)/step)%*%betalist[[j]]
          gross[j]=sum(abs(betalist[[j]])) 
          # if(card[j]>p*0.15){init.lambda=init.lambda[-(1:j)];break}
        }
      }
      
      vec.re.lambda[f] = r.test%*%(betalist[[1]]/ sum(betalist[[1]]))
      
    } ## END OF FOR LOOP FOLD
    
    # Calculate SR 
    #SR.test = mean(vec.re.lambda)/sd(vec.re.lambda)
    SR.test = mean(vec.re.lambda) - 0.5*var(vec.re.lambda)

    return (-SR.test)
  }
  
  # Number of iterations in Bayes Opt
  n.ite = 15
  
  # Fit initial GP model with initial observations. Later will be updated inside the optimizer
  if(exists("model")) rm(model)
  try (model <-  km( formula=~1, design= as.data.frame(init_observations[,1]), response= as.data.frame(init_observations[,2]),
                     covtype="matern5_2", coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                     noise.var = rep(noise.var,nrow(init_observations)),
                     estim.method="MLE",optim.method = "BFGS", control = list(trace = FALSE))
  )
  
  if(!exists("model"))
  { cat("Error occured during model update")
    ind.best = which.min(init_observations[,2])
    best.x = init_observations[ind.best,1]
    ybestvec = rep(0, n.ite)
    kmmeanvec = rep(0, n.ite)
    history.x = t(as.data.frame(init_observations[,1]))
    history.y = as.data.frame(init_observations[,2])
    coef.trend = coef.trend
    coef.cov = coef.cov
    coef.var = coef.var
    noise.var = noise.var
  }
  else
  {
    # Run Bayes Opt
    cat("BO starts at ",  as.character(Sys.time()),"\n")
    sink('null')
    optim.param = list()
    #optim.param$plugin.type = "min_km_mean" #Note: plugin.type is needed for EI criterion only
    optim.result = noisy.optimizer_withStopping(optim.crit = "EIwithVar", optim.param = optim.param,
                                                model = model,n.ite= n.ite,
                                                funnoise = l1fast_sr,noise.var = noise.var,
                                                lower = lambda.min, upper = lambda.max,
                                                CovReEstimate = TRUE, NoiseReEstimate = TRUE
    )
    sink()
    cat("BO ends at ",  as.character(Sys.time()),"\n")
    
    cat("OPTIM LAMBDA IS ", optim.result$best.x,"\n")

    best.x= optim.result$best.x
    ybestvec= optim.result$ybestvec
    kmmeanvec = optim.result$kmmeanvec
    
    # Get visited points
    history.x = optim.result$history.x
    history.y = optim.result$history.y
    
    # Get final GP hyperparameter settings
    km.coef = coef(optim.result$model)
    coef.trend = km.coef$trend
    coef.cov = km.coef$range
    coef.var = km.coef$sd2
    noise.var = optim.result$history.noise.var[length(optim.result$history.noise.var)]
  }
  
  return (list(best.x= best.x, ybestvec=ybestvec, kmmeanvec = kmmeanvec,
               history.x = history.x, history.y = history.y,
               coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var, noise.var =noise.var))
}