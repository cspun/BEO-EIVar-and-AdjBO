############################################################
# This file contains functions to evaluate out-of-sample performance of CV
# when tuning LPO with Sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate.
############################################################


# Evaluate out of sample performance of LPO with sample estimate (tuned with CV)
negSR.eval.l1fast <- function(r,n.train,n.test,delta,step,normalized=FALSE,
                              perturb=FALSE,lambda=NULL,logspaced=TRUE,coef=0.8,
                              nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train, n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  minval.list = vector("list", n.test+1)
  
  for (j in 1: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    
    # Get optimal lambda using CV
    lambda.result = cv.negSR.l1fast(r.train,delta,step,fold=4,normalized = normalized, perturb=perturb,
                                    lambda=lambda,logspaced=logspaced, coef=0.8,
                                    nlambda=30)
    optlambda = lambda.result$lambdaopt
    lambda.list[j] = optlambda
    minval.list[[j]] = lambda.result$vec.minval
    #cat("opt lambda", optlambda, "\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast(r.train,delta,step,normalized=normalized,
                              perturb=perturb,lambda=optlambda,coef=coef) 
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j <= n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%w.optlambda
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list = lambda.list, time.list = time.list, minval.list=minval.list))
}



# Evaluate out of sample performance of LPO with Bayes-Stein estimate (tuned with CV)
negSR.eval.l1fastbayesstein <- function(r,mu_0=NULL, tau=NULL,
                                        n.train,n.test,delta,step,shrinkage = FALSE, normalized=FALSE,
                                        lambda=NULL,logspaced=TRUE,coef=0.8,
                                        nlambda =ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train, n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  minval.list = vector("list", n.test+1)
  
  if (is.null(mu_0) |  is.null(tau)) #Prior values not specified -> estimate using data directly via Empirical Bayes
  {
    flag = TRUE
  }
  else
  {
    flag= FALSE
  }
  
  for (j in 1:(n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point

    if (flag)
    {
      # Get prior parameters from data (empirical bayes)
      BSprior = getbayessteinprior(r.train, shrinkage=shrinkage) 
      mu_0 = BSprior$mu_0
      tau = BSprior$tau
    }

    # Get optimal lambda using CV
    lambda.result = cv.negSR.l1fastbayesstein(r.train, mu_0=mu_0, tau=tau,
                                              delta,step,fold=4,shrinkage = shrinkage,
                                              normalized=normalized,
                                              lambda=lambda, nlambda=30)
    optlambda = lambda.result$lambdaopt
    lambda.list[j] = optlambda
    minval.list[[j]] = lambda.result$vec.minval
    cat("opt lambda", optlambda, "\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast.bayesstein(r.train,mu_0=mu_0,tau=tau,delta=delta,step=step,
                                         shrinkage = shrinkage, normalized=normalized,lambda=optlambda)
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j <= n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%w.optlambda
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list = lambda.list, time.list=time.list, minval.list=minval.list))
}


# Evaluate out of sample performance of LPO with Bayes estimate based on asset pricing model (tuned with CV)
negSR.eval.l1fastbayesassetpricing <- function(r, f, sigma_alpha, nu,
                                               n.train,n.test,delta,step,normalized=FALSE,
                                               lambda=NULL,logspaced=TRUE,coef=0.8,
                                               nlambda=ifelse(is.null(lambda),100,length(lambda)),
                                               lambda.max=NULL)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train, n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  minval.list = vector("list", n.test+1)
  
  for (j in 1:(n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data for asset returns
    f.train = f[indtrain,] #training data for factor returns
    indtest = rollwin.ind$test[,j] # get index of test point

    # Get optimal lambda using CV
    lambda.result = cv.negSR.l1fastbayesassetpricing(r.train,f.train,sigma_alpha = sigma_alpha,nu=nu,
                                                     delta,step,fold=4,normalized=normalized,
                                                     lambda=lambda, nlambda=30)
    optlambda = lambda.result$lambdaopt
    lambda.list[j] = optlambda
    minval.list[[j]] = lambda.result$vec.minval
    cat("opt lambda", optlambda, "\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast.bayesassetpricing(r.train,f.train, sigma_alpha = sigma_alpha,nu=nu,
                                                delta=delta,step=step,normalized=normalized,
                                                lambda=optlambda)
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j <= n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%w.optlambda
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list = lambda.list, time.list = time.list, minval.list=minval.list))
}


# Evaluate out of sample performance of LPO with Bayes Averaging estimate (tuned with CV)
negSR.eval.l1fastbayesaveraging <- function(r, random.subset=FALSE,
                                            n.test,delta,step,normalized=FALSE,
                                            lambda=NULL,logspaced=TRUE,coef=0.8,
                                            nlambda
                                            =ifelse(is.null(lambda),100,length(lambda)),
                                            lambda.max=NULL )
{
  cat("start at", as.character(Sys.time()),"\n")
  p = ncol(r)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  minval.list = vector("list", n.test+1)
  
  # At time t=1, the first model Model 1 is born, but since there's no historical data, set fixed pre-specified value
  # Later use burn-in period so the pre-specified values have negligible effects
  k_1_1 = 2
  delta_1_1 = 2
  mu_1_1 = r[1,] / k_1_1
  sigma_1_1 = (k_1_1*0.0001*diag(p) + r[1,]%*%t(r[1,]))/(delta_1_1*k_1_1)
  p_1_1 = 1
  
  # Now run from t=2
  prob_old = c(p_1_1)
  mu_old = list(mu_1_1)
  sigma_old = list(sigma_1_1)
  k_old = c(k_1_1)
  delta_old = c(delta_1_1)
  
  for (j in 2:(n.test+1))
  {
    start.time = Sys.time()
    indtrain = 1:j  #get index of training data 
    r.train = r[indtrain,] # training data for asset returns
    indtest = j+1 # get index of test point

    if (j<=19) #Burn-in period (results from this period are ignored)
    {
      lambda.list[j] = 0.8
      #cat("opt lambda", lambda.list[j], "\n")
    }
    else
    {
      # Get optimal lambda using CV
      lambda.result =cv.negSR.l1fastbayesaveraging(r.train, prob_old = prob_old, 
                                                   mu_old = mu_old, sigma_old =
                                                     sigma_old, 
                                                   k_old = k_old, delta_old = delta_old,
                                                   random.subset = random.subset,
                                                   delta=delta,step=step, fold=4,
                                                   normalized = normalized,lambda=lambda, nlambda=30)
      lambda.list[j] = lambda.result$lambdaopt
      minval.list[[j]] = lambda.result$vec.minval
      cat("opt lambda", lambda.list[j], "\n")
    }
    
    # Compute eta associated with optimal lambda
    l1fast.optlambda = l1fast.bayesaveraging(r.train, prob_old = prob_old, mu_old = mu_old,
                                             sigma_old =sigma_old, k_old = k_old, 
                                             delta_old = delta_old,
                                             random.subset = random.subset,seedval = j,
                                             delta=delta, step=step,
                                             normalized = normalized,lambda=lambda.list[j])
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Update 
    prob_old = l1fast.optlambda$prob_old
    mu_old = l1fast.optlambda$mu_old
    sigma_old = l1fast.optlambda$sigma_old
    k_old = l1fast.optlambda$k_old
    delta_old = l1fast.optlambda$delta_old
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j <= n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = tryCatch({r.test%*%w.optlambda}, error = function(c){return (NA)})
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list = lambda.list, time.list = time.list, minval.list=minval.list))
}


# Evaluate out of sample performance of 1/N strategy (equally weighted)
negSR.eval.EW <- function(r, n.train, n.test)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  p = ncol(r)
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train, n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  time.list = rep(0, n.test+1)
  
  for (j in 1: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    
    w.list[[j]] = rep(1/p, p)
    
    # Compute portolio return at test point
    if (j<=n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%rep(1/p, p) 
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, time.list=time.list))
}


# Evaluate out of sample performance of Plugin strategy
eval.plugin <- function(r, n.train, n.test, step)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  p = ncol(r)
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train, n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  time.list = rep(0, n.test+1)
  
  for (j in 1: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    
    # Compute weights 
    Sigmahat=cov(r.train)*(1-1/n.train)/step
    delta = rep(1,p)
    temp<-try(as.vector(ginv(Sigmahat)%*%delta),TRUE)
    if(is.vector(temp)){eta<-temp} 
    w.list[[j]] = eta/sum(eta)
    
    # Compute portolio return at test point
    if (j<=n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = crossprod(r.test, w.list[[j]])
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, time.list=time.list))
}