############################################################
# This file contains functions to find true optimal lambda 
# by incorporating future asset prices.
############################################################


# Calculate optimal lambda in backtest assuming future prices are known for LPO with sample estimates 
TrueOptimalLambda.lposample <- function(r, n.train, n.test,delta,step,normalized=FALSE,
                                        perturb=FALSE,lambda=NULL,logspaced=TRUE,coef=0.8,
                                        nlambda = ifelse(is.null(lambda),30,length(lambda)),lambda.max=NULL)
{
  # check whether lambda vector supplied is null, if null generate lambda
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  rollwin.ind = rollwindows.part(n.train, n.test)
  returns.mat = matrix(NA, nrow = n.test, ncol= nlambda)
  optlambda.list = rep(0,n.test) 
  
  for (j in 1:n.test)
  {
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    r.test = r[indtest,] #test point
    
    # Run l1 fast with whole training data 
    l1fast.init = l1fast(r.train, delta,step,normalized = normalized, perturb=perturb,
                         lambda=lambda, shrink=TRUE)
    
    for (l in 1:nlambda)
    {
      # compute return using (assuming) known r.test
      if (is.null(l1fast.init$betalist[[l]])) #case if lambda is not feasible
      {
        returns.mat[j,l] = NA
      }
      else
      {
        eta.lambda = l1fast.init$betalist[[l]]
        returns.mat[j,l] = r.test%*% (eta.lambda/sum(eta.lambda))
      }
    }
    
    if (j>1)
    {
      #Compute rolling sd
      roll.sd = apply(returns.mat,2, sd, na.rm=TRUE) #ignore NA values when computing sd
      #Compute rolling mean
      roll.mean = apply(returns.mat,2, mean, na.rm=TRUE) #ignore NA values when computing var
      #Find lambda that has max SR
      ind<-which.max(roll.mean/roll.sd) 
      lambdaopt<-lambda[ind]
    }
    else
    {
      #Find lambda that has max return when there is only 1 return 
      ind = which.max(returns.mat[1,])
      lambdaopt = lambda[ind]
    }
    #cat("lambda opt backtest ", lambdaopt, "\n")
    optlambda.list[j] = lambdaopt
  }
  
  return (list(optlambda.list=optlambda.list, returns.mat=returns.mat))
}




# Calculate optimal lambda in backtest assuming future prices are known for LPO with Bayes Stein estimates
TrueOptimalLambda.lpobs <- function(r, mu_0=NULL, tau=NULL,n.train, n.test,delta,step,shrinkage = FALSE,
                                    normalized=FALSE,
                                    lambda=NULL,logspaced=TRUE,coef=0.8,
                                    nlambda = ifelse(is.null(lambda),30,length(lambda)),lambda.max=NULL)
{
  # check whether lambda vector supplied is null, if null generate lambda
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  rollwin.ind = rollwindows.part(n.train, n.test)
  returns.mat = matrix(NA, nrow = n.test, ncol= nlambda)
  optlambda.list = rep(0,n.test) 
  
  if (is.null(mu_0) |  is.null(tau))
  {
    flag = TRUE
  }
  else
  {
    flag= FALSE
  }
  
  for (j in 1:n.test)
  {
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    r.test = r[indtest,] #test point
    if (flag)
    {
      #get prior parameters from data (empirical bayes)
      BSprior = getbayessteinprior(r.train, shrinkage=shrinkage) 
      mu_0 = BSprior$mu_0
      tau = BSprior$tau
    }
    
    # Run l1 fast with whole training data 
    l1fast.init = l1fast.bayesstein(r.train,mu_0=mu_0,tau=tau,delta=delta,step=step,shrinkage = shrinkage,
                                    normalized=normalized,lambda=lambda)
    for (l in 1:nlambda)
    {
      # compute return using (assuming) known r.test
      if (is.null(l1fast.init$betalist[[l]])) #case if lambda is not feasible
      {
        returns.mat[j,l] = NA
      }
      else
      {
        eta.lambda = l1fast.init$betalist[[l]]
        returns.mat[j,l] = r.test%*% (eta.lambda/sum(eta.lambda))
      }
    }
    
    if (j>1)
    {
      #Compute rolling sd
      #roll.sd = apply(returns.mat,2, sd, na.rm=TRUE) #ignore NA values when computing sd
      roll.var = apply(returns.mat,2,var, na.rm = TRUE)
      #Compute rolling mean
      roll.mean = apply(returns.mat,2, mean, na.rm=TRUE) #ignore NA values when computing var
      #Find lambda that has max SR
      #ind<-which.max(roll.mean/roll.sd) 
      ind = which.max(roll.mean - 0.5*roll.var)
      lambdaopt<-lambda[ind]
    }
    else
    {
      #Find lambda that has max return when there is only 1 return 
      ind = which.max(returns.mat[1,])
      lambdaopt = lambda[ind]
    }
    optlambda.list[j] = lambdaopt
    cat("test point",j,"done at",as.character(Sys.time()),"\n")
  }
  return (list(optlambda.list=optlambda.list, returns.mat=returns.mat))
}




# Calculate optimal lambda in backtest assuming future prices are known for LPO with Bayes asset pricing estimates
TrueOptimalLambda.lpoassetp <- function(r,f, sigma_alpha, nu,n.train, n.test,delta,step,
                                        normalized=FALSE,lambda=NULL,
                                        logspaced=TRUE,coef=0.8,
                                        nlambda = ifelse(is.null(lambda),30,length(lambda)),lambda.max=NULL)
{
  # check whether lambda vector supplied is null, if null generate lambda
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  rollwin.ind = rollwindows.part(n.train, n.test)
  returns.mat = matrix(NA, nrow = n.test, ncol= nlambda)
  optlambda.list = rep(0,n.test) 
  
  for (j in 1:n.test)
  {
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data for asset returns
    f.train = f[indtrain,] #training data for factor returns
    indtest = rollwin.ind$test[,j] # get index of test point
    r.test = r[indtest,] #test point
    
    # Run l1 fast with whole training data first
    l1fast.init = l1fast.bayesassetpricing(r.train,f.train, sigma_alpha = sigma_alpha,nu=nu,
                                           delta=delta,step=step,normalized=normalized,
                                           lambda=lambda)
    for (l in 1:nlambda)
    {
      # compute return using (assuming) known r.test
      if (is.null(l1fast.init$betalist[[l]])) #case if lambda is not feasible
      {
        returns.mat[j,l] = NA
      }
      else
      {
        eta.lambda = l1fast.init$betalist[[l]]
        returns.mat[j,l] = r.test%*% (eta.lambda/sum(eta.lambda))
      }
    }
    
    if (j>1)
    {
      #Compute rolling sd
      #roll.sd = apply(returns.mat,2, sd, na.rm=TRUE) #ignore NA values when computing sd
      roll.var = apply(returns.mat,2,var, na.rm = TRUE)
      #Compute rolling mean
      roll.mean = apply(returns.mat,2, mean, na.rm=TRUE) #ignore NA values when computing var
      #Find lambda that has max SR
      #ind<-which.max(roll.mean/roll.sd) 
      ind = which.max(roll.mean - 0.5*roll.var)
      lambdaopt<-lambda[ind]
    }
    else
    {
      #Find lambda that has max return when there is only 1 return 
      ind = which.max(returns.mat[1,])
      lambdaopt = lambda[ind]
    }
    optlambda.list[j] = lambdaopt
    cat("test point",j,"done at",as.character(Sys.time()),"\n")
  }
  return (list(optlambda.list=optlambda.list, returns.mat=returns.mat))
}




# Calculate optimal lambda in backtest assuming future prices are known for LPO with Bayes averaging estimates
TrueOptimalLambda.lpobayesave <- function(r,random.subset, n.test,delta,step,
                                          normalized=FALSE,lambda=NULL,
                                          logspaced=TRUE,coef=0.8,
                                          nlambda = ifelse(is.null(lambda),30,length(lambda)),lambda.max=NULL)
{
  # check whether lambda vector supplied is null, if null generate lambda
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  returns.mat = matrix(NA, nrow = n.test, ncol= nlambda)
  optlambda.list = rep(0,n.test)
  p = ncol(r)
  
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
  
  for (j in 2:n.test)
  {
    indtrain = 1:j  #get index of training data 
    r.train = r[indtrain,] # training data for asset returns
    indtest = j+1 # get index of test point
    r.test = r[indtest,] #test point
    l1fast.init = l1fast.bayesaveraging(r.train, prob_old = prob_old, mu_old = mu_old,
                                        sigma_old =sigma_old, k_old = k_old, 
                                        delta_old = delta_old,random.subset = random.subset,
                                        seedval= j, delta=delta,step=step,
                                        normalized = normalized,lambda=lambda)
    # Update 
    prob_old = l1fast.init$prob_old
    mu_old = l1fast.init$mu_old
    sigma_old = l1fast.init$sigma_old
    k_old = l1fast.init$k_old
    delta_old = l1fast.init$delta_old
    
    for (l in 1:nlambda)
    {
      # compute return using (assuming) known r.test
      if (is.null(l1fast.init$betalist[[l]])) #case if lambda is not feasible
      {
        returns.mat[j,l] = NA
      }
      else
      {
        eta.lambda = l1fast.init$betalist[[l]]
        returns.mat[j,l] = r.test%*% (eta.lambda/sum(eta.lambda))
      }
    }
    
    if (j>2)
    {
      #Compute rolling sd
      #roll.sd = apply(returns.mat,2, sd, na.rm=TRUE) #ignore NA values when computing sd
      roll.var = apply(returns.mat,2,var, na.rm = TRUE)
      #Compute rolling mean
      roll.mean = apply(returns.mat,2, mean, na.rm=TRUE) #ignore NA values when computing var
      #Find lambda that has max SR
      #ind<-which.max(roll.mean/roll.sd) 
      ind = which.max(roll.mean - 0.5*roll.var)
      lambdaopt<-lambda[ind]
    }
    else #for burn-in period (will be ignored)
    {
      lambdaopt = 0.8
    }
    optlambda.list[j] = lambdaopt
    cat("test point",j,"done at",as.character(Sys.time()),"\n")
  }
  return (list(optlambda.list=optlambda.list, returns.mat=returns.mat))
}