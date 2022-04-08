############################################################
# This file contains functions to evaluate out-of-sample performance of Bayes Opt (BO)
# when tuning LPO with Sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate.

# Bayes Opt being considered here use EIVar as acquisition function and with time series adjustment (AdjBO-EIVar).
############################################################


# Evaluate out of sample performance of LPO with sample estimate (tuned with BO)
bayesopt.eval.l1fast <- function(r, random.start = FALSE, n.train, n.test, delta,step,normalized=FALSE,
                                 perturb=FALSE,lambda=NULL,logspaced=TRUE,coef=0.8,
                                 nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train= n.train, n.test=n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test + 1)
  time.list = rep(0, n.test +1)
  ybest.vec = vector("list", n.test+1)
  kmmean.vec = vector("list", n.test+1)
  
  for (j in 1: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    
    # Get optimal lambda using Bayes Opt
    if (j==1) #at first time step, no previous results are available. 
    {
      optlambda = BO_noisy_optimizer_l1fast(r.train,  use.prevlambda = FALSE, prevlambda = NULL,
                                            coef.trend = NULL, coef.cov = NULL, coef.var = NULL, noise.var = 1e-10,
                                            random.start = random.start,delta=delta, step=step, fold=4,
                                            normalized = normalized)
      # Get previously visited points to pass to next time step
      prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
      # Get final GP hyperparameter settings to pass to next time step 
      coef.trend = optlambda$coef.trend
      coef.cov = optlambda$coef.cov
      coef.var = optlambda$coef.var
      noise.var = optlambda$noise.var
    }
    else #at subsequent time steps, use previous time step's optimization results as priors. 
    {
      optlambda = BO_noisy_optimizer_l1fast(r.train,use.prevlambda= TRUE, prevlambda = prevlambda,
                                            coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var, noise.var = noise.var,
                                            random.start = random.start,delta=delta, step=step, fold=4, 
                                            normalized = normalized)
      # Get previously visited points to pass to next time step
      prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
      # Get final GP hyperparameter settings to pass to next time step 
      coef.trend = optlambda$coef.trend
      coef.cov = optlambda$coef.cov
      coef.var = optlambda$coef.var
      noise.var = optlambda$noise.var
    }
    
    lambda.list[j] = optlambda$best.x
    ybest.vec[[j]] = optlambda$ybestvec
    kmmean.vec[[j]] = optlambda$kmmeanvec

    #cat("opt lambda", lambda.list[j], "\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast(r.train,delta,step,normalized=normalized,
                              perturb=perturb,lambda=lambda.list[j],coef=coef, shrink=TRUE) 
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j<=n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%w.optlambda
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list = lambda.list, time.list = time.list, ybest.vec=ybest.vec, kmmean.vec=kmmean.vec))
}



# Evaluate out of sample performance of LPO with Bayes Stein estimate (tuned with BO)
bayesopt.eval.l1fastbayesstein <- function(r, random.start= FALSE,mu_0=NULL, tau=NULL,n.train, n.test,
                                           delta,step, shrinkage = FALSE, normalized=FALSE,
                                           lambda=NULL,logspaced=TRUE,coef=0.8,
                                           nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train= n.train, n.test=n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  ybest.vec = vector("list", n.test+1)
  kmmean.vec = vector("list", n.test+1)
  
  if (is.null(mu_0) |  is.null(tau)) # Prior values not specified -> estimate based on data directly (empirical Bayes)
  {
    flag = TRUE
  }
  else
  {
    flag= FALSE
  }
  
  for (j in 1: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    indtest = rollwin.ind$test[,j] # get index of test point
    
    if (flag)
    {
      #get prior parameters from data (empirical bayes)
      BSprior = getbayessteinprior(r.train, shrinkage=shrinkage) 
      mu_0 = BSprior$mu_0
      tau = BSprior$tau
    }

    if (j==1) #at first time step, no previous results are available. 
    {
      optlambda = BO_noisy_optimizer_l1fastbs(r.train, mu_0 = mu_0, tau=tau,
                                              use.prevlambda = FALSE, prevlambda=NULL,
                                              coef.trend = NULL, coef.cov = NULL, coef.var = NULL, noise.var = 1e-10,
                                              random.start = random.start,
                                              delta=delta,
                                              step = step, fold=4, shrinkage=shrinkage,
                                              normalized=normalized)
      # Get previously visited points to pass to next time step
      prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
      # Get final GP hyperparameter settings to pass to next time step
      coef.trend = optlambda$coef.trend
      coef.cov = optlambda$coef.cov
      coef.var = optlambda$coef.var
      noise.var = optlambda$noise.var
    }
    else #at subsequent time steps, use previous time step's optimization results as priors. 
    {
      optlambda = BO_noisy_optimizer_l1fastbs(r.train, mu_0 = mu_0, tau=tau,
                                              use.prevlambda = TRUE, prevlambda = prevlambda,
                                              random.start = random.start,
                                              coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var, noise.var = noise.var,
                                              delta=delta,
                                              step = step, fold=4, shrinkage=shrinkage,
                                              normalized=normalized)
      # Get previously visited points to pass to next time step
      prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
      # Get final GP hyperparameter settings to pass to next time step
      coef.trend = optlambda$coef.trend
      coef.cov = optlambda$coef.cov
      coef.var = optlambda$coef.var
      noise.var = optlambda$noise.var
    }
 
    lambda.list[j] = optlambda$best.x
    ybest.vec[[j]] = optlambda$ybestvec
    kmmean.vec[[j]] = optlambda$kmmeanvec
    #cat("opt lambda", lambda.list[j], "\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast.bayesstein(r.train,mu_0=mu_0,tau=tau,delta=delta,step=step,
                                         shrinkage = shrinkage,
                                         normalized=normalized,lambda=lambda.list[j])
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j<=n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%w.optlambda
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list = lambda.list, time.list=time.list, ybest.vec=ybest.vec, kmmean.vec=kmmean.vec))
}



# Evaluate out of sample performance of LPO with Bayes estimate based on asset pricing model (tuned with BO)
bayesopt.eval.l1fastbayesassetp <- function(r,f, sigma_alpha, nu,random.start= FALSE, n.train, n.test,
                                            delta,step,normalized=FALSE,
                                            lambda=NULL,logspaced=TRUE,coef=0.8,
                                            nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL)
{
  cat("start at", as.character(Sys.time()),"\n")
  
  # get indices of rolling window (for matrix of returns)
  rollwin.ind = rollwindows.part(n.train= n.train, n.test=n.test+1)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  ybest.vec = vector("list", n.test+1)
  kmmean.vec = vector("list", n.test+1)
  
  for (j in 1: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = rollwin.ind$train[,j] #get index of training data 
    r.train = r[indtrain,] # training data
    f.train = f[indtrain,] #training data for factor returns
    indtest = rollwin.ind$test[,j] # get index of test point
    
    if (j==1) #at first time step, no previous results are available. 
    {
      optlambda = BO_noisy_optimizer_l1fastassetp(r.train, f.train, sigma_alpha = sigma_alpha,
                                                  nu=nu, use.prevlambda = FALSE, prevlambda= NULL,
                                                  coef.trend = NULL, coef.cov = NULL, coef.var = NULL, noise.var = 1e-10,
                                                  random.start=random.start,
                                                  delta=delta, step=step,
                                                  fold=4,normalized = normalized)
      # Get previously visited points to pass to next time step
      prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y)
      # Get final GP hyperparameter settings to pass to next time step
      coef.trend = optlambda$coef.trend
      coef.cov = optlambda$coef.cov
      coef.var = optlambda$coef.var
      noise.var = optlambda$noise.var
    }
    else #at subsequent time steps, use previous time step's optimization results as priors. 
    {
      optlambda = BO_noisy_optimizer_l1fastassetp(r.train, f.train, sigma_alpha = sigma_alpha,
                                                  nu=nu, use.prevlambda = FALSE, prevlambda= prevlambda,
                                                  coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var, noise.var = noise.var,
                                                  random.start=random.start,
                                                  delta=delta, step=step,
                                                  fold=4,normalized = normalized)
      # Get previously visited points to pass to next time step
      prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
      # Get final GP hyperparameter settings to pass to next time step
      coef.trend = optlambda$coef.trend
      coef.cov = optlambda$coef.cov
      coef.var = optlambda$coef.var
      noise.var = optlambda$noise.var
    }

    lambda.list[j] = optlambda$best.x
    ybest.vec[[j]] = optlambda$ybestvec
    kmmean.vec[[j]] = optlambda$kmmeanvec
    cat("opt lambda", lambda.list[j], "\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast.bayesassetpricing(r.train,f.train, 
                                                sigma_alpha = sigma_alpha,nu=nu,
                                                delta=delta,step=step,normalized=normalized,
                                                lambda=lambda.list[j]) 
    eta.optlambda = l1fast.optlambda$betalist[[1]]
    
    # Compute weights
    w.optlambda = eta.optlambda/ sum(eta.optlambda)
    w.list[[j]] = w.optlambda
    
    # Compute portolio return at test point
    if (j <=n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = r.test%*%w.optlambda
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list=lambda.list, time.list=time.list, ybest.vec=ybest.vec, kmmean.vec=kmmean.vec))
}



# Evaluate out of sample performance of LPO with Bayes Averaging estimate (tuned with BO)
bayesopt.eval.l1fastbayesave <- function(r, random.start= FALSE, random.subset=FALSE, n.test,
                                         delta,step,normalized=FALSE,
                                         lambda=NULL,logspaced=TRUE,coef=0.8)
{
  cat("start at", as.character(Sys.time()),"\n")
  p = ncol(r)
  
  w.list = vector("list", n.test+1)
  rp.list = rep(NA, n.test)
  lambda.list = rep(0, n.test+1)
  time.list = rep(0, n.test+1)
  ybest.vec = vector("list", n.test+1)
  kmmean.vec = vector("list", n.test+1)
  
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
  
  for (j in 2: (n.test+1))
  {
    start.time = Sys.time()
    indtrain = 1:j  #get index of training data 
    r.train = r[indtrain,] # training data for asset returns
    indtest = j+1 # get index of test point
    
    if (j<=19) #burn-in period (results on this period are ignored)
    {
      lambda.list[j] = 0.8
    }
    else
    {
      # Get optimal lambda using Bayes Opt
      if (j==20) #at first time step, no previous results are available. 
      {
        optlambda = BO_noisy_optimizer_l1fastbayesave(r.train, prob_old = prob_old, mu_old = mu_old,
                                                      sigma_old =sigma_old, k_old = k_old,
                                                      delta_old = delta_old,
                                                      random.subset = random.subset, random.start= random.start,
                                                      use.prevlambda = FALSE, prevlambda= NULL,
                                                      coef.trend = NULL, coef.cov = NULL, coef.var = NULL,noise.var = 1e-10,
                                                      delta=delta,
                                                      step=step, fold=4,
                                                      normalized = normalized,lambda=NULL)
        # Get previously visited points to pass to next time step
        prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
        # Get final GP hyperparameter settings to pass to next time step
        coef.trend = optlambda$coef.trend
        coef.cov = optlambda$coef.cov
        coef.var = optlambda$coef.var
        noise.var = optlambda$noise.var
      }
      else #at subsequent time steps, use previous time step's optimization results as priors. 
      {
        optlambda = BO_noisy_optimizer_l1fastbayesave(r.train, prob_old = prob_old, mu_old = mu_old,
                                                      sigma_old =sigma_old, k_old = k_old,
                                                      delta_old = delta_old,
                                                      random.subset = random.subset, random.start= random.start,
                                                      use.prevlambda = TRUE, prevlambda= prevlambda,
                                                      coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var, noise.var = noise.var,
                                                      delta=delta,
                                                      step=step, fold=4,
                                                      normalized = normalized,lambda=NULL)
        # Get previously visited points to pass to next time step
        prevlambda = data.frame(t(optlambda$history.x), optlambda$history.y) 
        # Get final GP hyperparameter settings to pass to next time step
        coef.trend = optlambda$coef.trend
        coef.cov = optlambda$coef.cov
        coef.var = optlambda$coef.var
        noise.var = optlambda$noise.var
      }
      lambda.list[j] = optlambda$best.x
      ybest.vec[[j]] = optlambda$ybestvec
      kmmean.vec[[j]] = optlambda$kmmeanvec
    }

    #cat("optimal lambda: ",lambda.list[j],"\n")
    
    # Compute eta associated with optimal lambda
    # Re fit first with whole training data
    l1fast.optlambda = l1fast.bayesaveraging(r.train, prob_old = prob_old, mu_old = mu_old,
                                             sigma_old =sigma_old, k_old = k_old, 
                                             delta_old = delta_old,
                                             random.subset = random.subset,seedval=j, 
                                             delta=delta,
                                             step=step,
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
    if ( j<=n.test)
    {
      r.test = r[indtest,] #test point
      rp.list[j] = tryCatch({r.test%*%w.optlambda}, error = function(c){return (NA)})
    }
    
    end.time = Sys.time()
    time.list[j] = end.time - start.time
    cat("test point",j,"done at",as.character(end.time),"\n")
  }
  
  return(list(w.list = w.list, rp.list = rp.list, lambda.list=lambda.list, time.list=time.list, ybest.vec=ybest.vec, kmmean.vec=kmmean.vec))
}