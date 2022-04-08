############################################################
# This file contains functions to find optimal lambda using CV 
# for LPO with sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate
############################################################


#Loss function in choosing lambda: minimize negative sharpe ratio
# negsharperatio<-function(r,eta,gamma=1){
#   # if eta is almost zero or null, set loss to 100000
#   if(sum(eta^2)<1e-06 | is.null(eta)){nsr=100000} else {w=eta/sum(eta);dailyr=r%*%w;nsr=-mean(dailyr)/sd(dailyr)}
#   nsr
# }

negCEQ<-function(r,eta,gamma=1){
  if(sum(eta^2)<1e-06 | is.null(eta)){nceq=100000} else {w=eta/sum(eta);dailyr=r%*%w;nceq=-(mean(dailyr)-gamma/2*var(dailyr))}
  nceq
}

# Get indexes of partition to train and test for k-fold CV
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

# Get indexes of partition for rolling windows
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


# Find optimal lambda using CV for l1fast with sample estimates
cv.negSR.l1fast<-function(r,delta,step,fold=6,normalized = FALSE,shrink=TRUE,
                          perturb=FALSE,lambda=NULL,logspaced=TRUE,coef=0.8,
                          nlambda = ifelse(is.null(lambda),50,length(lambda)),lambda.max=NULL)
{
  p<-ncol(r);n<-nrow(r)
  
  if (is.null(lambda)) { #generate lambdas if not supplied
    coef=ifelse(coef>1,1,coef)
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(lambda.max/50, lambda.max, length.out = nlambda)}
  }
  
  part.list<-cv.part(n,fold) #index of training and test set for each of the k fold
  loss.re<-matrix(Inf,nrow=fold,ncol=nlambda)
  minval = Inf
  vec.minval = rep(0,nlambda)
  
  for (jl in nlambda:1)
  {
    for (j in 1:fold)
    {
      indtrain<-part.list$trainMat[,j] #get index of training data for the current fold
      r.train<-r[indtrain,] # training data
      l1fastcv<-l1fast(r.train,delta=delta,step=step,normalized=normalized,shrink=shrink,
                       lambda= lambda[jl])
      indtest<-part.list$testMat[,j] #get index of test data for the current fold
      r.test <- r[indtest,]
      loss.re[j,jl]<- negCEQ(r.test, l1fastcv$betalist[[1]])
    }
    val = max(loss.re[,jl])
    if (val == 100000)
    {
      vec.minval[(1:jl)] = minval
      break
    }
    if (val < minval)
    {
      minval = val
    }
    vec.minval[jl] = minval
  }
  
  # Use robust SR as CV function
  loss.mean<-apply(loss.re, 2, max)
  ind<-which.min(loss.mean)
  lambdaopt<-lambda[ind]
  outlist<-list(lambdaopt=lambdaopt, lambdaindex=ind, vec.minval=rev(vec.minval))
  class(outlist)<-c("cv.clime")
  
  return(outlist)
}



# Find optimal lambda using CV for l1fast with Bayes Stein estimates
cv.negSR.l1fastbayesstein<-function (r,mu_0, tau,delta,step,fold=6,shrinkage = FALSE, normalized = FALSE,
                                     lambda=NULL,logspaced=TRUE,coef=0.8,
                                     nlambda = ifelse(is.null(lambda),50,length(lambda)),lambda.max=NULL)
{
  p<-ncol(r);n<-nrow(r)
  
  if (is.null(lambda)) { #generate lambdas if not supplied
    coef=ifelse(coef>1,1,coef)
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(lambda.max/50, lambda.max, length.out = nlambda)}
  }
  
  part.list<-cv.part(n,fold) #index of training and test set for each of the k fold
  loss.re<-matrix(Inf,nrow=fold,ncol=nlambda)
  minval = Inf
  vec.minval = rep(0,nlambda)
  
  for (jl in nlambda:1)
  {
    for (j in 1:fold)
    {
      indtrain<-part.list$trainMat[,j] #get index of training data for the current fold
      r.train<-r[indtrain,] # training data
      l1fastcv<-l1fast.bayesstein(r.train,mu_0=mu_0,tau=tau,delta=delta,
                                  step=step,shrinkage = shrinkage,
                                  normalized=normalized,lambda=lambda[jl])
      indtest<-part.list$testMat[,j] #get index of test data for the current fold
      r.test <- r[indtest,]
      loss.re[j,jl]<- negCEQ(r.test, l1fastcv$betalist[[1]])
    }
    val = max(loss.re[,jl])
    if (val == 100000)
    {
      vec.minval[(1:jl)] = minval
      break
    }
    if (val < minval)
    {
      minval = val
    }
    vec.minval[jl] = minval
  }
  
  # Use robust SR as CV function
  loss.mean<-apply(loss.re, 2, mean)
  ind<-which.min(loss.mean)
  lambdaopt<-lambda[ind]
  outlist<-list(lambdaopt=lambdaopt, lambdaindex=ind, vec.minval=rev(vec.minval))
  class(outlist)<-c("cv.clime")
  
  return(outlist)
}


# Find optimal lambda using CV for l1fast with Bayes Asset Pricing estimates
cv.negSR.l1fastbayesassetpricing<-function (r,f,sigma_alpha,nu,delta,step,fold=6,normalized = FALSE,
                                            lambda=NULL,logspaced=TRUE,coef=0.8,
                                            nlambda = ifelse(is.null(lambda),50,length(lambda)),lambda.max=NULL)
{
  f = as.matrix(f)
  p<-ncol(r);n<-nrow(r)
  
  if (is.null(lambda)) { #generate lambdas if not supplied
    coef=ifelse(coef>1,1,coef)
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(lambda.max/50, lambda.max, length.out = nlambda)}
  }
  
  part.list<-cv.part(n,fold) #index of training and test set for each of the k fold
  loss.re<-matrix(Inf,nrow=fold,ncol=nlambda)
  minval = Inf
  vec.minval = rep(0,nlambda)
  
  for (jl in nlambda:1)
  {
    for (j in 1:fold)
    {
      indtrain<-part.list$trainMat[,j] #get index of training data for the current fold
      r.train<-r[indtrain,] # training data
      f.train = f[indtrain,] #training data for factor returns
      l1fastcv<-l1fast.bayesassetpricing(r.train,f.train, sigma_alpha = sigma_alpha,
                                         nu=nu, delta=delta, step=step,
                                         normalized=normalized,lambda=lambda[jl])
      indtest<-part.list$testMat[,j] #get index of test data for the current fold
      r.test <- r[indtest,]
      loss.re[j,jl]<- negCEQ(r.test, l1fastcv$betalist[[1]])
    }
    val = max(loss.re[,jl])
    if (val == 100000)
    {
      vec.minval[(1:jl)] = minval
      break
    }
    if (val < minval)
    {
      minval = val
    }
    vec.minval[jl] = minval
  }
  
  # Use robust SR as CV function
  loss.mean<-apply(loss.re, 2, mean)
  ind<-which.min(loss.mean)
  lambdaopt<-lambda[ind]
  outlist<-list(lambdaopt=lambdaopt, lambdaindex=ind, vec.minval=rev(vec.minval))
  class(outlist)<-c("cv.clime")
  
  return(outlist)
}


# Find optimal lambda using CV for l1fast with Bayes Averaging estimates
# Note: For Bayes Averaging, only 1 pseudo-test set is used (taken to be last f=4 points). Since must preserve time dependency.
cv.negSR.l1fastbayesaveraging <- function(r, prob_old, mu_old, sigma_old, k_old, delta_old,
                                          random.subset = FALSE, delta, step,fold,
                                          normalized =FALSE,
                                          lambda=NULL,logspaced=TRUE,coef=0.8,
                                          nlambda=ifelse(is.null(lambda),50,length(lambda)),
                                          lambda.max=NULL)
{
  p<-ncol(r);t<-nrow(r)
  sample_cov = cov(r)
  
  return.vec = rep(0, nlambda)
  vec.minval=rep(0,nlambda)
  
  if (is.null(lambda)) {  # generate lambdas if not supplied
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  lambda.temp = lambda
  
  # 4 fold CV
  if (t>6)
  {
    # Take last 4 points as pseudo test set
    fold = fold
    returns.mat <- matrix(0,nrow=fold,ncol=nlambda)
    
    for (f in fold:1) 
    {
      # Use data from time 1 to t-f as training
      indtrain = 1:(t-f)
      n.train = length(indtrain)
      r.train = r[indtrain,]
      
      # Use data at time t-f + 1 as test point to get best lambda 
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
      betalist<-vector("list", nlambda) 
      # betalist is initialized with NULL and will be updated to store optimal eta associated with certain lambda.
      # if a lambda yields no feasible sol, the entry in betalist will still be the default NULL
      nzlist<-vector("list", nlambda)
      card<-rep(p,nlambda)
      fea<-rep(NA,nlambda)
      status<-rep(0,nlambda)
      betaeta<-rep(0,nlambda)
      diff<-rep(0,nlambda)
      gross<-rep(0,nlambda)
      if(t<p){
        A<-cbind(sigma_hat_t,-sigma_hat_t);A<-rbind(A,-A);b=c(delta,-delta);c=rep(-1,2*p);bbar=rep(1,2*p)
        for (j in nlambda:1) { # start with maximum lambda and decrease 
          temp<-tryCatch({ fastlpnocat(c,A,b+bbar*lambda[j])}, error = function(c){return (list(status=1))})
            
          status[j]=temp$status
          if(status[j]==1){lambda=lambda[-(1:j)];break} # if no feasible sol found, stop and output only lambda from previous iterations
          betalist[[j]]<-temp$sol[1:p]-temp$sol[-(1:p)] #optimal eta for this lambda
          fea[j]=max(abs(sigma_hat_t%*%betalist[[j]]-delta))-lambda[j] # check constraint
          if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}  # revert back from previous normalizations
          if(fea[j]>1e-03){lambda=lambda[-(1:j)];break} # if violated by > 1e-3 stop and output only lambda from previous iterations
          nzlist[[j]]=which(betalist[[j]]!=0) #indices where value nonzero for this eta
          card[j]=length(nzlist[[j]]) #count number of nonzero = sparsity of this eta (cardinality of this portfolio)
          betaeta[j]=t(delta)%*%betalist[[j]]
          diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/t)/step)%*%betalist[[j]]
          gross[j]=sum(abs(betalist[[j]])) #l1 norm of eta
          # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
        }
      } else {
        for (j in nlambda:1) {
          temp<-lpsimplex(sigma_hat_t,delta,lambda[j])
          status[j]=temp$status
          if(status[j]==2){lambda=lambda[-(1:j)];break} 
          betalist[[j]]<-temp$sol #optimal eta?
          fea[j]=max(abs(sigma_hat_t%*%betalist[[j]]-delta))-lambda[j]
          if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}
          if(fea[j]>1e-03){lambda=lambda[-(1:j)];break}
          nzlist[[j]]=which(betalist[[j]]!=0)
          card[j]=length(nzlist[[j]])
          betaeta[j]=t(delta)%*%betalist[[j]]
          diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/t)/step)%*%betalist[[j]]
          gross[j]=sum(abs(betalist[[j]])) 
          # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
        }
      }
      
      cat("nlambda before fold",f,"is",nlambda,"\n")
      cat("length lambda after fold",f,"is",length(lambda),"\n")
      # 
      # if (length(lambda) ==0)
      # {
      #   cat("no feasible lambda \n")
      #   break
      # }
      # 
      # if ((nlambda-length(lambda)) > 0)
      # {
      #   returns.mat = as.matrix(returns.mat[,-(1:(nlambda-length(lambda)))])
      #   betalist = betalist[-(1:(nlambda-length(lambda)))]
      # }
      # #cat("lossre ncol", ncol(returns.mat),"\n")
      # #cat("betalist length",length(betalist),"\n")
      # nlambda = length(lambda)
      # #cat("nlambda after fold",f,"is",nlambda,"\n")
      
      for (jl in 1:nlambda) 
      {
        returns.mat[f,jl]<- tryCatch({ r.test%*%(betalist[[jl]] / sum(betalist[[jl]]))},
                                     error = function(c){return (NA)})
          
      }
    } #end of for loop fold
    
    #Get optimal lambda
    returns.mean<-apply(returns.mat, 2, mean, na.rm=TRUE)
    # returns.sd <- apply(returns.mat, 2, sd)
    # returns.sr <- returns.mean/returns.sd
    
    returns.var = apply(returns.mat, 2, var, na.rm=TRUE)
    returns.sr = returns.mean - 0.5*returns.var
    cat(returns.sr,"\n")
    
    if (sum(is.na(returns.sr)) == nlambda)
    {
      cat("no feasible lambda \n")
      lambdaopt = lambda.max
    }
    else
    {
      cat("found feas lambda \n")
      ind<-which.max(returns.sr)
      cat("ind", ind,"\n")
      lambdaopt<-lambda.temp[ind]
    }
  }
  else 
  {
    lambdaopt<- 0.8
  }
  
  outlist<-list(lambdaopt=lambdaopt, vec.minval=vec.minval)
  return(outlist)
}
