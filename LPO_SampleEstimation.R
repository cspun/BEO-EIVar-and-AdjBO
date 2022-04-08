############################################################
# This file contains functions to find portfolio weights (l1fast)
# with Sample estimate
############################################################

library(fastclime)
library(lpSolve)
library(tawny)

lpsimplex<-function (Sigma, e, lambda) 
{
  p <- nrow(Sigma)
  if (p != ncol(Sigma)) 
    stop("Sigma should be a square matrix!")
  f.obj <- rep(1, 2 * p)
  con1 <- cbind(-Sigma, +Sigma)
  b1 <- lambda - e
  b2 <- lambda + e
  f.con <- rbind(-diag(2 * p), con1, -con1)
  f.dir <- rep("<=", 4 * p)
  f.rhs <- c(rep(0, 2 * p), b1, b2)
  lp.out <- lp("min", f.obj, f.con, f.dir, f.rhs)
  beta <- lp.out$solution[1:p] - lp.out$solution[(p + 1):(2 * 
                                                            p)]
  if (lp.out$status == 2) warning("No feasible solution!  Try a larger lambda maybe!")
  return(list(sol=beta,status=lp.out$status))
}

fastlpnocat<-function (obj, mat, rhs, lambda = 0) 
{
  m <- length(rhs)
  n <- length(obj)
  m0 <- dim(mat)[1]
  n0 <- dim(mat)[2]
  opt <- rep(0, n)
  status <- 0
  error <- 0
  if (m != m0 || n != n0) {
    cat("Dimensions do not match! \n")
    error <- 1
  }
  if (error == 0) {
    # fastlp(obj, mat, rhs, lambda=0) is used to solve a general linear programming 
    # in standard inequality form: "maximize obj*x, subject to: mat*x<=rhs, x>=0"
    
    str = .C("fastlp", as.double(obj), as.double(t(mat)), 
             as.double(rhs), as.integer(m0), as.integer(n0), as.double(opt), 
             as.integer(status), as.double(lambda), PACKAGE = "fastclime")
    opt <- unlist(str[6])
    status <- unlist(str[7])
    if (status == 0) {
      # cat("optimal solution found! \n")
      return(list(sol=opt,status=status))
    }
    else if (status == 1) {
      # cat("The problem is infeasible! \n")
      return(list(status=status))
    }
    else if (status == 2) {
      # cat("The problem is unbounded! \n")
      return(list(sol=rep(0,2*m0),status=status))
    }
  }
}


l1fast<-function(x,delta,step,normalized=FALSE,perturb=FALSE,lambda=NULL,shrink=FALSE,
                 nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL,
                 coef=0.8, logspaced=TRUE)
{
  n <- nrow(x);p <- ncol(x)
  Sigma<-cov(x)*(1-1/n)/step 
  
  if(shrink){Sigmap<-cov.shrink(x)/step} else{ 
    if (is.logical(perturb)){ 
      if (perturb) {
        eigvals<-eigen(Sigma, only.values = T)$values
        perturb<-max(max(eigvals)-p*min(eigvals),0)/(p-1)
      } else {perturb<-0}
    }
    Sigmap<-Sigma+diag(p)*perturb
  }
  
  if(normalized){
    d=diag(Sigmap)^(-1/2);D=diag(d)
    # normalize delta
    delta=d*delta;maxdelta=max(abs(delta));delta=delta/maxdelta 
    # normalize sigma: becomes correlation matrix
    Sigmap=D%*%Sigmap%*%D
  }
  
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  
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
  if(n<p){
    A<-cbind(Sigmap,-Sigmap);A<-rbind(A,-A);b=c(delta,-delta);c=rep(-1,2*p);bbar=rep(1,2*p)
    for (j in nlambda:1) { # start with maximum lambda and decrease 
      temp<-fastlpnocat(c,A,b+bbar*lambda[j]) 
      status[j]=temp$status
      if(status[j]==1){lambda=lambda[-(1:j)];break} # if no feasible sol found, stop and output only lambda from previous iterations
      betalist[[j]]<-temp$sol[1:p]-temp$sol[-(1:p)] #optimal eta for this lambda
      fea[j]=max(abs(Sigmap%*%betalist[[j]]-delta))-lambda[j] # check constraint
      if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}  # revert back from previous normalization
      if(fea[j]>1e-03){lambda=lambda[-(1:j)];break} # if violated by > 1e-3 stop and output only lambda from previous iterations
      nzlist[[j]]=which(betalist[[j]]!=0) #indices where value nonzero for this eta
      card[j]=length(nzlist[[j]]) #count number of nonzero = sparsity of this eta (cardinality of this portfolio)
      betaeta[j]=t(delta)%*%betalist[[j]]
      diff[j]=betaeta[j]-t(betalist[[j]])%*%Sigma%*%betalist[[j]]
      gross[j]=sum(abs(betalist[[j]])) #l1 norm of eta
      # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
    }
  } else {
    for (j in nlambda:1) {
      temp<-lpsimplex(Sigmap,delta,lambda[j])
      status[j]=temp$status
      if(status[j]==2){lambda=lambda[-(1:j)];break} 
      betalist[[j]]<-temp$sol #optimal eta?
      fea[j]=max(abs(Sigmap%*%betalist[[j]]-delta))-lambda[j]
      if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}
      if(fea[j]>1e-03){lambda=lambda[-(1:j)];break}
      nzlist[[j]]=which(betalist[[j]]!=0)
      card[j]=length(nzlist[[j]])
      betaeta[j]=t(delta)%*%betalist[[j]]
      diff[j]=betaeta[j]-t(betalist[[j]])%*%Sigma%*%betalist[[j]]
      gross[j]=sum(abs(betalist[[j]])) 
      # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
    }
  }
  return(list(gross=gross,diff=diff,betaeta=betaeta,betalist=betalist, lambda=lambda, nzlist=nzlist, cardinality=card, status=status, fea=fea))
}