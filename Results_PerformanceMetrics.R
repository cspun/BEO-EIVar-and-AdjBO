getSharpeRatio <- function(returnsVec)
{
  return (mean(returnsVec)/sd(returnsVec))
}


getReturnsVariance <- function(returnsVec)
{
  return (var(returnsVec))
}


getCEQ <- function(returnsVec, riskAversion)
{
  return (mean(returnsVec) - 0.5*riskAversion*var(returnsVec))
}


getWeightMatrix <- function(wlist)
{
  w.total = length(wlist)
  p = length(wlist[[1]])
  wmat = matrix(0,w.total,p) 
  for (i in 1:w.total)
  {
    wmat[i,] = wlist[[i]]
  }
  return (wmat)
}


getTurnover<-function(wlist,rtest){
  w = getWeightMatrix(wlist)
  ntest=nrow(rtest)
  wtemp<-w[1:ntest,]*(1+rtest)
  wtemp<-wtemp/apply(wtemp,1,sum)
  to=mean(apply(abs(w[2:(ntest+1),]-wtemp),1,sum))
  return (to)
}


# l1 norm differences with true optimal lambda
vsoptlambda.l1norm <- function(lambdalist, optlambdalist)
{
  return (sum(abs(lambdalist-optlambdalist)))
}


getTotalTime <- function(timelist) #in minutes
{
  return (sum(timelist))
}


getAverageTime <- function(timelist) #in minutes
{
  return (mean(timelist))
}