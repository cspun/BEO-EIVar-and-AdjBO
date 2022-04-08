############################################################
# This file contains functions to find portfolio weights (l1fast)
# with Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate
############################################################


## l1fast with Bayes-Stein estimate

# Additional function parameters:
# r - matrix of asset returns
# mu_0 - shrinkage target/prior grand mean 
# tau - prior precision of mu_0
l1fast.bayesstein<-function(r, mu_0,tau,delta,step, shrinkage = FALSE, normalized=FALSE,lambda=NULL,
                            nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL,
                            coef=0.8, logspaced=TRUE)
{
  n <- nrow(r);p <- ncol(r)
  
  # Sample mean and sample covariance matrix
  sample_mean = apply(r, 2, mean)
  if (shrinkage) # Apply shrinkage to prevent singularity in high-dimensional case
  {
    sample_cov = cov_shrink(r) 
  }
  else
  {
    sample_cov = cov(r)
  }
  
  if (tau != 0)
  {
    inv_sample_cov = solve(sample_cov) #will throw error if not invertible
  }
  else
  {
    cat("Prior precision must be specified. \n")
  }
  
  # Predictive mean and covariance matrix 
  w = (tau)/(n+tau)
  mu_n = (1-w)*sample_mean + w*rep(1,p)*mu_0
  sigma_n = sample_cov*(1+ 1/(n+tau)) + (tau/ (n*(n+1+tau)))*matrix(1,p,p)*(1/sum(inv_sample_cov%*%rep(1,p)))
  
  sigma_n = sigma_n / step 
  
  if(normalized){
    d=diag(sigma_n)^(-1/2);D=diag(d)
    # normalize delta
    delta=d*delta;maxdelta=max(abs(delta));delta=delta/maxdelta 
    # normalize sigma: becomes correlation matrix
    sigma_n=D%*%sigma_n%*%D
  }
  
  if (is.null(lambda)) { # generate lambda values
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    if (logspaced) {
      lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  
  #cat("inside l1fast bayes length lambda: ",nlambda,"\n")
  
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
    A<-cbind(sigma_n,-sigma_n);A<-rbind(A,-A);b=c(delta,-delta);c=rep(-1,2*p);bbar=rep(1,2*p)
    for (j in nlambda:1) { # start with maximum lambda and decrease 
      temp<-fastlpnocat(c,A,b+bbar*lambda[j]) 
      status[j]=temp$status
      if(status[j]==1){lambda=lambda[-(1:j)];break} # if no feasible sol found, stop and output only lambda from previous iterations
      betalist[[j]]<-temp$sol[1:p]-temp$sol[-(1:p)] #optimal eta for this lambda
      fea[j]=max(abs(sigma_n%*%betalist[[j]]-delta))-lambda[j] # check constraint
      if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}  # revert back from previous normalizations
      if(fea[j]>1e-03){lambda=lambda[-(1:j)];break} # if violated by > 1e-3 stop and output only lambda from previous iterations
      nzlist[[j]]=which(betalist[[j]]!=0) #indices where value nonzero for this eta
      card[j]=length(nzlist[[j]]) #count number of nonzero = sparsity of this eta (cardinality of this portfolio)
      betaeta[j]=t(delta)%*%betalist[[j]]
      diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/n)/step)%*%betalist[[j]]
      gross[j]=sum(abs(betalist[[j]])) #l1 norm of eta
      # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
    }
  } else {
    for (j in nlambda:1) {
      temp<-lpsimplex(sigma_n,delta,lambda[j])
      status[j]=temp$status
      if(status[j]==2){lambda=lambda[-(1:j)];break} 
      betalist[[j]]<-temp$sol #optimal eta?
      fea[j]=max(abs(sigma_n%*%betalist[[j]]-delta))-lambda[j]
      if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}
      if(fea[j]>1e-03){lambda=lambda[-(1:j)];break}
      nzlist[[j]]=which(betalist[[j]]!=0)
      card[j]=length(nzlist[[j]])
      betaeta[j]=t(delta)%*%betalist[[j]]
      diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/n)/step)%*%betalist[[j]]
      gross[j]=sum(abs(betalist[[j]])) 
      # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
    }
  }
  return(list(gross=gross,diff=diff,betaeta=betaeta,betalist=betalist, lambda=lambda, nzlist=nzlist, cardinality=card, status=status, fea=fea))
}


# Estimate prior values for Bayes Stein based on data directly (empirical Bayes)
# Return: prior parameters mu_0 and tau.
# Prior parameter mu_0 set to be mean past return of GMV portfolio.
# Prior parameter tau estimated based on formula derived in original paper.
getbayessteinprior <- function(r, shrinkage = FALSE)
{
  n <- nrow(r);p <- ncol(r)
  # Sample mean and sample covariance matrix
  sample_mean = apply(r, 2, mean)
  if (shrinkage)
  {
    sample_cov = cov_shrink(r) # use shrinkage estimate to address singularity of original sample cov matrix
  }
  else
  {
    sample_cov = cov(r) 
  }
  
  # Estimate of mu_0: Return of GMV portfolio
  w_MV = solve(sample_cov)%*%rep(1,p) #will throw error if not invertible
  w_MV = w_MV/sum(w_MV)
  mu_0 = sum(w_MV*sample_mean)
  
  # Estimate of tau (prior precision)
  d = t(sample_mean - rep(1,p)*mu_0)%*%solve(sample_cov)%*%(sample_mean - rep(1,p)*mu_0)
  tau = (p+2)/d
  
  return(list(mu_0 = mu_0, tau=tau[1,1]))
}



# l1fast with Bayesian "data-and-model" approach: Bayes estimate based on belief in asset pricing model.
# An investor has belief in certain asset pricing model and is captured by a prior 
# on the extent of model mispricing.

# For formula of predictive mean and covariance, refer to the Appendix of original paper. 

# Additional function parameters:
# r - matrix of asset returns
# F - matrix of factor returns 
# sigma_alpha - investor's degree of skepticism towards the model, i.e degree of mispricing. 0 means believe completely in the model. Higher value means more uncertain about the pricing model's accuracy.
l1fast.bayesassetpricing<-function(r, f, sigma_alpha, nu, delta,step,normalized=FALSE,lambda=NULL,
                                   nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL,
                                   coef=0.8, logspaced=TRUE)
{
  
  #cat("l1fast bayes start: ", as.character(Sys.time()),"\n")
  f = as.matrix(f)
  n_obs = nrow(r)
  n_assets = ncol(r)
  n_factor = ncol(f)
  p = n_assets
  n = n_obs
  sample_cov = cov(r)
  
  # Assumes r = ZA + U, where vec(U) ~ N(0, sigma x I_n_obs)
  # Calculate sample estimates
  Z = cbind(rep(1,n_obs), f)
  A = solve(t(Z)%*%Z) %*% t(Z) %*% r #A contains estimates for both alpha and B
  A_vec = c(A) # vectorization of matrix A
  res = r - Z%*%A
  sigma = (t(res)%*%res) / n_obs
  s2 = mean(diag(sigma))
  E_f = (t(f)%*%rep(1,n_obs)) / n_obs
  diff = f - rep(1,n_obs)%*%t(E_f)
  V_f = (t(diff)%*%diff) / n_obs
  
  # Posterior estimates of unknown parameters
  E_f.post = E_f
  V_f.post = (n_obs / (n_obs - n_factor - 2))*V_f
  
  D = matrix(0, n_factor+1, n_factor+1)
  D[1,1] = s2/(sigma_alpha^2)
  temp = D + t(Z)%*%Z
  temp.inv = solve(temp)
  A_vec.post = ( diag(n_assets) %x% (temp.inv%*%t(Z)%*%Z) ) %*% A_vec
  
  H = s2 * (nu - n_assets - 1) * diag(n_assets)
  Q = t(Z) %*% (diag(n_obs) - Z%*%temp.inv%*%t(Z)) %*% Z
  sigma.post = (H + n_obs*sigma + t(A)%*%Q%*%A)/(n_obs + nu - n_assets - n_factor - 1)
  
  covA_vec = sigma.post %x% temp.inv
  V_f_star = V_f.post + V_f/(n_obs - n_factor -2)
  
  # Construct covariance matrix of predictive distribution
  cov_n = matrix(0, n_assets, n_assets)
  for (i in 1:n_assets)
  {
    for (j in i:n_assets)
    {
      ind_bi = (i*(n_factor+1)-n_factor+1):(i*(n_factor+1))
      ind_bj = (j*(n_factor+1)-n_factor+1):(j*(n_factor+1))
      bi.post = A_vec.post[ind_bi]
      bj.post = A_vec.post[ind_bj]
      cov_bi_bj = covA_vec[ind_bi, ind_bj]
      first_term = t(bi.post)%*%V_f_star%*%bj.post + sum(diag(V_f_star%*%cov_bi_bj)) + sigma.post[i,j]
      
      ind_ai = (i*(n_factor+1) - n_factor):(i*(n_factor+1))
      ind_aj = (j*(n_factor+1) - n_factor):(j*(n_factor+1))
      cov_ai_aj = covA_vec[ind_ai, ind_aj]
      second_term = t(c(1, E_f.post))%*%cov_ai_aj%*%(c(1, E_f.post))
      
      cov_n[i,j] = first_term + second_term
    }
  }
  diagval = diag(cov_n)
  cov_n = cov_n + t(cov_n) - diag(diagval)
  
  cov_n = cov_n / step
  
  if(normalized){
    d=diag(cov_n)^(-1/2);D=diag(d)
    # normalize delta
    # delta = d*delta gives excess return per std dev (sharpe ratio) then divide by its max value to convert actual values capped at 1
    delta=d*delta;maxdelta=max(abs(delta));delta=delta/maxdelta 
    # normalize sigma: becomes correlation matrix
    cov_n=D%*%cov_n%*%D
  }
  
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  
  #cat("inside l1fast bayes length lambda: ",nlambda,"\n")
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
    A<-cbind(cov_n,-cov_n);A<-rbind(A,-A);b=c(delta,-delta);c=rep(-1,2*p);bbar=rep(1,2*p)
    for (j in nlambda:1) { # start with maximum lambda and decrease 
      temp<-fastlpnocat(c,A,b+bbar*lambda[j]) 
      status[j]=temp$status
      if(status[j]==1){lambda=lambda[-(1:j)];break} # if no feasible sol found, stop and output only lambda from previous iterations
      betalist[[j]]<-temp$sol[1:p]-temp$sol[-(1:p)] #optimal eta for this lambda
      fea[j]=max(abs(cov_n%*%betalist[[j]]-delta))-lambda[j] # check constraint
      if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}  # revert back from previous normalizations
      if(fea[j]>1e-03){lambda=lambda[-(1:j)];break} # if violated by > 1e-3 stop and output only lambda from previous iterations
      nzlist[[j]]=which(betalist[[j]]!=0) #indices where value nonzero for this eta
      card[j]=length(nzlist[[j]]) #count number of nonzero = sparsity of this eta (cardinality of this portfolio)
      betaeta[j]=t(delta)%*%betalist[[j]]
      diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/n)/step)%*%betalist[[j]]
      gross[j]=sum(abs(betalist[[j]])) #l1 norm of eta
      # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
    }
  } else {
    for (j in nlambda:1) {
      temp<-lpsimplex(cov_n,delta,lambda[j])
      status[j]=temp$status
      if(status[j]==2){lambda=lambda[-(1:j)];break} 
      betalist[[j]]<-temp$sol #optimal eta?
      fea[j]=max(abs(cov_n%*%betalist[[j]]-delta))-lambda[j]
      if(normalized){betalist[[j]]=betalist[[j]]*d*maxdelta}
      if(fea[j]>1e-03){lambda=lambda[-(1:j)];break}
      nzlist[[j]]=which(betalist[[j]]!=0)
      card[j]=length(nzlist[[j]])
      betaeta[j]=t(delta)%*%betalist[[j]]
      diff[j]=betaeta[j]-t(betalist[[j]])%*%(sample_cov*(1-1/n)/step)%*%betalist[[j]]
      gross[j]=sum(abs(betalist[[j]])) 
      # if(card[j]>p*0.15){lambda=lambda[-(1:j)];break}
    }
  }
  
  #cat("l1fast bayes end: ", as.character(Sys.time()),"\n")
  return(list(gross=gross,diff=diff,betaeta=betaeta,betalist=betalist, lambda=lambda, nzlist=nzlist, cardinality=card, status=status, fea=fea))
}




# l1fast with Bayesian Averaging estimate

# At time t, there exists t possible models of asset returns. This consists of t-1 existing models and one new model that is born at time t.
# Notation: When there are s models, P_s(m|F_t-1) is probability that model m is correct and no older model q<m is correct, conditioned on information available at time t-1 and that at least one of the s models is correct. 

# Existing models at t-1: Each model m assumes rt~N(mu_m, sigma_m)
# Prior: mu_m|sigma_m ~ N(mu_m_tminus1, sigma_m/ k_m_tminus1)
# Sigma_m ~ IW(sigma_m_tminus1 * delta_m_tminus1, delta_m_tminus1+n+1)
# m_tminus1 indicates beliefs in model m at time t-1

# Function parameters:
# r : matrix of asset returns from time 1 up to time t
# prob_old: vector of length t-1 containing probabilities of t-1 existing/older models, i.e. P_t-1(m|F_t-1)
# mu_old: vector of length t-1 containing beliefs on mean from each of the t-1 existing/older models
# sigma_old: vector of length t-1 containing beliefs on sigma from each of the t-1 existing/older models
# k_old: vector of length t-1 containing precision on mean beliefs from each of the t-1 existing/older models
# delta_old: vector of length t-1 containing precision on cov matrix beliefs from each of the t-1 existing/older models

l1fast.bayesaveraging <- function(r, prob_old, mu_old, sigma_old, k_old, delta_old, 
                                  random.subset = FALSE, seedval= seed, delta, step, normalized=FALSE,lambda=NULL,
                                  nlambda = ifelse(is.null(lambda),100,length(lambda)),lambda.max=NULL,
                                  coef=0.8, logspaced=TRUE)
{
  t = nrow(r)
  p = ncol(r)
  sample_cov = cov(r)
  
  if (random.subset)
  {
    # Take only a subset of the total assets
    p1 = 0.03*p
    # Sample p1 random indices 
    set.seed(seedval)
    ind_random = sample(x=ncol(r), size=p1, replace=FALSE)
  }
  else
  {
    p1 = p
    ind_random = 1:ncol(r)
  }
  
  # Step 1: A new model is born at time t
  # Here we assume returns at time t has not been observed yet
  # New model assumes rt~N(mu_t, sigma_t)
  # Prior for mu_t and sigma_t is Normal-Inverse Wishart
  # Prior beliefs of covariance betw any 2 assets is zero
  # Prior beliefs of mean excess returns and variances are identical for each asset
  r_prev = r[1:(t-1),]
  # Compute prior for newly born model t and append to existing collection of models
  if (t>2)
  {
    sample_mean = apply(r_prev,2,mean)
  }
  else
  {
    sample_mean = r_prev
  }
  mu_t_tminus1 = rep(1,p)* (sum(sample_mean)/p) #prior mean
  mu_old = append(mu_old, list(mu_t_tminus1))
  k_t_tminus1 = 1 #precision of beliefs on mean
  k_old = append(k_old, k_t_tminus1)
  delta_t_tminus1 = 1 #precision of beliefs on cov
  delta_old = append(delta_old, delta_t_tminus1)
  sqMat = (r_prev - matrix(sample_mean, nrow=t-1, ncol=p, byrow=TRUE))^2
  if (t>2)
  {
    lambda_vec_tminus1 = apply(sqMat, 2, sum)/(t-2)
    sigma_t_tminus1 = diag(p) * (sum(lambda_vec_tminus1)/p) #prior cov matrix
  }else
  {
    sigma_t_tminus1 = diag(p) * 0.0001 #prior cov matrix
  }
  sigma_old = append(sigma_old, list(sigma_t_tminus1))
  
  
  # Step 2: Probabilities for all models are adjusted for inclusion of newly born model t
  # Here we assume returns at time t has not been observed yet
  # Update probabilities from P_t-1(m|F_t-1) to P_t(m|F_t-1)
  # prob_afternewmodel now has length t 
  prob_afternewmodel = updateProb_SharingPrior(alpha=1, prob_old) #alpha = 1 is perfect sharing
  
  
  # Step 3: Excess returns at time t are observed
  r_t = r[t,]
  vec_likelihood = rep(0,t)
  for (m in 1:t)
  {
    # Step 3.1: Compute new parameters of each of the t models
    k_m_t = k_old[m] + 1
    delta_m_t = delta_old[m] + 1
    mu_m_t = (k_old[m]*mu_old[[m]] + r_t)/k_m_t
    sigma_m_t = ( (delta_old[m] * k_m_t * sigma_old[[m]]) +
                    (k_old[m]*(r_t-mu_old[[m]])%*%t(r_t-mu_old[[m]])) )/(delta_m_t*k_m_t)
    # Compute likelihood of observing r_t given model m and info at time t-1
    mvgammaquotient = 1
    for (i in 1:p1)
    {
      mvgammaquotient = mvgammaquotient * sqrt((delta_m_t+p1+1-i)/2)
    }
    #Get subset of matrix according to the p1 random indices
    sigma_m_t_reduced = sigma_m_t[,ind_random]
    sigma_m_t_reduced = sigma_m_t_reduced[ind_random,]
    sigma_m_tminus1_reduced = sigma_old[[m]][,ind_random]
    sigma_m_tminus1_reduced = sigma_m_tminus1_reduced[ind_random,]
    # Likelihood function
    a = sqrt( det(solve(delta_old[m]*sigma_m_tminus1_reduced)) )
    b = (det( ((delta_old[m]+1)/delta_old[m])*sigma_m_t_reduced*solve(sigma_m_tminus1_reduced) ))^((delta_old[m] + p1 +2)/2)
    detquotient = a/b
    vec_likelihood[m] = (k_old[m]/(pi*k_m_t))^(p1/2) * detquotient  * mvgammaquotient
    # Update parameters of each of the t models
    k_old[m] = k_m_t
    delta_old[m] = delta_m_t
    mu_old[[m]] = mu_m_t
    sigma_old[[m]] = sigma_m_t
  }
  
  # Step 3.2: Update probabilities using new return info at time t, i.e. from P_t(m|F_t-1) to P_t(m|F_t)
  prob_afterreturnobsv = rep(0, t)
  sum_prob = sum(vec_likelihood*prob_afternewmodel)
  for (m in 1:t)
  {
    prob_afterreturnobsv[m] = (vec_likelihood[m]*prob_afternewmodel[m])/sum_prob
  }
  # Update prob_old 
  prob_old = prob_afterreturnobsv

  
  # Step 4: Expectations of returns at time t+1 are formed
  # Compute estimates for E(r_t+1|F_t) and V(r_t+1|F_t)
  # Estimate for E(r_t+1|F_t)
  mu_hat_t = rep(0,p)
  for (m in 1:t)
  {
    mu_hat_t = mu_hat_t + mu_old[[m]]*prob_afterreturnobsv[m]
  }
  # Estimate for V(r_t+1|F_t)
  sigma_hat_t = matrix(0,p,p) 
  for (m in 1:t)
  {
    sigma_hat_t = sigma_hat_t + ( ((1+ k_old[m])/k_old[m])*sigma_old[[m]] 
                                  + mu_old[[m]]%*%t(mu_old[[m]]) )*prob_afterreturnobsv[m] 
  }
  sigma_hat_t = sigma_hat_t - mu_hat_t%*%t(mu_hat_t)

  # Step 5: Compute LPO weights based on estimates from step 4
  sigma_hat_t = sigma_hat_t / step 
  if(normalized){
    d=diag(sigma_hat_t)^(-1/2);D=diag(d)
    # normalize delta
    delta=d*delta;maxdelta=max(abs(delta));delta=delta/maxdelta 
    # normalize sigma: becomes correlation matrix
    sigma_hat_t=D%*%sigma_hat_t%*%D
  }
  
  if (is.null(lambda)) {
    coef=ifelse(coef>1,1,coef) 
    if (is.null(lambda.max)) {lambda.max=coef*max(abs(delta))}
    # generate lambda values
    if (logspaced) {lambda<-10^(seq(log10(lambda.max/50), log10(lambda.max), length.out=nlambda))} else {lambda<-seq(0, lambda.max, length.out = nlambda)}
  }
  
  cat("inside l1fast bayes length lambda: ",nlambda,"\n")
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
      temp<-fastlpnocat(c,A,b+bbar*lambda[j]) 
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
  
  cat("l1fast bayes end: ", as.character(Sys.time()),"\n")
  return(list(gross=gross,diff=diff,betaeta=betaeta,betalist=betalist, lambda=lambda, nzlist=nzlist,
              cardinality=card, status=status, fea=fea, 
              prob_old=prob_old, mu_old=mu_old, sigma_old=sigma_old, k_old=k_old, delta_old=delta_old, vec_likelihood=vec_likelihood))
  
}

# Function to update probabilities from P_t-1(m|F_t-1) to P_t(m|F_t-1) using sharing prior
# Function parameters:
# alpha: fraction of its previous probability that each model shares equally with itself and newer models
# Each model retains (1-alpha) of its previous probabilities
# prob_prev: vector containing previous probabilities before inclusion of newly born model t, i.e. P_t-1(m|F_t-1)
updateProb_SharingPrior <- function(alpha, prob_prev)
{
  # Length of prob_prev is t-1 since previously there are t-1 models
  t = length(prob_prev) + 1
  prob_new = rep(0,t)
  
  for (i in 1:(t-1))
  {
    terms1 = (1/( (t-1+1):(t-i+1) )) * prob_prev[1:i] 
    prob_new[i] = (1-alpha)*prob_prev[i] + alpha*(sum(terms1))
  }
  terms2 = (1/( (t-1+1):(t-(t-1)+1) )) * prob_prev
  prob_new[t] = alpha*(sum(terms2))
  
  return (prob_new)
}
