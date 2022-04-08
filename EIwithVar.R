############################################################
# This file contains functions to implement Expected Improvement with Variance acquisition function.
############################################################

library('stats')
library('utils')
library('data.table')
library('magrittr')
library('foreach')


# Generate initial values
Matrix_runif <- function(n, lower, upper) {
  foreach(i = seq_along(lower), .combine = "cbind") %do% {
    runif(n, min = lower[i], max = upper[i]) %>%
      pmin(., upper[i] - sqrt(.Machine$double.eps)) %>%
      pmax(., lower[i] + sqrt(.Machine$double.eps))
  } %>%
    matrix(., nrow = n, ncol = length(lower))
}


# EIVar acquisition function
EIwithVar <- function(new_val,GP_model, type, y.best, gamma=1)
{
  # Prediction for new x
  GP_pred = predict(object=GP_model, newdata=new_val, type= type, checkNames=FALSE)
  GP_mean = GP_pred$mean
  GP_sd = GP_pred$sd
  
  # Compute Exp Improvement
  z = (y.best - GP_mean)/GP_sd
  exp_imprv = (y.best-GP_mean)*pnorm(z) + GP_sd*dnorm(z)
  
  # Compute Var Improvement
  var_impv = GP_sd^2 * ( (z^2 + 1)*pnorm(z) + z*dnorm(z) ) - exp_imprv^2
  
  # Compute EIVar
  acq_val = exp_imprv - 0.5*gamma*var_impv
  
  return (list(acq_val = acq_val, exp_imprv=exp_imprv))
}


# Function to optimize (use negative since optim() will do minimization)
negEIwithVar <- function(new_val,GP_model, type, y.best, gamma=1)
{
  # Prediction for new x
  GP_pred = predict(object=GP_model, newdata=new_val, type= type, checkNames=FALSE)
  GP_mean = GP_pred$mean
  GP_sd = GP_pred$sd
  
  # Compute Exp Improvement
  z = (y.best - GP_mean)/GP_sd
  exp_imprv = (y.best-GP_mean)*pnorm(z) + GP_sd*dnorm(z)
  # Compute Var Improvement
  var_impv = GP_sd^2 * ( (z^2 + 1)*pnorm(z) + z*dnorm(z) ) - exp_imprv^2
  # Compute EIVar
  acq_val = exp_imprv - 0.5*gamma*var_impv
  
  return (-acq_val) # return negative since optim() will do minimization
}


# Maximization of EIVar
max_EIwithVar <- function(lower, upper, GP_model,type, y.best, gamma=1)
{
  # Generate random initial values
  mat_tries = Matrix_runif(100,lower,upper)
  
  # NegEIwithVar minimization
  mat_optim = foreach(i = 1:nrow(mat_tries), .combine = "rbind") %do%{
    optim_result <- optim(par = mat_tries[i,],
                          fn = negEIwithVar,
                          GP_model = GP_model, type=type, y.best=y.best,gamma=gamma,
                          method = "L-BFGS-B",
                          lower = lower,
                          upper = upper,
                          control = list(maxit = 100,
                                         factr = 5e11))
    c(optim_result$par, optim_result$value)
  }
  
  # Return best set
  ind.argmax = which.min(mat_optim[,ncol(mat_optim)])
  argmax = mat_optim[ind.argmax,1]
  
  return(argmax)
}