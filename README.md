# BEO-EIVar-and-AdjBO
Data and source codes for the paper ``Bayesian Estimation and Optimization for Learning Sequential Regularized Portfolios" by Godeliva Petrina Marisu and Chi Seng Pun (Nanyang Technological University, Singapore), where we developed a new acquisition function EIVar that targets high expected improvement with low improvement variance and we proposed adjustments to the Bayesian optimization to account for the sequential feature of portfolio constructions.

Main runner class:
1. BayesLPO_BayesOpt_prelimempstudy.R :
Preliminary empirical studies (3 months test points) to compare CV and different BO variations (with/without time series adjustment, EI acquisition function with/without variance). Used to tune LPO with Sample, Bayes-Stein, Bayes Asset Pricing (FF3) and Bayes Averaging estimate.
2. BayesLPO_BayesOpt_1yearempstudy.R:
1 year empirical study for comparing best BO variation from prelim empirical study in #1 (AdjBO-EIVar) with CV. Used to tune LPO with Sample, Bayes-Stein, Bayes Asset Pricing (FF3) and Bayes Averaging estimate.

Main Dependencies:
1. LPO_SampleEstimation.R : Find LPO weights with sample estimate of covar matrix.
2. LPO_BayesEstimation.R : Find LPO weights with Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate of covar matrix.
3. LPO_CV.R : Find optimal lambda using CV for LPO with sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate.
4. LPO_CV_TestEvaluation.R (depends on #3): Evaluate out-of-sample performance of CV when tuning LPO with Sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate.
5. LPO_BO_NoisyOptimizer.R (depends on #7) : Find optimal lambda using Bayes Opt for LPO with sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate.
6. LPO_BO_TestEvaluation.R (depends on #5) : Evaluate out-of-sample performance of Bayes Opt (BO) when tuning LPO with Sample, Bayes-Stein, Bayes Asset Pricing and Bayes Averaging estimate.
7. noisy_optimizer_withStopping.R : Extension of noisy.optimizer function of DiceOptim package. Extensions made: Add new acquisition function EI with Variance (EIVar), add stopping condition.
7. LPO_TrueLambdaBacktest.R: Find true optimal lambda by incorporating future asset prices.
