# Unbiased Log Gradients for Inverse Problems
# Deng Lu
# version 10
# 20,August,2019
# introduce parameter theta : the observation noise covariance \Tau=\theta*I

rm(list=ls())

library(MASS)
library(parallel)
library(doParallel)
library(lubridate)
# library(Rcpp)
library(smccbase15)

source("ModelPoiEqu.r")
source("smc_body.r")
source("smc_base.r")
source("gene_run.r")
source("visulise.r")
source("compare_model.r")
source("gene_run_sgd.r")
# sourceCpp("sm_Cbase.cpp")

# set parameters
Options = list(
  newData = FALSE ,  # (bool) generate new data
  trueVal = FALSE ,  # (bool) calculate true value via mlsmc
  run_unbia_smc = FALSE ,  # (bool) run unbiased estimation
  unbiased_index = 4 , # 1:unbiased; 2:unbiased with infinity Lmax; 3:unbiased coupled sum; 4:unbiased coupled sum with infinity Lmax
  run_compare_mlsmc = FALSE ,  # (bool) compare run mlsmc
  mlsmc_l_upper = 4 ,  ## l_upper when calculate mlsmc
  run_compare_smc = FALSE ,  # (bool) compare run smc
  smc_l_upper = 1,  ## l_upper when calculate smc
  check_mlsmc = FALSE ,  # (bool) check mlsmc rate
  run_test = FALSE  # (bool) run test mlsmc vs smc
)


Options_unbiased = list(
  unbia_smc_Lmax = 2 ,  # Lmax when do unbiased algorithm
  unbia_smc_Pmax_fixed = TRUE , # (bool) set Pmax as fixed when run unbiased estimation with infinity Lmax
  unbia_smc_Pmax = 2 ,  # Pmax when do unbiased algorithm with infinity Lmax, set Pmax as fixed
  C_max = 10,  # if Pmax is not fixed, Pmax = C-(beta-1)l
  unbia_smc_N = 1*10^0  # number of samples when do unbiased algorithm
)


Options_SGD = list(
  run_SGD = FALSE , # (bool) run SGD
  unbiased_index = 4 , # 1:unbiased; 2:unbiased with infinity Lmax; 3:unbiased coupled sum; 4:unbiased coupled sum with infinity Lmax
  theta_ini = 1 ,  # start point of SGD
  SGD_stepsize = 0.1,  # step size alpha in SGD
  SGD_stepsize_cons = FALSE ,  # # (bool) set stepsize as constant or not
  SGD_length = 20,  # total steps K in SGD
  run_SGD_MLSMC = TRUE , # (bool) run SGD
  epsilon_mlsmc_sgd = sqrt(2^(-c(12,13,14,15,16,17,18,19,20,22))),  # epsilon when do MLSMC
  mlsmc_l_upper_sgd = 2  # l_upper when calculate mlsmc
)

Param = list(
  cores = 4 ,  # parall cores
  ini_sigma = 0.3 ,  # sigma of mcmc kernel when sample from eta_0
  mcmc_sigma = 0.1 ,   # sigma of mcmc kernel in smc
  ini_mcmc = 50 ,  # initial mcmc runs
  B_burn = 3 ,  # burning size of number of samples
  cons = 4 , 
  l_level_trueVal = 10 ,  # l_upper when calculate true value
  mcmc_step = 3 ,  # mcmc steps in each MCMC kernel
  epsilon_mlsmc = sqrt(2^(-c(12,13,14,15,16,17,18,19,20,22))),  # epsilon when test mlsmc vs smc
  epsilon = sqrt(2^(-c(12,13,14,15,16)))
)

Mod_param = list(
  g_func_type = 2 ,  # 1: g function in paper, 2: g function in inverse problem
  big_k = 2 ,  # number of u_i
  small_k = 2 ,  # parameter in step-size
  theta_0 = 2 ,  # the observation noise covariance \Tau=\theta^{-1}*I
  big_M = 2 ,  # dimension of data
  # gp_x = seq(1/51, 1-1/51, by=1/51) ,  # G(u)=[p(x_1;u), p(x_2;u), p(x_M;u)]
  gp_x = c(0.25, 0.75) ,  # G(u)=[p(x_1;u), p(x_2;u), p(x_M;u)]
  theta_mu = 0,  #  theta is log-normal(mu, sigma^2)
  theta_sigma = 1, 
  x_g = 0.5 ,  # g(u) = p(0.5;u)
  var_pow = 4 ,   # variance power
  cost_pow = 1 ,  # cost power
  bias_pow = 2 ,  # power of bias
  l_g = 10,  # level l in g(u)
  l_simu = 10  # level l in simulate data
)

##########################################################################
# run unbiased algorithm
if (Options$run_unbia_smc){  
  algo_result = generate_and_run_unbiased(Param, Mod_param, Options, Options_unbiased)
}

if(Options_SGD$run_SGD){  # run SGD with unbiased estimator
  results_SGD = generate_and_run_SGD_unbiased(Param, Mod_param, Options, Options_SGD, Options_unbiased)
}

if(Options_SGD$run_SGD_MLSMC){  # run SGD with MLSMC
  results_SGD_mlsmc = generate_and_run_SGD_mlsmc(Param, Mod_param, Options, Options_SGD)
}

# if (Options$run_mlsmc){  # run mlsmc
#   runmlsmc(Param, Mod_param, Options)
# }

# if (Options$run_test){  # test mlsmc vs smc
#   test_mlsmc(Param, Mod_param, Options)
# # visulize_test( Param , Mod_param , Options )
# }

if (Options$trueVal){  # calculate truevalue
  generate_trueval(Param, Mod_param, Options)
}

if (Options$run_compare_mlsmc | Options$run_compare_smc){  # compare mlsmc vs smc
  compare_mlsmc_smc2(Param, Mod_param, Options)
}

if (Options$check_mlsmc){  # check mlsmc variance rate
  check_mlsmc_rate(Param, Mod_param, Options)
}

save.image( file = "myFinalRenv.rda" )

