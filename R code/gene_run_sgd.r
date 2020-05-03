# Unbiased Log Gradients for Inverse Problems
# SGD
# Deng Lu
# 5,Dec,2019


unbiased_est_withtheta <- function(param, mod_param, option, option_unbiased, option_sgd, mytheta, yt, Lmax, N_epsilon){
  # unbiased estimator using coupled samples
  yhat_est = c()  # store unbiased estimators
  L_record = c()  # store l
  N_record = c()  # store n
  cost_record1 = c()  # store cost time
  cost_record2 = c()  # store cost operation
  
  L_prob = rep(0 , Lmax + 1)  # probability of l
  var_pow = mod_param$var_pow  # rate of variance
  cost_pow = mod_param$cost_pow  # rate of cost
  for (l in 0 : Lmax){  # P_L(l) \propto 2^{-l(beta+gamma)/2}
    L_prob[l + 1] = 2^(-l*(var_pow + cost_pow)/2)
  }
  L_prob = L_prob / sum(L_prob)  # normalized probability
  
  # unbiase_record <- foreach(x=1:N_epsilon, .packages=c("smccbase10",'lubridate')) %dopar% {source("smc_body.r");source("smc_base.r");unbiase_generate_coupled_withtheta(param, mod_param, option, yt, Lmax, L_prob, mytheta)}
  if(option_sgd$unbiased_index == 1 | option_sgd$unbiased_index == 2){  # run unbiased estimator or with infinity Lmax
    unbiase_record <- foreach(x=1:N_epsilon, .packages=c("smccbase14",'lubridate')) %dopar% {source("smc_body.r");source("smc_base.r");unbiase_generate_withtheta(param, mod_param, option, option_unbiased, yt, Lmax, L_prob, mytheta)} 
  }else if(option_sgd$unbiased_index == 3 | option_sgd$unbiased_index == 4){  # run unbiased coupled sum estimator or with infinity Lmax
    unbiase_record <- foreach(x=1:N_epsilon, .packages=c("smccbase14",'lubridate')) %dopar% {source("smc_body.r");source("smc_base.r");unbiase_generate_coupled_withtheta(param, mod_param, option, option_unbiased, yt, Lmax, L_prob, mytheta)}
  }
  
  unbia_est_rec = c()
  unbia_cost2_rec = c()
  for(imc in 1:N_epsilon){
    yhat_est[imc] = unbiase_record[[imc]]$est  # record estmator
    L_record[imc] = unbiase_record[[imc]]$L_record  # record l_i
    N_record[imc] = unbiase_record[[imc]]$N_record  # record p_i
    cost_record1[imc] = unbiase_record[[imc]]$cost_record1  # record cost in seconds
    cost_record2[imc] = unbiase_record[[imc]]$cost_record2  # record cost in operation
  }
  
  unbia_cost2_rec <- mysum_c(cost_record2)  # calculate accumulative cost
  unbia_est_rec <- mymean_c(yhat_est)  # calculate accumulative estimates of unbiased
  
  return(list(est=unbia_est_rec, L_record=L_record, N_record=N_record, cost_record1=cost_record1,
              cost_record2=unbia_cost2_rec))
}


generate_and_run_SGD_unbiased <- function(param, mod_param, option, option_sgd, option_unbiased){
  # SGD
  # use unbiased estimator
  # theta = exp(gamma)
  if(option_sgd$unbiased_index == 1 | option_sgd$unbiased_index == 3){
    Lmax = option_unbiased$unbia_smc_Lmax  # upper bound of l when calculate unbiased estimator
  }else{
    Lmax = 20  # upper bound of l when calculate unbiased estimator for infinity Lmax
  }
  N_epsilon = option_unbiased$unbia_smc_N  # number of samples of unbiased estimator
  theta_ini = option_sgd$theta_ini  # start point of SGD
  sgd_stepsize = option_sgd$SGD_stepsize  # step size alpha in SGD
  sgd_length = option_sgd$SGD_length
  if(option_sgd$SGD_stepsize_cons){
    sgd_stepsize = rep(option_sgd$SGD_stepsize, sgd_length)   # step size alpha in SGD
  }else{  # if stepsize is not constant, set alpha_k = 1/k
    sgd_stepsize = rep(1, sgd_length)
    for(k in 1:sgd_length){
      sgd_stepsize[k] = 0.1 * 1/k
    }
  }
  
  # get work space
  workdir = getwd()
  # workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
  # iMC = as.numeric(workdir1[length(workdir1)])
  # workdir2 = ''
  # for (i in 1:(length(workdir1) - 1)){
    # workdir2 = paste0(workdir2, workdir1[i], '/')
  # }
  # workdir2 = paste0(workdir2, 'MC')
  workdir2 = paste0(workdir, 'MC')
  
  # generate or load data depending on options['newData']
  if (option$newData){  # generate data and store
    Yt = generate_data_example(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{  # load data
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
  
  coress = param$cores  # parallel computing
  cl <- makeCluster(coress)
  registerDoParallel(cl)
  
  time1 <- Sys.time()  # starting time
  theta_record = rep(0, sgd_length+1)  # record theta_k = exp(gamma_k)
  gamma_record = rep(0, sgd_length+1)  # record gamma_k
  cost_operation <- c()  # record cost in operation
  cost_time <- c()  # record cost in time
  theta_record[1] = theta_ini
  gamma_record[1] = log(theta_ini)
  for(k in 1:sgd_length){
    results_k = unbiased_est_withtheta(param, mod_param, option, option_unbiased, option_sgd, theta_record[k], Yt, Lmax, N_epsilon)
    gamma_record[k+1] = gamma_record[k] + sgd_stepsize[k] * results_k$est[length(results_k$est)]*exp(gamma_record[k])
    theta_record[k+1] = exp(gamma_record[k+1])
    if(k == 1){
      cost_operation[k] = results_k$cost_record2[length(results_k$cost_record2)]
    }else{
      cost_operation[k] = cost_operation[k-1] + results_k$cost_record2[length(results_k$cost_record2)]
    }
    time2 = Sys.time()
    runningtime = time2 - time1
    cost_time[k] = time_length(runningtime, unit='second')  # record cost time
    # cost_operation <- c(cost_operation, results_k$cost_record2[length(results_k$cost_record2)])
  }
  time2 <- Sys.time()
  runningtime <- time2 - time1
  
  write.table(theta_record, file=paste0(workdir2, '/unbiased_theta_record_run_', iMC, '.txt'))
  write.table(gamma_record, file=paste0(workdir2, '/unbiased_gamma_record_run_', iMC, '.txt'))
  write.table(cost_time, file=paste0(workdir2, '/unbiased_cost_time_run_', iMC, '.txt'))
  write.table(cost_operation, file=paste0(workdir2, '/unbiased_cost_operation_run_', iMC, '.txt'))
  
  stopCluster( cl )  # stop parallel
  return(list(theta_record=theta_record, gamma_record=gamma_record, cost_time=cost_time, cost_operation=cost_operation))
}


mlsmc_withtheta <- function(param, mod_param, option, option_sgd, yt, mytheta){
  # run MLSMC with theta
  if (mod_param$g_func_type == 1){  # g function in paper
    g_func = g_func_pap_c
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    g_func = g_func_inv_c
  }
  
  epsilon_ind = option_sgd$epsilon_mlsmc_sgd  # epsilon when calculate number of samples
  mlsmc_l_upper = option_sgd$mlsmc_l_upper_sgd  # upper bound L when do MLSMC
  time_imc_1 <- Sys.time()
  imc_temp = mlsmc_body_withtheta(param, mod_param, option, yt, mlsmc_l_upper, 1, TRUE, epsilon_ind[mlsmc_l_upper], mytheta)
  results_sample = imc_temp$smc_sample
  results_weights = imc_temp$smc_weights
  # calculate estimate
  samples_0 = results_sample[[1]]
  weights_0 = results_weights[[1]]
  yhat_est = mlsmc_est_smc_0_withtheta(samples_0, weights_0, g_func, param, mod_param, 1, yt, mytheta)
  if (mlsmc_l_upper > 1){
    for(ll in 2:mlsmc_l_upper){
      samples_ll = results_sample[[ll]]
      weights_ll = results_weights[[ll]]
      yhat_est = yhat_est + mlsmc_est_withtheta(samples_ll, weights_ll, g_func, param, mod_param, ll, yt, mytheta)
    }
  }
  mlsmc_record = yhat_est
  time_imc_2 <- Sys.time()
  time_diff = time_imc_2 - time_imc_1
  cost_mlsmc_record = time_length(time_diff, unit = 'second') 
  
  if(mlsmc_l_upper==0){
    cost_record2 = sum(imc_temp$number_of_sample[1])
  }else{
    cost_record2 = 0
    for(jj in 1:mlsmc_l_upper){
      cost_record2 = cost_record2+ imc_temp$number_of_sample[jj]*2^(jj-1)
    }
  }
  
  return(list(est=mlsmc_record, cost_time=cost_mlsmc_record, cost_operation=cost_record2))
}


generate_and_run_SGD_mlsmc <- function(param, mod_param, option, option_sgd){
  # SGD
  # use MLSMC
  theta_ini = option_sgd$theta_ini  # start point of SGD
  sgd_stepsize = option_sgd$SGD_stepsize  # step size alpha in SGD
  sgd_length = option_sgd$SGD_length
  if(option_sgd$SGD_stepsize_cons){
    sgd_stepsize = rep(option_sgd$SGD_stepsize, sgd_length)   # step size alpha in SGD
  }else{  # if stepsize is not constant, set alpha_k = 1/k
    sgd_stepsize = rep(1, sgd_length)
    for(k in 1:sgd_length){
      sgd_stepsize[k] = 0.1 * 1/k
    }
  }
  
  # get work space
  workdir = getwd()
  # workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
  # iMC = as.numeric(workdir1[length(workdir1)])
  # workdir2 = ''
  # for (i in 1:(length(workdir1) - 1)){
    # workdir2 = paste0(workdir2, workdir1[i], '/')
  # }
  # workdir2 = paste0(workdir2, 'MC')
  workdir2 = paste0(workdir, 'MC')
  
  # generate or load data depending on options['newData']
  if (option$newData){  # generate data and store
    Yt = generate_data_example(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{  # load data
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
  
  coress = param$cores  # parallel computing
  cl <- makeCluster(coress)
  registerDoParallel(cl)
  
  time1 <- Sys.time()  # starting time
  theta_record = rep(0, sgd_length+1)  # record theta_k = exp(gamma_k)
  gamma_record = rep(0, sgd_length+1)  # record gamma_k
  cost_operation <- c()  # record cost in operation
  cost_time <- c()  # record cost in time
  theta_record[1] = theta_ini
  gamma_record[1] = log(theta_ini)
  for(k in 1:sgd_length){
    results_k = mlsmc_withtheta(param, mod_param, option, option_sgd, Yt, theta_record[k])
    gamma_record[k+1] = gamma_record[k] + sgd_stepsize[k] * results_k$est*exp(gamma_record[k])
    theta_record[k+1] = exp(gamma_record[k+1])
    if(k == 1){
      cost_operation[k] = results_k$cost_operation
    }else{
      cost_operation[k] = cost_operation[k-1] + results_k$cost_operation
    }
    time2 = Sys.time()
    runningtime = time2 - time1
    cost_time[k] = time_length(runningtime, unit='second')  # record cost time
  }
  
  write.table(theta_record, file=paste0(workdir2, '/MLSMC_theta_record_run', iMC, '.txt'))
  write.table(gamma_record, file=paste0(workdir2, '/MLSMC_gamma_record_run', iMC, '.txt'))
  write.table(cost_time, file=paste0(workdir2, '/MLSMC_cost_time_run', iMC, '.txt'))
  write.table(cost_operation, file=paste0(workdir2, '/MLSMC_cost_operation_run', iMC, '.txt'))
  
  stopCluster( cl )  # stop parallel
  return(list(theta_record=theta_record, gamma_record= gamma_record, cost_time=cost_time, cost_operation=cost_operation))
}


pdf_norm <- function(x){
  return (exp(-x^2))
}

error_func <- function(x){  
  # error function erf(x) = 2/\sqrt{\pi} * \int_{0}^{x}e^{-t^2}dt
  return( (2/sqrt(pi)) * integrate(pdf_norm, 0, x)$value )
}

SGD_r_theta_log <- function(mytheta, y, theta_mu, theta_sigma){
  # - log likelihood function
  big_M = length(y)  # dimension of data
  gp_x = seq(1/(big_M+1), 1-1/(big_M+1), by=1/(big_M+1))  # G(u)=[ p(x1;u),p(x2;u),\dots,p(x_M;u)]
  G_u = sapply(gp_x, function(z) {0.5*(z^2-z)})
  log_r = 0.5 * (big_M-3) * log(mytheta) - 0.5*mytheta*(sum(y^2)-sum(G_u*y)^2/sum(G_u^2)) - 0.5*(log(mytheta)-theta_mu)^2 / theta_sigma^2
  
  s1 = sqrt(0.5*mytheta) * sqrt(sum(G_u^2)) * (1-sum(G_u*y)/sum(G_u^2))
  s2 = sqrt(0.5*mytheta) * sqrt(sum(G_u^2)) * (-1-sum(G_u*y)/sum(G_u^2)) 
  
  return (-(log_r + log(error_func(s1)-error_func(s2)))) # return - log(\gamma(\theta))
}

SGD_example_reference1 <- function(option, option_sgd, mod_param){
  # calculate reference solution of SGD
  # minimize - log likelihood
  # get work space
  workdir = getwd()
  # workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
  # iMC = as.numeric(workdir1[length(workdir1)])
  # workdir2 = ''
  # for (i in 1:(length(workdir1) - 1)){
    # workdir2 = paste0(workdir2, workdir1[i], '/')
  # }
  # workdir2 = paste0(workdir2, 'MC')
  workdir2 = paste0(workdir, 'MC')
  
  # generate or load data depending on options['newData']
  if (option$newData){  # generate data and store
    Yt = generate_data_example(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{  # load data
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
  
  theta_ini = option_sgd$theta_ini 
  theta_mu = mod_param$theta_mu
  theta_sigma = mod_param$theta_sigma
  sgd_mle = optim(theta_ini, SGD_r_theta_log, y=Yt, theta_mu=theta_mu, theta_sigma=theta_sigma, method='BFGS')
  write.table(sgd_mle$par, paste0(workdir2, '/SGD_example_truevalue.txt'))
  return (sgd_mle)
}


SGD_r_theta_log_derivative <- function(mytheta, y, theta_mu, theta_sigma){
  # derivative of log likelihood
  big_M = length(y)  # dimension of data
  gp_x = seq(1/(big_M+1), 1-1/(big_M+1), by=1/(big_M+1))  # G(u)=[ p(x1;u),p(x2;u),\dots,p(x_M;u)]
  G_u = sapply(gp_x, function(z) {0.5*(z^2-z)})
  log_r = 0.5 * (big_M-3) / mytheta - 0.5 * (sum(y^2)-sum(G_u*y)^2/sum(G_u^2)) - (log(mytheta)-theta_mu)/(mytheta*theta_sigma^2)
  
  s1 = sqrt(0.5*mytheta) * sqrt(sum(G_u^2)) * (1-sum(G_u*y)/sum(G_u^2)) 
  s2 = sqrt(0.5*mytheta) * sqrt(sum(G_u^2)) * (-1-sum(G_u*y)/sum(G_u^2))
  
  temp1 = exp(-s1^2) * sqrt(sum(G_u^2)) * (1-sum(G_u*y)/sum(G_u^2)) / (2*sqrt(2*mytheta))
  temp2 = exp(-s2^2) * sqrt(sum(G_u^2)) * (-1-sum(G_u*y)/sum(G_u^2)) / (2*sqrt(2*mytheta))
  log_r = log_r + (2/(error_func(s1)-error_func(s2))) * (1/sqrt(pi)) * (temp1 - temp2) 
  
  return (log_r)
}


SGD_example_reference2 <- function(option, option_sgd){
  # calculate reference solution of SGD
  # find zero point of derivative of log-likelihood
  # get work space
  workdir = getwd()
  # workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
  # iMC = as.numeric(workdir1[length(workdir1)])
  # workdir2 = ''
  # for (i in 1:(length(workdir1) - 1)){
    # workdir2 = paste0(workdir2, workdir1[i], '/')
  # }
  # workdir2 = paste0(workdir2, 'MC')
  workdir2 = paste0(workdir, 'MC')
  
  # generate or load data depending on options['newData']
  if (option$newData){  # generate data and store
    Yt = generate_data_example(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{  # load data
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
  
  theta_mu = mod_param$theta_mu
  theta_sigma = mod_param$theta_sigma
  sgd_mle2 = uniroot(SGD_r_theta_log_derivative, interval=c(0,100), y=Yt, theta_mu=theta_mu, theta_sigma=theta_sigma)
  write.table(sgd_mle2$root, paste0(workdir2, '/SGD_example_truevalue.txt'))
  return (sgd_mle)
}


# SGD_r_theta_u <- function(u, y ,mytheta){
#   # r(theta,u)
#   r_theta_u = rep(0, length(u))
#   for (i in 1:length(u)){
#     big_M = length(y)  # dimension of data
#     gp_x = seq(1/(big_M+1), 1-1/(big_M+1), by=1/(big_M+1))  # G(u)=[ p(x1;u),p(x2;u),\dots,p(x_M;u)]
#     G_u = sapply(gp_x, function(z) {0.5*(z^2-z)})
#     
#     r_theta_u[i] = mytheta^(0.5*big_M) * exp(-0.5*mytheta*sum((G_u*u[i]-y)^2))
#   }
#   
#   return (r_theta_u)
# }
# 
# SGD_r_theta <- function(mytheta, y, theta_mu, theta_sigma){
#   # r(theta)
#   inte = integrate(SGD_r_theta_u, lower=-1, upper=1, y=y, mytheta=mytheta)$value
#   inte = inte * exp(-(log(mytheta)-theta_mu)^2 / (2*theta_sigma^2)) / mytheta
#   return (-log(inte))
# }
# 
# 
# SGD_example_reference3 <- function(option, option_sgd, mod_param){
#   # calculate reference solution of SGD
#   # get work space
#   workdir = getwd()
#   workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
#   iMC = as.numeric(workdir1[length(workdir1)])
#   workdir2 = ''
#   for (i in 1:(length(workdir1) - 1)){
#     workdir2 = paste0(workdir2, workdir1[i], '/')
#   }
#   workdir2 = paste0(workdir2, 'MC')
# 
#   # generate or load data depending on options['newData']
#   if (option$newData){  # generate data and store
#     Yt = generate_data_example(param, mod_param)
#     write.table(Yt, paste0(workdir2, '/observation.txt'))
#   }else{  # load data
#     Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
#   }
#   
#   theta_ini = option_sgd$theta_ini 
#   theta_mu = mod_param$theta_mu
#   theta_sigma = mod_param$theta_sigma
#   sgd_mle3 = optim(theta_ini, SGD_r_theta, y=Yt, theta_mu=theta_mu, theta_sigma=theta_sigma, method='BFGS')
#   write.table(sgd_mle$root, paste0(workdir2, '/SGD_example_truevalue.txt'))
#   return (sgd_mle)
# }
