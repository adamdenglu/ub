# Unbiased Log Gradients for Inverse Problems
# Deng Lu
# 15,August,2019


generate_and_run_unbiased <- function(param, mod_param, option, option_unbiased){
  # unbiased estimator
  if(option$unbiased_index == 1 | option$unbiased_index == 3){
    Lmax = option_unbiased$unbia_smc_Lmax  # upper bound of l when calculate unbiased estimator
  }else{
    Lmax = 20  # upper bound of l when calculate unbiased estimator for infinity Lmax
  }
  mytheta = mod_param$theta_0
  
  # get work space
  workdir = getwd()
  # workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
  # iMC = as.numeric(workdir1[length(workdir1)])
  # workdir2 = ''
  # for (i in 1:(length(workdir1) - 1)){
    # workdir2 = paste0(workdir2, workdir1[i], '/')
  # }
  workdir3 = paste0(workdir, 'MC2')
  workdir2 = paste0(workdir, 'MC')
  
  # generate or load data depending on options['newData']
  if (option$newData){  # generate data and store
    Yt = generate_data(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{  # load data
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    # truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
    truevalue = read.table(paste0(workdir2, '/trueValue_inverse_example.txt'))[,1]
  }
  # truevalue_est = mean(truevalue[, 1])
  truevalue_est = truevalue
  
  coress = param$cores  # parallel computing
  cl <- makeCluster(coress)
  registerDoParallel(cl)
  
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
  
  N_epsilon = option_unbiased$unbia_smc_N  # number of samples of unbiased estimator
  
  if(option$unbiased_index == 1 | option$unbiased_index == 2){  # run unbiased estimator or with infinity Lmax
    unbiase_record <- foreach(x=1:N_epsilon, .packages=c("smccbase14",'lubridate')) %dopar% {source("smc_body.r");source("smc_base.r");unbiase_generate_withtheta(param, mod_param, option, option_unbiased, Yt, Lmax, L_prob, mytheta)} 
  }else if(option$unbiased_index == 3 | option$unbiased_index == 4){  # run unbiased coupled sum estimator or with infinity Lmax
    unbiase_record <- foreach(x=1:N_epsilon, .packages=c("smccbase14",'lubridate')) %dopar% {source("smc_body.r");source("smc_base.r");unbiase_generate_coupled_withtheta(param, mod_param, option, option_unbiased, Yt, Lmax, L_prob, mytheta)}
  }
 
  unbia_est_rec = c()
  unbia_mse_rec = c()
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
  unbia_mse_rec = (unbia_est_rec- truevalue_est) ^ 2  # calculate mse of unbiased
  
  if(option$unbiased_index == 1){
    write.table(unbia_est_rec, file=paste0(workdir3, '/unbiased_Lmax_', Lmax, 'est_run_', iMC, '.txt'))
    write.table(unbia_mse_rec, file=paste0(workdir2, '/unbiased_Lmax_', Lmax, 'mse_run_', iMC, '.txt'))
    write.table(unbia_cost2_rec, file=paste0(workdir2, '/unbiased_Lmax_', Lmax, 'cost2_run_', iMC, '.txt'))
  }else if(option$unbiased_index == 3){
    write.table(unbia_est_rec, file=paste0(workdir3, '/unbiased_coupled_Lmax_', Lmax, 'est_run_', iMC, '.txt'))
    write.table(unbia_mse_rec, file=paste0(workdir2, '/unbiased_coupled_Lmax_', Lmax, 'mse_run_', iMC, '.txt'))
    write.table(unbia_cost2_rec, file=paste0(workdir2, '/unbiased_coupled_Lmax_', Lmax, 'cost2_run_', iMC, '.txt'))
  }else if(option$unbiased_index == 2){
    Pmax = option_unbiased$unbia_smc_Pmax
    write.table(unbia_est_rec, file=paste0(workdir3, '/unbiased_infinityLmax_Pmax_', Pmax, 'est_run_', iMC, '.txt'))
    write.table(unbia_mse_rec, file=paste0(workdir2, '/unbiased_infinityLmax_Pmax_', Pmax, 'mse_run_', iMC, '.txt'))
    write.table(unbia_cost2_rec, file=paste0(workdir2, '/unbiased_infinityLmax_Pmax_', Pmax, 'cost2_run_', iMC, '.txt'))
  }else if(option$unbiased_index == 4){
    Pmax = option_unbiased$unbia_smc_Pmax
    write.table(unbia_est_rec, file=paste0(workdir3, '/unbiased_infinityLmax_coupled_Pmax_', Pmax, 'est_run_', iMC, '.txt'))
    write.table(unbia_mse_rec, file=paste0(workdir2, '/unbiased_infinityLmax_coupled_Pmax_', Pmax, 'mse_run_', iMC, '.txt'))
    write.table(unbia_cost2_rec, file=paste0(workdir2, '/unbiased_infinityLmax_coupled_Pmax_', Pmax, 'cost2_run_', iMC, '.txt'))
  }

  stopCluster( cl )  # stop parallel
  return(list(est=yhat_est, L_record=L_record, N_record=N_record, cost_record1=cost_record1,
              cost_record2=cost_record2))
}


generate_trueval <- function(param, mod_param, option){
  # calculate true value 
  workdir = getwd()
  true_mc = 10  # number of mc to calculate trueval
  
  # generate or load data depending on options['newData']
  if (option$newData){
    Yt = generate_data(param, mod_param)
    write.table(Yt, paste0(workdir, '/MC/observation.txt'))
  }else{
    Yt = read.table(paste0(workdir, '/MC/observation.txt'))[, 1]
  }

  coress = param$cores  # parallel computing
  cl <- makeCluster(coress)
  registerDoParallel(cl)
  
  l_level = param$l_level_trueVal  # l_upper when calculate true value
  eps = 0.0003  # epsilon when calculate true value
  trueVal_record = rep(0 , true_mc)
  for (iii in 1:true_mc){
    true_result = ml_smc(param, mod_param, option, Yt, l_level, 1, TRUE, eps)
    trueVal_record[iii] = true_result$trueVal
  }
  write.table(trueVal_record,  paste0(workdir,'/MC/trueValue.txt'))
  stopCluster(cl)
}


check_mlsmc_rate <- function(param, mod_param, option){
  # check variance rate
  # get work space
  workdir = getwd()
  # workdir1 = unlist(strsplit(workdir, split="/", fixed=T))
  # iMC = as.numeric(workdir1[length(workdir1)])  # number of MCs
  # workdir2 = ''
  # for (i in 1:(length(workdir1)-1)){
    # workdir2 = paste0(workdir2, workdir1[i], '/')
  # }
  workdir2 = paste0(workdir, 'MC')

  # generate or load data depending on options['newData']
  if (option$newData){
    Yt = generate_data(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }

	if (mod_param$g_func_type == 1){  # g function in paper
		g_func = g_func_pap_c
	}else if (mod_param$g_func_type == 2){  # g function in inverse problem
		g_func = g_func_inv_c
	}

  coress <- param$cores# parallel computing
  cl <- makeCluster( coress )
  registerDoParallel( cl )
  l_upper <- 7  # l_upper when do MLSMC
  epsilon_ind <- param$epsilon_mlsmc  # epsilons when calculate number of samples
  mytheta = mod_param$theta_0

  imc_results = mlsmc_body_withtheta(param, mod_param, option, Yt, l_upper, 1, TRUE, epsilon_ind[l_upper], mytheta)
  results_sample = imc_results$smc_sample
  results_weights = imc_results$smc_weights
  # calculate var
  est_re = rep(0, l_upper)
  for (l in 1:l_upper){
    samples_ll = results_sample[[l]]
    weights_ll = results_weights[[l]]
    est_re[l] = mlsmc_est_withtheta(samples_ll, weights_ll, g_func, param, mod_param, l, Yt, mytheta)
  }
  write.table(est_re, file=paste0(workdir2, '/checkrate_iMC', iMC, '.txt'))
  stopCluster(cl)
}


generate_trueval_example <- function(param, mod_param, option){
  # calculate true value of inverse problem (example one)
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
  
  mytheta = mod_param$'theta_0'
  theta_mu = mod_param$theta_mu
  theta_sigma = mod_param$theta_sigma
  big_M = length(Yt)  # dimension of data
  gp_x = seq(1/(big_M+1), 1-1/(big_M+1), by=1/(big_M+1))  # G(u)=[ p(x1;u),p(x2;u),\dots,p(x_M;u)]
  G_u = sapply(gp_x, function(z) {0.5*(z^2-z)})
  log_r = 0.5 * (big_M-3) / mytheta - 0.5 * (sum(Yt^2)-sum(G_u*Yt)^2/sum(G_u^2)) - (log(mytheta)-theta_mu)/(mytheta*theta_sigma^2)
  
  s1 = sqrt(0.5*mytheta) * sqrt(sum(G_u^2)) * (1-sum(G_u*Yt)/sum(G_u^2)) 
  s2 = sqrt(0.5*mytheta) * sqrt(sum(G_u^2)) * (-1-sum(G_u*Yt)/sum(G_u^2))
  
  temp1 = exp(-s1^2) * sqrt(sum(G_u^2)) * (1-sum(G_u*Yt)/sum(G_u^2)) 
  temp2 = exp(-s2^2) * sqrt(sum(G_u^2)) * (-1-sum(G_u*Yt)/sum(G_u^2)) 
  log_r = log_r + (1/(error_func(s1)-error_func(s2))) * (1/sqrt(2*pi*mytheta)) * (temp1 - temp2) 
  
  write.table(log_r,  paste0(workdir2,'/trueValue_inverse_example.txt'))
}
