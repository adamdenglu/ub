# Unbiased Log Gradients for Inverse Problems
# Deng Lu
# 7,August,2019
# smc algorithm


mlsmc_body_withtheta <- function(param, mod_param, option, yt, l_upper, n_num, setNofSamp, eps, mytheta){
  # main body of multilevel sequestial Monte Carlo
  # add parameter theta
  big_k = mod_param$big_k
  small_k = mod_param$small_k
  gp_x = mod_param$gp_x
  mcmc_sigma = param$mcmc_sigma
  mcmc_len = param$ini_mcmc
  ini_sigma = param$ini_sigma
  mcmc_step = param$mcmc_step
  
  if (setNofSamp){
    xNofSamp = rep(0, l_upper+1)
    var_pow = mod_param$var_pow
    cost_pow = mod_param$cost_pow
    K_l = 0
    for (ii in 1 : (l_upper + 1)){
      h_l = 2 ^ (-(ii - 1))
      K_l = K_l + h_l^(0.5*(var_pow - cost_pow))
    }
    h_l = 2 ^ (-l_upper)
    # nofSampIter1 = 1000 / ((eps^(-2))*(h_l^(0.5*(var_pow+cost_pow)))*K_l)
    nofSampIter1 = 1
    for ( ll in 1 : (l_upper + 1) ){
      h_l = 2 ^ ( - ( ll - 1 ) )
      xNofSamp[ ll ] = nofSampIter1 * (eps^(-2))*(h_l^(0.5*(var_pow+cost_pow))) * K_l
      xNofSamp[ ll ] = min( xNofSamp[ ll ] , 10^6 )
    }
    xNofSamp = ceiling( xNofSamp )
  }else{
    xNofSamp = n_num * rep(1, (l_upper + 1))
    xNofSamp = ceiling( xNofSamp )
  }
  
  sample_record = vector("list", l_upper + 1)  # record samples & weights
  # step 0, sample from eta_0
  # sample_0 <- foreach( x = 1 : xNofSamp[ 1 ] , .packages = c( "MASS" ) ) %dopar% {source("ModelPoiEqu.r"); source("smc_base.r") ; sample_eta_0( param , mod_param , yt )}
  sample_0 <- foreach(x=1:xNofSamp[1], .packages=c("smccbase14")) %dopar% {sample_eta_0_c(big_k,small_k, gp_x, ini_sigma, mcmc_len, mytheta, yt)}
  
  accept_record_level_0 = rep(0, big_k)
  for (lll in 1:xNofSamp[1]){
    accept_record_level_0 = accept_record_level_0 + sample_0[[lll]]$accept
  }
  accept_record_level_0 = accept_record_level_0 / xNofSamp[1]
  accept_record = matrix(0, nrow=big_k, ncol=max(l_upper,1))  # record accept_rate
  accept_record[, 1] = accept_record_level_0
  
  sample_record[[1]] = sample_0
  
  if (l_upper >= 2){
    for (ll in 1:(l_upper - 1)){
      sample_ll = sample_record[[ll]]
      weights_ll = rep(0, xNofSamp[ll])
      for (lll in 1:xNofSamp[ll]){
        weights_ll[lll] = sample_ll[[lll]]$prob
      }
      weights_ll_nor = exp(weights_ll) / sum(exp(weights_ll))  # normalize weights
      
      # resample step
      sample_ll_1_ind = sample(1:xNofSamp[ll], size=xNofSamp[ll + 1], replace=TRUE, prob=weights_ll_nor)
      sample_ll_1 = vector("list", xNofSamp[ll + 1])
      for (lll in 1:xNofSamp[ll + 1]){
        sample_ll_1[[lll]] = sample_ll[[sample_ll_1_ind[ lll ]]]
      }
      # mcmc step
      # sample_ll_1_1 = foreach( x = 1 : xNofSamp[ ll + 1 ] , .packages = c( "MASS" ) ) %dopar% {source("ModelPoiEqu.r") ; source("smc_base.r"); mcmc_smc( sample_ll_1[[ x ]]$samples , param , mod_param , ll , yt )}
      sample_ll_1_1 = foreach(x=1:xNofSamp[ll + 1], .packages=c("smccbase14")) %dopar% { mcmc_smc_c(sample_ll_1[[x]]$samples, ll, big_k, small_k, gp_x, mcmc_sigma, mytheta, yt, mcmc_step)}
      
      accept_record_level_ll = rep(0, big_k)
      for (lll in 1:xNofSamp[ll + 1]){
        accept_record_level_ll = accept_record_level_ll + sample_ll_1_1[[lll]]$accept_record_l
      }
      accept_record_level_ll = accept_record_level_ll / xNofSamp[ll + 1]
      accept_record[, ll + 1] = accept_record_level_ll
      
      sample_record[[ll + 1]] = sample_ll_1_1  # update samples
    }
  }
  ## the SMC stops
  
  # sort out results
  smc_sample = vector("list", l_upper)
  smc_weights = vector("list", l_upper)
  if (l_upper >= 1){
    for (ll in 1:l_upper){
      sample_level_ll = vector("list", xNofSamp[ll])
      weights_level_ll = rep(0, xNofSamp[ll])
      for (lll in 1:xNofSamp[ll]){
        sample_level_ll[[lll]] = sample_record[[ll]][[lll]]$samples
        weights_level_ll[lll] = sample_record[[ll]][[lll]]$prob
      }
      smc_sample[[ll]] = sample_level_ll
      smc_weights[[ll]] = weights_level_ll
    }
  }else{
    sample_level_ll = vector("list", xNofSamp[1])
    weights_level_ll = rep(0, xNofSamp[1])
    for (lll in 1:xNofSamp[1]){
      sample_level_ll[[lll]] = sample_record[[1]][[lll]]$samples
      weights_level_ll[lll] = sample_record[[1]][[lll]]$prob
    }
    smc_sample[[1]] = sample_level_ll
    smc_weights[[1]] = weights_level_ll
  }
  
  return (list(smc_sample=smc_sample, smc_weights=smc_weights, accept_rate=accept_record, number_of_sample=xNofSamp))
}


ml_smc <- function(param, mod_param, option, yt, l_level, N_num, setNofSamp, eps){
  # original multilevel smc algorithm
  if (mod_param$g_func_type == 1){  # g function in paper
	  g_func = g_func_pap_c
	}else if (mod_param$g_func_type == 2){  # g function in inverse problem
		g_func = g_func_inv_c
	}
  mytheta = mod_param$theta_0

  Results = mlsmc_body_withtheta(param, mod_param, option, yt, l_level, N_num, setNofSamp, eps, mytheta)

  results_sample = Results$smc_sample
  results_weights = Results$smc_weights
  # calculate estimate
  # timestart_es <- Sys.time()
  samples_0 = results_sample[[1]]
  weights_0 = results_weights[[1]]
  yhat_est = mlsmc_est_smc_0_withtheta(samples_0, weights_0, g_func, param, mod_param, 1, yt, mytheta)
  if (l_level > 1){
    for(ll in 2:l_level){
      samples_ll = results_sample[[ll]]
      weights_ll = results_weights[[ll]]
      yhat_est = yhat_est + mlsmc_est_withtheta(samples_ll, weights_ll, g_func, param, mod_param, ll, yt, mytheta)
    }
  }
  return(list(trueVal=yhat_est, results_smc=Results))
}


mlsmc_body_smc_withtheta <- function(param, mod_param, option, yt, l_upper, n_num, setNofSamp, eps, mytheta){
  # main body of multilevel sequestial Monte Carlo
  # only store samples at last level
  big_k = mod_param$big_k
  small_k = mod_param$small_k
  gp_x = mod_param$gp_x
  mcmc_sigma = param$mcmc_sigma
	mcmc_len = param$ini_mcmc
	ini_sigma = param$ini_sigma
  theta_0 = mytheta
  mcmc_step = param$mcmc_step

  xNofSamp = rep(0, l_upper + 1)  # number of samples
  if (setNofSamp){  # calculate number of samples according to paper
    var_pow = mod_param$var_pow
    cost_pow = mod_param$cost_pow
    K_l = 0
    for (ii in 1:(l_upper + 1)){
      h_l = 2 ^ (-(ii-1))
      K_l = K_l + h_l^(0.5*(var_pow - cost_pow))
    }
    h_l = 2 ^ (-l_upper)
    # ensure number of samples at last level is 1000
    nofSampIter1 = 1000 / ((eps^(-2))*(h_l^(0.5*(var_pow+cost_pow)))*K_l)
    # nofSampIter1 = 1
    for (ll in 1:(l_upper + 1)){
      h_l = 2 ^ (-(ll-1))
      xNofSamp[ll] = nofSampIter1 * (eps^(-2)) * (h_l^(0.5*(var_pow+cost_pow))) * K_l
		  xNofSamp[ll] = min(xNofSamp[ll], 10^7)  # ensure number is not too large
	  }
    xNofSamp = ceiling(xNofSamp)  # force the number to be integer
  }else{  # set number of samples to be the same
    xNofSamp = n_num * rep(1, (l_upper + 1))
    xNofSamp = ceiling(xNofSamp)
  }

  sample_record = vector("list", 2)  # record samples & weights

  # step 0, sample from eta_0
  # sample_0 <- foreach( x = 1 : xNofSamp[ 1 ] , .packages = c( "MASS" ) ) %dopar% {source("ModelPoiEqu.r"); source("smc_base.r") ; sample_eta_0( param , mod_param , yt )}
  sample_0 <- foreach(x=1:xNofSamp[1], .packages=c("smccbase14")) %dopar% {sample_eta_0_c(big_k,small_k, gp_x, ini_sigma, mcmc_len,theta_0, yt)}
  sample_record[[ 1 ]] = sample_0  # store samples at level 0

  accept_record_level_0 = rep(0, big_k)  # record accept rate at level 0
  for (lll in 1:xNofSamp[1]){
    accept_record_level_0 = accept_record_level_0 + sample_0[[lll]]$accept
  }
  accept_record_level_0 = accept_record_level_0 / xNofSamp[1]
  accept_record = matrix(0, nrow=big_k, ncol=max(l_upper,1))  # record accept_rate
  accept_record[, 1] = accept_record_level_0

  if (l_upper >= 2){  # if l_upper>2, run smc up to level l_upper-1
    for (ll in 1:(l_upper-1)){
      sample_ll = sample_record[[1]]  # samples at level ll-1
      weights_ll = rep(0, xNofSamp[ll])  # weights at level ll-1
      for (lll in 1:xNofSamp[ll]){
        weights_ll[lll] = sample_ll[[lll]]$prob
      }
      weights_ll_nor = exp(weights_ll) / sum(exp(weights_ll))  # normalize weights

      # resample step
      sample_ll_1_ind = sample(1:xNofSamp[ ll ], size=xNofSamp[ll+1], replace=TRUE, prob=weights_ll_nor)
      sample_ll_1 = vector("list",xNofSamp[ll+1])
      for (lll in 1:xNofSamp[ll + 1]){
        sample_ll_1[[lll]] = sample_ll[[sample_ll_1_ind[lll]]]
      }

      # mcmc step
      # sample_ll_1_1 = foreach( x = 1 : xNofSamp[ ll + 1 ] , .packages = c( "MASS" ) ) %dopar% {source("ModelPoiEqu.r") ; source("smc_base.r"); mcmc_smc( sample_ll_1[[ x ]]$samples , param , mod_param , ll , yt )}
      sample_ll_1_1 = foreach(x=1:xNofSamp[ll+1], .packages=c("smccbase14")) %dopar% { mcmc_smc_c(sample_ll_1[[x]]$samples, ll, big_k, small_k, gp_x, mcmc_sigma, theta_0, yt, mcmc_step)}

      accept_record_level_ll = rep(0, big_k)  # accept_rate at level ll
      for (lll in 1:xNofSamp[ll + 1]){
        accept_record_level_ll = accept_record_level_ll + sample_ll_1_1[[lll]]$accept_record_l
      }
      accept_record_level_ll = accept_record_level_ll / xNofSamp[ll+1]
      accept_record[, ll + 1] = accept_record_level_ll

      sample_record[[1]] = sample_ll_1_1  # update samples
    }
  }
  ## the SMC stops

  # sort out results
  smc_sample = vector("list", 1)  # store samples
  smc_weights = vector("list", 1)  # store weights
  if (l_upper >= 1){
    sample_level_ll = vector("list", xNofSamp[l_upper])
    weights_level_ll = rep(0, xNofSamp[l_upper])
    for (lll in 1:xNofSamp[l_upper]){
      sample_level_ll[[lll]] = sample_record[[1]][[lll]]$samples
      weights_level_ll[lll] = sample_record[[1]][[lll]]$prob
    }
    smc_sample[[1]] = sample_level_ll
    smc_weights[[1]] = weights_level_ll
  }else{
    sample_level_ll = vector("list", xNofSamp[1])
    weights_level_ll = rep(0, xNofSamp[1])
    for (lll in 1:xNofSamp[1]){
      sample_level_ll[[lll]] = sample_record[[1]][[lll]]$samples
      weights_level_ll[lll] = sample_record[[1]][[lll]]$prob
    }
    smc_sample[[1]] = sample_level_ll
    smc_weights[[1]] = weights_level_ll
  }

  return (list(smc_sample=smc_sample, smc_weights=smc_weights, accept_rate=accept_record))
}


mlsmc_body_withtheta_noparallel <- function(param, mod_param, option, yt, l_upper, n_num, setNofSamp, eps, mytheta){
  # main body of multilevel sequestial Monte Carlo
  # no parallel
  # add parameter theta
  big_k = mod_param$big_k
  small_k = mod_param$small_k
  gp_x = mod_param$gp_x
  mcmc_sigma = param$mcmc_sigma
  mcmc_len = param$ini_mcmc
  ini_sigma = param$ini_sigma
  mcmc_step = param$mcmc_step
  
  if (setNofSamp){
    xNofSamp = rep(0, l_upper+1)
    var_pow = mod_param$var_pow
    cost_pow = mod_param$cost_pow
    K_l = 0
    for ( ii in 1 : (l_upper + 1) ){
      h_l = 2 ^ (-(ii - 1))
      K_l = K_l + h_l^(0.5*(var_pow - cost_pow))
    }
    h_l = 2 ^ (-l_upper)
    nofSampIter1 = 1000 / ((eps^(-2))*(h_l^(0.5*(var_pow+cost_pow)))*K_l)
    for ( ll in 1 : (l_upper + 1) ){
      h_l = 2 ^ (-(ll - 1))
      xNofSamp[ll] = nofSampIter1 * (eps^(-2)) * (h_l^(0.5*(var_pow+cost_pow))) * K_l
      xNofSamp[ll] = min( xNofSamp[ll] , 10^7 )
    }
    xNofSamp = ceiling( xNofSamp )
  }else{
    xNofSamp = n_num * rep(1, (l_upper + 1))
    xNofSamp = ceiling( xNofSamp )
  }
  
  sample_record = vector("list", l_upper + 1)  # record samples & weights
  # step 0, sample from eta_0
  sample_0 = vector( "list" , xNofSamp[ 1 ] )
  for ( i in 1: xNofSamp[ 1 ]){
    sample_0[[ i ]] = sample_eta_0_c(big_k,small_k, gp_x, ini_sigma, mcmc_len, mytheta, yt)
  }
  
  accept_record_level_0 = rep(0, big_k)
  for (lll in 1:xNofSamp[1]){
    accept_record_level_0 = accept_record_level_0 + sample_0[[lll]]$accept
  }
  accept_record_level_0 = accept_record_level_0 / xNofSamp[1]
  accept_record = matrix(0, nrow=big_k, ncol=max(l_upper,1))  # record accept_rate
  accept_record[, 1] = accept_record_level_0
  
  sample_record[[1]] = sample_0
  
  if (l_upper >= 2){
    for (ll in 1:(l_upper - 1)){
      sample_ll = sample_record[[ll]]
      weights_ll = rep(0, xNofSamp[ll])
      for (lll in 1:xNofSamp[ll]){
        weights_ll[lll] = sample_ll[[lll]]$prob
      }
      weights_ll_nor = exp(weights_ll) / sum(exp(weights_ll))  # normalize weights
      
      # resample step
      sample_ll_1_ind = sample(1:xNofSamp[ll], size=xNofSamp[ll + 1], replace=TRUE, prob=weights_ll_nor)
      sample_ll_1 = vector("list", xNofSamp[ll + 1])
      for (lll in 1:xNofSamp[ll + 1]){
        sample_ll_1[[lll]] = sample_ll[[sample_ll_1_ind[ lll ]]]
      }
      # mcmc step
      sample_ll_1_1 = vector( "list" , xNofSamp[ ll + 1 ] )
      for ( i in 1: xNofSamp[ ll + 1 ]){
        sample_ll_1_1[[ i ]] = mcmc_smc_c( sample_ll_1[[ i ]]$samples, ll, big_k, small_k, gp_x,  mcmc_sigma, mytheta, yt, mcmc_step)
      }
      
      accept_record_level_ll = rep(0, big_k)
      for (lll in 1:xNofSamp[ll + 1]){
        accept_record_level_ll = accept_record_level_ll + sample_ll_1_1[[lll]]$accept_record_l
      }
      accept_record_level_ll = accept_record_level_ll / xNofSamp[ll + 1]
      accept_record[, ll + 1] = accept_record_level_ll
      
      sample_record[[ll + 1]] = sample_ll_1_1  # update samples
    }
  }
  ## the SMC stops
  
  # sort out results
  smc_sample = vector("list", l_upper)
  smc_weights = vector("list", l_upper)
  if (l_upper >= 1){
    for (ll in 1:l_upper){
      sample_level_ll = vector("list", xNofSamp[ll])
      weights_level_ll = rep(0, xNofSamp[ll])
      for (lll in 1:xNofSamp[ll]){
        sample_level_ll[[lll]] = sample_record[[ll]][[lll]]$samples
        weights_level_ll[lll] = sample_record[[ll]][[lll]]$prob
      }
      smc_sample[[ll]] = sample_level_ll
      smc_weights[[ll]] = weights_level_ll
    }
  }else{
    sample_level_ll = vector("list", xNofSamp[1])
    weights_level_ll = rep(0, xNofSamp[1])
    for (lll in 1:xNofSamp[1]){
      sample_level_ll[[lll]] = sample_record[[1]][[lll]]$samples
      weights_level_ll[lll] = sample_record[[1]][[lll]]$prob
    }
    smc_sample[[1]] = sample_level_ll
    smc_weights[[1]] = weights_level_ll
  }
  
  return (list(smc_sample=smc_sample, smc_weights=smc_weights, accept_rate=accept_record, number_of_sample=xNofSamp))
}


unbiase_generate_withtheta <- function(param, mod_param, option, option_unbiased, Yt, Lmax, L_prob, mytheta){
  var_pow = mod_param$var_pow  # rate of variance
  cons = param$cons
  B_burn = param$B_burn  # buruning size of N
  if (mod_param$g_func_type == 1){  # g function in paper
    g_func = g_func_pap_c
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    g_func = g_func_inv_c
  }

  time_imc_1 <- Sys.time()
  # step 1, sample L and N
  l_upper = sample(c(0 : Lmax), size=1, prob=L_prob)
  if(option$unbiased_index == 1){
    N_max = min(max(var_pow*(Lmax - l_upper)-B_burn-cons, 0), 15)  # Nmax = beta(Lmax-l_i)
  }else if(option$unbiased_index == 2){  # if run infinity Lmax
    N_max = option_unbiased$unbia_smc_Pmax
  }
  P_prob = rep(0, N_max + 1)  # probability of p_i
  for (p in 0:N_max){
    if (p < 4){  # $P_P(p) \propto 2^{-p+4}$ for p <4
      P_prob[p + 1] = 2^(-p+4)
    }else{  # $P_P(p) \propto 2^(-p)*p*[log2(p)]^{2}$ otherwise
      P_prob[p + 1] = 2^(-p)*p*(log2(p))^2
    }
  }
  P_prob = P_prob / sum(P_prob)
  N_num = sample(c(0 : N_max), size=1, prob=P_prob)

  # run smc with l level=l_upper, and number of samples=2^(N_num + B_burn)
  Results = mlsmc_body_withtheta_noparallel(param, mod_param, option, Yt, l_upper, 2^(N_num + B_burn), FALSE, 0.01, mytheta)
  results_sample = Results$smc_sample
  results_weights = Results$smc_weights
  # calculate estimate
  if (l_upper > 0){
    # est = eta^{N_p}(g^{l}G^{l-1})/eta^{N_p}(G^{l-1})-eta^{N_p}(g^{l-1})
    samples_l = results_sample[[l_upper]]
    weights_l = results_weights[[l_upper]]
    est_N = mlsmc_est_withtheta(samples_l, weights_l, g_func, param, mod_param, l_upper, Yt, mytheta)
  }else{
    # est =eta^{N_p}(g^{0})
    samples_0 = results_sample[[1]]
    weights_0 = results_weights[[1]]
    est_N = mlsmc_est_0_withtheta(samples_0, weights_0, g_func, param, mod_param, 0, Yt, mytheta)
  }

  if (N_num > 0){
    # run smc with l level=l_upper, and number of samples=2^(N_num-1+B_burn)
    Results_N1 = mlsmc_body_withtheta_noparallel(param, mod_param, option, Yt, l_upper, 2^(N_num+B_burn-1),
                             FALSE, 0.01, mytheta)
    results_sample_N1 = Results_N1$smc_sample
    results_weights_N1 = Results_N1$smc_weights
    if (l_upper > 0){
      samples_l_N1 = results_sample_N1[[l_upper]]
      weights_l_N1 = results_weights_N1[[l_upper]]
      est_N1 = mlsmc_est_withtheta(samples_l_N1, weights_l_N1, g_func, param, mod_param, l_upper, Yt, mytheta)
    }else{
      samples_0 = results_sample_N1[[1]]
      weights_0 = results_weights_N1[[1]]
      est_N1 = mlsmc_est_0_withtheta(samples_0, weights_0, g_func, param, mod_param, 0, Yt, mytheta)
    }
  }else{
    est_N1 = 0
  }

  unbiase_est = (est_N-est_N1) / (L_prob[l_upper + 1]*P_prob[N_num + 1])

  time_imc_2 <- Sys.time()
  time_diff = time_imc_2 - time_imc_1
  cost_record1 = time_length(time_diff, unit='second')  # record cost
  if(l_upper==0){
    cost_record2 = sum(Results$number_of_sample[1]) + ifelse(N_num>0, sum(Results_N1$number_of_sample[1]), 0)
  }else{
    cost_record2 = 0
    for(jj in 1:l_upper){
      cost_record2 = cost_record2+ Results$number_of_sample[jj]*2^(jj-1) + ifelse(N_num>0, Results_N1$number_of_sample[jj]*2^(jj-1), 0)
    }
  }

  return(list(est=unbiase_est, L_record=l_upper, N_record=N_num, cost_record1=cost_record1, cost_record2=cost_record2))
}


unbiase_generate_coupled_withtheta <- function(param, mod_param, option, option_unbiased, Yt, Lmax, L_prob, mytheta){
  # use coupled samples
  B_burn = param$B_burn  # buruning size of N
  cons = param$cons
  if (mod_param$g_func_type == 1){  # g function in paper
    g_func = g_func_pap_c
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    g_func = g_func_inv_c
  }
  var_pow = mod_param$var_pow  # rate of variance
  
  time_imc_1 <- Sys.time()
  # step 1, sample L and N
  l_upper = sample(c(0 : Lmax), size=1, prob=L_prob)
  if(option$unbiased_index == 3){
    N_max = min(max(var_pow*(Lmax - l_upper)-B_burn-cons, 0), 15)  # Nmax = beta(Lmax-l_i)-p0-C
  }else if(option$unbiased_index == 4){  # if run infinity Lmax
    N_max = option_unbiased$unbia_smc_Pmax
  }
  P_prob = rep(0, N_max + 1)  # probability of p_i
  for (p in 0:N_max){
    if (p < 4){  # $P_P(p) \propto 2^{-p+4}$ for p <4
      P_prob[p + 1] = 2^(-p+4)
    }else{  # $P_P(p) \propto 2^(-p)*p*[log2(p)]^{2}$ otherwise
      P_prob[p + 1] = 2^(-p)*p*(log2(p))^2
    }
  }
  P_prob = P_prob / sum(P_prob)
  N_num = sample(c(0 : N_max), size=1, prob=P_prob)
  fbar = c()
  for(k in 0:N_max){
    fbar[k+1] = sum(P_prob[c((k+1):(N_max+1))])
  }
  
  est = 0
  # run smc with l level=l_upper, and number of samples N0=2^(0+B_burn)
  result_record = vector("list", N_num+1)
  result_record[[1]] = mlsmc_body_withtheta_noparallel(param, mod_param, option, Yt, l_upper, 2^(B_burn), FALSE, 0.01, mytheta)
  # calculate estimate
  if (l_upper > 0){
    # est = eta^{N_p}(g^{l}G^{l-1})/eta^{N_p}(G^{l-1})-eta^{N_p}(g^{l-1})
    samples_l = result_record[[1]]$smc_sample[[l_upper]]
    weights_l = result_record[[1]]$smc_weights[[l_upper]]
    est_N0 = mlsmc_est_withtheta(samples_l, weights_l, g_func, param, mod_param, l_upper, Yt, mytheta)
  }else{
    # est =eta^{N_p}(g^{0})
    samples_l = result_record[[1]]$smc_sample[[1]]
    weights_l = result_record[[1]]$smc_weights[[1]]
    est_N0 = mlsmc_est_0_withtheta(samples_l, weights_l, g_func, param, mod_param, 0, Yt, mytheta)
  }
  est = est_N0 / fbar[1]
  # record cost
  if(l_upper==0){
    cost_record2 = sum(result_record[[1]]$number_of_sample[1]) 
  }else{
    cost_record2 = 0
    for(jj in 1:l_upper){
      cost_record2 = cost_record2 + result_record[[1]]$number_of_sample[jj]*2^(jj-1)
    }
  }
  
  if(N_num > 0){
    est_Nl_1 = est_N0
    for(pp in 1:N_num){
      result_record[[pp+1]] = mlsmc_body_withtheta_noparallel(param, mod_param, option, Yt, l_upper, 2^(B_burn+pp-1), FALSE, 0.01, mytheta)
      samples_l_1 = samples_l
      weights_l_1 = weights_l
      if (l_upper > 0){
        samples_l = c(samples_l_1, result_record[[pp+1]]$smc_sample[[l_upper]])
        weights_l = c(weights_l_1, result_record[[pp+1]]$smc_weights[[l_upper]])
        est_Nl = mlsmc_est_withtheta(samples_l, weights_l, g_func, param, mod_param, l_upper, Yt, mytheta)
      }else{
        samples_l = c(samples_l_1, result_record[[pp+1]]$smc_sample[[1]])
        weights_l = c(weights_l_1, result_record[[pp+1]]$smc_weights[[1]])
        est_Nl = mlsmc_est_0_withtheta(samples_l, weights_l, g_func, param, mod_param, 0, Yt, mytheta)
      }
      est = est + (est_Nl - est_Nl_1) / fbar[pp+1]
      est_Nl_1 = est_Nl
      # record cost
      if(l_upper==0){
        cost_record2 = cost_record2 + sum(result_record[[pp+1]]$number_of_sample[1]) 
      }else{
        for(jj in 1:l_upper){
          cost_record2 = cost_record2 + result_record[[pp+1]]$number_of_sample[jj]*2^(jj-1)
        }
      }
    }
  }
  
  unbiase_est = est / L_prob[l_upper + 1]
  
  time_imc_2 <- Sys.time()
  time_diff = time_imc_2 - time_imc_1
  cost_record1 = time_length(time_diff, unit='second')  # record cost
  
  return(list(est=unbiase_est, L_record=l_upper, N_record=N_num, cost_record1=cost_record1, cost_record2=cost_record2))
}

