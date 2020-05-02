compare_mlsmc_smc2 <- function(param, mod_param, option){
  # compare mlsmc vs smc
  # get work space
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split = "/", fixed = T)) 
  iMC = as.numeric(workdir1[length(workdir1)])  # number of MCs
  workdir2 = ''
  for (i in 1 : (length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC') 

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
  mytheta = mod_param$theta_0

  coress = param$cores  # parallel computing 
  cl <- makeCluster(coress)
  registerDoParallel(cl)

  # mlsmc
  if (option$run_compare_mlsmc){
    epsilon_ind = param$epsilon_mlsmc  # epsilon when calculate number of samples
    mlsmc_l_upper = option$mlsmc_l_upper
    time_imc_1 <- Sys.time()
    imc_temp = mlsmc_body_withtheta(param, mod_param, option, Yt, mlsmc_l_upper, 1, TRUE, 
    epsilon_ind[mlsmc_l_upper], mytheta)
    results_sample = imc_temp$smc_sample
    results_weights = imc_temp$smc_weights
    # calculate estimate
    samples_0 = results_sample[[1]]
    weights_0 = results_weights[[1]]
    yhat_est = mlsmc_est_smc_0_withtheta(samples_0, weights_0, g_func, param, mod_param, 1, Yt, mytheta)
    if (mlsmc_l_upper > 1){
      for(ll in 2:mlsmc_l_upper){
        samples_ll = results_sample[[ll]]
        weights_ll = results_weights[[ll]]
        yhat_est = yhat_est + mlsmc_est_withtheta(samples_ll, weights_ll, g_func, param, 
          mod_param, ll, Yt, mytheta)
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

  write.table(mlsmc_record, file=paste0(workdir2, '/mlsmc_level_', mlsmc_l_upper, 'imc', iMC, '.txt'))
  # write.table(cost_mlsmc_record, file=paste0(workdir2, '/cost_mlsmc_level_', 
       # mlsmc_l_upper, 'imc', iMC, '.txt'))
  write.table(cost_record2, file=paste0(workdir2, '/cost2_mlsmc_level_', mlsmc_l_upper, 'imc', iMC, '.txt'))
  }
	
  # smc
  if (option$run_compare_smc){
    var_pow = mod_param$var_pow  # variance rate
    smc_l_upper = option$smc_l_upper  # l upper when do SMC
    N0 = read.table(paste0(workdir2, '/N0.txt' ))[1,1]
    N_smc_test = ceiling(N0 * 2^(var_pow*smc_l_upper))  # number of samples when do SMC
    time_imc_1 <- Sys.time()
    # imc_temp = mlsmc_body(param , mod_param, option, Yt, smc_l_upper+1, N_smc, FALSE, 0.1)
    imc_temp = mlsmc_body_smc_withtheta(param , mod_param, option, Yt, smc_l_upper+1, N_smc_test, FALSE, 0.1, mytheta)
    results_sample = imc_temp$smc_sample
    results_weights = imc_temp$smc_weights
    # samples_l = results_sample[[smc_l_upper+1]]
    # weights_l = results_weights[[smc_l_upper+1]]
    samples_l = results_sample[[1]]
    weights_l = results_weights[[1]]
    smc_record = mlsmc_est_smc_0_withtheta(samples_l, weights_l, g_func, param, mod_param, 
    smc_l_upper, Yt, mytheta)  						 
    time_imc_2 <- Sys.time()
    time_diff = time_imc_2 - time_imc_1
    cost_smc_record = time_length(time_diff, unit='second')   
    write.table(smc_record , file=paste0(workdir2, '/smc_level_', smc_l_upper, 'imc', iMC, '.txt' ) )
    write.table(cost_smc_record , file=paste0(workdir2, '/cost_smc_level_', smc_l_upper, 'imc', iMC, '.txt'))
  }
  stopCluster( cl )
}


find_smc_N0 <- function(param, mod_param, option, smc_l_upper){
  # compare mlsmc vs smc
  # get work space
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split = "/", fixed = T)) 
  iMC = as.numeric(workdir1[length(workdir1)])  # number of MCs
  workdir2 = ''
  for (i in 1 : (length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC') 
    
  # generate or load data depending on options['newData']
  if (option$newData){
    Yt = generate_data(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
	
  if (mod_param$g_func_type == 1){  # g function in paper
    g_func = g_func_pap_c
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))  # truevalue
    # load mse of mlsmc
    # mse_of_mlsmc = read.table(paste0(workdir2, '/mse_of_mlsmc_paper.txt'))[, 1]  
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    g_func = g_func_inv_c  
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))  # truevalue
    # load mse of mlsmc
    # mse_of_mlsmc = read.table(paste0(workdir2, '/mse_of_mlsmc_inverse.txt'))[, 1]  
  }
  truevalue_est = mean(truevalue[, 1]) 
  mytheta = mod_param$theta_0

  coress = param$cores  # parallel computing 
  cl <- makeCluster(coress)
  registerDoParallel(cl)

  # smc 
  var_pow = mod_param$var_pow  # variance rate
  N_smc_test = 5000  # number of samples when do SMC
  # N_smc_test = N_smc_test_1
  mse_smc_record = rep(0, 50)
  for(j in 1:50){
    imc_temp = mlsmc_body_smc_withtheta(param , mod_param, option, Yt, smc_l_upper+1, N_smc_test, FALSE, 0.1, mytheta)
    results_sample = imc_temp$smc_sample
    results_weights = imc_temp$smc_weights
    samples_l = results_sample[[1]]
    weights_l = results_weights[[1]]
    smc_record = mlsmc_est_0_withtheta(samples_l, weights_l, g_func, param, mod_param, smc_l_upper, Yt, mytheta)
    mse_smc_record[j] = (smc_record - truevalue_est)^2
  }
  N_smc_test_1 = mean(mse_smc_record) / mse_of_mlsmc[smc_l_upper] * N_smc_test
  N0 = N_smc_test_1 / 2^(var_pow*smc_l_upper)
  
  write.table(N0 , file=paste0(workdir2, '/N0_level', smc_l_upper ,'.txt' ) )
  stopCluster( cl )
}


find_mlsmc_N0 <- function(param, mod_param, option, mlsmc_l_upper){
  # compare mlsmc vs smc
  # get work space
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split = "/", fixed = T)) 
  iMC = as.numeric(workdir1[length(workdir1)])  # number of MCs
  workdir2 = ''
  for (i in 1 : (length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC') 
  
  # generate or load data depending on options['newData']
  if (option$newData){
    Yt = generate_data(param, mod_param)
    write.table(Yt, paste0(workdir2, '/observation.txt'))
  }else{
    Yt = read.table(paste0(workdir2, '/observation.txt'))[, 1]
  }
	
  if (mod_param$g_func_type == 1){  # g function in paper
    g_func = g_func_pap_c
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))  # truevalue
    # load mse of mlsmc
    # mse_of_mlsmc = read.table(paste0(workdir2, '/mse_of_mlsmc_paper.txt'))[, 1]  
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    g_func = g_func_inv_c  
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))  # truevalue
    # load mse of mlsmc
    # mse_of_mlsmc = read.table(paste0(workdir2, '/mse_of_mlsmc_inverse.txt'))[, 1]  
  }
  truevalue_est = mean(truevalue[, 1]) 
  mytheta = mod_param$theta_0

  coress = param$cores  # parallel computing 
  cl <- makeCluster(coress)
  registerDoParallel(cl)
    
  # mlsmc 
  var_pow = mod_param$var_pow  # variance rate
  mse_of_mlsmc_l = 2^(-5)*2^(-var_pow*mlsmc_l_upper)
  N_mlsmc_test = 5000  # number of samples when do mlsmc
  # N_mlsmc_test = N_mlsmc_test_1
  mse_mlsmc_record = rep(0, 50)
  for(j in 1:50){
    imc_temp = ml_smc(param , mod_param, option, Yt, mlsmc_l_upper, N_mlsmc_test, FALSE, 0.1)
    mlsmc_record = imc_temp$trueVal
    mse_mlsmc_record[j] = (mlsmc_record - truevalue_est)^2
  }
  N_mlsmc_test_1 = mean(mse_mlsmc_record) / mse_of_mlsmc_l* N_mlsmc_test
  N0 = N_mlsmc_test_1 / 2^(var_pow*mlsmc_l_upper)
  
  write.table(N0 , file=paste0(workdir2, '/N0_level', mlsmc_l_upper ,'.txt' ) )
  stopCluster( cl )
}