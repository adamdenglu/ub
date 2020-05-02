# Unbiased Log Gradients for Inverse Problems
# Deng Lu
# 10,September,2019
# visulize results


visulize_interp2 <- function(param, mod_param, option, ind, n_interp, Lmax){
  # rm(list=ls())
  ind = c(1:10)
  n_interp = 10000
  Lmax = 1
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split="/", fixed=T)) 
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
  
  # load true value
  # if (mod_param$g_func_type == 1){  # g function in paper
  # truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  # }else if (mod_param$g_func_type == 2){  # g function in inverse problem
  # truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  # }
  # truevalue_est = mean(truevalue[, 1]) 
  
  cost_rec = list()
  mse_rec = list()
  cmin_temp = c()
  cmax_temp = c()
  for (i in 1:length(ind)){
    # temp_est = read.table(paste0(workdir2, '/Lmax', Lmax, 'est_run_', ind[i], '.txt'))
    temp_mse = read.table(paste0(workdir2, '/Lmax', Lmax, 'mse_run_', ind[i], '.txt'))
    temp_cost = read.table(paste0(workdir2, '/Lmax', Lmax, 'cost2_run_', ind[i], '.txt'))
    # temp_est = as.numeric(temp_est[, 1]) 
    temp_mse = as.numeric(temp_mse[, 1]) 
    temp_cost = as.numeric(temp_cost[, 1])
    
    cost_rec[[i]] = temp_cost
    mse_rec[[i]] = temp_mse
    
    cmin_temp = c(cmin_temp, cost_rec[[i]][1])
    cmax_temp = c(cmax_temp, cost_rec[[i]][length(temp_cost)])
  }
  
  # Interpolation
  mse_interp = list()
  cmin = max(cmin_temp)
  cmax = min(cmax_temp)
  mse_interp_mean = rep(0, n_interp)
  c_interp = seq(cmin, cmax, length.out=n_interp)
  for (i in 1:length(ind)){        
    temp = approx(cost_rec[[ i ]], y = mse_rec[[i]],  xout=c_interp)
    mse_interp[[i]] = temp[[2]]
    mse_interp_mean = mse_interp_mean + temp[[2]]
  }
  mse_interp_mean = mse_interp_mean / (length(ind))
  
  # # mlsmc
  # temp_est = read.table(paste0(workdir2, '/mlsmc_est_record.txt'))
  # temp_cost = read.table( paste0(workdir2, '/mlsmc_cost_record.txt'))
  # cost_mlsmc =  colMeans(temp_cost)
  # mse_mlsmc = colMeans((temp_est - truevalue_est)^2)
  
  pdf('mse vs cost_unbiased.pdf', width=10, height=10)
  par(mar=c(5, 6, 4, 2) )
  x_lab = seq(4, 26, by=2)  
  x_lab_1 <- c(expression(2^4), expression(2^6), expression(2^8), 
               expression(2^10), expression(2^12), expression(2^14), 
               expression(2^16), expression(2^18), expression(2^20), 
               expression(2^22), expression(2^24), expression(2^26))
  y_lab = seq(-16, 2, by=2) 
  y_lab_1 <- c(expression(2^-16),expression(2^-14), expression(2^-12), expression(2^-10), 
               expression(2^-8), expression(2^-6), expression(2^-4),
               expression(2^-2), expression(2^0), expression(2^2)) 
  
  plot(log2(c_interp), log2(mse_interp_mean), type="l", col="#E9967A",
       xlab='Cost', ylab='MSE \n', lwd=2, xaxt="n", yaxt='n', 
       cex.lab=1.5, xlim =c(x_lab[1], x_lab[length(x_lab)]), 
       ylim=c(y_lab[1], y_lab[length(y_lab)]))
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2)
  abline(v=x_lab, col="gray", lty="dotted", lwd=2)
  abline(h=y_lab, col="gray", lty="dotted", lwd=2)
  # lines(log2(mse_mlsmc), log2(cost_mlsmc), type="b", col="springgreen3", lwd=2)
  # legend('topright', inset=0.05, c('unbiased', 'mlsmc'), title='Type',
  # col=c("#E9967A", "springgreen3"), lt =c(1, 1),
  # lwd=c(2, 2))  
  dev.off()
  
  m1 <- lm(log2(mse_interp_mean) ~ log2(c_interp))
}


combind_unbiased <- function(param, mod_param, option){
  # rm(list=ls())
  # load data
  mse_lmax1 = read.table('unbiased_mse_Lmax1_intep.txt')[, 1]
  cost_lmax1 = read.table('unbiased_cost_Lmax1_intep.txt')[, 1]
  mse_lmax2 = read.table('unbiased_mse_Lmax2_intep.txt')[, 1]
  cost_lmax2 = read.table('unbiased_cost_Lmax2_intep.txt')[, 1]
  mse_lmax3 = read.table('unbiased_mse_Lmax3_intep.txt')[, 1]
  cost_lmax3 = read.table('unbiased_cost_Lmax3_intep.txt')[, 1]
  
  require("RColorBrewer")
  mycolor = brewer.pal(12, "Set2")
  pdf('mse vs cost of unbiased estimator.pdf', width=10, height=10)
  par(mar=c(5, 6, 4, 2) )
  x_lab = seq(2, 30, by=2)  
  x_lab_1 <-  c(expression(2^2), expression(2^4), expression(2^6),
                expression(2^8), expression(2^10), expression(2^12), 
                expression(2^14), expression(2^16), expression(2^18), 
                expression(2^20), expression(2^22), expression(2^24), 
                expression(2^26), expression(2^28), expression(2^30))
  y_lab = seq(-20, 4, by=2) 
  y_lab_1 <-  c(expression(2^-20), expression(2^-18),expression(2^-16),
                expression(2^-14), expression(2^-12), expression(2^-10), 
                expression(2^-8), expression(2^-6), expression(2^-4),
                expression(2^-2), expression(2^0), expression(2^2),
                expression(2^4)) 
  
  plot(cost_lmax1, mse_lmax1, type="l", col=mycolor[1],
       xlab='Cost', ylab='MSE \n', lwd=2, xaxt="n", yaxt='n', 
       cex.lab=1.5, xlim =c(x_lab[1], x_lab[length(x_lab)]), 
       ylim=c(y_lab[1], y_lab[length(y_lab)]))
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2)
  abline(v=x_lab, col="gray", lty="dotted", lwd=2)
  abline(h=y_lab, col="gray", lty="dotted", lwd=2)
  lines(cost_lmax2, mse_lmax2, type="l", col=mycolor[2], lwd=2)
  lines(cost_lmax3, mse_lmax3, type="l", col=mycolor[3], lwd=2)
  text1 = c('unbiased estimator with Lmax=3', 'unbiased estimator with Lmax=4', 
            'unbiased estimator with Lmax=5')
  legend('topright', inset=0.05, text1, 
        col=mycolor[1:3], lty =c(1, 1, 1), lwd=c(2, 2, 2))
  dev.off()
}


visulize_interp_mlsmc <- function(param, mod_param, option){
  # rm(list=ls())
  
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  }
  truevalue_est = mean(truevalue[, 1]) 
  
  # mlsmc
  temp_est = read.table('mlsmc_est_record.txt')
  temp_cost = read.table('mlsmc_cost_record.txt')
  cost_mlsmc =  colMeans(temp_cost)
  mse_mlsmc = colMeans((temp_est - truevalue_est)^2)
 
  # unbiased
  mse_lmax1 = read.table('unbiased_mse_Lmax1_intep.txt')[, 1]
  cost_lmax1 = read.table('unbiased_cost_Lmax1_intep.txt')[, 1]
  mse_lmax2 = read.table('unbiased_mse_Lmax2_intep.txt')[, 1]
  cost_lmax2 = read.table('unbiased_cost_Lmax2_intep.txt')[, 1]
  mse_lmax3 = read.table('unbiased_mse_Lmax3_intep.txt')[, 1]
  cost_lmax3 = read.table('unbiased_cost_Lmax3_intep.txt')[, 1]
  
  cost_unbiased = c()
  ind1 = which(mse_lmax1 < log2(mse_mlsmc)[1])[1]
  cost_unbiased[1] = cost_lmax1[ind1]
  ind2 = which(mse_lmax2 < log2(mse_mlsmc)[2])[1]
  cost_unbiased[2] = cost_lmax2[ind2]
  ind3 = which(mse_lmax3 < log2(mse_mlsmc)[3])[1]
  cost_unbiased[3] = cost_lmax3[ind3]
  
  pdf('mse vs cost.pdf', width=10, height=10)
  par(mar=c(5, 6, 4, 2) )
  y_lab = seq(12, 26, by=2)  
  y_lab_1 <-  c(expression(2^12), expression(2^14), expression(2^16), 
                expression(2^18), expression(2^20), expression(2^22),
                expression(2^24), expression(2^26))
  x_lab = seq(-18, -8, by=2) 
  x_lab_1 <-  c(expression(2^-18), expression(2^-16), expression(2^-14), 
                expression(2^-12), expression(2^-10), expression(2^-8)) 
  
  plot(log2(mse_mlsmc)[1:4], log2(cost_mlsmc)[1:4], type="b", col="#E9967A",
       xlab='MSE', ylab='Cost \n', lwd=2, xaxt="n", yaxt='n', 
       cex.lab=1.5, xlim =c(x_lab[1], x_lab[length(x_lab)]), 
       ylim=c(y_lab[1], y_lab[length(y_lab)]))
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2)
  abline(v=x_lab, col="gray", lty="dotted", lwd=2)
  abline(h=y_lab, col="gray", lty="dotted", lwd=2)
  lines(log2(mse_mlsmc)[c(1:length(cost_unbiased))], cost_unbiased, type="b", col="springgreen3", lwd=2)
  legend('topright', inset=0.05, c('MLSMC', 'unbiased algorithm'), title='Type',
  col=c("#E9967A", "springgreen3"), lt =c(1, 1),
  lwd=c(2, 2))
  dev.off()
  
  2^(cost_unbiased[1]) / 2^(log2(cost_mlsmc)[1])  # 2.343045
  2^(cost_unbiased[2]) / 2^(log2(cost_mlsmc)[2])  # 7.566411
  2^(cost_unbiased[3]) / 2^(log2(cost_mlsmc)[3])  # 66.57725
  
  lm(log2(cost_mlsmc)[1:4] ~ log2(mse_mlsmc)[1:4])
  # Call:
  #   lm(formula = log2(cost_mlsmc)[1:4] ~ log2(est_mlsmc)[1:4])
  # 
  # Coefficients:
  #   (Intercept)  log2(est_mlsmc)[1:4]  
  # 4.6941               -0.9097  
  lm(cost_unbiased ~ log2(mse_mlsmc)[c(1:length(cost_unbiased))])
  # Call:
  #   lm(formula = cost_unbiased ~ log2(mse_mlsmc)[c(1:length(cost_unbiased))])
  # 
  # Coefficients:
  #   (Intercept)  log2(mse_mlsmc)[c(1:length(cost_unbiased))]  
  # -2.390
}


visulize_cp_mse_cost_2 <- function(param , mod_param , option, ind_mlsmc, ind_smc){
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split="/", fixed=T)) 
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
  
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  }
  truevalue_est = mean(truevalue[, 1]) 
  
  imc_num = 50
  
  # results of mlsmc
  cost_rec = list()
  est_rec = list()
  cost_mlsmc = c()
  est_mlsmc = c()
  for (i in 1:length(ind_mlsmc)){
    cost_rec[[ i ]] = rep(0, imc_num)
    est_rec[[ i ]] = rep(0, imc_num)
    for (j in 1:imc_num ){
      temp_est = read.table(paste0(workdir2, '/mlsmc_level_', ind_mlsmc[i], 'imc', j, 
                                   '.txt'))
      temp_cost = read.table(paste0(workdir2, '/cost_mlsmc_level_', ind_mlsmc[i], 'imc', j, 
                                    '.txt'))
      cost_rec[[ i ]][j] = as.numeric(temp_cost[, 1])
      est_rec[[ i ]][j] = as.numeric(temp_est[, 1])
    }        
    cost_mlsmc[ i ] = mean(cost_rec[[i]])
    est_mlsmc[ i ] = mean((est_rec[[i]] - truevalue_est)^2)
  }
  # results of smc
  cost_rec_smc = list()
  est_rec_smc = list()
  cost_smc = c()
  est_smc = c()
  for (i in 1:length(ind_smc)){
    cost_rec_smc[[ i ]] = rep(0, imc_num)
    est_rec_smc[[ i ]] = rep(0, imc_num)
    for (j in 1:imc_num ){
      temp_est = read.table(paste0(workdir2, '/smc_level_', ind_smc[i], 'imc', j, '.txt'))
      temp_cost = read.table(paste0(workdir2, '/cost_smc_level_', ind_smc[i], 'imc', j, 
                                    '.txt'))
      cost_rec_smc[[ i ]][j] = as.numeric(temp_cost[, 1])
      est_rec_smc[[ i ]][j] = as.numeric(temp_est[, 1])
    }        
    cost_smc[ i ] = mean(cost_rec_smc[[i]])
    est_smc[ i ] = mean((est_rec_smc[[i]] - truevalue_est)^2)
  }
  
  pdf('compare mse vs cost.pdf',width=10,height=10)
  par(mar=c(5, 6, 4, 2) )
  x_lab = seq(-21, -9, by=3)  
  x_lab_1 <-  c(expression(2^-21), expression(2^-18), expression(2^-15), 
                expression(2^-12), expression(2^-9))
  y_lab = seq(0, 14, by=2) 
  y_lab_1 <-  c(expression(2^0),expression(2^2), expression(2^4), expression(2^6), 
                expression(2^8), expression(2^10), expression(2^12), 
                expression(2^14)) 
  plot(log2(est_mlsmc), log2(cost_mlsmc), type="b" , col="#E9967A", cex.lab=1.5,
       xlab='MSE', ylab='Cost (s)', lwd=2, xaxt="n", yaxt='n', 
       xlim=c(x_lab[1], x_lab[length(x_lab)]), ylim=c(y_lab[1], y_lab[length(y_lab)]))
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2)      
  lines(log2(est_smc), log2(cost_smc),  type="b", col="springgreen3", lwd=2)
  legend('topright', inset=0.05, c('mlsmc', 'smc') , title='Type',
         col=c("#E9967A", "springgreen3"), lty=c(1, 1),
         lwd=c(2, 2))  
  dev.off() 
  
  m1 <- lm(log2(cost_mlsmc) ~ log2(est_mlsmc))
  # Call:
  # lm(formula = log2(cost_mlsmc) ~ log2(est_mlsmc))
  
  # Coefficients:
  # (Intercept)  log2(est_mlsmc)  
  # -10.255           -1.034 
  m2 <- lm(log2(cost_smc) ~ log2(est_smc))
  # Call:
  # lm(formula = log2(cost_smc) ~ log2(est_smc))
  
  # Coefficients:
  # (Intercept)  log2(est_smc)  
  # -10.270         -1.122 
}


visulize_check_rate <- function(param, mod_param, option){
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split="/", fixed=T)) 
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
  
  est_rec = matrix(0, nrow=50, ncol=7)
  for(i in 1:nrow(est_rec)){
    est_rec[i, ] = read.table(paste0(workdir2, '/checkrate_iMC', i, '.txt'))[, 1]
  }
  var_rec = c()
  for(i in 1:ncol(est_rec)){
    var_rec = c(var_rec, var(est_rec[, i]))
  }
  small_k = param$small_k
  l_ind = c(0:-6) - small_k
  
  pdf('mlsmc variance rate estimate.pdf', width=10, height=10)
  par(mar=c(5, 6, 4, 2))
  x_lab = seq(-9, -3, by=1)  
  x_lab_1 <- c(expression(2^-9), expression(2^-8), expression(2^-7), 
               expression(2^-6), expression(2^-5), expression(2^-4),
               expression(2^-3))
  y_lab = seq(-38, -20, by=3) 
  y_lab_1 <- c(expression(2^-38), expression(2^-35), expression(2^-32), 
               expression(2^-29), expression(2^-26), expression(2^-23), 
               expression(2^-20)) 
  plot(l_ind, log2(var_rec), type="b", col="#E9967A", xaxt="n", yaxt='n', 
       xlab ='h_l', ylab='Variance \n', lwd=2, cex.lab=1.5)
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2)
  dev.off()
  
  l_upper = 7
  epsilon_ind = param$epsilon_mlsmc
  eps = epsilon_ind[l_upper]
  xNofSamp = rep(0, l_upper + 1)
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
  var_rec1 = var_rec * xNofSamp[1:l_upper]
  lm1 = lm(log2(var_rec1) ~ l_ind)
  beta_est = round(lm1$coefficients[2], 3)
  pdf('mlsmc variance rate estimate2.pdf', width=10, height=10)
  par(mar = c(8, 8, 4, 2), mgp=c(5, 1, 0))
  x_lab = seq(-9, -3, by=1)  
  x_lab_1 <- c(expression(2^-9), expression(2^-8), expression(2^-7), 
               expression(2^-6), expression(2^-5), expression(2^-4),
               expression(2^-3))
  y_lab = seq(-27, 0, by=3) 
  y_lab_1 <- c(expression(2^-27), expression(2^-24), expression(2^-21), 
               expression(2^-18), expression(2^-15), expression(2^-12), 
               expression(2^-9), expression(2^-6), expression(2^-3),
               expression(2^0)) 
  myxlab = expression(h[l]) 
  myylab = expression(paste(V[l]%~~%N[l],'var[', eta[l]^N[l], '(g', G[l], ')/',
                            eta[l]^N[l], '(', G[l], ')-', eta[l], '(g)]'))
  plot(l_ind, log2(var_rec1), type="b", col="#E9967A", xaxt="n", yaxt='n', 
       xlab=myxlab, ylab=myylab, lwd=2, cex.lab=1.5)
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2, mgp=c(5, 1, 0))
  text(-6, -18, expression(paste(beta, ' = 3.97')), cex=2)
  dev.off()
}


combine_result <- function(param, mod_param, option, ind_mlsmc){
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split="/", fixed=T)) 
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
  
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  }
  truevalue_est = mean(truevalue[, 1]) 
  
  imc_num = 50
  
  # results of mlsmc
  cost_rec = matrix(0, nrow=imc_num, ncol=length(ind_mlsmc))
  est_rec = matrix(0, nrow=imc_num, ncol=length(ind_mlsmc))
  for (i in 1:length(ind_mlsmc)){
    for (j in 1:imc_num ){
      temp_est = read.table(paste0(workdir2, '/mlsmc_level_', ind_mlsmc[i], 'imc', j, 
                                   '.txt'))
      temp_cost = read.table(paste0(workdir2, '/cost_mlsmc_level_', ind_mlsmc[i], 'imc', j, 
                                    '.txt'))
      cost_rec[j, i] = as.numeric(temp_cost[, 1])
      est_rec[j, i] = as.numeric(temp_est[, 1])
    }        
  }
  write.table(cost_rec , file='mlsmc_cost_record.txt' )
  write.table(est_rec , file='mlsmc_est_record.txt' )
}


visulize_test <- function(algo_result_test, param, mod_param, option){
  # visulize test results: mlsmc vs smc
  yhat_smc = algo_result_test$yhat_smc
  yhat_mlsmc_record = algo_result_test$yhat_mlsmc_record
  N_ind = algo_result_test$N_ind
  
  mse = (yhat_mlsmc_record - yhat_smc)^2
  pdf( 'mlsmc vs smc.pdf' ,width=10, height=10 )
  plot(N_ind, mse, type = "l")
  dev.off()
}


visulize_unbiased <- function(param, mod_param, option, est){
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split="/",fixed=T)) 
  # iMC = as.numeric( workdir1[ length( workdir1 ) ] )
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
     workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
  
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  }
  truevalue_est = mean(truevalue[, 1]) 
  
  est_record = c()
  for (i in 1:length(est)){
    est_record[i] = mean(est[1:i])
  }
  
  pdf('unbiased estimate.pdf', width=10, height=10)
  plot(c(1:length(est)), est_record, type="l", col = "#E9967A",
       xlab='sample size', ylab='estiamte', lwd=2)
  abline(h=truevalue_est, col="#00BFFF", lwd=2, lty=2)
  dev.off()
}


visulize_mlsmc_mse_cost <- function(param, mod_param, option, ind){
  workdir = getwd()
  workdir1 = unlist( strsplit(workdir, split="/", fixed=T)) 
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
  
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  }
  truevalue_est = mean(truevalue[, 1]) 
  
  cost_rec = list()
  est_rec = list()
  cost_mlsmc = c()
  est_mlsmc = c()
  for (i in 1:length(ind)){
    temp_est = read.table(paste0(workdir2, '/mlsmc_level_', ind[i], '.txt'))
    temp_cost = read.table(paste0(workdir2, '/cost_level_', ind[i], '.txt'))
    cost_rec[[i]] = as.numeric(temp_cost[, 1]) 
    est_rec[[i]] = as.numeric(temp_est[, 1]) 
    cost_mlsmc[i] = mean(cost_rec[[i]])
    est_mlsmc[i] = mean((est_rec[[i]] - truevalue_est)^2)
  }
    
  pdf('mlsmc mse vs cost.pdf', width=10, height=10)
  plot(est_mlsmc, cost_mlsmc, type="l", col="#E9967A",
       xlab='MSE', ylab='Cost (s)', lwd=2)
  dev.off()
}



visulize_interp <- function(param, mod_param, option, ind, n_interp, Lmax){
  workdir = getwd()
  workdir1 = unlist(strsplit(workdir, split="/", fixed=T)) 
  workdir2 = ''
  for (i in 1:(length(workdir1) - 1)){
    workdir2 = paste0(workdir2,  workdir1[i], '/') 
  }
  workdir2 = paste0(workdir2, 'MC')
    
  # load true value
  if (mod_param$g_func_type == 1){  # g function in paper
    truevalue = read.table(paste0(workdir2, '/trueValue_of_paper.txt'))
  }else if (mod_param$g_func_type == 2){  # g function in inverse problem
    truevalue = read.table(paste0(workdir2, '/trueValue_of_inverse.txt'))
  }
  truevalue_est = mean(truevalue[, 1]) 
    
  library(doParallel)
  coress = 4  # parallel computing  
  cl <- makeCluster(coress)
  registerDoParallel(cl)
  
  cost_rec = list()
  est_rec = list()
  mse_rec = list()
  cmin_temp = c()
  cmax_temp = c()
  for (i in 1:length(ind)){
    temp_est = read.table(paste0(workdir2, '/Lmax', Lmax, 'est_run_', ind[i], '.txt'))
    temp_cost = read.table(paste0(workdir2, '/Lmax', Lmax, 'cost2_run_', ind[i], '.txt'))
    temp_est = as.numeric(temp_est[, 1]) 
    temp_cost = as.numeric(temp_cost[, 1])
    cost_i = c()
    est_i = c()
    mse_i = c()
    
    cost_i <- foreach(x=1:length(temp_cost), .packages=c("smccbase7"), .combine='c') %dopar% {sum(temp_cost[1:x])}
    est_i <- foreach(x=1:length(temp_cost), .packages=c("smccbase7"), .combine='c') %dopar% {mean(temp_est[1:x])}
    mse_i <- (est_i - truevalue_est)^2
    
    # for (j in 1:length(temp_cost)){
        # cost_i[j] = sum(temp_cost[1:j])
        # est_i[j] = mean(temp_est[1:j])
        # mse_i[j] = (est_i[j] - truevalue_est)^2
    # }
    cost_rec[[i]] = cost_i
    est_rec[[i]] = est_i
    mse_rec[[i]] = mse_i
    
    cmin_temp = c(cmin_temp, cost_rec[[i]][1])
    cmax_temp = c(cmax_temp, cost_rec[[i]][length(temp_cost)])
  }
    
  # Interpolation
  mse_interp = list()
  cmin = max(cmin_temp)
  cmax = min(cmax_temp)
  mse_interp_mean = rep(0, n_interp)
  c_interp = seq(cmin, cmax, length.out=n_interp)
  for (i in 1:length(ind)){        
    temp = approx(cost_rec[[ i ]], y = mse_rec[[i]],  xout=c_interp)
    mse_interp[[i]] = temp[[2]]
    mse_interp_mean = mse_interp_mean + temp[[2]]
  }
  mse_interp_mean = mse_interp_mean / (length(ind))
    
  # # mlsmc
  # temp_est = read.table(paste0(workdir2, '/mlsmc_est_record.txt'))
  # temp_cost = read.table( paste0(workdir2, '/mlsmc_cost_record.txt'))
  # cost_mlsmc =  colMeans(temp_cost)
  # mse_mlsmc = colMeans((temp_est - truevalue_est)^2)
    
  pdf('mse vs cost_unbiased.pdf', width=10, height=10)
  par(mar=c(5, 6, 4, 2) )
  x_lab = seq(6, 24, by=2)  
  x_lab_1 <-  c(expression(2^6), expression(2^8), expression(2^10), 
                expression(2^12), expression(2^14), expression(2^16),
                expression(2^18), expression(2^20), expression(2^22),
                expression(2^24))
  y_lab = seq(-14, 2, by=2) 
  y_lab_1 <-  c(expression(2^-14), expression(2^-12), expression(2^-10), 
                expression(2^-8), expression(2^-6), expression(2^-4),
                expression(2^-2), expression(2^0), expression(2^2)) 
                 
  plot(log2(c_interp), log2(mse_interp_mean), type="l", col="#E9967A",
       xlab='Cost', ylab='MSE \n', lwd=2, xaxt="n", yaxt='n', 
       cex.lab=1.5, xlim =c(x_lab[1], x_lab[length(x_lab)]), 
       ylim=c(y_lab[1], y_lab[length(y_lab)]))
  axis(side=1, at=x_lab, labels=x_lab_1, cex.axis=1.5)
  axis(side=2, at=y_lab, labels=y_lab_1, cex.axis=1.5, las=2)
  abline(v=seq(6,24,2), col="gray", lty="dotted", lwd=2)
  abline(h=seq(-14, 2, by=2), col="gray", lty="dotted", lwd=2)
  # lines(log2(mse_mlsmc), log2(cost_mlsmc), type="b", col="springgreen3", lwd=2)
  # legend('topright', inset=0.05, c('unbiased', 'mlsmc'), title='Type',
          # col=c("#E9967A", "springgreen3"), lt =c(1, 1),
          # lwd=c(2, 2))  
  dev.off()
  
  m1 <- lm(log2(mse_interp_mean) ~ log2(c_interp))
  # Call:
  # lm(formula = log2(mse_interp_mean) ~ log2(c_interp))

  # Coefficients:
     # (Intercept)  log2(c_interp)  
           # 9.494          -0.999 
  unbiased_cost=log2(c_interp)[min(which(log2(mse_interp_mean)<log2(mse_mlsmc)[Lmax]))]
  write.table(unbiased_cost, file=paste0('unbiased_cost_Lmax', Lmax, '.txt'))

}





# visulize_cp_mse_cost <- function(param, mod_param, option, ind){
    # workdir = getwd()
    # workdir1 = unlist(strsplit(workdir, split="/", fixed=T)) 
    # workdir2 = ''
    # for (i in 1:(length(workdir1) - 1)){
       # workdir2 = paste0(workdir2,  workdir1[i], '/') 
    # }
    # workdir2 = paste0(workdir2, 'MC')
    
    # truevalue = read.table(paste0(workdir2, '/trueValue.txt'))
    # truevalue_est = mean(truevalue[, 1]) 
    
    # # results of mlsmc
    # cost_rec = list()
    # est_rec = list()
    # cost_mlsmc = c()
    # est_mlsmc = c()
    # for (i in 1:length(ind)){
        # temp_est = read.table(paste0(workdir2, '/mlsmc_level_', ind[i], '.txt'))
        # temp_cost = read.table(paste0(workdir2, '/cost_mlsmc_level_', ind[i], '.txt'))
        # cost_rec[[i]] = as.numeric(temp_cost[, 1]) 
        # est_rec[[i]] = as.numeric(temp_est[, 1]) 
        # cost_mlsmc[i] = mean(cost_rec[[i]])
        # est_mlsmc[i] = mean((est_rec[[i]] - truevalue_est)^2)
    # }
    # # results of smc
    # cost_rec_smc = list()
    # est_rec_smc = list()
    # cost_smc = c()
    # est_smc = c()
    # for ( i in 1 : length( ind ) ){
        # temp_est = read.table( paste0( workdir2 , '/smc_level_' , ind[i] , '.txt' ) )
        # temp_cost = read.table( paste0( workdir2 , '/cost_smc_level_' , ind[i] , '.txt' ) )
        # cost_rec_smc[[ i ]] = as.numeric( temp_cost[ , 1 ] ) 
        # est_rec_smc[[ i ]] = as.numeric( temp_est[ , 1 ] ) 
        # cost_smc[ i ] = mean( cost_rec_smc[[ i ]] )
        # est_smc[ i ] = mean( ( est_rec_smc[[ i ]] - truevalue_est ) ^ 2 )
    # }
    
    # pdf('compare mse vs cost.pdf',width=10,height=10)
    # x_lab = seq( -15, -8 , by = 1 )  
    # x_lab_1 <-  c( expression(2^-15) , expression(2^-14) , expression(2^-13) ,
                   # expression(2^-12) , expression(2^-11) , expression(2^-10) , 
                   # expression(2^-9) , expression(2^-8) )
    # y_lab = seq( 2 , 10 , by = 2 ) 
    # y_lab_1 <-  c( expression(2^2) , expression(2^4) ,  expression(2^6) , 
                   # expression(2^8) , expression(2^10) ) 
    # plot( log2( est_mlsmc ) , log2( cost_mlsmc ), type = "b" , col = "#E9967A" ,
          # xlab = 'MSE' , ylab = 'Cost (s)' , lwd = 2 , xaxt="n" , yaxt = 'n', 
          # xlim = c( -15, -8 ) , ylim = c( 2 , 10 ) )
    # axis( side = 1 , at = x_lab , labels = x_lab_1 , cex.axis = 1.5  )
    # axis( side = 2 , at = y_lab , labels = y_lab_1 , cex.axis = 1.5 , las = 2  )      
    # lines( log2( est_smc ) , log2( cost_smc ) ,  type = "b" , col = "springgreen3" ,lwd=2 )
    # legend( 'topright' , inset = 0.05 , c( 'mlsmc' , 'smc' ) , title = 'Type' ,
            # col = c( "#E9967A" , "springgreen3" ) , lty = c( 1 , 1 ) ,
            # lwd = c( 2,2) )  
    # dev.off() 
# }




