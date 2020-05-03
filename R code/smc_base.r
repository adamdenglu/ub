# Unbiased Log Gradients for Inverse Problems
# Deng Lu
# 7,August,2019
# smc_base algorithm


mlsmc_est_withtheta <- function(samples, weight, g_function, param, mod_param, l_level, y, mytheta){
  # estimate at level l_level: \eta_l(gG)/\eta_l(G) - \eta_l(g)
  # add parameter theta
  big_k = mod_param$big_k
  small_k = mod_param$small_k
  gp_x = mod_param$gp_x

  if(mod_param$g_func_type == 1){  # g function in paper
    x_g = mod_param$x_g
    l_g = mod_param$l_g
    eta_l_g = foreach(x=1:length(weight), .packages=c("smccbase14")) %dopar% {g_function(samples[[x]], big_k, small_k, l_g, x_g)}
    eta_l_g <- do.call('cbind',eta_l_g)
    eta_l_g_1 = eta_l_g
  }else if(mod_param$g_func_type == 2){  # g function in inverse problem
    eta_l_g = rep(0, length(weight))
    eta_l_g_1 = rep(0, length(weight))
    theta_mu = mod_param$theta_mu
    theta_sigma = mod_param$theta_sigma
    for (jjj in 1:length(eta_l_g)){
      eta_l_g[jjj] = g_function(samples[[jjj]], big_k, small_k, gp_x, mytheta, l_level, y, theta_mu, theta_sigma)
      eta_l_g_1[jjj] = g_function(samples[[jjj]], big_k, small_k, gp_x, mytheta, l_level-1, y, theta_mu, theta_sigma)
    }
  }
  
  eta_l_g_G = eta_l_g * exp(weight)
  est_l = mean(eta_l_g_G) / mean(exp(weight)) - mean(eta_l_g_1)
  return (est_l)
}


mlsmc_est_0_withtheta <- function(samples, weight, g_function, param, mod_param, l_level, y, mytheta){
  # estimate at level l_level: \eta_l(g)
  # add parameter theta
  big_k = mod_param$big_k
  small_k = mod_param$small_k
  gp_x = mod_param$gp_x
  x_g = mod_param$x_g
  l_g = mod_param$l_g
  
  if(mod_param$g_func_type == 1){  # g function in paper
    x_g = mod_param$x_g
    l_g = mod_param$l_g
    eta_l_g = foreach(x=1:length(weight), .packages=c("smccbase14")) %dopar% {g_function(samples[[x]], big_k, small_k, l_g, x_g)}
    eta_l_g <- do.call('cbind',eta_l_g)
  }else if(mod_param$g_func_type == 2){  # g function in inverse problem
    theta_mu = mod_param$theta_mu
    theta_sigma = mod_param$theta_sigma
    eta_l_g = rep(0, length(weight))
    for (jjj in 1:length(eta_l_g)){
      eta_l_g[jjj] = g_function(samples[[jjj]], big_k, small_k, gp_x, mytheta, l_level, y, theta_mu, theta_sigma)
    }
  }
  
  return (mean(eta_l_g))
}


mlsmc_est_smc_0_withtheta <- function(samples, weight, g_function, param, mod_param, l_level, y, mytheta){
  # estimate at level 0: \eta_0(gG)/\eta_0(G)
  # add parameter theta
  big_k = mod_param$big_k
  small_k = mod_param$small_k
  gp_x = mod_param$gp_x
  
  if(mod_param$g_func_type == 1){  # g function in paper
    x_g = mod_param$x_g
    l_g = mod_param$l_g
    eta_l_g = foreach(x=1:length(weight), .packages=c("smccbase14")) %dopar% {g_function(samples[[x]], big_k, small_k, l_g, x_g)}
    eta_l_g <- do.call('cbind',eta_l_g)
  }else if(mod_param$g_func_type == 2){  # g function in inverse problem
    theta_mu = mod_param$theta_mu
    theta_sigma = mod_param$theta_sigma
    eta_l_g = rep(0, length(weight))
    for (jjj in 1:length(weight)){
      eta_l_g[jjj] = g_function(samples[[jjj]], big_k, small_k, gp_x, mytheta, l_level, y, theta_mu, theta_sigma)
    }
  }
  
  eta_l_g_G = eta_l_g * exp(weight)
  est_l = mean(eta_l_g_G) / mean(exp(weight))
  return (est_l)
}

