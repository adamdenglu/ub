# Unbiased Log Gradients for Inverse Problems
# Deng Lu
# version 1
# 7,August,2019
# Model base functions
# Poisson Equation in one dimension
# -\delta * ( u_hat * \delta p ) = f


f_func <- function(x){
    # f(x) = 100*x
    return (100 * x) 
}


u_bar <- function(x){
    # \bar{u}(x) = 0.15
    return (0.15)
}


varphi <- function(x , k){
    if (k %% 2 == 0){  # if k is even
        return (cos(k * pi * x))
    }else if (k%%2 == 1){  # if k is odd
        return (sin(k * pi * x))
    }
}


psi_l <- function(x, xi, hl){
    # finite-element basis functions
    if ((xi - hl < x) & (x <= xi)){
        return ((x - xi + hl) / hl)
    }else if ((xi < x) & (x <= xi+hl)){
        return ((xi + hl - x) / hl)
    }else{
        return (0)
    }
}


u_hat <- function(x, u , sigma_k){
    # \hat{u}(x) = \bar{u}(x)+\sum_{k=1}^Ku_k*\sigma_k*\varphi_k(x)
    k_dim = length(u)
    u_hat_sum = u_bar(x) 
    for (ii in 1:k_dim){
        u_hat_sum = u_hat_sum + u[ii] * sigma_k[ii] * varphi(x, ii)
    }
    return (u_hat_sum)
}


inte_u_hat <- function(a, b, u, sigma_k){
    # integration of u_hat on [a,b]
    k_dim = length(u)
    inte = 0.15 * (b - a)
    for (ii in 1:k_dim){
        if (ii%%2 == 1){
            inte = inte - u[ii] * sigma_k[ii] * (cos(ii*pi*b) - cos(ii*pi*a)) / (ii * pi)
        }else{
            inte = inte + u[ii] * sigma_k[ii] * (sin(ii*pi*b) - sin(ii*pi*a)) / (ii * pi)
        }
    }
    return (inte)
}


# pde_solver_base_old <- function( u , sigma_k , l , small_k ){
    # # solve for p^l
    # dim_l = 2 ^ ( l + small_k ) - 1
    # h_l = 2 ^ ( - ( l + small_k ) )
    # xi = h_l * seq( 1 , dim_l )
    
    # # calculate vector f^l
    # f_l = rep( 0 , dim_l )
    # for ( i in 1 : dim_l ){
        # f_l[ i ] = ( 1 / h_l ) * ( 100 / 6 ) * ( ( xi[ i ] + h_l ) ^ 3 + ( xi[ i ] - h_l ) ^ 3 - 2 * xi[ i ] ^ 3 )
    # }
    
    # # # calculate matrix A^l
    # # A_l = matrix( 0 , nrow = dim_l , ncol = dim_l )
    # # for ( i in 1 : dim_l ){
        # # A_l[ i , i ] = ( 1 / ( h_l ^ 2 ) ) * ( inte_u_hat( xi[ i ] - h_l , xi [ i ] , u , sigma_k ) + inte_u_hat( xi[ i ] , xi[ i ] + h_l , u , sigma_k ) )
        # # if ( i > 1 ){
            # # A_l[ i - 1 , i ] = - ( 1 / ( h_l ^ 2 ) ) * inte_u_hat( xi[ i ] - h_l , xi [ i ] , u , sigma_k )
            # # A_l[ i , i - 1 ] = A_l[ i - 1 , i ]
        # # }
    # # }
    
    # # p_l = solve( A_l ) %*% f_l 
    
    # # use TDMA to solve equations A^l * p^l = f^l_g
    # # store tridiagonal of A^l
    # A_l_b = rep( 0 , dim_l )
    # A_l_a = rep( 0 , dim_l )
    # A_l_c = rep( 0 , dim_l )
    # for ( ii in 1 : dim_l ){
        # A_l_b[ ii ] = ( 1 / ( h_l ^ 2 ) ) * ( inte_u_hat( xi[ ii ] - h_l , xi [ ii ] , u , sigma_k ) + inte_u_hat( xi[ ii ] , xi[ ii ] + h_l , u , sigma_k ) )
        # if ( ii > 1 ){
            # A_l_a[ ii ] = - ( 1 / ( h_l ^ 2 ) ) * inte_u_hat( xi[ ii ] - h_l , xi [ ii ] , u , sigma_k )
            # A_l_c[ ii - 1 ] = A_l_a[ ii ]
        # }
    # }
    # # forward elimination
    # c_prime = rep( 0 , dim_l )
    # for ( ii in 1 : ( dim_l - 1 ) ){
        # if ( ii == 1 ){
            # c_prime[ ii ] = A_l_c[ ii ] / A_l_b[ ii ]
        # }else{
            # c_prime[ ii ] = A_l_c[ ii ] / ( A_l_b[ ii ] - A_l_a[ ii ] * c_prime[ ii - 1 ] )
        # }
    # }
    # d_prime = rep( 0 , dim_l )
    # for ( ii in 1 : dim_l ){
        # if ( ii == 1 ){
            # d_prime[ ii ] = f_l[ ii ] / A_l_b[ ii ]
        # }else{
            # d_prime[ ii ] = ( f_l[ ii ] - A_l_a[ ii ] * d_prime[ ii - 1 ] ) / ( A_l_b[ ii ] - A_l_a[ ii ] * c_prime[ ii - 1 ] )
        # }
    # }
    # # back substitution
    # p_l = rep( 0 , dim_l )
    # for ( ii in dim_l : 1 ){
        # if ( ii == dim_l ){
            # p_l[ ii ] = d_prime[ ii ]
        # }else{
            # p_l[ ii ] = d_prime[ ii ] - c_prime[ ii ] * p_l[ ii + 1 ]
        # }
    # }
       
    # return ( p_l )
# }


pde_solver_base <- function(u, sigma_k, l, small_k){
    # solve for p^l
    dim_l = 2^(l + small_k) - 1
    h_l = 2^(-(l + small_k))
    xi = h_l * seq(1, (dim_l + 1))
    
    # calculate vector f^l : f^l_i = < f , \psi^l_i >
    f_l = rep(0, dim_l)
    for (i in 1:dim_l){
        f_l[i] = (1/h_l) * (100/6) * ((xi[i] + h_l)^3 + (xi[i] - h_l)^3 - 2*xi[i]^3)
    }
    # f_l = ( 1 / h_l ) * ( 100 / 6 ) * sapply( 1 : dim_l , function(z) {( xi[ z ] + h_l ) ^ 3 + ( xi[ z ] - h_l ) ^ 3 - 2 * xi[ z ] ^ 3} )
    
    # # calculate matrix A^l
    # use TDMA to solve equations A^l * p^l = f^l_g
    # store tridiagonal of A^l
    inte_u_hat_store = rep(0, dim_l + 1)
    for (ii in 1:(dim_l + 1)){
        inte_u_hat_store[ii] = inte_u_hat(xi[ii] - h_l, xi[ii], u, sigma_k)
    }
    A_l_b = rep(0, dim_l)
    A_l_a = rep(0, dim_l)
    A_l_c = rep(0, dim_l)
    for (ii in 1:dim_l){
        A_l_b[ii] = (1/(h_l^2)) * (inte_u_hat_store[ii] + inte_u_hat_store[ii + 1])
        if (ii > 1){
            A_l_a[ii] = -(1/(h_l^2)) * inte_u_hat_store[ii]
            A_l_c[ii-1] = A_l_a[ii]
        }
    }
    # forward elimination
    c_prime = rep(0, dim_l)
    for (ii in 1:(dim_l - 1)){
        if (ii == 1){
            c_prime[ii] = A_l_c[ii] / A_l_b[ii]
        }else{
            c_prime[ii] = A_l_c[ii] / (A_l_b[ii] - A_l_a[ii] * c_prime[ii - 1])
        }
    }
    d_prime = rep(0, dim_l)
    for (ii in 1:dim_l){
        if (ii == 1){
            d_prime[ii] = f_l[ii] / A_l_b[ii]
        }else{
            d_prime[ii] = (f_l[ii] - A_l_a[ii] * d_prime[ii - 1]) / (A_l_b[ii] - A_l_a[ii] * c_prime[ii - 1])
        }
    }
    # back substitution
    p_l = rep(0, dim_l)
    for (ii in dim_l:1){
        if (ii == dim_l){
            p_l[ii] = d_prime[ii]
        }else{
            p_l[ii] = d_prime[ii] - c_prime[ii] * p_l[ii + 1]
        }
    }
       
    return (p_l)
}


pde_solver <- function(x, u, sigma_k, l, small_k, p_l){
    # p^l(x) = \sum_{i=1}^{2^(l+k)-1}p^l_i\Psi_i^l(x)
    p_sum = 0
    dim_l = 2^(l + small_k) - 1
    h_l = 2^(-(l + small_k)) 
    xi = h_l * seq(1, dim_l)
    # p_sum = sum( p_l * sapply( 1 : dim_l , function(z) { psi_l( x , xi[ z ] , h_l ) } ) )
    for (i in 1:dim_l){
        p_sum = p_sum + p_l[i] * psi_l(x, xi[i], h_l)
    }
    return (p_sum)
}


r_l <- function(u, sigma_k, l, small_k, gp_x, y, gamma_y){
    if (max(abs(u)) > 1){
        return (0)
    }else{
        p_l = pde_solver_base(u, sigma_k, l, small_k)
        gp = c(pde_solver(gp_x[1], u, sigma_k, l, small_k, p_l), pde_solver(gp_x[2], u, sigma_k, l, small_k, p_l))
        ss = exp(-0.5 * t(gp - y) %*% solve(gamma_y) %*% (gp - y))
        return (ss)
    }
}


r_l_log <- function(u, sigma_k, l, small_k, gp_x, y, gamma_y){
    if (max(abs(u)) > 1){
        return (-10^10) 
    }else{
        p_l = pde_solver_base(u, sigma_k, l, small_k)
        gp = c(pde_solver(gp_x[1], u, sigma_k, l, small_k, p_l), pde_solver(gp_x[2], u, sigma_k, l, small_k, p_l))
        ss = -0.5 * t(gp - y) %*% solve(gamma_y) %*% (gp - y)
        return (ss)
    }
}


g_func_inv <- function(u, param, mod_param, l_level, y){
    # g_function for inverse problem:\partile_\theta\log(r_theta^l(u))
    big_k = param$big_k  # number of u_i
    small_k = param$small_k  # parameter in step-size
    gp_x = param$gp_x  # G(u)=[ p(0.25;u),p(0.75;u) ]
    # x_g = mod_param$x_g
    # l_g = mod_param$l_g
    theta_0 = mod_param$theta_0
    y_dim = length(y)
    
    # sigma_k = (2/5)4^(-k)
    sigma_k = (2/5) * 4^(-seq(1, big_k))
    p_l = pde_solver_base(u, sigma_k, l_level, small_k)
    gp = c(pde_solver(gp_x[1], u, sigma_k, l_level, small_k, p_l), pde_solver(gp_x[2], u, sigma_k, l_level, small_k, p_l))
    ss =  0.5 * t(gp - y) %*% (gp - y) / (theta_0^2) - y_dim / (2*theta_0)
    return (ss)
}


# r_l_log <- function( u , sigma_k , l , small_k , gp_x , y , gamma_y ){
    # # test for normal
    # if ( max( abs( u ) ) > 1 ){
        # return ( -Inf ) 
    # }else{
        # ss = dnorm( u , mean = 0.5 + 2^(-l-1) , sd = 1 , log = TRUE )
        # return ( sum(ss) )
    # }
    # # ss = dnorm( u , mean = 1/(l+1) , sd = 1 , log = TRUE )
    # # return ( sum(ss) )
# }


generate_data <- function(param, mod_param){
  # generate data
  library(MASS)
  
  big_k = mod_param$big_k  # number of u_i
  small_k = mod_param$small_k  # parameter in step-size
  big_M = mod_param$big_M
  gp_x = seq(1/(big_M+1), 1-1/(big_M+1), by=1/(big_M+1))  # G(u)=[ p(x1;u),p(x2;u),\dots,p(x_M;u)]
  gamma_y = (mod_param$theta_0)^(-1) * diag(big_M)   # covariance matrix of noise in data
  
  l_simu = mod_param$l_simu  # h = 2^(-20), when generate the data
  
  # sigma_k = (2/5)4^(-k)
  sigma_k = (2/5) * 4^(- seq(1, big_k))
  # simulate u_i ~ Uniform(-1,1)
  u_simu = runif(big_k, min=-1, max=1)
  p_l_simu = pde_solver_base(u_simu, sigma_k, l_simu, small_k)
  gp_simu = sapply(gp_x, function(z) {pde_solver(z, u_simu, sigma_k, l_simu, small_k, p_l_simu)})
  # gp_simu = c(pde_solver(gp_x[1], u_simu, sigma_k, l_simu, small_k, p_l_simu), pde_solver(gp_x[2], u_simu, sigma_k, l_simu, small_k, p_l_simu))
  # generate noise
  noise_y = mvrnorm(n=1, mu=rep(0, nrow(gamma_y)), Sigma=gamma_y)
  y_data = gp_simu + noise_y
  return (y_data)
}


p_example <- function(x, u){
  # PDE solution in analytically calculable example
  return ((u/2) * (x^2-x))
}


generate_data_example <- function(param, mod_param){
  # generate data
  library(MASS)
  big_M = mod_param$big_M
  gp_x = mod_param$gp_x  # G(u)=[ p(x1;u),p(x2;u),\dots,p(x_M;u)]
  big_k = mod_param$big_k  # number of u_i
  # simulate u_i ~ Uniform(-1,1)
  u_simu = runif(big_k, min=-1, max=1)
  Gu = sapply(gp_x, function(z) {p_example(z, u_simu)})
  # generate noise
  gamma_y = (mod_param$theta_0)^(-1) * diag(big_M)   # covariance matrix of noise in data
  noise_y = mvrnorm(n=1, mu=rep(0, nrow(gamma_y)), Sigma=gamma_y)
  y_data = Gu + noise_y
  return (y_data)
}


g_func_pap <- function(u, param, mod_param){
    big_k = param$big_k
    small_k = param$small_k
    x_g = mod_param$x_g
    l_g = mod_param$l_g
    
    # sigma_k = (2/5)4^(-k)
    sigma_k = (2/5) * 4^(-seq(1, big_k))
    p_l = pde_solver_base(u, sigma_k, l_g, small_k)
    return (pde_solver(x_g, u, sigma_k, l_g, small_k, p_l))
}


# g_func <- function( u , param , mod_param , l_level, y ){  # test for normal
    # return ( sum(u) )
# }


my_sample <- function(x, num, prob){
    x_ret = c()
    c_prob = prob
    for(i in 1:length(prob)){
        c_prob[i] = sum(prob[1:i])
    }
    for(i in 1:num){
        runi = runif(1)
        i_ind = 1
        while(c_prob[i_ind] < runi){
            i_ind = i_ind+1
        }
        x_ret[i] = x[i_ind]
    }
    return(x_ret)
}




