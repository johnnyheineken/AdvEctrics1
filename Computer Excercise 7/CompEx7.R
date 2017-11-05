montecarlo <- function(gamma_1 = 0, rho = 0.5, N = 100, Rep =200){
  beta <- 1
  
  r <- 3
  pi_1 <- 1
  pi <- rep(0, r)
  pi[1] <- pi_1

  gamma <- rep(0, r)
  gamma[1] <- gamma_1
  

  Z <- matrix(rnorm(N*r,mean=0,sd=1), N, r) 
  sigma <- exp(Z %*% gamma)
  head(sigma)
  

  b_2SLS <- rep(0, Rep)
  b_GMM <- rep(0, Rep)
  b_EGMM <- rep(0, Rep)
  b_OLS <- rep(0, Rep)
  b_GMM <- rep(0, Rep)
  
  for (j in (1:Rep)){
    u <- rnorm(n = N, mean = 0, sd = sigma)
    v <- ((u/sigma)) * (rho + rnorm(N) * sqrt(1-rho^2))
    x <- Z %*% pi + v
    y <- x %*% beta + u
    b_OLS[j] <- solve(t(x) %*% x) %*% t(x) %*% y
    b_2SLS[j] <- (solve(t(x) %*% Z %*% solve(t(Z) %*% Z) %*%
                          t(Z) %*% x) %*% t(x) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% y)
    
    e = y - x %*% b_2SLS[j]
    
    
    S_hat <- matrix(rep(0, r*r),r, r)
    for (i in 1:N){
      u_sq <- (y[i] - x[i]*b_2SLS[j])^2
      S_hat <- S_hat + (u_sq * Z[i, ] %*% t(Z[i, ]))
    }
    S_hat <- S_hat * (1/N)
    
    b_GMM[j] <- (solve(t(x) %*% Z %*% solve(S_hat) %*% t(Z) %*% x) %*% 
      t(x) %*% Z %*% solve(S_hat) %*% t(Z)%*% y)
  }
  result <- list(b_OLS, b_2SLS, b_GMM)
  names(result) <- c("b_OLS", "b_2SLS", "b_GMM")
  return(result)
}

result_vector_creator <- function(beta_vec){
  return(c(mean(beta_vec) - beta, sd(beta_vec), sqrt(mean((beta_vec - beta)^2))))
}



results_creator <- function(beta_list){
  row_names<-names(beta_list)
  l_beta <- length(beta_list)
  results <- matrix(ncol = 3, nrow = l_beta)
  
  for (i in 1:l_beta){
    results[i, ] <- result_vector_creator(beta_list[[i]])
  }
  colnames(results) <- c("bias", "Variance", "RMSE")
  rownames(results) <- row_names
  return(results)
}



results_creator(montecarlo(gamma_1 = 0.8, rho = 0.5))
