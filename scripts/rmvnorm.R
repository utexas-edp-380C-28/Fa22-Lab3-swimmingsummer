#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
#' Discuss what this function does...
#'
#' @param n Describe this parameter here...
#' @param mu Describe this parameter here...
#' @param Sigma Describe this parameter here...
#' @return Describe what this returns...

rmvnorm <- function(n, mu, Sigma, rho){
  
  Z <- matrix(rnorm(length(mu) * n, 0, 1), ncol = length(mu))
  cor_mat <- matrix(rho, length(Sigma), length(Sigma))
  diag(cor_mat) <- 1
  sd_mat <- diag(length(Sigma))
  diag(sd_mat) <- Sigma
  
  cov_mat <- sd_mat %*% cor_mat %*% sd_mat
  X <- matrix(1, n) %*% t(mu) + Z %*% chol(cov_mat)
  return(X)
}
