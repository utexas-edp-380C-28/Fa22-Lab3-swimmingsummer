#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
# Name:  Suyoung Kim
#---------------------------------------------------#

## SETUP -------------------------------------------
# Source rmvnorm()
source("scripts/rmvnorm.R")

## Begin Code --------------------------------------

#1.a)

## Set seed / Defining parameters for x and y 

set.seed(17290)

p_x <- list(
  m_x = -2,
  sd_x = 3
)

p_y <- list(
  beta0 = 12,
  beta1 = 4,
  sd_e = 7
)

n <- 100

## Generating x and y

regfun <- function(n, p_x, p_y){
  x <- with(p_x, rnorm(n, m_x, sd_x))
  y <- with(p_y, beta0 + beta1*x + rnorm(n, 0, sd_e))
  
  return(data.frame(x = x, y = y))
}

## analyze function

anlfun <- function(mydata){
  
  m_x <- mean(mydata[, 1])
  m_y <- mean(mydata[, 2])
  s_x <- sd(mydata[, 1])
  s_y <- sd(mydata[, 2])
  cor <- cor(mydata[, 1], mydata[, 2])
  
  mdl <- lm(y ~ x, data = mydata)
  b0 <- mdl$coefficients[1]
  b1 <- mdl$coefficients[2]
  s_e <- sigma(mdl)
  
  return(c(m_x = m_x, m_y = m_y, s_x = s_x, s_y = s_y, 
           cor = cor, b0 = b0, b1 = b1, s_e = s_e))
}

## Generate 500 replication 

rep <- 500 
replicated_set <- lapply(1:rep, function(i) {
  one_rep <- regfun(n, p_x, p_y)})

results <- sapply(replicated_set, anlfun)

## Answer 
ans1a <- list(
  mean = rowMeans(results), 
  sd = apply(results, MARGIN = 1, sd)
)

# $mean
#  m_x            m_y            s_x            s_y            cor b0.(Intercept)           b1.x            s_e 
# -2.0061317      3.9810216      3.0034995     13.8596003      0.8631812     11.9780298      3.9873856      6.9749258 

# $sd
# m_x            m_y            s_x            s_y            cor b0.(Intercept)           b1.x            s_e 
# 0.31115637     1.41924855     0.20911389     1.03767609     0.02432744     0.86119359     0.24140973     0.49070350


#2.a)

OLS <- function(y, X){
  
  y <- as.matrix(y)
  X <- as.matrix(cbind(1, X))
  
  # Estimates
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  
  res <- y - (X %*% beta)
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  se <- sqrt(t(res) %*% res / df)
  
  beta_cov <- as.numeric(se^2) * (solve(t(X) %*% X))
  beta_se <- sqrt(diag(beta_cov))
  t_val <- beta / beta_se
  
  rsq <- (t(beta) %*% cov(X) %*% beta) / (t(beta) %*% cov(X) %*% beta + se^2)
  p_v <- round(2 * pt(abs(t_val), df, lower = FALSE), 6)  
  
  result <- matrix(NA, nrow = p + 2, ncol = 4)
  colnames(result) <- c("Estimate", "SE", "t value", "Pr(>|t|)")
  rownames(result) <- c(paste0('b', rep(0:(p-1))), "SD(e)", "R2")
  result[, 1] <- c(beta, se, rsq)
  result[, 2] <- c(beta_se, NA, NA)
  result[, 3] <- c(t_val, NA, NA)
  result[, 4] <- c(p_v, NA, NA)
  
  return(result)
}

# Function test 
ols_1 <- OLS(mtcars$mpg, cbind(mtcars$wt, mtcars$cyl, mtcars$gear))

## Comparing result
ols_2 <- summary(lm(mpg ~ wt + cyl + gear, data = mtcars))

## Answer
# My OLS function

# Estimate        SE    t value Pr(>|t|)
# b0    42.3863641 4.3789952  9.6794726 0.000000
# b1    -3.3920819 0.8208025 -4.1326409 0.000294
# b2    -1.5280010 0.4197533 -3.6402363 0.001093
# b3    -0.5228629 0.7788703 -0.6713094 0.507524
# SD(e)  2.5921848        NA         NA       NA
# R2     0.8182681        NA         NA       NA

# lm function

# lm(formula = mpg ~ wt + cyl + gear, data = mtcars)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -4.8443 -1.5455 -0.3932  1.4220  5.9416 

# Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  42.3864     4.3790   9.679 1.97e-10 ***
#  wt           -3.3921     0.8208  -4.133 0.000294 ***
#  cyl          -1.5280     0.4198  -3.640 0.001093 ** 
#  gear         -0.5229     0.7789  -0.671 0.507524    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 2.592 on 28 degrees of freedom
# Multiple R-squared:  0.8329,	Adjusted R-squared:  0.815 
# F-statistic: 46.53 on 3 and 28 DF,  p-value: 5.262e-11

# 3.a) 

## Set seed and parameters 

set.seed(21389)

n <- 100000
p_x1x2 <- list(
  rho = 0.3,
  mu = c(5, 10),
  Sigma = c(1, 2)
)

## Drawing X 

x1 <- with(p_x1x2, rmvnorm(n, mu, Sigma, rho))

# 3.b)

## Set seed and parameters 

set.seed(23921)
p_y1 <- list(
  beta = c(1, 1),
  r2 = 0.6, 
  mu = 10
)

## Drawing sigma_e 

eq22 <- with(p_y1, t(beta) %*% cov(x1) %*% beta)
eq23 <- with(p_y1, eq22 * (1/r2 -1))
sigma_e <- sqrt(eq23)
beta0 <- with(p_y1, mu - (t(p_x1x2$mu) %*% beta))

## Drawing y 

y1 <- with(p_y1, (matrix(1, n) %*% beta0)+ x1 %*% beta + rnorm(n, 0, sigma_e))

## OLS function

ols_3 <- OLS(y1, x1)

ans_3b <- rbind(b0 = cbind(beta0, ols_3[1 ,1], beta0 - ols_3[1, 1]),
                b1 = cbind(p_y1$beta[1], ols_3[2, 1], p_y1$beta[1] - ols_3[2, 1]),
                b2 = cbind(p_y1$beta[2], ols_3[3, 1], p_y1$beta[1] - ols_3[3, 1]),
                R2 = cbind(p_y1$r2, ols_3[5, 1], p_y1$r2 - ols_3[5, 1]),
                SE = cbind(sigma_e, ols_3[4, 1], sigma_e - ols_3[4, 1]),
                mu = cbind(p_y1$mu, mean(y1), p_y1$mu - mean(y1)))

colnames(ans_3b) <- c("Population", "Estimate", "Difference")
rownames(ans_3b) <- c("b0", "b1", "b2", "R2", "SD(e)", "mu")


## Answer for population difference

# Population   Estimate    Difference
# b0     -5.000000 -4.9370527 -0.0629472651
# b1      1.000000  0.9891609  0.0108391155
# b2      1.000000  0.9997359  0.0002641366
# R2      0.600000  0.5978317  0.0021682893
# SD(e)   2.033295  2.0364085 -0.0031133486
# mu     10.000000 10.0016495 -0.0016494593
# 3.c)

## Set seed and parameters

set.seed(123782)
p_y2 <- list(
  r = c(0.3, -0.4),
  mu = 10,
  sd = 5
)

## Cov_mat / Cor_mat 

x.cor_mat <- with(p_x1x2, matrix(c(1, rho, rho, 1), nrow = 2, ncol =2))
cor_mat <- with(p_y2, rbind(c(1, r), cbind(r, x.cor_mat)))
Sigma <- diag(c(p_y2$sd, p_x1x2$Sigma)) %*% cor_mat %*% diag(c(p_y2$sd, p_x1x2$Sigma))

## Beta 
beta_n <- solve(Sigma[2:3, 2:3]) %*% Sigma[2:3, 1]

## Residual Variance
eq26 <- with(p_y2, sd^2 - t(beta_n) %*% (cor_mat[2:3, 1]))
sen <- sqrt(eq26)

## Rsquared 
rsqn <- cor_mat[1, -1] %*% solve(cor_mat[-1, -1]) %*% cor_mat[-1, 1]


## Beta_0

beta0_n <- with(p_y2, mu - (t(p_x1x2$mu) %*% beta_n))

y2 <- with(p_y2, (matrix(1, n) %*% beta0_n)+ x1 %*% beta_n + rnorm(n, 0, sen))

ols_3c <- OLS(y2, x1)

ans_3c <- rbind(b0 = cbind(beta0_n, ols_3c[1, 1], beta0_n - ols_3c[1, 1]),
                b1 = cbind(beta_n[1, 1], ols_3c[2, 1], beta_n[1, 1] - ols_3c[2, 1]),
                b2 = cbind(beta_n[2, 1], ols_3c[3, 1],beta_n[2, 1] - ols_3c[3, 1]),
                SE = cbind(sen, ols_3c[4, 1], sen - ols_3c[4, 1]), 
                R2 = cbind(rsqn, ols_3c[5, 1], rsqn - ols_3c[5, 1]))


colnames(ans_3c) <- c("Population", "Estimate", "Difference")
rownames(ans_3c) <- c("b0", "b1", "b2", "SD(e)", "R2")

## Answer

# Population   Estimate   Difference
# b0    11.9230769 11.7478467  0.175230230
# b1     2.3076923  2.3204565 -0.012764167
# b2    -1.3461538 -1.3348942 -0.011259675
# SD(e)  4.8753698  4.8719020  0.003467843
# R2     0.3538462  0.2697996  0.084046551

## 3.d)

gen_m1 <- function(n, p_x, p_y){
  
  # Method 1
  
  ## Drawing X 
  
  p <- length(p_x$mu)
  x <- with(p_x, rmvnorm(n, mu, Sigma, rho))
  
  beta <- with(p_y, matrix(beta, nrow = p))
  eq23 <- with(p_y, t(beta) %*% cov(x) %*% beta * (1/r2 -1))
  sigma_e <- sqrt(eq23)
  beta0 <- with(p_y, mu - (t(p_x$mu) %*% beta))
  
  y1 <- with(p_y, (matrix(1, n) %*% beta0)+ x %*% beta + rnorm(n, 0, sigma_e))
  
  output <- cbind(y1,x)
  attr(output, 'par_xy') <- list(p_x, p_y)
  return(output) 
}

## Parameters 
set.seed(6972)
n <- 100000
p_x <- list(mu = rep(0, 5),
            Sigma = c(1, sqrt(2), sqrt(3), 2, sqrt(5)), 
            rho = 0.15)

p_y <- list(beta = rep(1, 5), mu = 25, r2 = 0.5)

m1 <- gen_m1(n, p_x, p_y)

m_3d1 <- OLS(m1[, 1], cbind(m1[,2:6]))

# Estimate          SE    t value Pr(>|t|)
# b0    25.0219149 0.015222210 1643.77672        0
# b1     0.9812633 0.015651423   62.69483        0
# b2     1.0199224 0.011081135   92.04132        0
# b3     0.9881525 0.009059712  109.07107        0
# b4     0.9852174 0.007842738  125.62161        0
# b5     1.0050576 0.007055326  142.45375        0
# SD(e)  4.8135563          NA         NA       NA
# R2     0.4991145          NA         NA       NA

#### Answer 
ans_3d1 <- rbind(b1 = cbind(p_y$beta[1], m_3d1[2, 1], p_y$beta[1] - m_3d1[2, 1]),
                 b2 = cbind(p_y$beta[2], m_3d1[3, 1],p_y$beta[2] - m_3d1[3, 1]),
                 b3 = cbind(p_y$beta[3], m_3d1[4, 1],p_y$beta[3] - m_3d1[4, 1]),
                 b4 = cbind(p_y$beta[4], m_3d1[5, 1],p_y$beta[4] - m_3d1[5, 1]),
                 b5 = cbind(p_y$beta[5], m_3d1[6, 1],p_y$beta[5] - m_3d1[6, 1]),
                 rsq = cbind(p_y$r2, m_3d1[8, 1], p_y$r2 - m_3d1[8, 1]),
                 SE = cbind(sigma_e, m_3d1[7, 1], sigma_e - m_3d1[7, 1]),
                 mu = cbind(p_y$mu, mean(m1[,1]),p_y$mu - mean(m1[,1])))

colnames(ans_3d1) <- c("Population", "Estimate", "Difference")
rownames(ans_3d1) <- c("b1", "b2", "b3", "b4", "b5", "R2", "SD(e)", "mu")


# Population   Estimate    Difference
# b1      1.000000  0.9812633  0.0187366927
# b2      1.000000  1.0199224 -0.0199223718
# b3      1.000000  0.9881525  0.0118474752
# b4      1.000000  0.9852174  0.0147825927
# b5      1.000000  1.0050576 -0.0050575912
# R2      0.500000  0.4991145  0.0008855123
# SD(e)   4.835865  4.8135563  0.0223090964
# mu     25.000000 25.0334707 -0.0334707122

## Generating Method 2 
gen_m2 <- function(n, p_x, p_y){
  
  # Method 2
  p <- length(p_x$mu)
  x <- with(p_x, rmvnorm(n, mu, Sigma, rho))
  
  x_cor_mat <- with(p_x, matrix(rho, length(Sigma), length(Sigma)))
  diag(x_cor_mat) <- 1
  cor_mat <- with(p_y, rbind(c(1, rho), cbind(rho, x_cor_mat)))
  Sigma <- diag(c(p_y$sd, p_x$Sigma)) %*% cor_mat %*% diag(c(p_y$sd, p_x$Sigma))
  
  ## Beta 
  beta <- solve(Sigma[2:(p+1), 2:(p+1)]) %*% Sigma[2: (p +1), 1]
  
  ## Residual Variance
  eq26 <- with(p_y, sd^2 - t(beta) %*% Sigma[2:(p + 1), 1])
  se <- sqrt(eq26)
  
  ## Rsquared 
  rsqn <- cor_mat[1, -1] %*% solve(cor_mat[-1, -1]) %*% cor_mat[-1, 1]
  
  
  ## Beta_0
  
  beta0 <- with(p_y, mu - (t(p_x$mu) %*% beta))
  
  y2 <- with(p_y, (matrix(1, n) %*% beta0)+ x %*% beta + rnorm(n, 0, se))
  
  output <- cbind(y2,x)
  attr(output, 'par_xy') <- c(p_x, p_y, beta, beta0, rsqn, se)
  return(output) 
}

## set seed and param
set.seed(1237)
p_x <- list(mu = rep(0, 5),
            Sigma = c(1, sqrt(2), sqrt(3), 2, sqrt(5)), 
            rho = 0.15)

p_y <- list(rho = c(-0.15, -0.5, 0.15, 0.3, 0.20), 
            mu = 10, 
            sd = 4)

## gen data using method 2 
m2 <- gen_m2(10000, p_x, p_y)
m_3d2 <- OLS(m2[,1], cbind(m2[,2:6]))

## Answer
ans_3d2 <- rbind(b0 = cbind(beta0, m_3d2[1, 1], beta0 - m_3d2[1, 1]),
                 b1 = cbind(beta[1], m_3d2[2, 1], beta[1] - m_3d2[2, 1]),
                 b2 = cbind(beta[2], m_3d2[3, 1], beta[2] - m_3d2[3, 1]),
                 b3 = cbind(beta[3], m_3d2[4, 1], beta[3]- m_3d2[4, 1]),
                 b4 = cbind(beta[4], m_3d2[5, 1], beta[4]- m_3d2[5, 1]),
                 b5 = cbind(beta[5], m_3d2[6, 1], beta[5]- m_3d2[6, 1]),
                 se = cbind(se, m_3d2[7, 1], se - m_3d2[7, 1]),
                 r2 = cbind(rsqn, m_3d2[8, 1], rsqn - m_3d2[8, 1])
                 
)

colnames(ans_3d2) <- c("Population", "Estimate", "Difference")
rownames(ans_3d2) <- c("b0", "b1", "b2", "b3", "b4", "b5", "SD(e)", "R2")


#### Answer

#Population   Estimate   Difference
# b0    10.0000000  9.9958695  0.004130497
# b1    -0.7058824 -0.7106237  0.004741309
# b2    -1.6637807 -1.6813177  0.017536996
# b3     0.4075414  0.3964946  0.011046756
# b4     0.7058824  0.7044516  0.001430782
# b5     0.4209069  0.4194456  0.001461357
# SD(e)  2.8284271  2.8030957  0.025331412
# R2     0.5000000  0.5030344 -0.003034398
