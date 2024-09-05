library(MASS)
library(optimbase)
library(ncvreg)
library(leaps)

set.seed(12)

n = 50 #Sample size (n = 50, 100, 200)
standard_deviation = 1 #Variance of model (sigma = 1, 3)

#Setting up covariance matrix for Design Matrix
mean_vector = rep(0, 8)
covariance_matrix_entries = c(1,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,
                             0.5,1,0.5,0.25,0.125,0.0625,0.03125,0.015625,
                             0.25,0.5,1,0.5,0.25,0.125,0.0625,0.03125,
                             0.125,0.25,0.5,1,0.5,0.25,0.125,0.0625,
                             0.0625,0.125,0.25,0.5,1,0.5,0.25,0.125,
                             0.03125,0.0625,0.125,0.25,0.5,1,0.5,0.25,
                             0.015625,0.03125,0.0625,0.125,0.25,0.5,1,0.5,
                             0.0078125,0.015625,0.03125,0.0625,0.125,0.25,0.5,1)
covariance_matrix = matrix(covariance_matrix_entries, nrow = 8, ncol = 8, byrow = TRUE)

number_of_experiments = 1000 #1000 simulations to consider

#Setting up vectors for MRME
GCV_RME_vector = rep(0, number_of_experiments)
BIC_RME_vector = rep(0, number_of_experiments)
oracle_RME_vector = rep(0, number_of_experiments)

#Setting up vectors for number of incorrect/correct zeros
count_zero_GCV_I = rep(0, number_of_experiments)
count_zero_BIC_I = rep(0, number_of_experiments)
count_zero_GCV_C = rep(0, number_of_experiments)
count_zero_BIC_C = rep(0, number_of_experiments)

#For loop that conducts 1000 simulations
for (z in (1:number_of_experiments)){
  x = mvrnorm(n, mu = mean_vector, Sigma = covariance_matrix, empirical = FALSE) #Design Matrix
  colnames(x) = paste0("x", 1:8)
  rownames(x) = paste0(1:n)
  
  beta = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1) #True coefficients
  
  #Obtaining y data through regression equation
  y = rep(0, n)
  for (i in (1:n)){
    y[i] = rnorm(1, mean = sum(x[i,] * beta), sd = standard_deviation)
  }
  
  y_hat_true = x %*% beta #True y_hat values (minus the noise)
  
  #Oracle Model (linear model fixing relevant predictors only)
  oracle_model = lm(y ~ subset(x, select = c(x1,x2,x5)) - 1)
  y_hat_oracle = predict(oracle_model)
  ME_oracle = mean((y_hat_oracle - y_hat_true)^2) #Model Error (Oracle)
  
  #Full Linear Model (all predictors considered)
  full_linear_model = lm(y ~ x - 1)
  y_hat_SF = predict(full_linear_model)
  ME_SF = mean((y_hat_SF - y_hat_true)^2) #Model Error (Unpenalised OLS)
  
  RME_oracle = ME_oracle/ME_SF #Relative Model Error (Oracle)
  
  #GCV Approach 
  scad_model = cv.ncvreg(x, y, family = c("gaussian"), penalty = "SCAD")
  scad_coef = coef(scad_model) #GCV coefficients
  
  #Calculating predicted values (GCV)
  y_hat_GCV = rep(0, n)
  for (i in (1:n)){
    y_hat_GCV[i] = sum(scad_coef[-1] * x[i,])
  }
  
  ME_GCV = mean((y_hat_GCV - y_hat_true)^2) #Model Error (GCV)
  RME_GCV = ME_GCV/ME_SF #Relative Model Error (GCV)
  
  #BIC Approach
  p_prime = function(t, l){
    return(l * ((t <= l) + (((3.7 * l - t) * ((3.7 * l - t) > 0))/(2.7 * l)) * (t > l)))
  } #Function for derivative of SCAD penalty
  
  BIC = rep(0, 1000) #Setting up vector containing all BIC values for 1000 lambdas considered
  
  #For loop calculating all BIC values 
  for (i in (1:1000)){
    lambda = i/1000
    lambda_model = ncvfit(x, y, penalty = "SCAD", lambda = lambda) #fitting model with lambda candidate
    b = lambda_model$beta #coefficients for lambda candidate
    y_hat_lambda = rep(0, n)
    for (j in (1:n)){
      y_hat_lambda[j] = sum(b * x[j,])
    } #predicted values for lambda candidate
    sigma_hat = mean((y_hat_lambda - y)^2) #training error
    #Solving for DF and accounting for division by 0 by reducing dimensions
    p_vector = rep(0, 8)
    for (k in (1:8)){
      p_vector[i] = p_prime(abs(b[k]), lambda)/abs(b[k])
    }
    vanish = rep(1,8)
    for (o in (1:8)){
      if (b[o] == 0){
        vanish[o] = 0
      }
    }
    x_new = NULL
    for (p in (1:8)){
      if (vanish[p] == 1){
        x_new = cbind(x_new, x[,p])
      }
    }
    p_vector = p_vector[!is.na(p_vector) & !is.infinite(p_vector)]
    diagonal_matrix = diag(p_vector, sum(vanish), sum(vanish))
    DF = sum(diag(x_new %*% (t(x_new) %*% x_new + n * diagonal_matrix) %*% t(x_new)))
    BIC[i] = log(sigma_hat) + DF * log(n) / n #BIC for lambda candidate
  }
  
  BIC_lambda = which.min(BIC)/1000 #lambda corresponding to minimum BIC value
  
  BIC_model = ncvfit(x, y, penalty = "SCAD", lambda = BIC_lambda)
  BIC_coeffs = BIC_model$beta #BIC coefficients
  
  #Calculating predicted values (BIC)
  y_hat_BIC = rep(0, n)
  for (i in (1:n)){
    y_hat_BIC[i] = sum(BIC_coeffs * x[i,])
  }
  
  ME_BIC = mean((y_hat_BIC - y_hat_true)^2) #Model Error (BIC)
  RME_BIC = ME_BIC/ME_SF #Relative Model Error (BIC)
  
  #Adding Relative Model Errors for this single simulation to vector
  GCV_RME_vector[z] = RME_GCV
  BIC_RME_vector[z] = RME_BIC
  oracle_RME_vector[z] = RME_oracle
  
  #Adding Number of Incorrect Zeros for this single simulation to vector
  count_zero_GCV_I[z] = sum(beta != 0 & scad_coef[-1] == 0)
  count_zero_BIC_I[z] = sum(beta != 0 & BIC_coeffs == 0)
  
  #Adding Number of Correct Zeros for this single simulation to vector
  count_zero_GCV_C[z] = sum(beta == 0 & scad_coef[-1] == 0)
  count_zero_BIC_C[z] = sum(beta == 0 & BIC_coeffs == 0)
}

#Number of Zeros (I)
mean(count_zero_GCV_I)
mean(count_zero_BIC_I)

#Number of Zeros (C)
mean(count_zero_GCV_C)
mean(count_zero_BIC_C)

#MRME (%)
median(GCV_RME_vector)
median(BIC_RME_vector)
median(oracle_RME_vector)


