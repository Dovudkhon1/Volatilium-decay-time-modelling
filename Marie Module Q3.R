#Q3: Marie Model

#Loading data
surv_data <- read.csv("surv.csv")

#Make some important data
init <- surv_data$init
fin <- surv_data$fin
times <- surv_data$times


############################################################
#Replicating the Q1 process in order to get phi1_mle, phi2_mle, sigma_mle
#This part of the code is the same as Q1

neg_log_likelihood_Henri <- function(params, log_fin, log_init, times) {
  phi1 <- params[1]
  phi2 <- params[2]
  log_sigma_sq <- params[3]  # Log-variance to enforce positivity
  sigma_sq <- exp(log_sigma_sq)
  
  residuals <- log_fin - log_init - phi1 - phi2 * times
  n <- length(log_fin)
  
  nll <- (n/2) * (log(2 * pi) + log_sigma_sq) + (1/(2 * sigma_sq)) * sum(residuals^2)
  return(nll)
}

#Calculate the logarithmic data
log_fin <- log(fin)
log_init <- log(init)

#The optim function will be used to find the optiming phi1, phi2 and sigma for the model given the data.
initial_params_Henri <- c(0, 0, 0)  #Simple initial value
result_Henri <- optim(
  par       = initial_params_Henri,
  fn        = neg_log_likelihood_Henri,
  log_fin   = log_fin,
  log_init  = log_init,
  times     = times,
  method    = "BFGS",
  hessian   = TRUE
)


phi1_mle <- result_Henri$par[1]
phi2_mle <- result_Henri$par[2]
sigma_sq_mle <- exp(result_Henri$par[3])
sigma_mle <- sqrt(sigma_sq_mle)

cat("MLE Estimates:\n",
    "φ1 =", phi1_mle, "\n",
    "φ2 =", phi2_mle, "\n",
    "σ =", sigma_mle, "\n")


############################################################
#Start with Marie Module in Q3

#Define negative log-likelihood function of the residuals. 
#Params is a numeric vector of length 3 that holds (theta1, theta2, log_s_sq).
#fin, init, and times are vectors of the same number of observations.
#Extract three variables from the parameter vector.
neg_log_likelihood_Marie <- function(params, fin, init, times) {
  theta1   <- params[1]      
  theta2   <- params[2]      
  log_s_sq <- params[3]     
  
  #Convert log_s_sq to s^2 by exponentiating
  s_sq <- exp(log_s_sq)      
  
  #Compute the predicted final count under the model
  pred <- init * (theta1 * exp(-theta2 * times))
  
  #Calculate residuals (observed - predicted)
  resid <- fin - pred
  
  #Compute the negative log-likelihood contribution.
  nll <- sum(
    0.5 * log(2 * pi) + 0.5 * (log_s_sq + 2 * log(init)) + 0.5 * (resid^2) / (s_sq * init^2)
  )
  
  return(nll) #Return the total negative log-likelihood
}

#Q1-based MLE to initialise Q3 parameters

theta1_init <- exp(phi1_mle)
theta2_init <- -phi2_mle
log_s_sq_init <- log(sigma_mle^2)

initial_params_Marie <- c(theta1_init, theta2_init, log_s_sq_init)

#Optimisation Marie Module
#The optim function will be used to find the optiming theta1, theta2 and log_s_sq for the model given the data.
result_Marie <- optim(
  par   = initial_params_Marie,
  fn    = neg_log_likelihood_Marie,
  fin   = fin,
  init  = init,
  times = times,
  method = "BFGS",
  hessian = TRUE
)

#Extract maximum likelihood estimates from the optimization output
theta1_mle <- result_Marie$par[1] #Extract the estimate for theta1.
theta2_mle <- result_Marie$par[2] #Extract the estimate for theta2.
s_sq_mle   <- exp(result_Marie$par[3]) #The estimate for s^2 is exp(log_s_sq).
s_mle      <- sqrt(s_sq_mle) #The standard deviation s is the square root of s^2.

#print the estimated data
cat("MLE Estimate Marie Module:\n")
cat("theta1 =", theta1_mle, "\n")
cat("theta2 =", theta2_mle, "\n")
cat("s      =", s_mle, "\n")



#Calculate the predicted final counts under Marie's model
predicted_fin_Marie <- init * (theta1_mle * exp(-theta2_mle * times))

#Compute the 95% confidence intervals
lower_bound_Marie <- predicted_fin_Marie - 1.96 * s_mle * init
upper_bound_Marie <- predicted_fin_Marie + 1.96 * s_mle * init

#Plot Observed vs. Predicted final counts, including a line of the best fit
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
plot(predicted_fin_Marie, fin,
     xlab = "Predicted Final Count", 
     ylab = "Observed Final Count",
     main = "Observed vs Predicted (Q3: Marie Module)",
     pch = 19, col = "#1f77b4", log = "xy")
abline(a = 0, b = 1, col = "#ff7f0e", lwd = 2)
arrows(predicted_fin_Marie, lower_bound_Marie, 
       predicted_fin_Marie, upper_bound_Marie,
       length = 0.05, angle = 90, code = 3, col = "darkgray")
#Draw vertical error bars from lower_bound to upper_bound at each predicted value.

legend("topleft",
       legend = c("Predicted", "1:1 Line", "95% CI"),
       col    = c("#1f77b4", "#ff7f0e", "darkgray"),
       pch    = c(19, NA, NA),
       lty    = c(NA, 1, 1))

#Residuals vs Fitted, their log values
residuals_Marie <- fin - predicted_fin_Marie
plot(predicted_fin_Marie, residuals_Marie,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs Fitted (Q3: Marie Module)",
     pch = 19, col = "#2ca02c")

#Draw a horizontal line at y=0, helping to see if residuals center around zero.
abline(h = 0, col = "#d62728", lty = 2, lwd = 2)

#comments for Marie Module
#The estimated variance of the residuals is around 0.015.
#As a result, the model may not accurately capture lower count behaviour
#but for higher counts it performs better.
