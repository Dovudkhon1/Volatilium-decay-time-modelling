library(stats)
library(ggplot2)

# Loading data
surv_data <- read.csv("/Users/dovudkhon/Desktop/Bristol/Inference/surv.csv")

# Extract some important data
init <- surv_data$init
fin <- surv_data$fin
times <- surv_data$times

log_fin <- log(surv_data$fin)
log_init <- log(surv_data$init)
 

# Question 1: Henri's model
# Negative log-likelihood function of the residuals, what is happening here is because the
# residuals follow a normal, the variable (log_fin - log_init - phi1 - phi2 * times) is 
# normally distributed, allowing us to find the likelihood function to then find phi1, phi2
# and sigma. Writing this function will be important as the optim function can be used to 
# find the needed estimates.
neg_log_likelihood_Henri <- function(params, log_fin, log_init, times) {
  phi1 <- params[1]
  phi2 <- params[2]
  log_sigma_sq <- params[3]  # Log-variance to make sure of positivity
  sigma_sq <- exp(log_sigma_sq)
  
  residuals <- log_fin - log_init - phi1 - phi2 * times
  n <- length(log_fin)
  
  nll <- (n/2) * (log(2 * pi) + log_sigma_sq) + (1/(2 * sigma_sq)) * sum(residuals^2)
  return(nll)
}

# The optim function will be used to find the optimizing phi1, phi2 and sigma for the model given the data.
initial_params <- c(0, 0, 0)
result_Henri <- optim(
  par = initial_params,
  fn = neg_log_likelihood_Henri,
  log_fin = log_fin,
  log_init = log_init,
  times = times,
  method = "BFGS",
  hessian = TRUE
)
print(result_Henri)
# Final estimates and printing the outcomes
phi1_mle <- result_Henri$par[1]
phi2_mle <- result_Henri$par[2]
sigma_sq_mle <- exp(result_Henri$par[3])
sigma_mle <- sqrt(sigma_sq_mle)  # Standard deviation
print(sigma_mle)
cat("MLE Estimates for Henri's model:\n",
    "φ1 =", phi1_mle, "\n",
    "φ2 =", phi2_mle, "\n",
    "σ =", sigma_mle, "\n")

# Calculate predictions, for later usage
predicted_log_fin_Henri <- log_init + phi1_mle + phi2_mle * times
predicted_fin_Henri <- exp(predicted_log_fin_Henri)

# Calculate 95% confidence intervals
se <- sigma_mle  # Standard error of the residuals
lower_bound_Henri <- exp(predicted_log_fin_Henri - 1.96 * se)  # Lower bound (original scale)
upper_bound_Henri <- exp(predicted_log_fin_Henri + 1.96 * se)  # Upper bound (original scale)

# Creating needed plots
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))

# Observed vs Predicted final counts, including a line of the best fit
plot(fin, predicted_fin_Henri, 
     xlab = "Predicted Final Count", ylab = "Observed Final Count",
     main = "Observed vs Predicted Counts(Henri's model)",
     pch = 19, col = "#1f77b4", log="xy" )
abline(a = 0, b = 1, col = "#ff7f0e", lwd = 2)

# Add 95% confidence intervals
arrows(fin, lower_bound_Henri, fin, upper_bound_Henri, 
       length = 0.05, angle = 90, code = 3, col = "darkgray")

# Add legend
legend("topleft", 
       legend = c("Predicted", "1:1 Line", "95% CI"), 
       col = c("#1f77b4", "#ff7f0e", "darkgray"), 
       pch = c(19, NA, NA), 
       lty = c(NA, 1, 1))

# Residuals vs Fitted, their log values
residuals_log_Henri <- log_fin - predicted_log_fin_Henri
plot(predicted_log_fin_Henri, residuals_log_Henri,
     xlab = "Fitted Values(log)", ylab = "Residuals(log)",
     main = "Residuals vs Fitted(Henri's model)",
     pch = 19, col = "#2ca02c")
abline(h = 0, col = "#d62728", lty = 2, lwd = 2)

# QQ-plot of residuals
qqnorm(residuals_log_Henri, main = "QQ Plot of Residuals(Henri's model)")
qqline(residuals_log_Henri, col = "red")

#Comments about the model


# The estimated variance of the residuals is around 0.09, so 6/70 of the residuals of the lower
# observed counts are outliers as they differ by at least twice the variance from 0, which could suggest in most
# test this model may not be so adequate when discribing models of lower count, but when it comes
# to higher predictions, it is mostly accurate.

# Question 2: Hypothesis test for Pierre's claim
# H0: phi1 = 0, H1 = phi1!= 0.
# Extract standard error of phi1 from the inverse Hessian matrix

se_phi1 <- sqrt(solve(result_Henri$hessian)[1, 1])  # First diagonal element is variance of phi1

# Compute Wald test statistic
Z_phi1 <- phi1_mle / se_phi1

# Compute two-tailed p-value
p_value <- 2 * (1 - pnorm(abs(Z_phi1)))

# Print results
cat("Wald Test Statistic for H0: phi1 = 0:", Z_phi1, "\n")
cat("p-value:", p_value, "\n")

# Decision rule
if (p_value < 0.05) {
  cat("Reject H0: phi1 is significantly different from 0. Pierre's claim is not supported.\n")
} else {
  cat("Fail to reject H0: No significant evidence against Pierre's claim. His hypothesis holds.\n")
}

#Question 3: Marie's Model

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
# Define the negative log-likelihood function for Marie's model
# Use Q1-based MLE to initialise Q3 parameters

theta1_init <- exp(phi1_mle)
theta2_init <- -phi2_mle
log_s_sq_init <- log(sigma_mle^2)

initial_params_Marie <- c(theta1_init, theta2_init, log_s_sq_init)

# Optimising Marie's model
# The optim function will be used to find the optimizing 
# theta1, theta2 and log_s_sq for the model given the data.
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
cat("MLE Estimate for Marie's Module:\n",
    "theta1 =", theta1_mle, "\n",
    "theta2 =", theta2_mle, "\n",
    "s      =", s_mle, "\n")

#Calculate the predicted final counts under Marie's model
predicted_fin_Marie <- init * (theta1_mle * exp(-theta2_mle * times))

#Compute the 95% confidence intervals
lower_bound_Marie <- predicted_fin_Marie - 1.96 * s_mle * init
upper_bound_Marie <- predicted_fin_Marie + 1.96 * s_mle * init

#Plot Observed vs. Predicted final counts, including a line of the best fit
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
plot(fin, predicted_fin_Marie,
     xlab = "Predicted Final Count", 
     ylab = "Observed Final Count",
     main = "Observed vs Predicted (Marie's model)",
     pch = 19, col = "#1f77b4", log = "xy")
abline(a = 0, b = 1, col = "#ff7f0e", lwd = 2)
arrows(fin, lower_bound_Marie, 
       fin, upper_bound_Marie,
       length = 0.05, angle = 90, code = 3, col = "darkgray")
#Draw vertical error bars from lower_bound to upper_bound at each predicted value.

legend("topleft",
       legend = c("Predicted", "1:1 Line", "95% CI"),
       col    = c("#1f77b4", "#ff7f0e", "darkgray"),
       pch    = c(19, NA, NA),
       lty    = c(NA, 1, 1))

#Residuals vs Fitted
residuals_Marie = fin - predicted_fin_Marie
plot(predicted_fin_Marie, residuals_Marie,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs Fitted (Marie's model)",
     pch = 19, col = "#2ca02c")

#Draw a horizontal line at y=0, helping to see if residuals center around zero.
abline(h = 0, col = "#d62728", lty = 2, lwd = 2)

# QQ-plot of residuals
qqnorm(residuals_Marie, main = "QQ Plot of Residuals(Marie's model)")
qqline(residuals_Marie, col = "red")

#comments for Marie's model
#The estimated variance of the residuals is around 0.015.
#As a result, the model may not accurately capture lower count behaviour
#but for higher counts it performs better.

AIC_Henri <- 2 * length(result_Henri$par) + 2 * result_Henri$value
AIC_Marie <- 2 * length(result_Marie$par) + 2 * result_Marie$value

cat("AIC (Henri):", AIC_Henri, "\n")
cat("AIC (Marie):", AIC_Marie, "\n")
