#Loading data
data <- read.csv("surv.csv")

#Making some important data
log_fin <- log(data$fin)
log_init <- log(data$init)
times <- data$times
y <- log_fin - log_init 

# Negative log-likelihood function of the residuals, what is happening here is because the
# residuals follow a normal, the variable (log_fin - log_init - phi1 - phi2 * times) is 
# normally distributed, allowing us to find the likelihood function to then find phi1, phi2
# and sigma. Writing this function will be important as the optim function can be used to 
# find the needed estimates.
neg_log_likelihood <- function(params, log_fin, log_init, times) {
  phi1 <- params[1]
  phi2 <- params[2]
  log_sigma_sq <- params[3]  # Log-variance to enforce positivity
  sigma_sq <- exp(log_sigma_sq)
  
  residuals <- log_fin - log_init - phi1 - phi2 * times
  n <- length(log_fin)
  
  nll <- (n/2) * (log(2 * pi) + log_sigma_sq) + (1/(2 * sigma_sq)) * sum(residuals^2)
  return(nll)
}

# The optim function will be used to find the optiming phi1, phi2 and sigma for the model given the data.
initial_params <- c(0, 0, 0)
result <- optim(
  par = initial_params,
  fn = neg_log_likelihood,
  log_fin = log_fin,
  log_init = log_init,
  times = times,
  method = "BFGS",
  hessian = TRUE
)

# Final estimates and printing the outcomes
phi1_mle <- result$par[1]
phi2_mle <- result$par[2]
sigma_sq_mle <- exp(result$par[3])  # Exponentiate to get variance
sigma_mle <- sqrt(sigma_sq_mle)     # Standard deviation

cat("MLE Estimates:\n",
    "φ1 =", phi1_mle, "\n",
    "φ2 =", phi2_mle, "\n",
    "σ =", sigma_mle, "\n")

# Calculate predictions, for later usage
predicted_log_fin <- log_init + phi1_mle + phi2_mle * times
predicted_fin <- exp(predicted_log_fin)

######################################### Calculate 95% confidence intervals
se <- sigma_mle  # Standard error of the residuals
lower_bound <- exp(predicted_log_fin - 1.96 * se)  # Lower bound (original scale)
upper_bound <- exp(predicted_log_fin + 1.96 * se)  # Upper bound (original scale)

# Creating needed plots
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))

# Observed vs Predicted final counts, including a line of the best fit
plot(predicted_fin, data$fin, 
     xlab = "Predicted Final Count", ylab = "Observed Final Count",
     main = "Observed vs Predicted Counts",
     pch = 19, col = "#1f77b4", log="xy" )
abline(a = 0, b = 1, col = "#ff7f0e", lwd = 2)

########################################## Add 95% confidence intervals
arrows(predicted_fin, lower_bound, predicted_fin, upper_bound, 
       length = 0.05, angle = 90, code = 3, col = "darkgray")

############################################### Add legend
legend("topleft", 
       legend = c("Predicted", "1:1 Line", "95% CI"), 
       col = c("#1f77b4", "#ff7f0e", "darkgray"), 
       pch = c(19, NA, NA), 
       lty = c(NA, 1, 1))

# Residuals vs Fitted, their log values
residuals_log <- log_fin - predicted_log_fin
plot(predicted_log_fin, residuals_log,
     xlab = "Fitted Values (log)", ylab = "Residuals (log)",
     main = "Residuals vs Fitted",
     pch = 19, col = "#2ca02c")
abline(h = 0, col = "#d62728", lty = 2, lwd = 2)

#Comments about the model


# The estimated variance of the residuals is around 0.09, so 6/70 of the residuals of the lower
# observed counts are outliers as they differ by at least twice the variance from 0, which could suggest in most
# test this model may not be so adequate when discribing models of lower count, but when it comes
# to higher predictions, it is mostly accurate.


