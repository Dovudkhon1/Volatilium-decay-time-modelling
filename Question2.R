# Extract standard error of phi1 from the inverse Hessian matrix
se_phi1 <- sqrt(solve(result$hessian)[1, 1])  # First diagonal element is variance of phi1

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
