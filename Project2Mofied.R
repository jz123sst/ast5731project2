# Import required libraries
# Load necessary libraries
library(ggplot2)
library(MASS)
library(dplyr)
library(coda)
library(mcmcse)

set.seed(5731)

# Read the data
df <- read.table('C:/Users/epicn/Downloads/SCPUnion2.1_mu_vs_z.txt', comment.char='#', header=FALSE)
data <- df[, c(2, 3)]
data[, 2] <- 10^((data[, 2] + 5) / 5 - 6) # convert distance modulus to distance [Mpc]

# Filter data with z < 0.1
local_data <- data[data[, 1] < 0.1, ]

# Scatter plot
ggplot(local_data, aes(x = local_data[, 1], y = local_data[, 2])) +
  geom_point() +
  labs(x = 'Redshift', y = 'Distance [Mpc]', title = 'Data Visualization')

# Linear regression model
X <- cbind(1, local_data[, 1])
Y <- local_data[, 2]
n <- nrow(X)
p <- ncol(X)

# MLE
beta_MLE <- solve(t(X) %*% X) %*% t(X) %*% Y
SSE <- sum((Y - X %*% beta_MLE)^2)
lambda_MLE <- n / SSE
MSE <- SSE / (n - p)



cat(sprintf('MLE for β: %f Mpc\n', beta_MLE))
cat(sprintf('MLE for λ: %f Mpc^{-2}\n', lambda_MLE))
cat(sprintf('MLE for H_0: %f km/s/Mpc\n', 3e5 / beta_MLE[2]))

# Calculate residuals
residuals <- Y - X %*% beta_MLE

# Create a data frame for residuals
residuals_df <- data.frame(residuals = as.vector(residuals))

# Generate QQ Plot using ggplot2
ggplot(residuals_df, aes(sample = residuals)) +
  stat_qq(color = "blue") +
  stat_qq_line(color = "red") +
  labs(
    title = "QQ Plot of Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()

# Gibbs sampler
n_samples <- 1e5

beta_mean <- t(t(c(0, 4283)))

# Prior parameters
m <- beta_mean
c <- c(1e-1, 1e-3)
a <- 1e-3 
b <- 1e-3

# Initial values
beta <- beta_MLE
lambda_ <- lambda_MLE

shape <- a + n / 2


beta_samples <- matrix(0, n_samples, p)
lambda_samples <- numeric(n_samples)

for (i in 1:n_samples) {
  # Construct the prior precision matrix using diag(c)
  prior_precision <- diag(c)  # Diagonal matrix with elements of c on the diagonal
  
  
  # beta | lambda, Y, X
  cov <- solve(prior_precision * diag(p) + lambda_ * t(X) %*% X)
  mean <- cov %*% (lambda_ * t(X) %*% Y + c * m)
  beta <- mvrnorm(1, mu = mean, Sigma = cov)
  
  # lambda | beta, Y, X
  scale <- 2 / (2 * b + t(Y - X %*% beta) %*% (Y - X %*% beta))
  lambda_ <- rgamma(1, shape = shape, scale = scale)
  
  beta_samples[i, ] <- beta
  lambda_samples[i] <- lambda_
}
# Define parameter names for Model 1
param_names <- c("beta_0", "beta_1")

# Loop over each parameter to calculate CI and generate plots
for (i in 1:p) {
  # Calculate 95% credible intervals
  beta_CI <- quantile(beta_samples[, i], probs = c(0.025, 0.975))
  
  # Print posterior mean and CI
  cat(sprintf('Posterior mean for %s: %f\n', param_names[i], beta_mean[i]))
  cat(sprintf('95%% CI for %s: [%f, %f]\n\n', param_names[i], beta_CI[1], beta_CI[2]))
  
  # Create histogram with density and credible interval lines
  plot <- ggplot(data.frame(beta_samples = beta_samples[, i]), aes(x = beta_samples)) +
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "white", color = "black") +
    geom_vline(xintercept = beta_mean[i], color = 'red', linetype = 'dashed', size = 1, 
               label = "Mean") +
    geom_vline(xintercept = beta_CI[1], color = 'blue', linetype = 'dashed', size = 1, 
               label = "2.5% CI") +
    geom_vline(xintercept = beta_CI[2], color = 'blue', linetype = 'dashed', size = 1, 
               label = "97.5% CI") +
    labs(
      title = sprintf('Marginal Posterior of %s', param_names[i]),
      x = sprintf('%s (Mpc)', param_names[i]),
      y = 'Density'
    ) +
    theme_minimal()
  
  # Print the plot
  print(plot)
}

# Results for lambda
lambda_mean <- mean(lambda_samples)
lambda_CI <- quantile(lambda_samples, c(0.025, 0.975))
cat(sprintf('Posterior mean: %f\n', lambda_mean))

ggplot(data.frame(lambda_samples = lambda_samples), aes(x = lambda_samples)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "white", color = "black") +
  geom_vline(xintercept = lambda_mean, color = 'red', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = lambda_CI[1], color = 'blue', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = lambda_CI[2], color = 'blue', linetype = 'dashed', size = 1) +
  labs(title = 'Marginal Posterior of λ', x = 'λ (Mpc^{-2})', y = 'Density') +
  theme_minimal()

# Estimate H0
H0_samples <- 3e5 / beta_samples[, 2]

H0_mean <- mean(H0_samples)
H0_CI <- quantile(H0_samples, c(0.025, 0.975))
cat(sprintf('Posterior mean: %f km/s/Mpc\n', H0_mean))
cat(sprintf('95%% credible interval: [%f, %f] km/s/Mpc\n', H0_CI[1], H0_CI[2]))

ggplot(data.frame(H0_samples = H0_samples), aes(x = H0_samples)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "white", color = "black") +
  geom_vline(xintercept = H0_mean, color = 'red', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = H0_CI[1], color = 'blue', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = H0_CI[2], color = 'blue', linetype = 'dashed', size = 1) +
  labs(title = 'Posterior of H0', x = 'H0 (km/s/Mpc)', y = 'Density') +
  theme_minimal()

# Convert samples to mcmc objects
model1_mcmc_beta <- mcmc(beta_samples)
model1_mcmc_lambda <- mcmc(lambda_samples)

#Find ESS
model1_mcmc <- mcmc(cbind(beta_samples, lambda_samples))
effectiveSize(model1_mcmc)
minESS(3)

# Select the last 5000 iterations for trace plots
last_5000_beta <- tail(model1_mcmc_beta, 5000)
last_5000_lambda <- tail(model1_mcmc_lambda, 5000)

# Individual trace plots for beta parameters
for (i in 1:p) {
  traceplot(last_5000_beta[, i], 
            main = paste("Trace Plot for beta_", i-1),
            xlab = "Iteration",
            ylab = expression(beta))
}

# Trace plot for lambda parameter
traceplot(last_5000_lambda, 
          main = "Trace Plot for lambda",
          xlab = "Iteration",
          ylab = expression(lambda))

#individual acf plots
for (i in 1:p) {
  acf_plot <- autocorr.plot(model1_mcmc_beta[, i], main = paste("ACF for beta_", i-1))
}

acf_plot <- autocorr.plot(model1_mcmc_lambda, main = "ACF for lambda")

# Number of posterior predictive samples (you can adjust this as needed)
n_pp_samples <- 5000
indices <- sample(1:n_samples, n_pp_samples)

# Initialize matrix to store predictive samples
Y_pred_samples <- matrix(0, nrow = nrow(X), ncol = n_pp_samples)

for (k in 1:n_pp_samples) {
  s <- indices[k]
  beta_s <- beta_samples[s, ]
  lambda_s <- lambda_samples[s]
  
  # Predicted mean
  mu_pred <- X %*% beta_s
  
  # Generate predictive Y from the posterior predictive distribution
  Y_pred_samples[, k] <- rnorm(nrow(X), mean = mu_pred, sd = sqrt(1 / lambda_s))
}

# Predictive mean and intervals
Y_pred_mean <- rowMeans(Y_pred_samples)
Y_pred_CI <- apply(Y_pred_samples, 1, quantile, probs = c(0.025, 0.975))

# Actual Y values
Y_actual <- Y

# Create data frame for plotting
plot_data <- data.frame(
  Y_actual = Y_actual,
  Y_pred_mean = Y_pred_mean,
  Y_pred_lower = Y_pred_CI[1, ],
  Y_pred_upper = Y_pred_CI[2, ]
)

# Scatterplot of actual Y vs predicted Y
ggplot(plot_data, aes(x = Y_actual, y = Y_pred_mean)) +
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymin = Y_pred_lower, ymax = Y_pred_upper), width = 0.2, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
  labs(
    title = 'Posterior Predictive Check: Actual vs Predicted Distannce (Mpc) (Model 1)',
    x = 'Actual Distance (Mpc)',
    y = 'Predicted Distance (Mpc)'
  ) +
  theme_minimal()



# Model 2: No intercept term
X <- as.matrix(local_data[, 1, drop = FALSE])
Y <- as.matrix(local_data[, 2, drop = FALSE])
n <- nrow(X)
p <- ncol(X)

# MLE
beta_MLE <- solve(t(X) %*% X) %*% t(X) %*% Y
SSE <- sum((Y - X %*% beta_MLE)^2)
lambda_MLE <- n / SSE
MSE <- SSE / (n - p)

cat(sprintf('MLE for β: %f Mpc\n', beta_MLE[1]))
cat(sprintf('MLE for λ: %f Mpc^{-2}\n', lambda_MLE))
cat(sprintf('MLE for H_0: %f km/s/Mpc\n', 3e5 / beta_MLE[1]))

# Gibbs sampler
n_samples <- 1e5

# Prior parameters
m <- 4283
c <- 1e-3
a <- 1e-3
b <- 1e-3

# Initial values
beta <- beta_MLE
lambda_ <- lambda_MLE

shape <- a + n / 2

beta_samples <- numeric(n_samples)
lambda_samples <- numeric(n_samples)

for (i in 1:n_samples) {
  # beta | lambda, Y, X
  cov <- solve(c + lambda_ * t(X) %*% X)
  mean <- cov %*% (lambda_ * t(X) %*% Y + c * m)
  beta <- mvrnorm(1, mu = mean, Sigma = cov)
  
  # lambda | beta, Y, X
  scale <- 2 / (2 * b + t(Y - X %*% beta) %*% (Y - X %*% beta))
  lambda_ <- rgamma(1, shape = shape, scale = scale)
  
  beta_samples[i] <- beta
  lambda_samples[i] <- lambda_
}

# Results for beta
beta_mean <- mean(beta_samples)
beta_CI <- quantile(beta_samples, c(0.025, 0.975))
cat(sprintf('Posterior mean: %f\n', beta_mean))

ggplot(data.frame(beta_samples = beta_samples), aes(x = beta_samples)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "white", color = "black") +
  geom_vline(xintercept = beta_mean, color = 'red', linetype = 'dashed', size=1) +
  geom_vline(xintercept = beta_CI[1], color = 'blue', linetype = 'dashed', size=1) +
  geom_vline(xintercept = beta_CI[2], color = 'blue', linetype = 'dashed', size=1) +
  labs(title = 'Marginal Posterior of β', x = 'β (Mpc)', y = 'Density') +
  theme_minimal()

# Results for lambda
lambda_mean <- mean(lambda_samples)
lambda_CI <- quantile(lambda_samples, c(0.025, 0.975))
cat(sprintf('Posterior mean: %f\n', lambda_mean))

ggplot(data.frame(lambda_samples = lambda_samples), aes(x = lambda_samples)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "white", color = "black") +
  geom_vline(xintercept = lambda_mean, color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = lambda_CI[1], color = 'blue', linetype = 'dashed', size=1) +
  geom_vline(xintercept = lambda_CI[2], color = 'blue', linetype = 'dashed', size=1) +
  labs(title = 'Marginal Posterior of λ', x = 'λ (Mpc^{-2})', y = 'Density',size=1) +
  theme_minimal()

# Estimate H0
H0_samples <- 3e5 / beta_samples

H0_mean <- mean(H0_samples)
H0_CI <- quantile(H0_samples, c(0.025, 0.975))
cat(sprintf('Posterior mean: %f km/s/Mpc\n', H0_mean))
cat(sprintf('95%% credible interval: [%f, %f] km/s/Mpc\n', H0_CI[1], H0_CI[2]))

ggplot(data.frame(H0_samples = H0_samples), aes(x = H0_samples)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "white", color = "black") +
  geom_vline(xintercept = H0_mean, color = 'red', linetype = 'dashed', size=1) +
  geom_vline(xintercept = H0_CI[1], color = 'blue', linetype = 'dashed', size=1) +
  geom_vline(xintercept = H0_CI[2], color = 'blue', linetype = 'dashed', size=1) +
  labs(title = 'Posterior of H0', x = 'H0 (km/s/Mpc)', y = 'Density', size=1) +
  theme_minimal() 



# Convert samples to mcmc objects
model2_mcmc_beta <- mcmc(beta_samples)
model2_mcmc_lambda <- mcmc(lambda_samples)

#Find ESS
model1_mcmc <- mcmc(cbind(beta_samples, lambda_samples))
effectiveSize(model1_mcmc)
minESS(2)


# Select the last 5000 iterations for trace plots
last_5000_beta <- tail(model1_mcmc_beta, 5000)
last_5000_lambda <- tail(model1_mcmc_lambda, 5000)

# Individual trace plots for beta parameters
for (i in 1:p) {
  traceplot(last_5000_beta[, i], 
            main = paste("Trace Plot for beta"),
            xlab = "Iteration",
            ylab = expression(beta))
}

# Trace plot for lambda parameter
traceplot(last_5000_lambda, 
          main = "Trace Plot for lambda",
          xlab = "Iteration",
          ylab = expression(lambda))


#  individual ACF plots
autocorr.plot(model2_mcmc_beta, main = "ACF for beta")
autocorr.plot(model2_mcmc_lambda, main = "ACF for lambda")

# Number of posterior predictive samples (you can adjust this as needed)
n_pp_samples <- 5000
indices <- sample(1:n_samples, n_pp_samples)

# Initialize matrix to store predictive samples
Y_pred_samples <- matrix(0, nrow = nrow(X), ncol = n_pp_samples)

for (k in 1:n_pp_samples) {
  s <- indices[k]
  beta_s <- beta_samples[s]
  lambda_s <- lambda_samples[s]
  
  # Predicted mean
  mu_pred <- X * beta_s
  
  # Generate predictive Y from the posterior predictive distribution
  Y_pred_samples[, k] <- rnorm(nrow(X), mean = mu_pred, sd = sqrt(1 / lambda_s))
}

# Predictive mean and intervals
Y_pred_mean <- rowMeans(Y_pred_samples)
Y_pred_CI <- apply(Y_pred_samples, 1, quantile, probs = c(0.025, 0.975))

# Actual Y values
Y_actual <- Y

# Create data frame for plotting
plot_data <- data.frame(
  Y_actual = Y_actual,
  Y_pred_mean = Y_pred_mean,
  Y_pred_lower = Y_pred_CI[1, ],
  Y_pred_upper = Y_pred_CI[2, ]
)

# Scatterplot of actual Y vs predicted Y
ggplot(plot_data, aes(x = Y_actual, y = Y_pred_mean)) +
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymin = Y_pred_lower, ymax = Y_pred_upper), width = 0.2, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
  labs(
    title = 'Posterior Predictive Check: Actual vs Predicted Distance (Mpc) (Model 2)',
    x = 'Actual Distance (Mpc)',
    y = 'Predicted Distance (Mpc)'
  ) +
  theme_minimal()
