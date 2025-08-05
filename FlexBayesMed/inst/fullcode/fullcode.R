# Here is the code to reproduce the results in
# Modeling Skewed and Heavy-Tailed Errors in Bayesian Mediation Analysis Using Standard Improper Priors
# The first part "MATPT (Section 2)" gives the code to reproduce the figures in Section 2
# The second part "Case Study (Section 6)" gives the code to reproduce all the results
# in Section 6 Case Study
# The third part gives a sample code to reproduce the results in Section 5 Simulation studies.
# One can make some minor changes to the sample code to reproduce all the results of the simulation studies.



##################################
###      MATPT (Section 2)     ###
##################################

###### Plot the densities of the centred two-piece distribution ######
######           This corresponds to Figure 1           ######
library(FlexBayesMed)

xseq <- seq(-6, 6, length.out = 2001)

# RIght skewed
plot(NA, xlim = range(xseq), ylim = c(0, 0.42),
     xlab = "x", ylab = "Density",
     main = "")
xseq   <- seq(-6, 6, length.out = 2001)
lines(xseq, dctpt(xseq,gamma=1.5,nu=3), col = "red", lwd = 2.5, lty=5)
lines(xseq, dctpt(xseq,gamma=3,nu=5), col = "darkgreen", lwd = 2.5, lty=2)
lines(xseq, dctpt(xseq,gamma=1.5, nu="normal"), col="royalblue",lwd=2.5,lty=3)
lines(xseq, dctpt(xseq,gamma=3, nu="normal"), col="orange",lwd=2.5,lty=4)
lines(xseq, dnorm(xseq), col="black",lwd=2,lty=1)

legend("topright", inset = 0.04, cex   = 1.14,
       legend = c(expression(paste(gamma, " = 1.5, ", nu, " = 3")),
                  expression(paste(gamma, " = 3,   ", nu, " = 5")),
                  expression(paste(gamma, " = 1.5, ", nu, " = ", infinity)),
                  expression(paste(gamma, " = 3,   ", nu, " = ", infinity)),
                  expression(N(0,1))),
       col = c("red", "darkgreen", "royalblue", "orange", "black"),
       lty = c(5,2,3,4,1),
       lwd = 1,
       bty = "n")

# Left Skewed
plot(NA, xlim = range(xseq), ylim = c(0, 0.42),
     xlab = "x", ylab = "Density",
     main = "")
xseq   <- seq(-6, 6, length.out = 2001)
lines(xseq, dctpt(xseq,gamma=2/3,nu=3), col = "red", lwd = 2.5, lty=5)
lines(xseq, dctpt(xseq,gamma=1/3,nu=5), col = "darkgreen", lwd = 2.5, lty=2)
lines(xseq, dctpt(xseq,gamma=2/3, nu="normal"), col="royalblue",lwd=2.5,lty=3)
lines(xseq, dctpt(xseq,gamma=1/3, nu="normal"), col="orange",lwd=2.5,lty=4)
lines(xseq, dnorm(xseq), col="black",lwd=2,lty=1)

legend("topleft", inset = 0.001, cex   = 1.14,
       legend = c(expression(paste(gamma, " = 2/3, ", nu, " = 3")),
                  expression(paste(gamma, " = 1/3,   ", nu, " = 5")),
                  expression(paste(gamma, " = 2/3, ", nu, " = ", infinity)),
                  expression(paste(gamma, " = 1/3,   ", nu, " = ", infinity)),
                  expression(N(0,1))),
       col = c("red", "darkgreen", "royalblue", "orange", "black"),
       lty = c(5,2,3,4,1),
       lwd = 1,
       bty = "n")


######          Plot the Skewness of the CTPT          ######
######           This corresponds to Figure 2           ######
library(FlexBayesMed)

# Fisher
gamma_seq <- seq(0.25, 3, length.out = 200)

## Skewness for nu=5
sk_nu5   <- sapply(gamma_seq, function(g) skewness_ctpt(g, 5))

## Skewness for nu=10
sk_nu10  <- sapply(gamma_seq, function(g) skewness_ctpt(g, 10))

## Skewness for two-piece normal
sk_infty <- sapply(gamma_seq, function(g) skewness_ctpt(g, "normal"))

plot(gamma_seq, sk_nu5,
     type = "l",
     col = "blue", lwd = 2.5, lty=1,
     ylim = range(c(sk_nu5, sk_nu10, sk_infty)),
     xlab = expression(gamma),
     ylab = expression(paste("SK(X|",gamma,",",nu,")")),
     main = "")

lines(gamma_seq, sk_nu10, col = "red",  lty=2, lwd = 2.5)
lines(gamma_seq, sk_infty, col = "darkgreen", lty=3, lwd = 2.5)

legend("bottomright",cex=1.2, inset=0.04,
       legend = c(expression(nu==5), expression(nu==10), expression(nu==infinity)),
       col = c("blue", "red", "darkgreen"), lty = c(1,2,3), lwd = 2.5,
       bty = "n")

##################################
###   Case Study (Section 6)   ###
##################################

library(e1071)
# Load the data
data <- read.table(file="C:/Users/Li/Downloads/rdata.txt", header=TRUE)

# Check the OLS residual
x <- data$ME
m <- data$HE
y <- data$math

reg1 <- lm(m~x)
summary(reg1)
kurtosis(reg1$residuals)
skewness(reg1$residuals)
qqnorm(reg1$residuals)
qqline(reg1$residuals,col="red")

reg2 <- lm(y~x+m)
summary(reg2)
kurtosis(reg2$residuals)
skewness(reg2$residuals)
qqnorm(reg2$residuals)
qqline(reg2$residuals,col="red")


library(FlexBayesMed)

# FIt the Full Model
fit_full <- flex_med(X=x, M=m, Y=y, model="full", iter = 30000, refresh = 3000,
                     warmup = 6000, chains = 1, seed = 123)

fit_full_H1_alpha <- fit_full$X2M_with_alpha
fit_full_H1_beta <- fit_full$XM2Y_with_beta
fit_full_H0_alpha <- fit_full$X2M_without_alpha
fit_full_H0_beta <- fit_full$XM2Y_without_beta

# FIt the Gamma-Only Model
fit_gamma <- flex_med(X=x, M=m, Y=y, model="gamma-only", iter = 30000, refresh = 3000,
                      warmup = 6000, chains = 1, seed = 123)

fit_gamma_H1_alpha <- fit_gamma$X2M_with_alpha
fit_gamma_H1_beta <- fit_gamma$XM2Y_with_beta
fit_gamma_H0_alpha <- fit_gamma$X2M_without_alpha
fit_gamma_H0_beta <- fit_gamma$XM2Y_without_beta

# FIt the Nu-Only Model
fit_nu <- flex_med(X=x, M=m, Y=y, model="nu-only", iter = 30000, refresh = 3000,
                   warmup = 6000, chains = 1, seed = 123)

fit_nu_H1_alpha <- fit_nu$X2M_with_alpha
fit_nu_H1_beta <- fit_nu$XM2Y_with_beta
fit_nu_H0_alpha <- fit_nu$X2M_without_alpha
fit_nu_H0_beta <- fit_nu$XM2Y_without_beta

# FIt the Normal Model
fit_normal <- flex_med(X=x, M=m, Y=y, model="normal", iter = 30000, refresh = 3000,
                       warmup = 6000, chains = 1, seed = 123)

fit_normal_H1_alpha <- fit_normal$X2M_with_alpha
fit_normal_H1_beta <- fit_normal$XM2Y_with_beta
fit_normal_H0_alpha <- fit_normal$X2M_without_alpha
fit_normal_H0_beta <- fit_normal$XM2Y_without_beta


# Compare four models for the X to M path
set.seed(123)
bridge_full_H1_alpha <- bridge_sampler(fit_full_H1_alpha)
bridge_gamma_H1_alpha <- bridge_sampler(fit_gamma_H1_alpha)
bridge_nu_H1_alpha <- bridge_sampler(fit_nu_H1_alpha)
bridge_normal_H1_alpha <- bridge_sampler(fit_normal_H1_alpha)

bf_full_gamma <- as.numeric(bf(bridge_full_H1_alpha, bridge_gamma_H1_alpha))[1]
log(bf_full_gamma)
bf_full_nu <- as.numeric(bf(bridge_full_H1_alpha, bridge_nu_H1_alpha))[1]
log(bf_full_nu)
bf_full_normal <- as.numeric(bf(bridge_full_H1_alpha, bridge_normal_H1_alpha))[1]
log(bf_full_normal)

bf_gamma_nu <- as.numeric(bf(bridge_gamma_H1_alpha, bridge_nu_H1_alpha))[1]
log(bf_gamma_nu)
bf_gamma_normal <- as.numeric(bf(bridge_gamma_H1_alpha, bridge_normal_H1_alpha))[1]
log(bf_gamma_normal)

bf_nu_normal <- as.numeric(bf(bridge_nu_H1_alpha, bridge_normal_H1_alpha))[1]
log(bf_nu_normal)

# Check the output of the Full Model for path: X to M
print(fit_full_H1_alpha,digits_summary = 3)
print(fit_gamma_H1_alpha,digits_summary = 3) # Compare with the gamma-only model
samples_full_H1_alpha <- extract(fit_full_H1_alpha)
samples_gamma_H1_alpha <- extract(fit_gamma_H1_alpha)
samples_nu_H1_alpha <- extract(fit_nu_H1_alpha)
samples_normal_H1_alpha <- extract(fit_normal_H1_alpha)

alpha_full <- samples_full_H1_alpha$alpha
alpha_gamma <- samples_gamma_H1_alpha$alpha
alpha_nu <- samples_nu_H1_alpha$alpha
alpha_normal <- samples_normal_H1_alpha$alpha

gamma_full <- samples_full_H1_alpha$gamma
nu_full <- samples_full_H1_alpha$nu
sigma_full <- samples_full_H1_alpha$sigma

density_est_gamma <- density(gamma_full)
plot(density_est_gamma, main="", xlab=expression(gamma), ylab="Density")
polygon(density_est_gamma, col="grey", border="black",lty=1,lwd=2)

density_est_nu <- density(nu_full,from = 2)
plot(density_est_nu, main="", xlab=expression(nu), ylab="Density")
polygon(density_est_nu, col="grey", border="black",lty=1,lwd=2)

density_est_sigma <- density(sigma_full)
plot(density_est_sigma, main="", xlab=expression(sigma), ylab="Density")
polygon(density_est_sigma, col="grey", border="black",lty=1,lwd=2)

density_est_alpha_full <- density(alpha_full)
density_est_alpha_gamma <- density(alpha_gamma)
density_est_alpha_nu <- density(alpha_nu)
density_est_alpha_normal <- density(alpha_normal)

# Plot the estimated posterior density of alpha
plot(density_est_alpha_full, main="", xlab=expression(alpha), ylab="Density")
lines(density_est_alpha_full, col="red",lty=1,lwd=4)
lines(density_est_alpha_gamma, col="blue",lty=2,lwd=4)
lines(density_est_alpha_nu, col="green4",lty=3,lwd=4)
lines(density_est_alpha_normal, col="black",lty=4,lwd=4)
legend("topleft", legend = c("Full Model", expression(paste(gamma,"-Only Model")),
                              expression(paste(nu,"-Only Model")), "Normal Model"),
       col = c("red","blue","green4","black"),
       lty = c(1, 2, 5, 4), lwd = 3, cex = 1.15)

# Plot the estimated posterior of gamma with its prior
norm_const <- pgamma(20, shape = 2, rate = 2) - pgamma(0.05, shape = 2, rate = 2)
prior_seq <- seq(0.051,19.99,length.out=2000)
prior_y <- dgamma(prior_seq, shape = 2, rate = 2) / norm_const

plot(density_est_gamma, main = "", xlab = expression(gamma), ylab = "Density", xlim=c(0.05,4))
lines(density_est_gamma, col = "red", lty = 1, lwd = 3)
lines(prior_seq, prior_y, col = "blue", lwd = 3, lty = 2)  # prior

legend("topright", legend = c("Posterior", "Prior"),col = c("red", "blue"),lty = c(1, 2),
       lwd = c(2, 2), bty = "n", cex=1.2)

# Plot the estimated posterior of nu with its prior
upper <- max(density_est_nu$x)
d <- 0.01
lower <- 2
prior_y <- ifelse(density_est_nu$x > lower, d*exp( -d*(density_est_nu$x-lower) ),0)

plot(density_est_nu, main = "", xlab = expression(nu), ylab = "Density", xlim = c(lower, upper), ylim = c(0,0.01))

lines(density_est_nu, col = "red",  lty = 1, lwd = 3)
lines(density_est_nu$x, prior_y, col = "blue", lwd = 3, lty = 2)  # prior

legend("topright",legend = c("Posterior", "Prior"), cex=1.2, col= c("red", "blue"), lty = c(1,2), lwd = c(2, 2), bty = "n")


# Compare four models for the XM to Y path
set.seed(123)
bridge_full_H1_beta <- bridge_sampler(fit_full_H1_beta)
bridge_gamma_H1_beta <- bridge_sampler(fit_gamma_H1_beta)
bridge_nu_H1_beta <- bridge_sampler(fit_nu_H1_beta)
bridge_normal_H1_beta <- bridge_sampler(fit_normal_H1_beta)

bf_full_gamma <- as.numeric(bf(bridge_full_H1_beta, bridge_gamma_H1_beta))[1]
log(bf_full_gamma)
bf_full_nu <- as.numeric(bf(bridge_full_H1_beta, bridge_nu_H1_beta))[1]
log(bf_full_nu)
bf_full_normal <- as.numeric(bf(bridge_full_H1_beta, bridge_normal_H1_beta))[1]
log(bf_full_normal)

bf_gamma_nu <- as.numeric(bf(bridge_gamma_H1_beta, bridge_nu_H1_beta))[1]
log(bf_gamma_nu)
bf_gamma_normal <- as.numeric(bf(bridge_gamma_H1_beta, bridge_normal_H1_beta))[1]
log(bf_gamma_normal)

bf_nu_normal <- as.numeric(bf(bridge_nu_H1_beta, bridge_normal_H1_beta))[1]
log(bf_nu_normal)

# Check the output of the Full Model for path: XM to Y
print(fit_full_H1_beta,digits_summary = 3)
print(fit_nu_H1_beta,digits_summary = 3) # Compare with the nu-only model
samples_full_H1_beta <- extract(fit_full_H1_beta)
samples_nu_H1_beta <- extract(fit_nu_H1_beta)
samples_gamma_H1_beta <- extract(fit_gamma_H1_beta)
samples_normal_H1_beta <- extract(fit_normal_H1_beta)

beta_full <- samples_full_H1_beta$beta
beta_gamma <- samples_gamma_H1_beta$beta
beta_nu <- samples_nu_H1_beta$beta
beta_normal <- samples_normal_H1_beta$beta

gamma_full <- samples_full_H1_beta$gamma
nu_full <- samples_full_H1_beta$nu
sigma_full <- samples_full_H1_beta$sigma

density_est_gamma <- density(gamma_full)
plot(density_est_gamma, main="", xlab=expression(gamma), ylab="Density")
polygon(density_est_gamma, col="grey", border="black",lty=1,lwd=2)

density_est_nu <- density(nu_full)
plot(density_est_nu, main="", xlab=expression(nu), ylab="Density")
polygon(density_est_nu, col="grey", border="black",lty=1,lwd=2)

density_est_sigma <- density(sigma_full)
plot(density_est_sigma, main="", xlab=expression(sigma), ylab="Density")
polygon(density_est_sigma, col="grey", border="black",lty=1,lwd=2)

density_est_beta_full <- density(beta_full)
density_est_beta_nu <- density(beta_nu)
density_est_beta_gamma <- density(beta_gamma)
density_est_beta_normal <- density(beta_normal)

# Plot the estimated posterior density of beta
plot(density_est_beta_full, main="", xlab=expression(beta), ylab="Density")
lines(density_est_beta_full, col="red",lty=1,lwd=4)
lines(density_est_beta_gamma, col="blue",lty=2,lwd=4)
lines(density_est_beta_nu, col="green4",lty=3,lwd=4)
lines(density_est_beta_normal, col="black",lty=4,lwd=4)
legend("topright", legend = c("Full Model", expression(paste(gamma,"-Only Model")),
                             expression(paste(nu,"-Only Model")), "Normal Model"),
       col = c("red","blue","green4","black"),
       lty = c(1, 2, 3, 4), lwd = 3, cex = 1.15)

# Plot the estimated posterior of gamma with its prior
norm_const <- pgamma(20, shape = 2, rate = 2) - pgamma(0.05, shape = 2, rate = 2)
prior_seq <- seq(0.051,19.99,length.out=2000)
prior_y <- dgamma(prior_seq, shape = 2, rate = 2) / norm_const

plot(density_est_gamma, main = "", xlab = expression(gamma), ylab = "Density", xlim=c(0.05,4))
lines(density_est_gamma, col = "red", lty = 1, lwd = 3)
lines(prior_seq, prior_y, col = "blue", lwd = 3, lty = 2)  # prior

legend("topright", legend = c("Posterior", "Prior"),col = c("red", "blue"),lty = c(1, 2),
       lwd = c(2, 2), bty = "n", cex=1.2)

# Plot the estimated posterior of nu with its prior
upper <- max(density_est_nu$x)
d <- 0.01
lower <- 2
prior_y <- ifelse(density_est_nu$x > lower, d*exp( -d*(density_est_nu$x-lower) ),0)

plot(density_est_nu, main = "", xlab = expression(nu), ylab = "Density", xlim = c(lower, upper), ylim = c(0,0.3))

lines(density_est_nu, col = "red",  lty = 1, lwd = 3)
lines(density_est_nu$x, prior_y, col = "blue", lwd = 3, lty = 2)  # prior

legend("topright",legend = c("Posterior", "Prior"), cex=1.2, col= c("red", "blue"), lty = c(1,2), lwd = c(2, 2), bty = "n")



# Plot the estimated density of mediation effect
mean(fit_full$samples_med)
samples_med_sorted <- sort(fit_full$samples_med)
c(samples_med_sorted[600],samples_med_sorted[23400]) # 95% credible interval
density_est_med_full <- density(fit_full$samples_med)
density_est_med_gamma <- density(fit_gamma$samples_med)
density_est_med_nu <- density(fit_nu$samples_med)
density_est_med_normal <- density(fit_normal$samples_med)

# Plot the estimated posterior density of alpha
plot(density_est_med_full, main="", xlab=expression(paste(alpha,beta)), ylab="Density")
lines(density_est_med_full, col="red",lty=1,lwd=4)
lines(density_est_med_gamma, col="blue",lty=2,lwd=4)
lines(density_est_med_nu, col="green4",lty=3,lwd=4)
lines(density_est_med_normal, col="black",lty=4,lwd=4)
legend("topright", legend = c("Full Model", expression(paste(gamma,"-Only Model")),
                             expression(paste(nu,"-Only Model")), "Normal Model"),
       col = c("red","blue","green4","black"),
       lty = c(1, 2, 3, 4), lwd = 3, cex = 1.15)


# Bayes factor of the mediation effect
fit_full$bf_med



##################################
###   Simulation (Section 5)   ###
##################################

# Load the FlexBayesMed package
library(FlexBayesMed)
set.seed(123)

# The following five functions convert standardised coefficients into unstandardised coefficients
# Please refer to Section B in the Supplementary Materials
library(Rmpfr)
f <- function(ga,nu){
  nu_mpfr <- mpfr(nu, precBits=512)
  log_numerator <- log(2) + log(nu_mpfr) + lgamma((nu_mpfr + 1)/2)
  log_denominator <- 0.5 * (log(pi) + log(nu_mpfr)) + log(nu_mpfr - 1) + lgamma(nu_mpfr/2)
  log_ratio <- log_numerator - log_denominator
  return(as.numeric(exp(log_ratio))*(ga-1/ga))
}

H <- function(gamma, nu){
  return((nu/(nu-2)) * (gamma^2 - 1 + (1/(gamma^2))) - f(gamma, nu)^2)
}

coef_1 <- function(stand_alpha,gamma_1,nu_1){ # Convert unstandardised alpha
  H1 <- H(gamma_1,nu_1)
  s2 <- H1/((1/(stand_alpha^2))-1)
  return(sqrt(s2))
}

coef_2 <- function(stand_beta,stand_alpha,stand_tau,gamma_1,nu_1,gamma_2,nu_2){ # Convert unstandardised beta
  H1 <- H(gamma_1,nu_1)
  H2 <- H(gamma_2,nu_2)
  alpha <- coef_1(stand_alpha,gamma_1,nu_1)
  K <- alpha^2 + H1
  l1 <- K*(1-stand_tau^2)/(H2*(stand_beta^2))
  l2 <- K/H2
  s2 <- 1/(l1-l2)
  return(sqrt(s2))
}

coef_tau <- function(stand_beta,stand_tau,gamma_2,nu_2){ # Convert unstandardised tau
  H2 <- H(gamma_2,nu_2)
  tau <- stand_tau * sqrt(H2/(1 - (stand_tau)^2 - (stand_beta)^2))
  return(tau)
}

# Here we only give an example under one simulation setting
# One can change the following parameters for other settings easily
N <- c(50) ## Sample size for each dataset
stand_alpha <- c(0.14) ## Standardised coefficient alpha
stand_beta <- c(0.14) ## Standardised coefficient beta
stand_tau <- c(0.14) ## Standardised coefficient tau
gamma <- c(0.33)
nu <- c(3) # nu = 3, 10; to set normal tail, use "nu = "normal

# Loop parameters
exp_num <- 1000 # The number of simulated datasets


# Generate data under H_1: mediation effect exists (we fix sigma=1)
## Output directory
output_dir <- "../sample_simulation_data_H1"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

for (i in 1:length(N)) {
  size <- N[i]
  for (j in 1:length(stand_alpha)) {
    for (k in 1:length(stand_beta)) {
      for (l in 1:length(gamma)) {
        for (m in 1:length(nu)) {
          for (p in 1:length(stand_tau)){
            alpha_val <- coef_1(stand_alpha[j], gamma[l], nu[m])
            beta_val <- coef_2(stand_beta[k], stand_alpha[j], stand_tau[p], gamma[l], nu[m], gamma[l], nu[m])
            tau_val <- coef_tau(stand_beta[k], stand_tau[p], gamma[l], nu[m])

            X <- rnorm(size * exp_num)
            e1 <- rctpt(size * exp_num, gamma[l], nu[m]) # Generate CTPT errors for path X to M
            e2 <- rmatpt(size * exp_num, gamma[l], nu[m]) # Generate CTPT errors for path XM to Y

            M <- alpha_val * X + e1
            Y <- tau_val * X + beta_val * M + e2

            exp_data_X <- matrix(X, ncol = exp_num)
            exp_data_M <- matrix(M, ncol = exp_num)
            exp_data_Y <- matrix(Y, ncol = exp_num)

            file_name_X <- paste0("X_N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], ".csv")
            file_name_M <- paste0("M_N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], ".csv")
            file_name_Y <- paste0("Y_N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], ".csv")

            output_subdir <- paste0(output_dir,"/N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], "/")

            if (!dir.exists(output_subdir)) {
              dir.create(output_subdir)
            }

            write.csv(exp_data_X, file.path(output_subdir, file_name_X))
            write.csv(exp_data_M, file.path(output_subdir, file_name_M))
            write.csv(exp_data_Y, file.path(output_subdir, file_name_Y))
          }
        }
      }
    }
  }
}

# Generate data under H_0: mediation effect does not exist (we fix sigma=1)
## Output directory
output_dir <- "../sample_simulation_data_H0"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Nested loops to generate and save data
for (i in 1:length(N)) {
  size <- N[i]
  for (j in 1:length(stand_alpha)) {
    for (k in 1:length(stand_beta)) {
      for (l in 1:length(gamma)) {
        for (m in 1:length(nu)) {
          for (p in 1:length(stand_tau)){
            X <- rnorm(size * exp_num)
            M <- rep(NA, size * exp_num)
            Y <- rep(NA, size * exp_num)
            e1 <- rctpt(size * exp_num, gamma[l], nu[m]) # Generate CTPT errors for path X to M
            e2 <- rctpt(size * exp_num, gamma[l], nu[m]) # Generate CTPT errors for path XM to Y

            idx <- sample(c(1,2,3),size = size * exp_num, replace = TRUE) ## We assume equal conditional prior probability for the three cases under H0

            group_00 <- which(idx==1)

            alpha_val <- coef_1(0, gamma[l], nu[m])
            beta_val <- coef_2(0, 0, stand_tau[p], gamma[l], nu[m], gamma[l], nu[m])
            tau_val <- coef_tau(0, stand_tau[p], gamma[l], nu[m])

            M[group_00] <- alpha_val * X[group_00] + e1[group_00]
            Y[group_00] <- tau_val * X[group_00] + beta_val * M[group_00] + e2[group_00]

            group_01 <- which(idx==2)

            alpha_val <- coef_1(0, gamma[l], nu[m])
            beta_val <- coef_2(stand_beta[k], 0, stand_tau[p], gamma[l], nu[m], gamma[l], nu[m])
            tau_val <- coef_tau(stand_beta[k], stand_tau[p], gamma[l], nu[m])

            M[group_01] <- alpha_val * X[group_01] + e1[group_01]
            Y[group_01] <- tau_val * X[group_01] + beta_val * M[group_01] + e2[group_01]

            group_10 <- which(idx==3)

            alpha_val <- coef_1(stand_alpha[j], gamma[l], nu[m])
            beta_val <- coef_2(0, stand_alpha[j], stand_tau[p], gamma[l], nu[m], gamma[l], nu[m])
            tau_val <- coef_tau(0, stand_tau[p], gamma[l], nu[m])

            M[group_10] <- alpha_val * X[group_10] + e1[group_10]
            Y[group_10] <- tau_val * X[group_10] + beta_val * M[group_10] + e2[group_10]

            exp_data_X <- matrix(X, ncol = exp_num)
            exp_data_M <- matrix(M, ncol = exp_num)
            exp_data_Y <- matrix(Y, ncol = exp_num)

            # Create unique file names for each matrix
            file_name_X <- paste0("X_N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], ".csv")
            file_name_M <- paste0("M_N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], ".csv")
            file_name_Y <- paste0("Y_N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], ".csv")

            # Output directory
            output_subdir <- paste0(output_dir,"/N_", size, "_alpha_", stand_alpha[j], "_beta_", stand_beta[k], "_tau_", stand_tau[p], "_gamma_", gamma[l], "_nu_", nu[m], "/")

            # Create the output directory if it doesn't exist
            if (!dir.exists(output_subdir)) {
              dir.create(output_subdir)
            }

            # Save matrices as CSV files
            write.csv(exp_data_X, file.path(output_subdir, file_name_X))
            write.csv(exp_data_M, file.path(output_subdir, file_name_M))
            write.csv(exp_data_Y, file.path(output_subdir, file_name_Y))

          }
        }
      }
    }
  }
}


library(rstan)
library(bridgesampling)
X_matrix <- read.csv(file.path("../sample_simulation_data_H0/N_50_alpha_0.14_beta_0.14_tau_0.14_gamma_0.33_v_3/X_N_50_alpha_0.14_beta_0.14_tau_0.14_gamma_0.33_v_3.csv"),header=TRUE)
M_matrix <- read.csv(file.path("../sample_simulation_data_H0/N_50_alpha_0.14_beta_0.14_tau_0.14_gamma_0.33_v_3/M_N_50_alpha_0.14_beta_0.14_tau_0.14_gamma_0.33_v_3.csv"),header=TRUE)
Y_matrix <- read.csv(file.path("../sample_simulation_data_H0/N_50_alpha_0.14_beta_0.14_tau_0.14_gamma_0.33_v_3/Y_N_50_alpha_0.14_beta_0.14_tau_0.14_gamma_0.33_v_3.csv"),header=TRUE)

X_matrix <- as.matrix(X_matrix[1:50,2:1001])
M_matrix <- as.matrix(M_matrix[1:50,2:1001])
Y_matrix <- as.matrix(Y_matrix[1:50,2:1001])

# mode, mean, sd, 2.5, 25, 50, 75, 97.5
alpha <- matrix(rep(NA,1000*7),ncol=7)
beta <- matrix(rep(NA,1000*7),ncol=7)
med <- matrix(rep(NA,1000*7),ncol=7)
tau <- matrix(rep(NA,1000*7),ncol=7)
gamma_H1_alpha <- matrix(rep(NA,1000*7),ncol=7)
gamma_H1_beta <- matrix(rep(NA,1000*7),ncol=7)
v_H1_alpha <- matrix(rep(NA,1000*7),ncol=7)
v_H1_beta <- matrix(rep(NA,1000*7),ncol=7)
gamma_H0_alpha <- matrix(rep(NA,1000*7),ncol=7)
gamma_H0_beta <- matrix(rep(NA,1000*7),ncol=7)
sigma_alpha_H1 <- matrix(rep(NA,1000*7),ncol=7)
sigma_beta_H1 <- matrix(rep(NA,1000*7),ncol=7)
sigma_alpha_H0 <- matrix(rep(NA,1000*7),ncol=7)
sigma_beta_H0 <- matrix(rep(NA,1000*7),ncol=7)
v_H0_alpha <- matrix(rep(NA,1000*7),ncol=7)
v_H0_beta <- matrix(rep(NA,1000*7),ncol=7)
bf_alpha <- rep(NA,1000)
bf_beta <- rep(NA,1000)
bf_med <- rep(NA,1000)



