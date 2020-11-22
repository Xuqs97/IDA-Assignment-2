#Question 2b
log_like_norm = function(dataex2,mu) {
  r = dataex2[,2]
  y = dataex2[,1]
  sum(r*dnorm(y, mean = mu, sd = 1.5, log = TRUE) + (1 - r) * pnorm(y,mean = mu, sd = 1.5, log.p = TRUE))
}
require(maxLik)
mle = maxLik(logLik = log_like_norm, data=dataex2,start = (mu = 10))
summary(mle)
#Question 4
observed = dataex4[complete.cases(dataex4),]
missed = dataex4[!complete.cases(dataex4),]
Q_beta_t = function(beta) {
  a = sum(observed$Y * (beta[1] + observed$X * beta[2]))
  b = sum((beta[1] + missed$X * beta[2]) * exp(beta_old[1] + missed$X * beta_old[2])/(1 + exp(beta_old[1] + missed$X * beta_old[2])))
  c = sum(log(1+exp(beta[1] + dataex4$X * beta[2])))
  return(a + b - c)
}
diff = 1
eps = 1e-8
beta = c(1,1)
while (diff > eps) {
  beta_old = beta
  #M step
  beta = coef(maxLik(logLik = Q_beta_t, start = c(0,0)))
  diff = max(abs(beta - beta_old))
}
beta
#Question 5b
theta_mle = function(p0,mu0,sigma0squared,lamda0) {
  diff =1
  eps = 1e-8
  p = p0
  mu = mu0
  sigmasquared = sigma0squared
  lamda = lamda0
  n = length(dataex5)
  while (diff > eps) {
    p_old = p
    mu_old = mu
    sigmasquared_old = sigmasquared
    lamda_old = lamda
    #E step
    pt = p * dlnorm(dataex5,mu_old,sqrt(sigmasquared_old))/(p * dlnorm(dataex5,mu_old,sqrt(sigmasquared_old)) + (1 - p) * dexp(dataex5,lamda_old))
    #M step
    p = sum(pt)/n
    mu = sum(pt * log(dataex5))/sum(pt)
    sigmasquared = sum(pt * (log(dataex5) - mu)**2)/sum(pt)
    lamda = sum(1 - pt)/sum((1 - pt) * dataex5)
    diff_p = abs(p - p_old)
    diff_mu = abs(mu - mu_old)
    diff_sigmasquared = abs(sigmasquared - sigmasquared_old)
    diff_lamda = abs(lamda - lamda_old)
    diff = max(diff_p, diff_mu, diff_sigmasquared, diff_lamda)
  }
  return(c(p, mu, sigmasquared, lamda))
}
theta = theta_mle(p0 = 0.1, mu0 = 1, sigma0squared = 0.25, lamda0 = 2)
p = theta[1]
mu = theta[2]
sigmasquared = theta[3]
lamda = theta[4]

f_mixture = function(data) {
  p * dlnorm(data, mu, sqrt(sigmasquared)) + (1-p) * dexp(data, lamda)
}

library(ggplot2)
data_ex5 = data.frame(dataex5)
ggplot(
  data_ex5,
  aes(
    x = dataex5,
  )
) + 
  geom_histogram(aes(y = stat(density)),bins = 100) +
  geom_function(fun = f_mixture) + 
  labs(title = "Histogram of the data and the estimated density")