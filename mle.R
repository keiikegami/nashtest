robust_norm_censored_4 <- read_csv("~/Desktop/2017_Aterm/卒論/code/robust_norm_censored_4.csv")
data <- robust_norm_censored_4
library(stats4)
library(bbmle)
# mle for BR model
# theta = [beta1, beta2, delta]
br1 <- function(theta){
  xb1 <- data[["Pop"]]*theta[1]+data[["Dist1"]]*theta[2]
  xb2 <- data[["Pop"]]*theta[1]+data[["Dist2"]]*theta[2]
  logl <- data[["num0"]]*(log(pnorm(-xb1)*pnorm(-xb2)) - log(1-pnorm(-xb1)*pnorm(-xb2)-pnorm(xb1+theta[3])*pnorm(xb2+theta[3])))
  + data[["num2"]]*(log(pnorm(xb1+theta[3])*pnorm(xb2+theta[3])) - log(1-pnorm(-xb1)*pnorm(-xb2)-pnorm(xb1+theta[3])*pnorm(xb2+theta[3])))
  + 1000*log(1-pnorm(-xb1)*pnorm(-xb2)-pnorm(xb1+theta[3])*pnorm(xb2+theta[3]))
  return(-sum(logl))
}

br2 <- function(beta1, beta2, delta){
  xb1 <- data[["Pop"]]*beta1+data[["Dist1"]]*beta2
  xb2 <- data[["Pop"]]*beta1+data[["Dist2"]]*beta2
  logl <- data[["num0"]]*(log(pnorm(-xb1)*pnorm(-xb2)) - log(1-pnorm(-xb1)*pnorm(-xb2)-pnorm(xb1+delta)*pnorm(xb2+delta)))
  + data[["num2"]]*(log(pnorm(xb1+delta)*pnorm(xb2+delta)) - log(1-pnorm(-xb1)*pnorm(-xb2)-pnorm(xb1+delta)*pnorm(xb2+delta)))
  + 1000*log(1-pnorm(-xb1)*pnorm(-xb2)-pnorm(xb1+delta)*pnorm(xb2+delta))
  return(-sum(logl))
}

# optimization
oppai <- optim(c(1,-1,-0.5),br1,hessian = F, method = "Nelder-Mead")
result2 <- mle2(br2, start = list(beta1 = 1, beta2 = -1, delta = -0.5),method = "Nelder-Mead")
