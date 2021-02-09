# hbic_int function

#library(tidyverse)
#library(devtools)
#install_github('seonjoo/smm')
#library(smm)

# Data

#data = read_csv("data.csv") %>%
#  dplyr::select(-np_speed_attention, -np_memory, -np_vocab) %>%
#  na.omit()
#
#N = nrow(data) %>% as.numeric()
#V = 45
#
#X = as.numeric(scale(data$age))
#Y = as.numeric(scale(data$np_reasoning))
#M = as.matrix(scale(data[,4:48]))

# Fit the model

## Fitted model
#fit =
#  sparse.txtmedint.sgrlasso.largep_omega(
#    X, M, Y,
#    lambda1 = exp(seq(-1, -6, length = 20)),
#    lambda2 = 0.1,
#    alpha = c(0.5, 0.75, 0.9), non.zeros.stop = 500, Omega.out = TRUE)
#
## Null model
#fit.n =
#  sparse.txtmedint.sgrlasso.largep_omega(
#    X, M, Y,
#    lambda1 = exp(2),
#    lambda2 = 0.1,
#    alpha = 0.9, non.zeros.stop = 500, Omega.out = TRUE)


# HBIC performance

#' Function to do HBIC computation for exposure by mediation interaction model
#'
#' @param fit
#' @param fit.n
#'
#' @return
#' @export
#'
#' @examples
hbic_int <- function(fit, fit.n){

  sum.fit = NULL
  sum.null = NULL
  sum.diff = NULL
  lik.diff = NULL
  hbic = NULL

  for (k in 1:length(fit$alpha)) {

    if (k < length(fit$Omega[[1]]) + 1){
      j = k
      #print(j)

      # calculate last term (sum)
      sum1 = 0
      sum2 = 0
      for (i in 1:N){

        s1 =
          as.numeric(
            t(M[i,] - fit$hata[,k] * X[i]) %*%
              data.matrix(data.frame(fit$Omega[[1]][j])) %*%
              (M[i,] - fit$hata[,k] * X[i]))
        sum1 = sum1 + s1

        s2 =
          as.numeric(
            t(M[i,] - fit.n$hata * X[i]) %*%
              data.matrix(data.frame(fit.n$Omega)) %*%
              (M[i,] - fit.n$hata * X[i]) )
        sum2 = sum2 + s2

      }

      sum.fit[k] <- sum1
      sum.null[k] <- sum2
      sum.diff[k] <- sum.fit[k] - sum.null[k]

      # calculate difference in 2*log-likelihood (l(fit.n) - l(fit))
      lik.diff[k] =
        - N*log(fit.n$sigmasq) + N*log(fit$sigmasq[k])
      + N*log(det(data.matrix(data.frame(fit.n$Omega)))) - N*log(det(data.matrix(data.frame(fit$Omega[[1]][j]))))
      + 1/fit.n$sigmasq * rowsum(
        ((Y - X*fit.n$c - M %*% fit.n$hatb1 - (X*M) %*% fit.n$hatb2))^2, rep(1,N))
      - 1/fit$sigmasq[k] * rowsum(
        ((Y - X*fit$c[k] - M %*% fit$hatb1[,k] - (X*M) %*% fit$hatb2[,k]))^2, rep(1,N))
      + sum.diff[k]

      # calculate hbic
      hbic[k] = lik.diff[k] - (fit.n$nump - fit$nump[k])*log(N/(2*pi))

    } else if (k > length(fit$Omega[[1]]) + length(fit$Omega[[2]])) {
      j = k - length(fit$Omega[[1]]) - length(fit$Omega[[2]])

      # calculate last term (sum)
      sum1 = 0
      sum2 = 0
      for (i in 1:N){
        s1 = as.numeric(
          t(M[i,] - fit$hata[,k] * X[i]) %*%
            data.matrix(data.frame(fit$Omega[[3]][j])) %*%
            (M[i,] - fit$hata[,k] * X[i]) )
        sum1 = sum1 + s1
        s2 = as.numeric(
          t(M[i,] - fit.n$hata * X[i]) %*%
            data.matrix(data.frame(fit.n$Omega)) %*%
            (M[i,] - fit.n$hata * X[i]) )
        sum2 = sum2 + s2
      }

      sum.fit[k] <- sum1
      sum.sat[k] <- sum2
      sum.diff[k] <- sum.fit[k]-sum.sat[k]

      # calculate difference in 2*log-likelihood (l(fit.n) - l(fit))
      lik.diff[k] =
        - N*log(fit.n$sigmasq) + N*log(fit$sigmasq[k])
      + N*log(det(data.matrix(data.frame(fit.n$Omega)))) - N*log(det(data.matrix(data.frame(fit$Omega[[3]][j]))))
      - 1/fit.n$sigmasq * rowsum(
        ((Y - X*fit.n$c - M %*% fit.n$hatb1 - (X*M) %*% fit.n$hatb2))^2, rep(1,N))
      + 1/fit$sigmasq[k] * rowsum(
        ((Y - X*fit$c[k] - M %*% fit$hatb1[,k] - (X*M) %*% fit$hatb2[,k]))^2, rep(1,N))
      + sum.diff[k]

      # calculate hbic
      hbic[k] = lik.diff[k] - (fit.n$nump - fit$nump[k])*log(N/(2*pi))

    } else{
      j = k - length(fit$Omega[[1]])

      # calculate last term (sum)
      sum1 = 0
      sum2 = 0
      for (i in 1:N){
        s1 =
          as.numeric(
            t(M[i,] - fit$hata[,k] * X[i]) %*%
              data.matrix(data.frame(fit$Omega[[2]][j])) %*%
              (M[i,] - fit$hata[,k] * X[i]) )
        sum1 = sum1 + s1

        s2 =
          as.numeric(
            t(M[i,] - fit.n$hata * X[i]) %*%
              data.matrix(data.frame(fit.n$Omega)) %*%
              (M[i,] - fit.n$hata * X[i]) )
        sum2 = sum2 + s2
      }

      sum.fit[k] <- sum1
      sum.sat[k] <- sum2
      sum.diff[k] <- sum.fit[k] - sum.sat[k]

      # calculate difference in 2*log-likelihood (l(fit.n) - l(fit))
      lik.diff[k] =
        - N*log(fit.n$sigmasq) + N*log(fit$sigmasq[k])
      + N*log(det(data.matrix(data.frame(fit$Omega[[2]][j])))) -
        N*log(det(data.matrix(data.frame(fit.n$Omega))))
      - 1/fit.n$sigmasq * rowsum(
        ((Y - X*fit.n$c - M %*%
            fit.n$hatb1 - (X*M) %*%
            fit.n$hatb2))^2, rep(1,N)) +
        1/fit$sigmasq[k] * rowsum(
          ((Y - X*fit$c[k] - M %*%
              fit$hatb1[,k] - (X*M) %*%
              fit$hatb2[,k]))^2, rep(1,N))
      + sum.diff[k]

      # calculate hbic
      hbic[k] = lik.diff[k] - (fit.n$nump - fit$nump[k])*log(N/(2*pi))
    }

  }

  print(hbic)

}

#hbic_int = hbic_int(fit, fit.n)
