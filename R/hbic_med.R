# hbic_med function

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
#  sparse.mediation.sgrplasso.largep_omega(
#    X, M, Y,
#    lambda1 = exp(seq(-1, -6, length = 20)),
#    lambda2 = 0.1,
#    alpha = c(0.5, 0.75, 0.9), non.zeros.stop = 500, Omega.out = TRUE)
#
## Null model
#fit.n =
#  sparse.mediation.sgrplasso.largep_omega(
#    X, M, Y,
#    lambda1 = exp(2),
#    lambda2 = 0.1,
#    alpha = 0.9, non.zeros.stop = 500, Omega.out = TRUE)


# HBIC performance

hbic_med <- function(fit, fit.n){

  sum.fit = NULL
  sum.null = NULL
  sum.diff = NULL
  lik.diff = NULL
  hbic = NULL

  for (k in 1:length(fit$alpha)) {

    # calculate last term (sum)
    sum1 = 0
    sum2 = 0

    for (i in 1:N) {
      s1 =
        as.numeric(
          t(M[i,] - fit$hata[,k] * X[i]) %*% #1x45
            data.matrix(data.frame(fit$Omega[[k]])) %*% # 45x45
            (M[i,] - fit$hata[,k] * X[i]) #1x45
        )
      sum1 = sum1 + s1

      s2 =
        as.numeric(
          t(M[i,] - fit.n$hata * X[i]) %*%
            data.matrix(data.frame(fit.n$Omega)) %*%
            (M[i,] - fit.n$hata * X[i])
        )
      sum2 = sum2 + s2

    }
    #print(sum1)
    sum.fit[k] <- sum1
    sum.null[k] <- sum2
    sum.diff[k] <- sum.fit[k] - sum.null[k]

    # calculate difference in 2*log-likelihood (l(fit.n) - l(fit))
    lik.diff[k] =
      - N*log(fit.n$sigma_y_sq) + N*log(fit$sigma_y_sq[k]) # nlog(sigma^2)
    + N*log(det(data.matrix(data.frame(fit.n$Omega)))) - N*log(det(data.matrix(data.frame(fit$Omega[[k]])))) # nlog|omega|
    - 1/fit.n$sigma_y_sq * rowsum(((Y - X*fit.n$c - M %*% fit.n$hatb))^2, rep(1,N)) + 1/fit$sigma_y_sq[k] * rowsum(((Y - X*fit$c[k] - M %*% fit$hatb[,k]))^2, rep(1,N)) # 1/sigma^2 sum (Yi - c'Xi - Mib)^2
    + sum.diff[k] # the last term (sum)

    # calculate hbic
    hbic[k] = lik.diff[k] - (fit.n$nump - fit$nump[k])*log(N/(2*pi))
  }
  print(hbic)
}

#hbic_med = hbic_med(fit, fit.n)
