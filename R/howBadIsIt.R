#' \code{howBadIsIt}
#'
#' Given a dataset, this function calculates the correlation between sample membership and the variable of interest, the effective sample size, and more.
#' @param dat string. A data.frame containing two columns: R, which is denotes sample membership and must take values {0,1};
#' and Y, which is the variable of interest. The function will throw an error is colnames(dat) != c("R","Y").
#' @return A data.frame with four columns: d1, which is the empirical correlation between sample membership and the variable of interest; estimated_d1,
#' which is an estimate of d1 obtained by solving Meng's formula for estimation error for d1 (this value reflects selection bias plus other sources of error in the data);
#' d2, which is sqrt((1-f)/f), i.e. the data quantity term; d3, which is the problem difficulty term, sd(Y); error, which
#' is the difference between the sample average and population average, given by d1*d2*d3; and finally, nEff, which is the effective size of the sample, that is, the
#' size of the random sample that would yield an estimate of the population average with same error.
#' @details The function returns various quantities of interest defined in Meng, X. L. (2018). Statistical paradises and paradoxes in big data (I): Law of large populations, big data paradox, and the 2016 us presidential election. Annals of Applied Statistics, 12(2), 685-726. https://doi.org/10.1214/18-AOAS1161SF
#' @export

howBadIsIt <- function(dat) {

  if(!is.data.frame(dat)) stop("dat must be a data.frame")

  if(!colnames(dat)[1] == "R") stop("first column name must be R")

  if(!colnames(dat)[2] == "Y") stop("second column name must be Y")

  # define variables
  n <- nrow(dat[dat$R == 1, ])

  n1 <- nrow(dat[dat$R == 1 & dat$Y == 1, ])

  ybar_sample <- n1/n

  N <- nrow(dat)

  N1 <- nrow(dat[dat$Y == 1, ])

  ybar_pop <- N1/N

  f <- n/N

  # define terms
  d1 <- cor(dat$R,
            dat$Y,
            use = "pairwise.complete.obs",
            method = "pearson")

  #d2 <- sqrt((N-n)/n)

  d2 <- sqrt((1-f)/f)

  d3 <- sd(dat$Y)

  # calculate error
  error <- d1*d2*d3

  # calculate effective sample size
  # first estimate d1

  est_d1 <- error / d2 / d3

  # then calc nEff
  nEff <- (f/(1-f)) * (1/(est_d1^2))

  return(data.frame(d1 = d1,
                    estimated_d1 = est_d1,
                    d2 = d2,
                    d3 = d3,
                    n = n,
                    N = N,
                    ybar_pop = ybar_pop,
                    ybar_sample = ybar_sample,
                    error = error,
                    nEff = nEff))

}

