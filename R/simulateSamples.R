#' \code{simulateSamples}
#'
#' Simulates many simple random samples (SRSs) of size n_eff and larger, based samples of size n. n_eff, n and other arguments may be taken from the outputs of \code{howBadIsIt}. For each sample,
#' the function returns various statistics (see below).
#' @param dat string. A data.frame containing two columns: R, which is denotes sample membership and must take values {0,1};
#' and Y, which is the variable of interest. The function will throw an error is colnames(dat) != c("R","Y").
#' @param N Numeric. Population size.
#' @param n Numeric. Sample size.
#' @param n_eff Numeric. Effective sample size. Might be taken from the outputs of \code{howBadIsIt}.
#' @param ddc Numeric. Data defect correlation. Might be taken from the outputs of \code{howBadIsIt}.
#' @param nSim Numeric. Number of samples to simulate. nSim SRSs and large, biased samples will be simulated.
#' @param truth Numeric. Population mean.
#' @return A data.frame with eight columns: smallSRSEstimate = Estimate of the population mean from the SRS of size n_eff; largeBiasedEstimate = estimate of population mean from the biased sample of size n;
#' smallSRSCoverage = % of (Wald) confidence intervals from the simulated SRSs covering the true population mean;
#' largeBiasedCoverage =  % of (Wald) confidence intervals from the simulated biased samples of size n covering the true population mean;
#' smallSRSError = sample mean - population mean for the SRS samples; largeBiasedError = sample mean - population mean for the biased samples of size n;
#' smallSRSSquaredError = smallSRSError squared; largeBiasedSquaredError = largeBiasedError squared; sim = simulation number. The data.frame has nSim rows, one for each simulation.
#' @details The function returns various quantities of interest defined in Meng, X. L. (2018). Statistical paradises and paradoxes in big data (I): Law of large populations, big data paradox, and the 2016 us presidential election. Annals of Applied Statistics, 12(2), 685-726. https://doi.org/10.1214/18-AOAS1161SF
#' @export

simulateSamples <- function(dat, N, n, n_eff, ddc, nSim, truth) {

  stats <- lapply(1:nSim,
                  function(x) {

                    ## SRS of size n_eff

                    samp <- sample(x = dat$Y,
                                   size = n_eff,
                                   replace = F)

                    av <- mean(samp)

                    halfWidth <- sqrt(av*((1-av)/length(samp)))

                    #halfWidth <- sd(samp) / sqrt(n_eff)

                    upper <- av + (1.96 * halfWidth)

                    lower <- av - (1.96 * halfWidth)

                    error <- av - truth

                    ## sample of size n

                    corr.mat <- matrix(c(1, ddc,
                                         ddc, 1), nrow = 2)

                    samp2 <- simstudy::genCorGen(N, nvars = 2,
                                                 params1 = c(mean(dat$R), mean(dat$Y)),
                                                 corMatrix = corr.mat,
                                                 dist = "binary",
                                                 method = "ep", wide = TRUE)

                    samp2 <- samp2$V2[samp2$V1 == 1] # true values at sampled locations

                    av2 <- mean(samp2)

                    halfWidth2 <- sqrt(av2*((1-av2)/length(samp2)))

                    #halfWidth2 <- sd(samp2) / sqrt(n)

                    upper2 <- av2 + (1.96 * halfWidth2)

                    lower2 <- av2 - (1.96 * halfWidth2)

                    error2 <- av2 - truth

                    data.frame(SRS.mean = av,
                               SRS.upper = upper,
                               SRS.lower = lower,
                               SRS.error = error,
                               mean = av2,
                               upper = upper2,
                               lower = lower2,
                               error = error2)

                  }
  )

  stats <- do.call("rbind", stats)

  ## coverage of SRS CIs

  covered <- ifelse(stats$SRS.upper > truth & stats$SRS.lower < truth, 1, 0)

  coverage <- length(covered[covered == 1]) / length(covered)

  ## coverage of biased sample CIs

  covered2 <- ifelse(stats$upper > truth & stats$lower < truth, 1, 0)

  coverage2 <- length(covered2[covered2 == 1]) / length(covered2)

  data.frame(smallSRSEstimate = stats$SRS.mean,
             largeBiasedEstimate = stats$mean,
             smallSRSCoverage = coverage,
             largeBiasedCoverage = coverage2,
             smallSRSError = stats$SRS.error,
             largeBiasedError = stats$error,
             smallSRSSquaredError = stats$SRS.error^2,
             largeBiasedSquaredError = stats$error^2,
             sim = 1:nrow(stats))

}
