#' \code{IPW}
#'
#' Estimate a population average using inverse probability weighted sample average.
#' @param dat string. A data.frame containing at least three columns: Y, which is the variable of interest; PIP, which is the pseudo-inclusion probabilities; and period, which denotes the time period of interest. 
#' Pseudo-inclusion probabilities may be specific to one period or not. 
#' @param truncate Locial. Whether or not to replace zeros in the pseudo-inclusion probabilities with the smallest non-zero value in the focal time period.
#' @param reduceVar String. Method used to reduce the variance in the inverse probability weights. Defaults to FALSE, in which case weights are simply 1/PIP. 
#' @return A data.frame with four columns: d1, which is the empirical correlation between sample membership and the variable of interest; estimated_d1,
#' which is an estimate of d1 obtained by solving Meng's formula for estimation error for d1 (this value reflects selection bias plus other sources of error in the data);
#' d2, which is sqrt((1-f)/f), i.e. the data quantity term; d3, which is the problem difficulty term, sd(Y); error, which
#' is the difference between the sample average and population average, given by d1*d2*d3; and finally, nEff, which is the effective size of the sample, that is, the
#' size of the random sample that would yield an estimate of the population average with same error.
#' @details The function returns various quantities of interest defined in Meng, X. L. (2018). Statistical paradises and paradoxes in big data (I): Law of large populations, big data paradox, and the 2016 us presidential election. Annals of Applied Statistics, 12(2), 685-726. https://doi.org/10.1214/18-AOAS1161SF
#' @export

IPW <- function(dat, 
                truncate = TRUE, 
                reduceVar = FALSE) {
  
  if(! "PIP" %in% colnames(dat)) stop("dat must contain a column for pseudo-inclusion probabilities called PIP.")
  if(! "Y" %in% colnames(dat)) stop("dat must contain a column for the response variable of interest called Y.")
  if (!is.data.frame(dat)) stop("dat must be a dataframe.")
  
  stats <- lapply(unique(dat$period),
                  function(x) {
                    
                    if (truncate) {
                      
                      d <- dat[dat$period == x, ]
                      
                      minPIP <- min(d$PIP[d$PIP != 0])
                      
                      print(paste0("Replacing zeros in PIP for period ", x, " with the smallest non-zero value in that period, ", minPIP))
                      
                      d[d$PIP == 0] <- minPIP
                      
                    }
                    
                    w <- 1 / d$PIP
                    
                    wAv <- sum(w * d$Y) / sum(w)
                    
                    sampleAv <- mean(d$Y)
                    
                    data.frame(period = x,
                               sample_av = sampleAv,
                               weighted_av = wAv)
                    
                  }
  )
  
  do.call("rbind", stats)
  
}