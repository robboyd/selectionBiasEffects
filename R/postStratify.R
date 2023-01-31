#' \code{postStratify}
#'
#' Estimate a population average using post-stratified sample mean.
#' @param dat string. A data.frame containing columns for sample membership (1 or 0), the variable of interest (can be NA for rows where sample membership = 0) and the variables on which to stratify. 
#' @param Z String or vector. Names of the variables on which to poststratify. These must match the relevant column names in dat.
#' @param discretize String or vector. Any variables in Z that are measured on a continuous scale and need to be discretized. 
#' @param breaks Numeric. Number of bins into which the variables in discretize will be split. Bins are of equal length.
#' @param plotP Logical. Whether or not to display a mosaic plot showing the proportion of the population in each poststratification cell. Only plots if the number of cells is less than 100.
#' @return The poststraified sample mean.
#' @details At present, the function ignores cells that are not in the population.
#' @export

postStratify <- function(dat, Z, discretize, breaks, plotP) {
  
  for (i in 1:ncol(dat)) {
    
    if(names(dat)[i] %in% discretize) {
      
      dat[,i] <- cut(dat[,i], 
                     breaks = breaks,
                     labels = FALSE,
                     include.lowest = TRUE)
      
      dat[,i] <- as.numeric(dat[,i])
      
    }
    
  }
  
  popTab <- table(dat[,Z])
  
  propTab <- popTab/sum(popTab)
  
  if(plotP & length(propTab) <= 100) plot(propTab)
  
  samp <- dat[dat$sampled_units == 1, ]
  
  tab <- propTab
  
  tab[] <- NA
  
  comb <- data.frame(expand.grid(rownames(tab), colnames(tab)))
  
  strataMeans <- sapply(1:nrow(comb),
                        function(x) {
                          mean(samp$heather_true_dist[samp$X_Annual.Mean.Temperature == comb$Var1[x] &
                                                        samp$cs_envzones == comb$Var2[x]])
                        }
  )
  
  comb$mean <- strataMeans
  
  Ps <- sapply(1:nrow(comb),
               function(x) {
                 propTab[comb$Var1[x], comb$Var2[x]]
               }
  )
  
  comb$P <- Ps
  
  weighted_mean <- comb$mean * comb$P
  
  sum(weighted_mean, na.rm = T)
  
}