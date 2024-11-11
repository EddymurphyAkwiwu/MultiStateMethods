compute_Weibull_cdf <- function(shape_hat, intercept, coeff1, coeff2, time_range, covariates) {
  
  
  # Initialize vectors to store results
  CDF_values <- rep(NA, length(time_range))
  
  error_margins <- rep(NA, length(time_range))
  
  
  # Loop over time values
  for (i in 1:length(time_range)) {
    time <- time_range[i]
    
    # Initialize vectors to store intermediate results
    CDF_vals <- rep(NA, nrow(covariates))
 
    # Loop over covariates
    for (j in 1:nrow(covariates)) {
      covariate1 <- covariates[j, 1]
      covariate2 <- covariates[j, 2]
      
      # Compute the CDF using estimated parameters
      scale_hat <- exp(intercept + coeff1 * covariate1 + coeff2 * covariate2)
      CDF_hat <- 1 - exp(-((time / scale_hat) ^ shape_hat)) # or use pweibull(q=time,shape =shape_hat,  scale= scale_hat)
     
      

      
      # Store results in vectors
      CDF_vals[j] <- CDF_hat

    }
    
    # Store averaged results in matrices
    CDF_values[i] <-  mean(CDF_vals)
 
  }
  
 
  return(CDF_values)
}
