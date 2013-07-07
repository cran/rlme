lmr <- function(f, data, se=FALSE, method='L-BFGS-B') {
  #for large data, use method='L-BFGS-B'
  
  response_name = as.character(as.formula(f)[[2]])
  covariate_names = attributes(terms(f))$term.labels
  
  x = as.matrix(data[covariate_names])
  x = apply(x, 2, function(x) {
    x - mean(x)
  })
  
  y = as.matrix(data[[response_name]])
  
  # Fit with minimize_dispersion
  
  estimate = minimize_dispersion(x, y, se = se, method=method)
  #for large data, use method='L-BFGS-B'
  
  fit = list()
  fit$num.obs = nrow(y)
  fit$y = y
  fit$formula = f
  fit$method = "lmr"
  fit$ehat = estimate$ehat
  
  if(se == TRUE) {
    intercept.se = taustar(estimate$ehat, ncol(x)) / sqrt(nrow(x))
    
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), 
                                   Estimate = estimate$theta,
                                   StdError = c(intercept.se, estimate$standard.error))
  } else {
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), 
                                   Estimate = estimate$theta)
  }
  
  class(fit) = "rlme"
  
  return(fit)
}