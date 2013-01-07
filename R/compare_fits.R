compare.fits = function(x, fit1, fit2) {
  
  cov = as.matrix(data.frame(y = rep(1, length(x[,1])), x))
  
  # Extract beta estimates
  fit1.beta = fit1$fixed.effects$Estimate
  fit2.beta = fit2$fixed.effects$Estimate
  
  # Extract beta variance matrix
  var.b = fit1$var.b
  
  tdbeta = fitdvcov(cov, fit1.beta, fit2.beta, var.b)$tdbeta
  cfits = fitdvcov(cov, fit1.beta, fit2.beta, var.b)$cfits
  
  return(tdbeta)
  return(cfits) #added 1/25/2013
}