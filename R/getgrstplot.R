#Obtain st residual plots( st res vs. yhat and QQ of ehat) in GR
#4 graphs
getgrstplot=function(rlme.fit){
  
#dev.off()
  y=rlme.fit$y
  xstar=rlme.fit$xstar
  ystar=rlme.fit$ystar
  residstar.gr=rlme.fit$ehats #resid after the star
  
  ystarhat.gr=(ystar-residstar.gr) #
  yhat.gr=y-rlme.fit$effect.err #y-marginal error yields fitted(y) as in REML
  
  temp <- stanresidgr(x=xstar,y=ystar,resid=residstar.gr,delta=.80,param=2,conf=.95)
  standr.gr <- temp$stanr
  
  trim=2
  par(mfrow = c(1, 2),font.main=1)
  
  plot(standr.gr ~ yhat.gr, pch = "o", xlab="Fit",ylab="Standardized Residual", main="Stand. Residuals vs. Fits in GR")
  abline(h = c(-trim, trim), col = "red")
  
  qqnorm(standr.gr, pch = "o", main = "Normal Q-Q Plot")
  qqline(standr.gr, col = "red")
  
  par(mfrow=c(1, 1))
  #in lme
#   standr.lme=residuals(lme1, type="pearson") #or
#   standr.lme=as.vector(residuals(lme1, type="pearson"))
  
}

#Obtain st residual plots( st res vs. yhat and QQ of ehat) in GR
#4 graphs
getlmestplot=function(rlme.fit){
  
  #dev.off()
  standr.lme=rlme.fit$standr.lme
  y=rlme.fit$y
  yhat.lme=y-rlme.fit$effect.err

  trim=2
  par(mfrow = c(1, 2),font.main=1)
  
  plot(standr.lme ~ yhat.lme, pch = "o", xlab="Fit",ylab="Standardized Residual", main="Stand. Residuals vs. Fits in GR")
  abline(h = c(-trim, trim), col = "red")
  
  qqnorm(standr.lme, pch = "o", main = "Normal Q-Q Plot")
  qqline(standr.lme, col = "red")
  
  par(mfrow=c(1, 1))
  
}
