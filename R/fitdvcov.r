fitdvcov= function(x1,beta1,beta2,vcw){
#     By McKean.
#     x1 is the design with intercept in.   The
#     number of colums is the same as the length of beta1.
#
#     Obtains tdbetas for the fits based on estimates beta1 and bets2
#     Standardization is by Wilcoxon var-cov matrix which is vcw.

#     TDBETAS
  n=dim(x1)[1]
  p=dim(x1)[2]
  bd = beta1 - beta2
  tdbeta=t(bd)%*%solve(vcw)%*%bd
  bmtd=(4*(p+1)^2)/n
#     CFITS
  fit1 = x1 %*% beta1
  fit2 = x1 %*% beta2
  xv = x1%*%vcw%*%t(x1) 
  cfits = rep(0,n)
  for(i in 1:n){
    cfits[i] = (fit1[i] - fit2[i])/sqrt(xv[i,i])
  } 
  bmcf=2*sqrt((p+1)/n)
  list(tdbeta=c(tdbeta),bmtd=bmtd,cfits=c(cfits),bmcf=bmcf)
}
