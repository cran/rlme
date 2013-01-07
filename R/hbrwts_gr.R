hbrwts_gr=function(xmat, y,
                   robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
                   percent=0.95, intest=myltsreg(xmat,y)$coef) {
  # By J. W. McKean from the ww pack
  # HBR weights as a vector for gee.
  # This is modified hbr in ww package.
  
  xmat=as.matrix(xmat)
  y=as.matrix(y)
  n=dim(xmat)[1]
  p=dim(xmat)[2]
  cut=qchisq(percent,p)
  resids=y-intest[1]-xmat%*%as.matrix(intest[2:(p+1)])
  sigma=mad(resids)
  m=psi(cut/robdis2)
  a=resids/(sigma*m)
  c=(median(a)+3*mad(a))^2
  h=sqrt(c)/a
  #tmp=pairup(h)
  #ans=psi(abs(tmp[,1]*tmp[,2]))
  ans=psi(abs(h))
  ans 
}
