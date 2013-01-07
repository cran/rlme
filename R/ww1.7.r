# This is the ww package functions. See the authors Terpstra and McKean (2005). This has wil and hbr scores.
# What I changed is in solve(). i replaced with ginv(). Generalized inverse is used when singular.
# Rfit would be used for general scores.
# R Functions for Weighted Wilcoxon Estimation and Inference  
#
# Version:      1.7
# Last Update:  12/04/09
# Requirements: R 2.8.1 or higher
#               'quantreg' package
#               'MASS' package (R 1.9.0 or higher)
#               'lqs' package (R 1.8.1 or lower)
# 
# Authors:      Jeff Terpstra                     Joe McKean
# 
# Address:      Department of Statistics          Department of Statistics
#               Western Michigan University       Western Michigan University
#               5504 Everett Tower                5506 Everett Tower
#               Kalamazoo, MI  49008              Kalamazoo, MI  49008
#
# Email:        jeffrey.terpstra@wmich.edu        joseph.mckean@wmich.edu
#
# Phone:        269-387-0611                      269-387-4541
#
# Fax:          269-387-1419                      269-387-1419                     
#
# Notes:        These functions require the 'quantreg' and 'MASS'/'lqs' packages
#               which are available from http://cran.r-project.org/
#               For example,
#                 -install.packages("quantreg",lib="/destdir") 
#                 -library("quantreg",lib.loc="./") loads the "quantreg" package
#                  on a UNIX machine
#                 -library("MASS") loads the "MASS" package
#
############################################################################### 
#                                                                             #
# Functions that calculate WW-estimates (and inference)                       #
#                                                                             #
###############################################################################
#
# wwfit returns the estimates, residuals, and the weight matrix for an 
# arbitrary set of weights
# 
wwfit=function(x, y, bij=wilwts(as.matrix(x)), center=F) {
  if(!any(search()=="package:quantreg"))
    stop("wwfit:  The 'quantreg' package is not loaded.")
  x=as.matrix(x)
  n=dim(x)[1]
  p=dim(x)[2]
  if(center) {
    xbar=apply(x,2,mean)
    x=apply(x,2,function(x){x-mean(x)})
  }
  ypairs=pairup(y)
  yi=ypairs[,1]
  yj=ypairs[,2]
  xpairs=pairup(x)
  xi=xpairs[,1:p]
  xj=xpairs[,(p+1):(2*p)]
  newy=bij*(yi-yj)
  newx=bij*(xi-xj)
  #
  # "br" corresponds to l1fit in S (recommended when n << 5,000, p << 20).
  # "fnb" is designed for fairly large problems.  See help(rq).
  #
  if (((n*(n-1)/2) < 5000) & (p < 20))  
    tmp=rq.fit.br(newx,newy,tau=0.5,ci=F)
  else
    tmp=rq.fit.fnb(cbind(newx),cbind(newy),tau=0.5)
  est=tmp$coefficients
  #
  # Using the median of the (initial) residuals as an estimate of intercept.
  #
  int=median(y - (x %*% as.matrix(est)))
  resid=as.vector(y - int - (x %*% as.matrix(est)))
  if(center) {
    int=int-(t(as.matrix(est))%*%as.matrix(xbar)) 
  }
  wts=matrix(0,n,n)
  index=pairup(1:n)
  wts[index]=bij
  wts[index[,2:1]]=bij
  ans=list(coefficients=c(int,est), residuals=resid, weights=wts)
  ans 
}
#
# This function performs a Weighted Wilcoxon analysis using (the default)
# WIL, THEIL, GR, HBR, or BL weights.  If bij is numeric then GR (i.e. non-random)
# weights are assumed.
# 
wwest=function(x, y, bij="WIL", center=F, print.tbl=T) {
  if(is.character(bij)) {
    type=bij
    bij=switch(bij,
      WIL = wilwts(x),
    THEIL = theilwts(x),
       GR = grwts(x),
      HBR = hbrwts(x, y),
       BL = blwts(x,y),
      stop("wwest:  The weight type should be WIL, THEIL, GR, HBR, or BL"))
  }
  else {
    type="GR"
  }
  tmp1=wwfit(x, y, bij, center)
  n=length(y)
  p=length(tmp1$coef)-1
  ans=cbind(tmp1$coef)
  tmp2=switch(type,
    WIL = varcov.gr(x,tmp1$weights,tmp1$residuals),
  THEIL = varcov.gr(x,tmp1$weights,tmp1$residuals),
     GR = varcov.gr(x,tmp1$weights,tmp1$residuals),
    HBR = varcov.hbr(x,tmp1$weights,tmp1$residuals),
     BL = varcov.hbr(x,tmp1$weights,tmp1$residuals))
#  ans=cbind(ans,sqrt(diag(tmp2$varcov))) #I added to omit the NANs in GR
bb <- diag(tmp2$varcov)
bb[bb<0] <-0
  ans=cbind(ans,sqrt(bb))
  ans=cbind(ans,ans[,1]/ans[,2])
  #
  # note: The t distribution is used for p-values.  See, for
  #       example, page 147 of HM.
  #
  ans=cbind(ans,2*pt(abs(ans[,3]),n-p-1,lower.tail=FALSE))
  ans=round(ans,4)
  dimnames(ans)=list(paste("BETA",0:p,sep=""),c("EST","SE","TVAL","PVAL"))
  if (print.tbl) {
  tmp3=wald(tmp1$coef,tmp2$varcov,amat=cbind(rep(0,p),diag(p)),true=rep(0,p),n=n)
  BETA=""
  for(i in 1:p) {BETA=paste(BETA,"BETA",i,"=",sep="")}
  cat("\n")
  cat(paste("Wald Test of H0: ",BETA,"0\n",sep=""))
  cat(paste("TS:",round(tmp3[1],4),"PVAL:",round(tmp3[2],4),"\n"))
  cat("\n")
  if (type=="WIL") {
    tmp4=regrtest(x,y,print.tbl=F)
    cat(paste("Drop Test of H0: ",BETA,"0\n",sep=""))
    cat(paste("TS:",round(tmp4$fr,4),"PVAL:",round(tmp4$pval,4),"\n"))
    cat("\n")
  }
  prmatrix(ans,na.print="")
  repeat {
    cat("\n")
    cat("Would you like to see residual plots (y/n)?","\n")
    yn=as.character(readline())
    if(yn=="y" | yn=="Y" | yn=="yes") {
      studres=switch(type,
        WIL = studres.gr(x,tmp1$weights,tmp1$residuals),
      THEIL = studres.gr(x,tmp1$weights,tmp1$residuals),
         GR = studres.gr(x,tmp1$weights,tmp1$residuals),
        HBR = studres.hbr(x,tmp1$weights,tmp1$residuals),
         BL = studres.hbr(x,tmp1$weights,tmp1$residuals))
      yhat=y-tmp1$residuals
      par(mfrow=c(2,2))
      plot(yhat,tmp1$residuals,xlab="Fit",ylab="Residual",
           main="Residuals vs. Fits")
      hist(tmp1$residuals,freq=FALSE,main="Histogram of Residuals",
           xlab="Residual")
      plot(studres,xlab="Case",ylab="Studentized Residual",
           main="Case Plot of\nStudentized Residuals")
      abline(h=c(-2,2))
      qqnorm(tmp1$residuals,main="Normal Q-Q Plot of Residuals")
      qqline(tmp1$residuals)
      #boxplot(tmp1$residuals,
      #         horizontal=T,main="Boxplot of Residuals")
      break
    }
    if(yn=="n" | yn=="N" | yn=="no" ) break 
  }
  }
  invisible(list(tmp1=tmp1,tmp2=tmp2,ans=ans))
}
#
# Function that sets up pairwise comparisons
#
pairup=function(x,type="less") {
  x=as.matrix(x)
  n=dim(x)[1]
  i=rep(1:n,rep(n,n))
  j=rep(1:n,n)
  c1=apply(x,2,function(y){rep(y,rep(length(y),length(y)))})
  c2=apply(x,2,function(y){rep(y,length(y))})
  ans=cbind(c1,c2)
  ans=switch(type, less=ans[(i<j), ], leq=ans[i<=j, ], neq=ans)
  ans
}
###############################################################################
#                                                                             #
# Functions that calculate weights                                            #
#                                                                             #
###############################################################################
#
# Wilcoxon weights
# 
wilwts=function(xmat) {
  xmat=as.matrix(xmat)
  n=dim(xmat)[1]
  ans=rep(1,n*(n-1)/2)
  ans
}
#
# Theil weights
#
theilwts=function(xmat) {
  xmat=as.matrix(xmat)
  p=dim(xmat)[2]
  xpairs=pairup(xmat)
  xi=xpairs[,1:p]
  xj=xpairs[,(p+1):(2*p)] 
  diff=as.matrix(xi-xj)
  ans=apply(diff,1,function(y){sqrt(sum(y*y))})
  ans=1/ans
  ans[ans==Inf]=0
  ans
}
#
# GR weights 
# 
grwts=function(xmat, robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
               percent=0.95, k=2) {
  xmat=as.matrix(xmat)
  n=dim(xmat)[1]
  p=dim(xmat)[2]
  cut=qchisq(percent,p)
  h=pmin(1,((cut/robdis2)^(k/2)))
  tmp=pairup(h)
  ans=tmp[,1]*tmp[,2]
  ans
} 
#
# HBR weights
# 
hbrwts=function(xmat, y,
                robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
                percent=0.95, intest=myltsreg(xmat,y)$coef) {
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
  tmp=pairup(h)
  ans=psi(abs(tmp[,1]*tmp[,2]))
  ans 
}
#
# psi function in HBR weights
# 
psi=function(x)
{
x[x==-Inf]=-100
x[x==Inf]=100
ans=-1*(x<=-1) + x*(-1<x & x<1) + 1*(x>=1)
ans
}
#
# Factored HBR weights that only downweight bad leverage points
# 
blwts=function(xmat, y,
               robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
               percent=0.95, k=2, intest=myltsreg(xmat,y)$coef) {
  xmat=as.matrix(xmat)
  y=as.matrix(y)
  n=dim(xmat)[1]
  p=dim(xmat)[2]
  cut1=qchisq(percent,1)
  cutp=qchisq(percent,p)
  resids=y-intest[1]-xmat%*%as.matrix(intest[2:(p+1)])
  sigma=mad(resids)
  ind1=as.numeric(abs(resids)>sigma*sqrt(cut1))
  ind2=as.numeric(robdis2>cutp)
  tmp=(cutp/robdis2)^(k/2)
  h=1-(ind1*ind2*(1-tmp))
  tmp2=pairup(h)
  ans=tmp2[,1]*tmp2[,2]
  ans
}
#
# wrapper function that can be used to get WIL, THEIL, GR, HBR, and BL weights
# 
wts=function(xmat, y, type="WIL", percent=0.95, k=2,
             robdis2=if(type!="WIL") mycov.rob(as.matrix(xmat),
               method="mcd")$robdis2 else NULL,
             intest=if(type=="HBR" | type=="BL") myltsreg(xmat,y)$coef else NULL) {
  xmat=as.matrix(xmat)
  y=as.matrix(y)
  switch(type,
    WIL = wilwts(xmat),
  THEIL = theilwts(xmat),
     GR = grwts(xmat, robdis2, percent, k),
    HBR = hbrwts(xmat, y, robdis2, percent, intest),
     BL = blwts(xmat, y, robdis2, percent, k, intest),
    stop("wts:  TYPE should be WIL, THEIL, GR, HBR or BL"))
}
###############################################################################
#                                                                             #
# Functions that calculate (robust) measures of location, dispersion, and     #
# distance                                                                    #
#                                                                             #
###############################################################################
#
# New version of mahalanobis
#
# Changed the matrix statement in the third row from matrix(x, ncol = length(x))
# to matrix(x, nrow = length(x)) so that the function works for p=1. 
#
mymahalanobis=function(x, center, cov, inverted = FALSE, tol.inv = 1e-07) {
  x <- if (is.vector(x))
      matrix(x, nrow = length(x))
  else as.matrix(x)
  x <- sweep(x, 2, center)
  if (!inverted)
      cov <- solve(cov, tol = tol.inv)
  retval <- rowSums((x %*% cov) * x)
  names(retval) <- rownames(x)
  retval
}
#
# R-Based MCD/MVE/Classical measures of location and dispersion
#
# Changes to cov.rob:
#  - eliminate the use of "seed" and set the seed to the same value every time.
#  - use mymahalanobis instead of mahalanobis.
#  - send back robdis2.
# 
mycov.rob=function (x, cor = FALSE, quantile.used = floor((n + p + 1)/2), 
    method = c("mve", "mcd", "classical"), nsamp = "best") {
    if(v1.9.0()) {
      if(!any(search()=="package:MASS"))
        stop("mycov.rob:  The 'MASS' package is not loaded.")
      PACK="MASS"}
    else {
      if(!any(search()=="package:lqs"))
        stop("mycov.rob:  The 'lqs' package is not loaded.")
      PACK="lqs"}
    method <- match.arg(method)
    x <- as.matrix(x)
    xcopy=x
    if (any(is.na(x)) || any(is.infinite(x))) 
        stop("mycov.rob:  missing or infinite values are not allowed")
    n <- nrow(x)
    p <- ncol(x)
    if (n < p + 1) 
        stop(paste("mycov.rob:  At least", p + 1, "cases are needed"))
    if (method == "classical") {
        center=colMeans(x)
        cov=var(x)
        robdis2=mymahalanobis(xcopy,center,cov)
        ans <- list(center = colMeans(x), cov = var(x), robdis2 = robdis2)
    }
    else {
        if (quantile.used < p + 1) 
            stop(paste("mycov.rob:  quantile must be at least", p + 1))
        divisor <- apply(x, 2, IQR)
        if (any(divisor == 0)) 
            stop("mycov.rob:  at least one column has IQR 0")
        x <- x/rep(divisor, rep(n, p))
        qn <- quantile.used
        ps <- p + 1
        nexact <- choose(n, ps)
        if (is.character(nsamp) && nsamp == "best") 
            nsamp <- if (nexact < 5000) 
                "exact"
            else "sample"
        if (is.numeric(nsamp) && nsamp > nexact) {
            warning(paste("only", nexact, "sets, so all sets will be tried"))
            nsamp <- "exact"
        }
        samp <- nsamp != "exact"
        if (samp) {
            if (nsamp == "sample") 
                nsamp <- min(500 * ps, 3000)
        }
        else nsamp <- nexact
        if (exists(".Random.seed", envir = .GlobalEnv)) {
            save.seed <- .Random.seed
            on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
        }
        set.seed(123)
        z <- .C("mve_fitlots", as.double(x), as.integer(n), as.integer(p), 
            as.integer(qn), as.integer(method == "mcd"), as.integer(samp), 
            as.integer(ps), as.integer(nsamp), crit = double(1), 
            sing = integer(1), bestone = integer(n), PACKAGE = PACK)
        z$sing <- paste(z$sing, "singular samples of size", ps, 
            "out of", nsamp)
        crit <- z$crit + 2 * sum(log(divisor)) + if (method == 
            "mcd") 
            -p * log(qn - 1)
        else 0
        best <- seq(n)[z$bestone != 0]
        if (!length(best)) 
            stop("mycov.rob:  x is probably collinear")
        means <- colMeans(x[best, , drop = FALSE])
        rcov <- var(x[best, , drop = FALSE]) * (1 + 15/(n - p))^2
        dist <- mymahalanobis(x, means, rcov)
        cut <- qchisq(0.975, p) * quantile(dist, qn/n)/qchisq(qn/n, 
            p)
        center=colMeans(x[dist < cut, , drop = FALSE]) * divisor
        cov <- divisor * var(x[dist < cut, , drop = FALSE]) * 
            rep(divisor, rep(p, p))
        robdis2=mymahalanobis(xcopy,center,cov)
        attr(cov, "names") <- NULL
        ans <- list(center = center, cov = cov, robdis2 = robdis2,
                    msg = z$sing, crit = crit, best = best)
    }
    if (cor) {
        sd <- sqrt(diag(ans$cov))
        ans <- c(ans, list(cor = (ans$cov/sd)/rep(sd, rep(p, 
            p))))
    }
    ans$n.obs <- n
    ans
}
###############################################################################
#                                                                             #
# Functions that calculate high breakdown regression estimates                #
#                                                                             #
###############################################################################
#
# R-Based LTS regression 
#
myltsreg=function(xmat, y) 
{
if(v1.9.0()) {
  if(!any(search()=="package:MASS"))
    stop("mycov.rob:  The 'MASS' package is not loaded.")}
else {
  if(!any(search()=="package:lqs"))
    stop("mycov.rob:  The 'lqs' package is not loaded.")}
xmat=as.matrix(xmat)
if (exists(".Random.seed", envir = .GlobalEnv)) {
   save.seed <- .Random.seed
   on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
}
set.seed(123)
tmp=ltsreg(xmat,y,intercept=T)
ans=list(coefficients=tmp$coefficients, residuals=tmp$residuals)
ans
}
#
# R-Based LMS regression 
#
mylmsreg=function(xmat, y) 
{
if(v1.9.0()) {
  if(!any(search()=="package:MASS"))
    stop("mycov.rob:  The 'MASS' package is not loaded.")}
else {
  if(!any(search()=="package:lqs"))
    stop("mycov.rob:  The 'lqs' package is not loaded.")}
xmat=as.matrix(xmat)
if (exists(".Random.seed", envir = .GlobalEnv)) {
   save.seed <- .Random.seed
   on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
}
set.seed(123)
tmp=lmsreg(xmat,y,intercept=T)
ans=list(coefficients=tmp$coefficients, residuals=tmp$residuals)
ans
}
###############################################################################
#                                                                             #
# Functions that calculate scale parameter estimates                          #
#                                                                             #
###############################################################################
#
# Density estimate of wilcoxontau = 1/(sqrt(12)*\int f^2(x)dx).
# Discussed in Section 3.7.1 of HM.
#
# note:  Use delta=0.80 if n/p > 5 and delta=0.95 otherwise (HM (p.184)).
# note:  This is Joe's version (from joe3.r) which contains a degrees of
#        freedom correction, Huber's correction (Huber, p.174), and a
#        standardization.
#
#  resd: is the vector of residuals.
#     p: is the number of regression coefficients (without the intercept), but
#        we use p+1 for the DF/Huber corrections.
# delta: Bandwidth parameter for estimate of tau.   Values  should be 
#        between .8 and .99.  Larger values result in a more conservative test.
#        If n/p is less than 5, recommend delta = .95.   Default is .8.
# param: This is a parameter for Huber's degrees of freedom correcton
#        as discussed in Huber (1981, p.174) and it is used in the estimator
#        of tau.  Smaller values lead to more conservative estimates of tau.
#        Default is 2.   It is the benchmark for labeling a standardized
#        residual as a potential outlier.
#
wilcoxontau<-function(resd,p,delta=if((length(resd)/p) > 5) 0.80 else 0.95,param=2){
  eps<-.000001
  n<-length(resd)
  temp<-pairup(resd)
  dresd<-sort(abs(temp[,1]-temp[,2]))
  dresd = dresd[(p+1):choose(n,2)]
  tdeltan<-quantile(dresd,delta)/sqrt(n)
  w<-rep(0,length(dresd))
  w[dresd <= tdeltan]<-1
  cn<-2/(n*(n-1))
  scores = sqrt(12)*((1:n)/(n+1) - .5)
  mn = mean(scores)
  con = sqrt(sum((scores - mn)^2)/(n+1))
  scores = (scores - mn)/con
  dn = scores[n] - scores[1]
  wilcoxontau<-sqrt(n/(n-p-1))*((2*tdeltan)/(dn*sum(w)*cn))
  w<-rep(0,n)
  stan<-(resd-median(resd))/mad(resd)
  w[abs(stan) < param]<-1
  hubcor<-sum(w)/n
  if(hubcor < eps){hubcor<-eps}
  fincor<-1 + (((p+1)/n)*((1-hubcor)/hubcor))
  wilcoxontau<-fincor*wilcoxontau
  names(wilcoxontau)<-NULL
  wilcoxontau
}
#
# Confidence interval estimate of taustar = 1/(2*f(0)).
# See, for example, HM (p.7-8 and p.25-26).
#
# resid: full model residuals.
#     p: is the number of regression coefficients (without the intercept), but
#        we use p+1 for the DF correction.
#  conf: confidence level of CI used.  Default set at 0.95.
#
taustar = function(resid,p,conf=.95){
n = length(resid)
zc = qnorm((1+conf)/2)
c1 = (n/2) - ((sqrt(n)*zc)/2) - .5
ic1 = floor(c1)
if(ic1 < 0){ic1=0}
z = sort(resid)
l = z[ic1+1]
u = z[n-ic1]
df = sqrt(n)/sqrt(n-p-1)
taustar = df*((sqrt(n)*(u-l))/(2*zc))
taustar 
}
###############################################################################
#                                                                             #
# Functions that calculate variance-covariance matrices                       #
#                                                                             #
###############################################################################
#
# var-cov matrix for GR-estimate.
# Given on page 288 of HM.  It does both GR and Nonrandom weights cases.
#
varcov.gr=function(x,bmat,res,delta=0.80) {
  x=as.matrix(x)
  xbar=as.matrix(apply(x,2,mean))
  bmat=as.matrix(bmat)
  res=as.vector(res)
  n=dim(x)[1]
  p=dim(x)[2]
  # Calculate W matrix
  diag(bmat)=rep(0,n)
  w=-1*bmat
  diag(w)=bmat%*%as.matrix(rep(1,n))
  w=(1/n)*w
  # Calculate inverse of C matrix
  cmat=(1/n)*t(x)%*%w%*%x
# cinv=solve(cmat)
  cinv=ginv(cmat) #generalized inverse when singular

  # Calculate V matrix
  vmat=(1/n)*t(x)%*%w%*%w%*%x
  # Calculate estimates of wilcoxontau and taustar 
  tau=wilcoxontau(res,p,delta)
  tau1=taustar(res,p)
  # Calculate variance-covariance matrix for GR
  varcov22=(tau^2/n)*cinv%*%vmat%*%cinv
  varcov12=-1*(t(xbar))%*%varcov22
  varcov11=(tau1^2/n) + (t(xbar))%*%varcov22%*%xbar
  varcov=cbind(rbind(varcov11,t(varcov12)),rbind(varcov12,varcov22))
  attr(varcov,"names")=NULL
  ans=list(varcov=varcov, tau1=tau1, tau=tau, wmat=w, cmat=cmat, vmat=vmat)
  ans
}
# 
# var-cov matrix for HBR-estimate.
# Uses Joe's method.  See, for example, Chang et al. (1999). High
# Breakdown Rank Regression.  JASA, 94(445), p.205-219.
#
varcov.hbr=function(x,bmat,res,delta=0.80) {
  x=as.matrix(x)
  xbar=as.matrix(apply(x,2,mean))
  bmat=as.matrix(bmat)
  res=as.vector(res)
  n=dim(x)[1]
  p=dim(x)[2]
  # Calculate estimates of tau and tau1
  tau=wilcoxontau(res,p,delta)
  tau1=taustar(res,p)
  # Calculate W matrix
  diag(bmat)=rep(0,n)
  w=-1*bmat
  diag(w)=bmat%*%as.matrix(rep(1,n))
  w=(1/(sqrt(12)*tau))*w
  # Calculate inverse of C matrix
  cmat=(1/n^2)*t(x)%*%w%*%x
  cinv=solve(cmat)
  # Calculate V matrix
  u=(1/n)*(bmat-diag(c(bmat%*%cbind(rep(1,n)))))%*%x
  u=u*(1-2*rank(res)/n)
  vmat=var(u)
  # Calculate variance-covariance matrix for HBR 
  varcov22=(1/(4*n))*cinv%*%vmat%*%cinv
  varcov12=-1*(t(xbar))%*%varcov22
  varcov11=(tau1^2/n) + (t(xbar))%*%varcov22%*%xbar
  varcov=cbind(rbind(varcov11,t(varcov12)),rbind(varcov12,varcov22))
  attr(varcov,"names")=NULL
  ans=list(varcov=varcov, tau1=tau1, tau=tau, wmat=w, cmat=cmat, vmat=vmat)
  ans
}
###############################################################################
#                                                                             #
# Functions that calculate Wald statistics and p-values                       #
#                                                                             #
# note:  The F distribution is used for the p-value.  See, for example,       #
#        p.151 and p.173 of HM.                                               #
#                                                                             #
# note:  'est' contains the intercept term so amat should as well.            # 
#                                                                             #
###############################################################################
#
# Wald Statistic and P-value
#
wald=function(est, varcov, amat, true, n) {
  true=as.matrix(true)
  est=as.matrix(est)
  amat=as.matrix(amat)
  p=dim(est)[1]-1
  q=dim(amat)[1]   
  # Calculate the wald test statistic and p-value
  temp1=as.matrix(amat%*%est - true)
  temp2=as.matrix(amat%*%varcov%*%t(amat))
  T2=t(temp1)%*%solve(temp2)%*%temp1
  T2=T2/q
  pvalue=1-pf(T2,q,n-p-1)
  #pvalue=1-pchisq(T2,q)
  c(T2,pvalue)
}
#
#         This routine obtains the pseudo-observations.
#         It assumes the Wilcoxon fit of y = x*beta.   These observations
#         are discussed on page 244 of Hettmansperger and McKean (1998).
#         They can be fed into a LS package as the observations.
#         Resulting contrast tests will be Wald type tests.
#
#         x          = Design martix without intercept.  The routine centers x.
#         y          = Response vector
#         delta      = parameter for estimate of tau
#         param      = Huber's correction for df correction for estimate of tau
#
wilcoxonpseudo = function(x,y,delta=.80,param=2){
  x = as.matrix(x)
  n = length(x[,1])
  p = length(x[1,])
  one = matrix(rep(1,n),ncol=1)
  x = x - (one%*%t(one)/n)%*%x
  tempw = wwest(x,y,"WIL")
  residw = tempw$tmp1$residuals
  fitw = y - residw
  arr = order(residw)
  jr = rep(0,n)
  for(i in 1:n){jr[arr[i]] = i}
  sc = sqrt(12)*((jr/(n+1)) - .5)
  zeta = sqrt((n-p-1)/sum(sc^2))
  tau = wilcoxontau(residw,p,delta,param)
  wilcoxonpseudo = fitw + tau*zeta*sc
  wilcoxonpseudo
}
###############################################################################
#                                                                             #
# Functions that perform (Wilcoxon based) testing for general linear          #
# hypotheses.                                                                 #
#                                                                             #
###############################################################################
#
# Function that performs the Drop in Dispersion test for Wilcoxon scores.
#
#   xmat    = Design matrix where intercept column has been stripped off.
#   y       = vector of responses.
#
#             Assumed model is  y = xmat%*%beta + e
#
#   amat    = Hypothesis matrix.  The hypothesis is
#
#             H_0:  amat%*%beta = 0  vs  H_A:  amat%*%beta neq 0
#
#             where beta is the vector of regression coefficients.
#
#   Assumes that amat has the same number of columns as xmat.
#
#   note:  Currently (12/9/03) this is only used with Wilcoxon
#          weights.
#
droptest<-function(xmat,y,amat,delta=.80,param=2,print.tbl=T){
  xmat = as.matrix(xmat)
  amat = rbind(amat)
  p = length(xmat[1,])
  pp1 = p+1
  n = length(xmat[,1])
  q = length(amat[,1])
  if (p!=dim(amat)[2])
    stop("droptest:  The number of columns in amat and xmat are different.")
  #
  #   Full model (Wilcoxon) fit.
  #
  full = wwfit(xmat,y)
  dfull = wildisp(full$residuals)
  tauhat = wilcoxontau(full$residuals,p,delta,param)
  #
  #  Reduced model (Wilcoxon) fit.
  #
  if(q < p) {
    xuse = xmat
    ause = amat
    xred = redmod(xuse,ause)
    red = wwfit(xred,y)
    dred = wildisp(red$residuals)
    rd = dred-dfull
    mrd = rd/q
    fr = mrd/(tauhat/2)}
  else {
    warning("droptest:  H_0: beta=0 is being tested since q>=p.",call.=F)
    q=p
    dred = wildisp(y)
    rd = dred-dfull
    mrd = rd/q
    fr = mrd/(tauhat/2)}
  df2 = n - p - 1
  ts2 = tauhat/2
  pval=1-pf(fr,q,df2)
  if (print.tbl) {
    cnames=c("RD","DF","MRD","TS","PVAL")
    rnames=c("H0","Error")
    ans=cbind(c(rd,NA),c(q,df2),c(mrd,ts2),c(fr,NA),c(pval,NA))
    ans=round(ans,4)
    dimnames(ans)=list(rnames,cnames)
    cat("\n")
    prmatrix(ans,na.print="")
    cat("\n")
  }
  invisible(list(full=full,dred=dred,dfull=dfull,tauhat=tauhat,q=q,fr=fr,pval=pval))
}
#   
#   Reduced model design matrix.
#
#   This function is called by the routine droptest.
#   Obtains the reduced model design matrix
#   for full model design matrix xmat and 
#   hypothesis matrix amat.  Checks to see
#   if amat has full row rank.
#   Uses Theorem 3.7.2 of H&M (1998, p. 186)
#
redmod<-function(xmat,amat){
  xmat=as.matrix(xmat)
  amat=rbind(amat)
  q<-length(amat[,1])
  p<-length(xmat[1,])
  temp<-qr(t(amat))
  if(temp$rank != q)
    stop("redmod:  The hypothesis matrix is not full row rank.")
  else
    {zed<-qr.qty(temp,t(xmat))
    redmod<-rbind(zed[(q+1):p,])}
  t(redmod)
}
#
#   Returns the Wilcoxon dispersion function at resid.
#
wildisp<-function(resid){
  n = length(resid)
  sresid = sort(resid)
  scores = sqrt(12)*((1:n)/(n+1) - .5)
  mn = mean(scores)
  con = sqrt(sum((scores - mn)^2)/(n+1))
  scores = (scores - mn)/con
  sum(scores*sresid)
}
#
#   (Wilcoxon) test of regression significance.
#
#   xmat    = Design matrix.  Assumes that xmat's intercept 
#             column has been stripped off.
#   y       = Vector of responses.
#
#   Assumed model is:
#                 y =  xmat%*%beta  + e*
#   where e* = alpha + e  (alpha is intercept).
#   Hypotheses are H_0: beta = 0 vs H_A:  beta neq 0
#
regrtest = function(xmat,y,delta=0.80,param=2,print.tbl=T) {
  xmat=as.matrix(xmat)
  p=dim(xmat)[2]
  ans=suppressWarnings(droptest(xmat,y,diag(rep(1,p)),delta,param,print.tbl))
  invisible(ans)
}
#
#      Returns the (Wilcoxon) drop in dispersion test results and a robust
#      ANOVA table for cell means model.
#
#      y      = vector of responses
#      levels = vector of indicators of levels
#               which are assumed to be 1 thru some positive
#               integer k.  Need not be sorted.
#               Levels are used to create a cell mean design matrix.
#      amat   = Hypotheses matrix corresponding to the
#               cell mean design matrix (rows must be contrast)
#      q      = number of constraints, (i.e. the number of rows of amat).
#
#      note:    for building the contrast matrix:  The vector of means is 
#               assumed to be (mu_1, mu_2, ... , mu_k).
#
cellmntest = function(y,levels,
                      amat=cbind(rep(1,max(levels)-1),-1*diag(max(levels)-1)),
                      delta=.80,param=2,print.tbl=T){
  amat=rbind(amat)
  xcell = cellmnxy(levels)
  p = length(xcell[1,])
  xmat = xcell[,2:p]
  amat = amat[,2:p]
  ans=suppressWarnings(droptest(xmat,y,amat,delta,param,print.tbl))
  ans$full$coef=ans$full$coef+c(0,rep(ans$full$coef[1],p-1))
  invisible(ans)
}
#
#   cell mean design matrix.
#
#   Column i corresponds to treatment i and
#   has a 1 where ever level=i. 
#
#   This routine assumes the levels are
#   1 through k.  They need not be sorted.
#
cellmnxy=function(levels){
  k = max(levels)
  n = length(levels)
  cellmnxy = matrix(rep(0,n*k),ncol=k)
  for(i in 1:k){cellmnxy[,i][levels==i]=1}
  cellmnxy
}
#
# This function conducts pairwise Wilcoxon-based drop in dispersion
# tests for the one-way analysis of variance model.
#
#      y: vector of response variables 
# levels: integer vector indicating group membership (must be 1,...,p)
#  delta: tuning parameter for scale estimate
#  param: tuning parameter for huber correction of scale
#
# For a given design the output consists of a matrix with one column.  The
# column represents the p-value for the drop in dispersion test of mu_i versus
# mu_j.  The order of the rows is given as (1,2), (1,3), ..., (1,p), (2,3),
# (2,4), ..., (2,p), ..., (p-1,p) where p represents the number of groups.
#
pwcomp<-function(y,levels,delta=.80,param=2){
  p<-max(levels)
  m<-pairup(1:p)
  rnames<-NULL
  pval<-NULL
  for(i in 1:dim(m)[1]) {
    a<-rep(0,p)
    a[m[i, 1]] <- 1
    a[m[i, 2]] <- -1
    rnames[i]<-paste("G",m[i,1],"-","G",m[i,2],sep="")
    pval[i]<-cellmntest(y,levels,a,delta=delta,param=param,
                        print.tbl=F)$pval
  }
  pval<-cbind(round(pval,4))
  dimnames(pval)<-list(rnames,"PVAL")
  pval
}
###############################################################################
#                                                                             #
# Functions that calculate studentized residuals                              #
#                                                                             #
###############################################################################
#
# This calculates the studentized residuals for the 
# GR estimate as given on page 300 of HM.  However,
# there are three typos on page 300 of HM.
#
#   1.  J = (1/n)11' not 11'
#   2.  K4 = sqrt(12)*tau*xi not sqrt(12)*xi/tau
#   3.  K5 = sqrt(12)*tau*tau_s*delta5 not sqrt(12)*tau_s*delta5/tau
#
# Question:  Does this match up with Joe's 'stanresid' when Wilcoxon weights 
#            are used?
#   Answer:  Very close when p=1, but some differences when p > 1.  Why?
#
studres.gr=function(x,bmat,res,delta=0.80,center=T) {
  x=as.matrix(x)
  if(center) {
    x=apply(x,2,function(x){x-mean(x)})
  }
  bmat=as.matrix(bmat)
  res=as.vector(res)
  n=dim(x)[1]
  p=dim(x)[2]
  # Calculate matrices
  # W matrix
  diag(bmat)=rep(0,n)
  w=-1*bmat
  diag(w)=bmat%*%as.matrix(rep(1,n))
  w=(1/n)*w
  # Kw, H, I, and J matrices
  Kw=x%*%solve(t(x)%*%w%*%x)%*%t(x)%*%w
  H=x%*%solve(t(x)%*%x)%*%t(x)
  I=diag(n)
  J=matrix(1/n,n,n)
  # Calculate estimates of constants
  sigma2=(mad(res))^2
  tau1=taustar(res,p)
  tau=wilcoxontau(res,p,delta)
  delta.s=(n/(n-p-1))*mean(abs(res))
  K3=2*tau1*delta.s-(tau1)^2
  tmp=pairup(res,"neq")
  xi=mean(tmp[,1]*sign(tmp[,1]-tmp[,2]))
  K4=sqrt(12)*tau*xi
  delta5=mean(sign(tmp[,1])*sign(tmp[,1]-tmp[,2]))
  K5=sqrt(12)*tau*tau1*delta5
  v=sigma2*I - K3*J - (K4*I - K5*J)%*%t(Kw) + (tau^2)*Kw%*%t(Kw)
  # In case any variances are negative...
  diag(v)[diag(v)<=0]=sigma2*diag(I-(1/n + H))[diag(v)<=0]
  as.vector(res/sqrt(diag(v)))
}
#
# This calculates the studentized residuals for the 
# HBR estimate as given on page 210 of Change et al.  
#
# Note: In (19) of Chang (1/n)*sum should be 2*(1/n)*sum.
#       See K2 below.
#
studres.hbr=function(x,bmat,res,delta=0.80,center=T) {
  x=as.matrix(x)
  if(center) {
    x=apply(x,2,function(x){x-mean(x)})
  }
  bmat=as.matrix(bmat)
  res=as.vector(res)
  n=dim(x)[1]
  p=dim(x)[2]
  # Calculate estimates of constants
  sigma2=(mad(res))^2
  tau1=taustar(res,p)
  tau=wilcoxontau(res,p,delta)
  K1=(n/(n-p-1))*mean(abs(res))
  K2=2*mean((rank(res)/(n+1) - 0.5)*res)
  # Calculate matrices
  H=x%*%solve(t(x)%*%x)%*%t(x)
  I=diag(n)
  J=matrix(1/n,n,n)
  # W matrix
  diag(bmat)=rep(0,n)
  w=-1*bmat
  diag(w)=bmat%*%as.matrix(rep(1,n))
  w=(1/(sqrt(12)*tau))*w
  # Calculate inverse of C matrix
  cmat=(1/n^2)*t(x)%*%w%*%x
  cinv=solve(cmat)
  # Calculate V matrix
  u=(1/n)*(bmat-diag(c(bmat%*%cbind(rep(1,n)))))%*%x
  u=u*(1-2*rank(res)/n)
  vmat=var(u)
  v=sigma2*I + tau1^2*J +
    (1/4)*x%*%((1/n^2)*cinv)%*%vmat%*%((1/n^2)*cinv)%*%t(x) -
    2*tau1*K1*J - sqrt(12)*tau*K2*(w%*%x%*%((1/n^2)*cinv)%*%t(x) +
    x%*%((1/n^2)*cinv)%*%t(x)%*%w)
  # In case any variances are negative...
  diag(v)[diag(v)<=0]=sigma2*diag(I-(1/n + H))[diag(v)<=0]
  as.vector(res/sqrt(diag(v)))
}
#
#    standardized residuals from the Wilcoxon fit
#    Discussed in McKean, Sheather and Hettmansperger (1990, JASA).
#
#    x       = (non-centered) design matrix
#    y       = responses
#    delta   = parameter for estimate of tau
#    param   = Huber's correction for df correction for estimate of tau
#    conf    = Confidence coefficient used in CI estimate of taustar
#
#    stanresid = Standardized residuals
#    ind     = Indicates if a standardization was negative with a 1.
#              Otherwise it is 0. If negative [mad^2*(1-h_i)]^.5
#              is used.
#
# Question: Does this match up with studres.gr (with Wilcoxon weights)?
#   Answer: Very close when p=1, but some differences when p > 1.  Why?
#     Note: Currently not being used.
#
stanresid = function(x,y,delta=.80,param=2,conf=.95){
# center x
xc = as.matrix(centerx(x))
n = length(y)
p = length(xc[1,])
pp1 = p+1
tempw = wwest(x,y,"WIL",print.tbl=F)
resid = tempw$tmp1$residuals
hc = diag(xc%*%solve(t(xc)%*%xc)%*%t(xc))
#
#    get taus
#
tau = wilcoxontau(resid,p,delta=.80,param=2)
taus = taustar(resid,p,conf=.95)
deltas = sum(abs(resid))/(n-pp1)
delta = wildisp(resid)/(n-pp1)
sig = mad(resid)
k1 = (taus^2/sig^2)*(((2*deltas)/taus)-1)
k2 = (tau^2/sig^2)*(((2*delta)/tau)-1)
s1 = sig^2*(1-(k1/n)-k2*hc)
s2 = s1
s2[s1 <= 0] = sig^2*(1-(1/n)-hc[s1 <= 0])
ind = rep(0,n)
ind[s1 <= 0] = 1
stanresid = resid/sqrt(s2)
list(stanr = stanresid,ind=ind,rawresids=resid,betaw=tempw$tmp1$coef,
tau=tau,taustar=taus)
}
###############################################################################
#                                                                             #
# Functions that perform diagnostics between different fits.                  #
#                                                                             #
###############################################################################
#
# Diagnostics (TDBETA and CFIT) that detect differences between fits.
#
# This function returns the diganostics (TDBETA and CFIT) between
# any two comparisons between LS, WIL, GR, and HBR.  All comparisons
# are standardized using the WIL fit.
# See, for example, McKean, Naranjo and Sheather in the articles: 
# (1996, Computational Stat., 223-243) and
# (1996,Comm. in Stat-Theory, 2575-2595) or section 5.5 of HM.
#
# x       = (non-centered)design matrix
# y       = responses
# est     = vector containing any two of "LS", "WIL", "GR", or "HBR"
# delta   = parameter for estimate of tau
# param   = Huber's correction for df correction for estimate of tau
# conf    = Confidence coefficient used in CI estimate of taustar
#
# tdbeta  = total difference in fit (standardized at WIL)
# bmtd    = Benchmark for tdbeta
# cfit    = Caswise difference diagnostics
# bmcf    = Benchmark for cfit
#
fitdiag=function(x,y,est=c("WIL","GR"),delta=0.80,param=2,conf=0.95){
  x=as.matrix(centerx(x))
  n=dim(x)[1]
  p=dim(x)[2]
  # Wilcoxon estimate
  tempw=wwest(x,y,"WIL",print.tbl=F)
  residw=tempw$tmp1$residuals
  tempvc=varcov.gr(centerx(x),tempw$tmp1$weights,tempw$tmp1$residuals)
  vcw=as.matrix(tempvc$varcov)
  # Other estimates, TDBETAs, and CFITs
  tempgr=NULL
  temphbr=NULL
  templs=NULL
  if (any("WIL"==est) & any("GR"==est)) {
    tempgr=wwest(x,y,"GR",print.tbl=F)
    diff=tempw$tmp1$coef - tempgr$tmp1$coef
  }
  if (any("WIL"==est) & any("HBR"==est)) {
    temphbr=wwest(x,y,"HBR",print.tbl=F)
    diff=tempw$tmp1$coef - temphbr$tmp1$coef
  }  
  if (any("GR"==est) & any("HBR"==est)) {
    tempgr=wwest(x,y,"GR",print.tbl=F)
    temphbr=wwest(x,y,"HBR",print.tbl=F)
    diff=tempgr$tmp1$coef - temphbr$tmp1$coef
  }
  if (any("WIL"==est) & any("LS"==est)) {
    templs=lsfit(x,y)
    diff=tempw$tmp1$coef - templs$coef
  }
  if (any("GR"==est) & any("LS"==est)) {
    tempgr=wwest(x,y,"GR",print.tbl=F)
    templs=lsfit(x,y)
    diff=tempgr$tmp1$coef - templs$coef
  }
  if (any("HBR"==est) & any("LS"==est)) {
    temphbr=wwest(x,y,"HBR",print.tbl=F)
    templs=lsfit(x,y)
    diff=temphbr$tmp1$coef - templs$coef
  }
  tdbeta=t(cbind(diff))%*%solve(vcw)%*%cbind(diff)
  bmtd=(4*(p+1)^2)/n
  xmat=cbind(rep(1,n),x)
  diffc=xmat%*%diff
  diffvc=xmat%*%vcw%*%t(xmat)
  cfit=diffc/(sqrt(diag(diffvc)))
  bmcf=2*sqrt((p+1)/n)
  se=sqrt(diag(vcw))  
  list(tdbeta=c(tdbeta),bmtd=bmtd,cfit=c(cfit),bmcf=bmcf,est=est,
  betaw=tempw$tmp1$coef,betagr=tempgr$tmp1$coef,betahbr=temphbr$tmp1$coef,
  betals=templs$coef,
  vcw=vcw,tau=tempvc$tau,taus=tempvc$tau1,se=se)
}
#
# plotdiagfit
#   plots the results from 'diagfit'.
#
plotfitdiag=function(result) {
  n=length(result$cfit)
  main1=paste("CFITS for",result$est[1],"and",result$est[2])
  main2=paste("TDBETA:",round(result$tdbeta,2),"Benchmark:",
            round(result$bmtd,2))
  plot(c(1,n),c(min(result$cfit,-1*result$bmcf),
    max(result$cfit,result$bmcf)),type="n",
    main=paste(main1,"\n",main2),xlab="CASE",ylab="CFIT")
  points(1:n,result$cfit)
  abline(h=c(-1*result$bmcf,result$bmcf))
}
#
#    Returns the diganostics between the Wilcoxon and LS fits
#    See, for example, McKean, Naranjo and Sheather (1999,
#    J. Nonparametrics, 161-188).
#
#    x       = design matrix
#    y       = responses
#    delta   = parameter for estimate of tau
#    param   = Huber's correction for df correction for estimate of tau
#    conf    = Confidence coefficient used in CI estimate of taustar
#
#    tdbeta  = total difference in fit (standardized)
#    bmtd    = Benchmark for tdbeta
#    tdint   = Difference in intercepts (standardized).
#    cd      = Caswise differences diagnostics.
#    bmcd    = Benchmark for casewise 
#    
diffwls = function(x,y,delta=.80,param=2,conf=.95){
x=as.matrix(centerx(x))
n = length(x[,1])
p = length(x[1,])
# Wilcoxon estimate
tempw = wwest(x,y,"WIL",print.tbl=F)
residw = tempw$tmp1$residuals
tempvc=varcov.gr(centerx(x),tempw$tmp1$weights,tempw$tmp1$residuals)
vcw = as.matrix(tempvc$varcov)
# Least Squares estimate
templs = lsfit(x,y)
diff = tempw$tmp1$coef - templs$coef
vcwint = vcw[1,1]
pp1 = length(x[1,]) + 1
vcwbeta = vcw[2:pp1,2:pp1]
tdbeta = t(diff)%*%solve(vcw)%*%diff
tdint = diff[1]^2/vcwint
bmtd = (4*pp1^2)/n
xmat = cbind(rep(1,n),x)
diffc = xmat%*%diff[1:pp1]
diffvc = xmat%*%vcw%*%t(xmat)
cd = diffc/(sqrt(diag(diffvc)))
bmcd = 2*sqrt(pp1/n)
se = sqrt(diag(vcw))
list(tdbeta=tdbeta,tdint=tdint,bmtd=bmtd,cfit=cd,bmcd=bmcd,est=c("WIL","LS"),
betaw=tempw$tmp1$coef,betals=templs$coef,vcw=vcw,tau=tempvc$tau,taus=tempvc$tau1,se=se)
}
###############################################################################
#                                                                             #
# Miscellaneous functions                                                     #
#                                                                             #
###############################################################################
#
# centerx
#
centerx = function(x){
  x = as.matrix(x)
  n = length(x[,1])
  one = matrix(rep(1,n),ncol=1)
  x - (one%*%t(one)/n)%*%x 
}
# 
# Function used to check for correct version number(s)
#
v1.9.0=function() {
  major=version$major
  minor=version$minor
  n=as.numeric(paste(major,minor,sep=""))
  n>=19
}
#
# 'print.matrix' command for R 2.4.0 (supperseeded by 'prmatrix')
#
#print.matrix = function (x, digits = NULL, quote = TRUE, na.print = NULL, print.gap = NULL, 
#                         right = FALSE, max=NULL, ...) { 
#    noOpt = missing(digits) && missing(quote) && missing(na.print) && 
#            missing(print.gap) && missing(right) && length(list(...)) == 0
#    .Internal(print.default(x, digits, quote, na.print, print.gap, right, max, noOpt))
#}
