#    standardized residuals from the Wilcoxon fit for GR
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
#     Note: try!
#

stanresidgr = function(x,y,resid,delta=.80,param=2,conf=.95){
# center x
xc = as.matrix(centerx(x))
n = length(y)
p = length(xc[1,])
pp1 = p+1
#tempw = wwest(x,y,"WIL",print.tbl=F)
#resid = tempw$tmp1$residuals
hc = diag(xc%*%ginv(t(xc)%*%xc)%*%t(xc))	#use solve instead of ginv
#
#    get taus
#
tau = wilcoxontau(resid,p,delta=.80,param=2) #param ?
taus = taustar(resid,p,conf=.95)
deltas = sum(abs(resid))/(n-pp1)
delta = wildisp(resid)/(n-pp1)
sig = mad(resid) #try dispvar for sig
k1 = (taus^2/sig^2)*(((2*deltas)/taus)-1)
k2 = (tau^2/sig^2)*(((2*delta)/tau)-1)
s1 = sig^2*(1-(k1/n)-k2*hc)
s2 = s1
s2[s1 <= 0] = sig^2*(1-(1/n)-hc[s1 <= 0])
ind = rep(0,n)
ind[s1 <= 0] = 1
stanresid = resid/sqrt(s2)
list(stanr = stanresid,ind=ind,rawresids=resid,tau=tau,taustar=taus)
}
