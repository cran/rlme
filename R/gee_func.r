weightf <- function(ehats,ahats,med){
w=ahats/(ehats-med)	
w[is.infinite(w)] <- 0 
w[is.nan(w)] <- 0	
w[w==0] <- max(w) 
list(w=w)
}

# geels1 <- function(x,ehat,r){
# m0=solve(t(x)%*%r%*%x,tol=tol)
# vary=ehat%*%t(ehat)+ x%*%m0%*%t(x) #this is a bias correction bias=x%*%solve(crossprod(xstar))%*%t(x)and vary=vary+bias
# m1=(t(x)%*%r%*%vary%*%r%*%x) 
# varb=m0%*%m1%*%m0
# se1=sqrt(diag(varb))	#se with var(Y) #df=n-p-(I-1)-(9-J)=16 in mixed, n-p in lm so t dist will be affected..
# #dd[is.nan(dd)]=0
# list(se=se1)
# }
# 
# sisinp_sum <- function(ahats,school){
# xx<-as.matrix(cbind(school,ahats))
# aa<-unique(school)
# nc<-length(aa)
# for( j in 1:nc ) {
# 	x1<-as.matrix(xx[school==aa[j],])
# 	n1<-nrow(x1)
# 	xx[school==aa[j],]<-x1[,2]-matrix(rep(1,n1),ncol=1)%*%sum(x1[,2])
# }
# si=xx[,1]
# list(si=si)
# }
# 
# sisinp_med <- function(ahats,school){
# xx<-as.matrix(cbind(school,ahats))
# aa<-unique(school)
# nc<-length(aa)
# for( j in 1:nc ) {
# 	x1<-as.matrix(xx[school==aa[j],])
# 	n1<-nrow(x1)
# 	xx[school==aa[j],]<-x1[,2]-matrix(rep(1,n1),ncol=1)%*%median(x1[,2])
# }
# si=xx[,1]
# list(si=si)
# }
# 
# pcalc <- function(pval1){
# for(k in 1:pp){
# print(mean(pval1[,k]<.05,na.rm=T))
# }
# }
