LM_est=function(x, y, dat, method="REML"){s
  model = as.formula(paste("y ~ 1 + ", paste(colnames(x), collapse=' + ')))
  
  fit.lme = lme(model, data=dat, random=~1|school/section, method = method) #method = c("REML", "ML")
                                          
  theta0 <-extract.lme.cov2(fit.lme, dat, start.level=3)$V[1]	#err
  theta2 <-extract.lme.cov2(fit.lme, dat, start.level=2)$V[[1]][1,1]-theta0	#sect
  theta1 <-extract.lme.cov2(fit.lme, dat, start.level=1)$V[[1]][1,1]-(theta0+theta2) #school
  sigma.l <-c(theta1,theta2,theta0) 	#sigma school,sec,error
  theta <- as.vector(summary(fit.lme)$tTable[,1]) #Beta estimates
  intra_err.lm<-(theta0)/(sum(sigma.l))		#intra error
  intra_sch.lm<-(theta1)/(sum(sigma.l)) 	 	#intra school or rho2
  intra_sect.lm<-(theta1+theta2)/(sum(sigma.l)) 	#intra section or rho1
  ses <- as.vector(summary(fit.lme)$tTable[,2])	#df=61=72-(3-1+9-3+3)=n-(I-1+sum(Ji)-I+pp)=(sum(Ji)+pp-1)
  
  ehat <- as.vector(fit.lme$residuals[,1])  #get all with fit.lme$resid  #try as.vector(fit.lme$residuals[,1])
  varb=fit.lme$varFix 
                                          
effect_sch=random.effects(fit.lme,level=1)[,1] #school random effects
effect_sec=random.effects(fit.lme,level=2)[,1] #section in school random effects
effect_err = as.vector(fit.lme$residuals[,3])

  standr.lme=residuals(fit.lme, type="pearson") #st res
  standr.lme=as.vector(residuals(fit.lme, type="pearson"))
                                          
  list(theta=theta, ses=ses, varb=varb, sigma=sigma.l, ehat = ehat, effect_sch = effect_sch, effect_sec = effect_sec, effect_err = effect_err, standr.lme=standr.lme)
}

#############JR
JR_est=function(x,y,I,sec,mat,school,section, rprpair = 'hl-disp'){

pp <- dim(x)[2]+1	#number of X+intec, dim(x)[2]
wilfit = wilonestep(y,x)
theta = wilfit$theta[1:(pp)]  #fixed effects
ehat=wilfit$ehat 
ahat =scorewil(ehat)$scorewil
taus=wilfit$taus		#for intercept, L1
tauhat=wilfit$tauhat	#wilc for Beta dist

##sigmas and intras (rho1=intra_sect and rho2=intra_sch needed in stud res.)
scale.fit=rprmeddis(I,sec,mat,ehat,location,scale, rprpair = rprpair)
sigma=c(scale.fit$siga2,scale.fit$sigw2,scale.fit$sigmae2) #sigma sch, sect, err ?sure
effect_sch = scale.fit$frei	#school random eff
effect_sec = scale.fit$frew	#sec random eff
effect_err = scale.fit$free

#Kloke's: intercept's var matrix var_alpha
var_alpha=interc_se(x,ehat,school,section,taus)$var_alpha

#"WSM" rho's
rho1 <- rhosect(ahat,school,section)  				# rho matrix for phi for section (rho_phi_sect)    
rho1_est_9<- sum(sum( (apply(rho1$rho2,1,sum,na.rm=T)/apply(rho1$npair,1,sum) ) * t(apply(rho1$mat,1,sum))) )/ sum(sum(rho1$mat)) # # weighted Simple Moment:WSM (with sample size)
rho2 <- rhosch(ahat,school,section) 	# rho matrix for phi for school (rho_phi_sch)
rho2_est_5 <- sum( ( (rho2$rho2 / rho2$npair ) *  rho2$nvec )/ sum(rho2$nvec),na.rm=T )		# weighted Simple Moment:WSM
v1=1
v2=rho1_est_9
v3=rho2_est_5

#fixed parameters' var  matrix V. no intercept in V as cob(beta_est)
V=beta_var(x, school, tauhat, v1, v2, v3, section)$var #in JR, Fixed effect (Beta) est distribution's variance est
ses <- c(sqrt(var_alpha), sqrt(diag(V)) )

theta <- as.vector(theta) #fixed parameter estimate
sigma <-sigma

list(theta=theta, ses=ses, varb=V, sigma=sigma, ehat = ehat, effect_sch = effect_sch, effect_sec = effect_sec, effect_err = effect_err )
}



##GR estimate ?this is Iter GR. Get GR only with numstp = 2
GR_est=function(x, y, I, sec, mat, school, section, rprpair = 'hl-disp'){

#I=school, section, mat is size matrix for sections
init = T
sigmaa2=1
sigmaw2=1
sigmae2=1
thetaold=c(0)
numstp = 50
eps = .0001
iflag2 = 0
i = 0
is = 0
J = sum(sec)
n=length(y)

if(is.null(dim(x)[2])){
nopar_scale=1+1+3
} else {
nopar_scale=dim(x)[2]+1+3		#1 intrcp + 2 X +3 sigma parameters
}

coll = matrix(rep(0,nopar_scale),ncol=nopar_scale)	
collresch = matrix(rep(0,I),ncol=I)
collresec = matrix(rep(0,J),ncol=J)
collreserr = matrix(rep(0,n),ncol=n)

while(is == 0){
   i = i + 1
   if(i == 1){
      wilfit = wilstep(I,sec,mat,init,y,x,sigmaa2,sigmaw2,sigmae2,thetaold,eps,iflag2, rprpair = rprpair)
      coll = rbind(coll,c(wilfit$theta,wilfit$sigmaa2,wilfit$sigmaw2,wilfit$sigmae2))


effect_sch = wilfit$rea	#school random eff
effect_sec = wilfit$rew	#sec random eff
effect_err = wilfit$ree

      thetaold = wilfit$theta	#why here??
	ehat = wilfit$ehat
	theta=coll[dim(coll)[1],1:(dim(x)[2]+1)]			#fixed int, trtm, cov
	sigma=c(wilfit$sigmaa2,wilfit$sigmaw2,wilfit$sigmae2) 	#sigma sch, sect, err
   }

sigmaa2 = wilfit$sigmaa2
sigmaw2 = wilfit$sigmaw2
sigmae2 = wilfit$sigmae2

init = F

wilfit = wilstep(I,sec,mat,init,y,x,sigmaa2,sigmaw2,sigmae2,thetaold,eps,iflag2, rprpair = rprpair)
if((wilfit$iflag2==1) | (i > numstp)){is=1}
coll = rbind(coll,c(wilfit$theta,wilfit$sigmaa2,wilfit$sigmaw2,wilfit$sigmae2)) #this is storage for fixed and scale estimates

effect_sch = wilfit$rea	#school random eff
effect_sec = wilfit$rew	#sec random eff
effect_err = wilfit$ree

thetaold = wilfit$theta #to stop for the next iter, we need it.
ehat = wilfit$ehat
theta=coll[dim(coll)[1],1:(dim(x)[2]+1)]			#fixed int, trtm, cov
sigma=c(wilfit$sigmaa2,wilfit$sigmaw2,wilfit$sigmae2) 	#sigma sch, sect, err
i
}


#now se's
#fixed parameters' var  matrix V
sigmay = sigymake(I,sec,mat,sigma[1],sigma[2],sigma[3]) #sigmay$sigy2 is cov(e)
sigma12inv = matrix(sigmay$sigy12i,ncol=length(y))
#theta
xstar = sigma12inv%*%cbind(1,x)
ystar = sigma12inv%*%y  
ehats	= ystar-xstar%*%theta #this is ehat* Please confirm
ahats = scorewil(ehats)$scorewil

ehat	= y-cbind(1,x)%*%theta #this is ehat* Please confirm

taus=taustar(ehats,p=dim(x)[2],conf=.95)
tauhat=wilcoxontau(ehats,p=dim(x)[2]) #?? check if 2 or 3

XXinv <- solve(crossprod(xstar))

x_c=xstar-apply(xstar,2,mean)	#centering

aa=taus^2*(XXinv%*%t(xstar))%*%{rep(1,n)%*%solve(t(rep(1,n))%*%rep(1,n))%*%rep(1,n)}%*%(xstar%*%XXinv)
bb=tauhat^2*XXinv%*%t(xstar)%*%{x_c%*%solve(t(x_c)%*%x_c)%*%t(x_c)}%*%(xstar%*%XXinv)

varb=aa+bb   #cov(beta_est): (p+1)x(p+1)
ses <- sqrt(diag(varb)) 

list(theta=theta, ses=ses, sigma=sigma, varb=varb, ehat = ehat, ehats = ehats, effect_sch = effect_sch, effect_sec = effect_sec, effect_err = effect_err, iter=i, coll=coll, xstar=xstar, ystar=ystar)
}

############
#GEE
GEER_est=function(x,y,I,sec,mat,school,section, weight="wil", rprpair = 'hl-disp'){
  
weight = tolower(weight)
if(weight == "wil") {
  weight = 1
}
if(weight == "hbr") {
  weight = 2
}

n = sum(mat)
  
#when weight=2, it uses hbr weight in gee. when 1, it is same as the gee paper weight (phi/(e-med))
#initial b estimates
#k=0
  
fitw = wwest(x,y,print.tbl=F)
b0 = fitw$tmp1$coef

ehat0=y-cbind(1,x)%*%b0 

fitvc = rprmeddis(I,sec,mat,ehat=ehat0,location,scale, rprpair = rprpair) #??ask if ok rprpair.
sigmaa2 = fitvc$siga2
sigmaw2 = fitvc$sigw2
sigmae2 = fitvc$sigmae2

x=cbind(1,x)
collb <- b0 
collsigma <- c(sigmaa2,sigmaw2,sigmae2) #var estimates/weights

iter<-0
chk <-1
b2<-b0
max.iter<-2 #try up to 50.
ww <- 0

while( (chk > 0.0001) && (iter < max.iter)  )  {	

iter<-iter+1
b1<-b2	#b1 is from wwest initial

#k=1
sigmay = sigymake(I,sec,mat,sigmaa2,sigmaw2,sigmae2)
sigma12inv = matrix(sigmay$sigy12i,ncol=n) #check if different
siggma12<-matrix(sigma12(sigmay$sigy2),ncol=n) ##siggma12 <- sig*diag(1,n) ##
ystar=sigma12inv %*%y
xstar=sigma12inv %*%x
ehats=ystar-xstar%*%b2
med <- median(ehats)
ahats=scorewil(ehats)$scorewil

yss=y-siggma12%*%rep(1,n)*med
ys=ystar-siggma12%*%rep(1,n)*med

if(weight==1){
w <- weightf(ehats,ahats,med)$w	
}
if(weight==2){
w <- hbrwts_gr(xstar,ehats) #hbr
}

ww <- rbind(ww,t(w)) 
w <- diag(as.vector(w))
r <- sigma12inv %*% w %*% sigma12inv 
b2 <- solve(t(x)%*%r%*%x,tol=tol)%*%(t(x)%*%r%*%yss) 
ehat=y-x%*%b2

fitvc <- rprmeddis(I,sec,mat,ehat=ehat,location,scale, rprpair = rprpair)
sigmaa2 = fitvc$siga2
sigmaw2 = fitvc$sigw2
sigmae2 = fitvc$sigmae2
chk <- sum((b2 - b1)^2)/sum(b1^2)
collb <- rbind(collb,c(b2))
collsigma <- rbind(collsigma,c(sigmaa2,sigmaw2,sigmae2))
}

iter=iter+1
b <- collb[iter,]	#gee-wr beta
theta=b

ehat=y-x%*%b #lm(y~x-1)$coef
ahat=scorewil(ehat)$scorewil
tauhat=wilcoxontau(ehat,p=dim(x)[2])

fitvc <- rprmeddis(I,sec,mat,ehat=ehat,location,scale, rprpair = rprpair)
sigmaa2 = fitvc$siga2
sigmaw2 = fitvc$sigw2
sigmae2 = fitvc$sigmae2
sigma  <- c(sigmaa2,sigmaw2,sigmae2)

#prediction/estimate errors: e=sch+sec+err
effect_sch = fitvc$frei	#school random eff
effect_sec = fitvc$frew	#sec random eff
effect_err = fitvc$free #epsilon

sigmay = sigymake(I,sec,mat,sigmaa2,sigmaw2,sigmae2)
sigma12inv = matrix(sigmay$sigy12i,ncol=n) 
siggma12<-matrix(sigma12(sigmay$sigy2),ncol=n) ##siggma12 <- sig*diag(1,n) ##
if(weight==1){
w <- weightf(ehats,ahats,med)$w	#weights standardized
}
if(weight==2){
w <- hbrwts_gr(xstar,ehats) #hbr
}
w <- diag(as.vector(w))
r <- sigma12inv %*% w %*% sigma12inv  

ystar=sigma12inv %*%y
xstar=sigma12inv %*%x
ehats=ystar-xstar%*%b 
ahats=scorewil(ehats)$scorewil
tauhats=wilcoxontau(ehats,p=dim(x)[2])

##NP way to get ses
#se with NP1 p.10, as in gee paper: SiSi from ai-sum(ai),w
#si=sisinp_sum(ahats,school)$si
#sisi=si%*%t(si)
m01=solve(t(x)%*%r%*%x,tol=tol)
#m11=(t(x)%*%sigma12inv%*%sisi%*%sigma12inv%*%x)
#varb=m01%*%m11%*%m01
#se21=sqrt(diag(varb))	#se with Si

#se22 with NP2 p.10 SiSi from ai, w
#sisi2=ahats%*%t(ahats)
m03=solve(t(x)%*%sigma12inv%*%sigma12inv%*%x,tol=tol)
#m12=(t(x)%*%sigma12inv%*%sisi2%*%sigma12inv%*%x)
#varb=m01%*%m12%*%m01
#se22=sqrt(diag(varb))	#se with Si

#se with NP3 p.10: tauhats^2, SiSi from ai-sum(ai),w=I/t
#m13=(t(x)%*%sigma12inv%*%sisi%*%sigma12inv%*%x)
#varb=tauhats^2*m03%*%m13%*%m03
#se23=sqrt(diag(varb))	#se with Si

##AP ways to get se, use this one by default.
#se31 with AP1 p.10, gee paper: sisi=I, w from i/taus (from y is used?)
sisi3=diag(rep(1,n))
m15=(t(x)%*%sigma12inv%*%sisi3%*%sigma12inv%*%x)	#actually, paper replace weigths with 1/tauh 
varb=tauhats^2*m03%*%m15%*%m03 
se31=sqrt(diag(varb))	#se with Si, same as 23
###check m03=m15 so varb=tauhats^2*m03

#se with AP2 p.10, gee paper:w from calc, sisi=I (so might need ys when w is used?)
#varb=m01%*%m15%*%m01	#or tauhat^3 ??
#se32=sqrt(diag(varb))	#se with Si

#CS ways to get se(beta). 
#se with CS1 for var(phi) as in JR: SiSi=CS, w 
rho1=rhosect(ahats,school,section)  	# rho matrix for phi for section (rho_phi_sect)    
rho1_est_9=sum(sum( (apply(rho1$rho2,1,sum)/apply(rho1$npair,1,sum) ) * t(apply(rho1$mat,1,sum))) )/ sum(sum(rho1$mat)) # # weighted Simple Moment:WSM (with sample size)
rho2=rhosch(ahats,school,section) 	# rho matrix for phi for school (rho_phi_sch)
#fix this part, remove when school has one section
aaa=( (rho2$rho2 / rho2$npair ) *  rho2$nvec )
aaa=aaa[!is.na(aaa)]
rho2_est_5 = sum( aaa/ sum(rho2$nvec) )  	# weighted Simple Moment:WSM
#rho2_est_5 = sum( ( (rho2$rho2 / rho2$npair ) *  rho2$nvec )/ sum(rho2$nvec) )		# weighted Simple Moment:WSM
v1=1
v2=rho1_est_9
v3=rho2_est_5

#since nested levels affect, it is reflected in CS design
sisi4=0

for (i in unique(school)) {
	sisi4=adiag(sisi4, Bmat_sch(v1, v2, v3,section[school==i]))
		}
sisi4=sisi4[2:dim(sisi4)[1],2:dim(sisi4)[2]]

sisi5=matrix(sisi4,ncol=n)
m16=(t(x)%*%sigma12inv%*%sisi4%*%sigma12inv%*%x)

 if(weight==2){
 varb=m01%*%m16%*%m01
 se41=sqrt(diag(varb))	#se with Si. this is good with hbr in gee. I called CS2. Please modify for I-school case.
 }

#se with CS2 for var(phi) as in JR: SiSi = CS, w from I/tauhats (this should be *)
 if(weight==1){
varb=tauhats^2*m03%*%m16%*%m03
se41=sqrt(diag(varb))	#se with Si.i used this in thesis as CS.originally it was se42. i changed 1/1/2013.
 }

list(theta=theta, ses_AP=se31, ses_CS=se41, varb=varb, sigma=sigma, ehat = ehat, effect_sch = effect_sch, effect_sec = effect_sec, effect_err = effect_err, iter=iter, w=diag(w))

}
