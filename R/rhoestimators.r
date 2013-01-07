##############functions#########################################
## Three-level nested design functions
## B block matrices with rho1 and rho2 in school i
## Some functions and notations are modified and extended version of jfit by Kloke.

Bmat_sch <-
function (v1, v2, v3, section) 
{
#v1=rho_sect, v2=rho_sch for error of phi depending on the context
    vblock <- unique(section)
    m <- length(vblock)
    N <- length(section)
    B <- matrix(0, nrow = N, ncol = N)
    ind <- 1
    for (k in 1:m) {
        nk <- sum(section == vblock[k])
		if(nk==1){
     		B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- 1
		} else {
        	Bk <- diag(rep(v1 - v2, nk)) + v2
        	B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- Bk
		}
      	ind <- ind + nk
    }
    B[B==0] <- v3
    invisible(B)
}



##Fixed effect (Beta) est distribution's variance est
beta_var <- 
function (x, school, tauhat, v1, v2, v3, section) 
{
    x <- as.matrix(x)
    ublock <- unique(school)
    I <- length(ublock)
    p <- ncol(x)
    sig2 <- matrix(rep(0, p^2), ncol = p)
    for (i in 1:I) {
        x_i <- as.matrix(x[school == ublock[i], ])
        n_i <- nrow(x_i)

	  s_i <- Bmat_sch(v1, v2, v3, section[school == ublock[i]])
        
        sig2 <- sig2 + t(x_i) %*% s_i %*% x_i
    }

    xx_inv <- solve(crossprod(x))
    V <- tauhat^2 * xx_inv %*% sig2 %*% xx_inv
list(var=V)
}


##number of sections in each school for 3-level
sec_vec <-
function (school,section) {

    ublock <- unique(school)
    I <- length(ublock)
    vvec <- matrix(0, nrow=I)

    for (i in 1:I){

   	vblock <- unique(section[school==ublock[i]])
   	J <- length(vblock)
 	vvec[i] <- J
}
return(vvec)
#so this vector is called 'sec' in the simulations
}


##number of sections in each school for 3-level
sch_vec <-
function (school,section) {

    ublock <- unique(school)
    I <- length(ublock)
    vvec <- as.vector(sec_vec(school,section))
    sctn=rep(0, length(school))
    coll=0

    	for (i in 1:length(vvec)){
	vblock <- unique(section[school==ublock[i]])
      cols <- rep(0, vvec[i])

		for (j in 1:vvec[i]){
		zz=section[school == ublock[i]] == vblock[j]
		cols[j]=length(zz[zz==TRUE])
		}
	sctn=rep(1:vvec[i],cols)
	coll=cbind(coll,t(sctn))
	}
coll=coll[-1]		
return(coll)
#so this vector is called 'section' in the simulations
}


##section and school vector in 3-level data
sec_sch_vec <-
function (y,school1,section1) {

one= 1+0*y
schsize=aggregate(one, by=list(y), FUN=sum)[,2]	#this gives school sizes
school=factor(rep(1:length(schsize),schsize))	#this generates school nos as I did alternative

sss=aggregate(one, by=list(school1, section1), FUN=sum)[,3]	#as alternative to my ss in main.begin. try
section=factor(rep(1:length(sss),sss))
#section=factor(rep(1:sss[i],sss[i]))


    ublock <- unique(school)
    I <- length(ublock)
    vvec <- as.vector(sec_vec(school,section))
    sctn=rep(0, length(school))
    coll=0

    	for (i in 1:length(vvec)){
	vblock <- unique(section[school==ublock[i]])
      cols <- rep(0, vvec[i])

		for (j in 1:vvec[i]){
		zz=section[school == ublock[i]] == vblock[j]
		cols[j]=length(zz[zz==TRUE])
		}
	sctn=rep(1:vvec[i],cols)
	coll=cbind(coll,t(sctn))
	}
coll=coll[-1]		
return(coll)
#so this vector is called 'school' in the simulations
}


#mat matrix
mat_vec <-
function (school,section) 
{

    ublock <- unique(school)
    I <- length(ublock)
 
    mat <- matrix(0, nrow=I, ncol=max(as.numeric(section)))

    for (i in 1:I) {
	vblock <- unique(section[school==ublock[i]])

      a_i <- school[school == ublock[i]]	#score for school i
	a_i <-a_i[!is.na(a_i)]

 	for (j in 1:length(vblock)) {
         a_ij <- a_i[section[school == ublock[i]] == vblock[j]]	#score for section j in school i
	   a_ij <-a_ij[!is.na(a_ij)]
	   mat[i,j] <- length(a_ij)

	}
    }
return(mat)
}



##3-way nested design rho_phi for section calculation
rhosect <-
function (ahat,school,section) 
{
#Only rho2 is used. this is Moment Estimates.
  
    ublock <- unique(school)
    I <- length(ublock)
 
    nvec <- matrix(0, nrow=I, ncol=max(as.numeric(section)))
    rho1_vec <- rho2_vec <- rho3_vec <- rho4_vec <- nvec #rho1 is mine, rho2 is Kloke's way
    
    for (i in 1:I) {
	vblock <- unique(section[school==ublock[i]])

      a_i <- ahat[school == ublock[i]]	#score for school i
	a_i <-a_i[!is.na(a_i)]


 	for (j in 1:length(vblock)) {
         a_ij <- a_i[section[school == ublock[i]] == vblock[j]]	#score for section j in school i
	   a_ij <-a_ij[!is.na(a_ij)]

	   nvec[i,j] <- length(a_ij)
	   pair=as.matrix(pairup1(a_ij)$ans)

		 if(dim(pair)[2]==1){
	   rho1_vec[i,j]= 0  #it was NA.
	   rho2_vec[i,j]= 0
	   rho3_vec[i,j]= 0
	   rho4_vec[i,j]= 0
				} else {
	   rho1_vec[i,j]=as.numeric(cor(pair[,1],pair[,2])) 			#just correlation way
	   rho2_vec[i,j]= sum(pair[,1]*pair[,2])					#simple moment way 
#	   rho3_vec[i,j]= (1/length(pair[,1]))*(sum(pair[,1]*pair[,2])-length(pair[,1])*mean(pair[,1])*mean(pair[,2]) ) /(sd(pair[,1]-mean(pair[,1]))*sd(pair[,2]-mean(pair[,2])))  #subtract section average
	   rho3_vec[i,j]= (1/length(pair[,1]))*(sum(pair[,1]*pair[,2])-length(pair[,1])*mean(a_ij)^2 ) /( sd(pair[,1]-mean(a_ij))*sd(pair[,2]-mean(a_ij)) )  #subtract section average

#	   rho4_vec[i,j]= (1/length(pair[,1]))*(sum(pair[,1]*pair[,2])-length(pair[,1])*mean(pairup1(a_i)$ans[,1])*mean(pairup1(a_i)$ans[,2]) ) /(sd(pair[,1]-mean(pairup1(a_i)$ans[,1]))*sd(pair[,2]-mean(pairup1(a_i)$ans[,2])) )  #subtract school average
	   rho4_vec[i,j]= (1/length(pair[,1]))*(sum(pair[,1]*pair[,2])-length(pair[,1])*mean(a_i)^2 ) / ( sd(pair[,1]-mean(a_i) )*sd(pair[,2]-mean(a_i)) )  #subtract school average
	 }
	}
    }
list(rho1=rho1_vec, rho2=rho2_vec, rho3=rho3_vec, rho4=rho4_vec, mat=nvec, npair=choose(nvec,2))
}



#######rho_phi for school##
rhosch <-
function (ahat,school,section) 
{
  #Only rho2 is used. this is Moment Estimates.
  
    ublock <- unique(school)
    I <- length(ublock)
 
    nvec <- vector(I, mode = "numeric")
    rho1_vec <- rho2_vec <- rho3_vec <- rho4_vec <-npair <- nvec 

    for (i in 1:I) {
	vblock  <- unique(section[school==ublock[i]])
      m <- choose(length(vblock),2)  
#	rho1_vecc <- rho2_vecc <- vector(m, mode = "numeric") 

	nn <- 0
      rho1 <- rho2 <- rho3 <- rho4 <-0

      a_i <- ahat[school == ublock[i]]	#scores for school i

 	for (j in 1:(length(vblock)-1)) {
	   for (k in (j+1):length(vblock)) {
			
	     #fixed 1/25/2013
	     if(length(vblock)==1){   a_ij <- a_i[section[school == ublock[i]] == vblock[j]]  #scores for section j in school i
	                              a_ik <- a_ij
	                              pair=1 #as.matrix(pairup2(a_ij,a_ik)$ans)
	                              nn <- 0 #nn+pairup2(a_ij,a_ik)$n	#number of pairwise
	                              rho1 <- 0 #rho1+as.numeric(cor(pair[,1],pair[,2])) #correlation way
	                              rho2 <- 0 #rho2 + sum(pair[,1]*pair[,2]) #moment est way 
	                              rho3 <- 0 #rho3 + ( sum(pair[,1]*pair[,2])-length(pair[,1])*mean(a_ij)*mean(a_ik) ) / ( sd(pair[,1]-mean(a_ij))*sd(pair[,2]-mean(a_ik)) )	#subtract section average
	                              rho4 <- 0 #rho4 + ( sum(pair[,1]*pair[,2])-length(pair[,1])*mean(a_i)*mean(a_i)^2 ) / ( sd(pair[,1]-mean(a_i))*sd(pair[,2]-mean(a_i)) )	#subtract school average
	     } #I fixed when data=instruction, 1 section in 1 school.
	     
	     else {
	       
       a_ij <- a_i[section[school == ublock[i]] == vblock[j]]  #scores for section j in school i
	   a_ik <- a_i[section[school == ublock[i]] == vblock[k]]
	   pair=as.matrix(pairup2(a_ij,a_ik)$ans)
	   nn <- nn+pairup2(a_ij,a_ik)$n	#number of pairwise
         
	   rho1 <- rho1+as.numeric(cor(as.numeric(pair[,1]), as.numeric(pair[,2]))) #correlation way
	   rho2 <- rho2 + sum(as.numeric(pair[,1])*as.numeric(pair[,2])) #moment est way 
	   rho3 <- rho3 + ( sum(as.numeric(pair[,1])*as.numeric(pair[,2]))-length(pair[,1])*mean(a_ij)*mean(a_ik) ) / ( sd(pair[,1]-mean(a_ij))*sd(pair[,2]-mean(a_ik)) )	#subtract section average
	   rho4 <- rho4 + ( sum(pair[,1]*pair[,2])-length(pair[,1])*mean(a_i)*mean(a_i)^2 ) / ( sd(pair[,1]-mean(a_i))*sd(pair[,2]-mean(a_i)) )	#subtract school average

	          }
	   }
	 }
	rho1_vec[i]=rho1
	rho2_vec[i]=rho2
	rho3_vec[i]=rho3
	rho4_vec[i]=rho4
	nvec[i] <- length(a_i) 
	npair[i]=nn
	}
list(rho1=rho1_vec, rho2=rho2_vec, rho3=rho3_vec, rho4=rho4_vec, nvec=nvec, npair=npair)
}

##intercept se

interc_se = function(x,ehat,school,section,taus){

	rho1 <- rhosect(sign(ehat),school,section)
#	c1 <- sum(sum(rho1$rho2))/(sum(sum(rho1$npair))-dim(x)[2])	#Simple Moment
	c1 <- sum(sum( (apply(rho1$rho2,1,sum,na.rm=T)/apply(rho1$npair,1,sum) ) * t(apply(rho1$mat,1,sum)) )) / sum(sum(rho1$mat)) # S Moment's weighted with sample size

	rho2 <- rhosch(sign(ehat),school,section)
#	c2 <- sum(sum( rho2$rho2 )) / (sum(sum(rho2$npair)-dim(x)[2]))	#Simple Moment
	c2 <- sum( ( (rho2$rho2 / rho2$npair ) *  rho2$nvec )/ sum(rho2$nvec) ,na.rm=T)	# weighted S Moment

	N <- sum(sum(rho1$mat))
#	a <- N+2*sum(rho1$npair)+2*sum(rho2$npair)
	a <- N^2-dim(x)[2]
	var_alpha <- (taus^2/a ) * (N + 2*sum(sum(rho1$npair))*c1 + 2*sum(rho2$npair)*c2 ) #? divided by N
	
	list(var_alpha=var_alpha)

	}



##wilc score##
scorewil = function(ehat){
#u is resid
    N=length(ehat)
    u=rank(ehat,ties.method=c("random"))/(N+1)
    scorewil = sqrt(12)*(u-.5)
    list(scorewil=scorewil) 
}


##pairup for a vector with its components for rho_sect##
pairup1=function(x,type="less") {
  x=as.matrix(x)
  n=dim(x)[1]
  i=rep(1:n,rep(n,n))
  j=rep(1:n,n)
  c1=apply(x,2,function(y){rep(y,rep(length(y),length(y)))})
  c2=apply(x,2,function(y){rep(y,length(y))})
  ans=cbind(c1,c2)
  ans=ans[(i<j), ]
  n=dim(ans)[1]
  list(ans=ans,n=n)
}

##two vectors, pairup all combinations for rho_sch##
pairup2=function(x,y) {
  x=as.matrix(x);y=as.matrix(y)
  n1=dim(x)[1];  n2=dim(y)[1]
  c1=rep(x,rep(n2,n1))
  c2=rep(y,n1)
  ans=cbind(as.numeric(c1),as.numeric(c2))
  list(ans=ans,n=n1*n2)
}

# ##point-wise correlation formula for two vectors##
# correl=function(x,y) {
# a=x%*%x
# b=y%*%y
# c=x%*%y
# n=length(x)
# xb=mean(x)
# yb=mean(y)
# rho=(c-n*xb*yb)/sqrt((a-n*xb^2)*(b-n*yb^2))
# list(rho=rho)
# }



##Simple Mixed Model in I schools, no section effect!##
#Kloke's function. Not used in the rlme
jrfit2 <-
function (x, y, block) 
{
    x <- as.matrix(x)
    p <- ncol(x)
    fit <- rfit(y ~ x, symmetric = FALSE)
    betahat <- fit$coef[2:(p + 1)]
    alphahat <- fit$coef[1]
    ehat <- fit$resid
    ahat <- scorewil(ehat) #?getScore??
    ublock <- unique(block)
    m <- length(ublock)
    nvec <- vector(m, mode = "numeric")
    totals <- total <- 0
    for (k in 1:m) {
        ak <- ahat[block == ublock[k]]
        ek <- ehat[block == ublock[k]]
        nvec[k] <- length(ak)
        for (i in 1:(nvec[k] - 1)) {
            for (j in (i + 1):nvec[k]) {
                total <- total + ak[i] * ak[j]
                totals <- total + sign(ek[i]) * sign(ek[j])
            }
        }
    }
    M <- sum(sum(choose(nvec, 2))) - p
    rho <- total/M
    rhos <- totals/M
    taus <- taustar(ehat, p)
    nstar <- sum(sum(nvec * (nvec - 1)))/sum(sum(nvec))
    sigmastar <- 1 + nstar * rhos
    names(alphahat) <- "Intercept"
    coef <- c(alphahat, betahat)
    res <- list(coefficients = coef, residuals = ehat, fitted.values = y - 
        ehat, rhohat = rho, tauhat = fit$tauhat, rhohats = rhos, 
        tauhats = taus, sigmastar = sigmastar, x = x, y = y, 
        block = block) #scores = scores
    class(res) <- list("jrfit")
    res
}



##wil one step
wilonestep = function(y,x){
#wwest using wil scores, yields fixed estimates(theta) abd sigmas
#or, use rfit for

	   n = length(y) 

         fitw = wwest(x,y,print.tbl=F)
         theta = fitw$tmp1$coef
         ehat = fitw$tmp1$residuals
    taus=taustar(ehat,p=dim(x)[2],conf=.95)
    tauhat=wilcoxontau(ehat,p=dim(x)[2])


     list(theta=theta,ehat=ehat,tauhat=tauhat,taus=taus)
}


# #4 graphs
# makepng=function(stresid.lmer1,yhat.lmer1,stanresid.gr,yhat.gr){
# 
# #dev.off()
# png("example_data/00.png")
# 
# trim=2
# par(mfrow = c(2, 2),font.main=1)
# 
# #plot(density(stresid.lmer1), main = "REML via LMER")
# 
# plot(stresid.lmer1 ~ yhat.lmer1, pch = "o", xlab="Fit", ylab="Standardized Residual", main="Stand. Residuals vs. Fits in REML" )
# abline(h = c(-trim, trim), col = "red")
# 
# qqnorm(stresid.lmer1 , pch = "o", main = "Normal Q-Q Plot")
# qqline(stresid.lmer1 , col = "red")
# 
# #dffits = abs(resid(lmer1, "dffits"))
# #plot(dffits, type = "h")
# 
# #plot(density(stanresid.gr), main = "Residuals in GR")
# plot(stanresid.gr ~ yhat.gr, pch = "o", xlab="Fit",ylab="Standardized Residual", main="Stand. Residuals vs. Fits in GR")
# abline(h = c(-trim, trim), col = "red")
# 
# qqnorm(stanresid.gr, pch = "o", main = "Normal Q-Q Plot")
# qqline(stanresid.gr, col = "red")
# 
# #title(main="Original Data", outer=T) #  with Outlier Med-MAD HL-DISP
# dev.off()
# 
# }


